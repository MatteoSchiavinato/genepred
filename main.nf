#!/usr/bin/env nextflow



if (params.help) {
  println """

  ------------------------------------------------------------------------------

  Matteo Schiavinato
  matteo.schiavinato.90@gmail.com
  December 2020 - May 2021
  Gene Prediction pipeline (originally for Beta vulgaris)

  ------------------------------------------------------------------------------

  WARNING: AUGUSTUS uses a lot of RAM
  Make sure you have approximately 20GB RAM x num. threads available on your PC

  ------------------------------------------------------------------------------

  OPTIONS

  ### miscellaneous ###

  --help                      This help message
  --output_directory          Create this directory and put everything in there
  --threads                   Number of parallel threads to use in mapping of each file (16)

  ### input ###

  --genome_fasta              Reference FASTA genome
  --repeats_gff               GFF file containing the annotation of repeats
                              The annotations must have been performed on the same FASTA file
                              declared in --genome_fasta
                              The feature type (3rd column) must be "nonexonpart"
                              This feature must be included in the extrinsic.cfg file that you pass
                              with --extrinsic

  --bam_long                  BAM mapping files from pacbio/nanopore (minimap2 + filterBam + sorting)
                              Will be used to calculate anchor hints

  --bam                       BAM mapping files (hisat2/bowtie2 + filterBam + sorting)
  --psl                       PSL mapping files (BLAT + filterPSL.pl + sorting)

  ### augustus ###

  --gene_model                "complete" or "partial" (complete)
  --species_name              Parameters to use from the Augustus config/species directory
                              (choose from: "augustus*/config/species/...")  --noise_reduction_script    Path to script to perform noise reduction. Get it from
                              Minoche et al., 2015 (https://doi.org/10.1186/s13059-015-0729-7)
                              In the supplementary material

  --augustus_bin              Path to the scripts for BAM hints preparation
                              (usually: "/bin" or "/auxprogs")
  --augustus_scripts          Path to the scripts for PSL hints preparation
  --augustus_config_path      Path to the /config/ folder inside the Augustus directory. This
                              folder contains the /species/ directory, inside which there has
                              to be a folder named exactly like the declared [--species_name]
  --extrinsic                 "extrinsic.cfg" file for gene prediction parameters
                              (usually: "augustus*/config/extrinsic/...")

  ### prediction filtering ###

  --rep_overlapping_feature   This feature in the GFF file will be checked to see if any of them
                              overlap any annotated repeat. Suggestion: "CDS" or "exon"
                              (Default: exon)

  --min_peptide_length        Minimum length of a peptide in order to retain it after the prediction
                              (Default: 10)

  --min_evidence              Minimum evidence from RNASeq data in order to retain a gene model
                              (Default: 1)

  ------------------------------------------------------------------------------

  """
  exit 0
}


// -----------------------------------------------------------------------------

// setting output directory
WD="${params.output_directory}"

// -----------------------------------------------------------------------------

// ##############################################

if (params.bam_long) {

  Channel.fromPath("${params.bam_long}/*.bam")
    .map { it -> [it.baseName, it] }
    .set { Filt_bams_long }

  process bam2hints_long {

    executor = "local"
    cpus = 1
    maxForks = params.threads

    // executor = "${params.executor}"
    // module = []
    // clusterOptions = "-N 1 --ntasks-per-node 16 --account ${params.account} --partition ${params.partition} --qos ${params.qos}"

    publishDir "${WD}/process/evidence", mode: 'copy', pattern: '*.bam.long.hints.gff'

    input:
      tuple val(baseName), file(query) from Filt_bams_long

    output:
      file "${baseName}.bam.long.hints.gff" into Bam_Hints_long

    script:
      """
      ${params.augustus_bin}/bam2hints \
      --source=PB \
      --maxintronlen=${params.max_intron_size} \
      --priority=${params.long_read_priority} \
      --ep_cutoff=20 \
      --exonhints \
      --in=${query} \
      --out=${baseName}.bam.long.hints.gff
      """
  }
} else {
  Bam_Hints_long = Channel.empty()
}


if (params.psl) {

  Channel.fromPath("${params.psl}/*.psl")
    .map { it -> [it.baseName, it] }
    .set { Filt_psls_intron }

  Channel.fromPath("${params.psl}/*.psl")
    .map { it -> [it.baseName, it] }
    .set { Filt_psls_exon }

  process blat2hints {

    executor = "local"
    cpus = 1
    maxForks = params.threads

    // executor = "${params.executor}"
    // module = []
    // clusterOptions = "-N 1 --ntasks-per-node 16 --account ${params.account} --partition ${params.partition} --qos ${params.qos}"

    publishDir "${WD}/process/blat2hints", mode: 'copy', pattern: '*.psl.introns.gff'

    input:
      tuple val(baseName), file(psl) from Filt_psls_intron

    output:
      tuple val(baseName), file("${baseName}.psl.introns.gff") into Blat2Hints

    script:
      """
      ${PERL} \
      ${params.augustus_scripts}/blat2hints.pl \
      --intronsonly \
      --source=E \
      --priority=${params.psl_priority} \
      --maxintronlen=${params.max_intron_size} \
      --in=${psl} \
      --out=${baseName}.psl.introns.gff
      """
  }


  process aln2wig {

    executor = "local"
    cpus = 1
    maxForks = params.threads

    // executor = "${params.executor}"
    // module = []
    // clusterOptions = "-N 1 --ntasks-per-node 16 --account ${params.account} --partition ${params.partition} --qos ${params.qos}"

    publishDir "${WD}/process/aln2wig", mode: 'copy', pattern: '*.psl.exons.wig'

    input:
      tuple val(baseName), file(psl) from Filt_psls_exon

    output:
      tuple val(baseName), file("${baseName}.psl.exons.wig") into Aln2Wig

    script:
      """
      ${params.augustus_bin}/aln2wig \
      -f ${psl} \
      > ${baseName}.psl.exons.wig
      """
  }

  // join channels based on baseName
  // output values in channel are like this:
  // [ baseName , file_gff , file_wig ]
  Blat2Hints
    .join(Aln2Wig)
    .set { Psl_IntronWigs }


  process psl_noise_reduction {

    executor = "local"
    cpus = 1
    maxForks = params.threads

    // executor = "${params.executor}"
    // module = []
    // clusterOptions = "-N 1 --ntasks-per-node 16 --account ${params.account} --partition ${params.partition} --qos ${params.qos}"

    publishDir "${WD}/process/evidence", mode: 'copy', pattern: '*.psl.{introns,exons}.nr.{gff,wig}'

    input:
      tuple val(baseName), file(intron), file(wig) from Psl_IntronWigs

    output:
      file "${baseName}.psl.introns.nr.gff" into Psl_Nr_Introns
      tuple val(baseName), file("${baseName}.psl.exons.nr.wig") into Psl_Nr_Wigs

    script:
      """
      ${PERL} \
      ${params.noise_reduction_script} \
      ${intron} \
      ${wig} \
      ${baseName}.psl.introns.nr.gff \
      ${baseName}.psl.exons.nr.wig
      """
  }


  process psl_wig2hints {

    executor = "local"
    cpus = 1
    maxForks = params.threads

    // executor = "${params.executor}"
    // module = []
    // clusterOptions = "-N 1 --ntasks-per-node 16 --account ${params.account} --partition ${params.partition} --qos ${params.qos}"

    publishDir "${WD}/process/evidence", mode: 'copy', pattern: '*.psl.exons.nr.gff'

    input:
      tuple val(baseName), file(wig) from Psl_Nr_Wigs

    output:
      file "${baseName}.psl.exons.nr.gff" into Psl_Nr_Exons

    script:
    """
    cat ${wig} | \
    ${PERL} \
    ${params.augustus_scripts}/wig2hints.pl \
    --width=10 --margin=10 --minthresh=2 \
    --minscore=4 --prune=0.1 --src=W --type=ep \
    --radius=4.5 --pri=${params.psl_priority} \
    --strand="." \
    --UCSC=${baseName}.unstranded.track \
    > ${baseName}.psl.exons.nr.gff
    """
  }

} else {
  Psl_Nr_Exons = Channel.empty()
  Psl_Nr_Introns = Channel.empty()
}


// BAM FILES

if (params.bam) {

  Channel.fromPath("${params.bam}/*.bam")
    .map { it -> [it.baseName, it] }
    .set { Filt_bams_intron }

  Channel.fromPath("${params.bam}/*.bam")
    .map { it -> [it.baseName, it] }
    .set { Filt_bams_exon }

  process bam2hints {

    executor = "local"
    cpus = 1
    maxForks = params.threads

    // executor = "${params.executor}"
    // module = []
    // clusterOptions = "-N 1 --ntasks-per-node 16 --account ${params.account} --partition ${params.partition} --qos ${params.qos}"

    publishDir "${WD}/process/bam2hints", mode: 'copy', pattern: '*.bam.introns.gff'

    input:
      tuple val(baseName), file(bam) from Filt_bams_intron

    output:
      tuple val(baseName), file("${baseName}.bam.introns.gff") into Bam2Hints

    script:
      """
      ${params.augustus_bin}/bam2hints \
      --intronsonly \
      --source=E \
      --priority=${params.bam_priority} \
      --maxintronlen=${params.max_intron_size} \
      --in=${bam} \
      --out=${baseName}.bam.introns.gff
      """
  }


  process bam2wig {

    executor = "local"
    cpus = 1
    maxForks = params.threads

    // executor = "${params.executor}"
    // module = []
    // clusterOptions = "-N 1 --ntasks-per-node 16 --account ${params.account} --partition ${params.partition} --qos ${params.qos}"

    publishDir "${WD}/process/bam2wig", mode: 'copy', pattern: '*.bam.exons.wig'

    input:
      tuple val(baseName), file(bam) from Filt_bams_exon

    output:
      tuple val(baseName), file("${baseName}.bam.exons.wig") into Bam2Wig

    script:
      """
      ${params.augustus_bin}/bam2wig \
      ${bam} \
      > ${baseName}.bam.exons.wig
      """
  }


  // join channels based on baseName
  // output values in channel are like this:
  // [ baseName , file_gff , file_wig ]
  Bam2Hints
    .join(Bam2Wig)
    .set { Bam_IntronWigs }


  process bam_noise_reduction {

    executor = "local"
    cpus = 1
    maxForks = params.threads

    // executor = "${params.executor}"
    // module = []
    // clusterOptions = "-N 1 --ntasks-per-node 16 --account ${params.account} --partition ${params.partition} --qos ${params.qos}"

    publishDir "${WD}/process/evidence", mode: 'copy', pattern: '*.bam.{introns,exons}.nr.{gff,wig}'

    input:
      tuple val(baseName), file(intron), file(wig) from Bam_IntronWigs

    output:
      file "${baseName}.bam.introns.nr.gff" into Bam_Nr_Introns
      tuple val(baseName), file("${baseName}.bam.exons.nr.wig") into Bam_Nr_Wigs

    script:
      """
      ${PERL} \
      ${params.noise_reduction_script} \
      ${intron} \
      ${wig} \
      ${baseName}.bam.introns.nr.gff \
      ${baseName}.bam.exons.nr.wig
      """
  }

  process bam_wig2hints {

    executor = "local"
    cpus = 1
    maxForks = params.threads

    // executor = "${params.executor}"
    // module = []
    // clusterOptions = "-N 1 --ntasks-per-node 16 --account ${params.account} --partition ${params.partition} --qos ${params.qos}"

    publishDir "${WD}/process/evidence", mode: 'copy', pattern: '*.bam.exons.nr.gff'

    input:
      tuple val(baseName), file(wig) from Bam_Nr_Wigs

    output:
      file "${baseName}.bam.exons.nr.gff" into Bam_Nr_Exons

    script:
    """
    cat ${wig} | \
    ${PERL} \
    ${params.augustus_scripts}/wig2hints.pl \
    --width=10 --margin=10 --minthresh=2 \
    --minscore=4 --prune=0.1 --src=W --type=ep \
    --radius=4.5 --pri=${params.bam_priority} \
    --strand="." \
    --UCSC=${baseName}.unstranded.track \
    > ${baseName}.bam.exons.nr.gff
    """
  }
} else {
  Bam_Nr_Introns = Channel.empty()
  Bam_Nr_Exons = Channel.empty()
}


// repeats

if (params.repeats_gff) {
  process repeat_hints {

    executor = "local"
    cpus = 1
    maxForks = params.threads

    // executor = "${params.executor}"
    // module = []
    // clusterOptions = "-N 1 --ntasks-per-node 16 --account ${params.account} --partition ${params.partition} --qos ${params.qos}"

    publishDir "${WD}/process/evidence", mode: 'copy', pattern: "repeat_hints.gff"

    input:
      file params.repeats_gff

    output:
      file "repeat_hints.gff" into repeat_hints

    script:
      """
      cat ${params.repeats_gff} | \
      ${AWK} '{print \$1"\t"\$2"\tnonexonpart\t"\$4"\t"\$5"\t"\$6"\t"\$7"\t"\$8"\t"\$9}' | \
      cut -f 1-8 | \
      ${AWK} '{print \$0"\tsrc=RM;pri=${params.rep_priority};"}' \
      > repeat_hints.gff
      """
  }
}

Bam_Nr_Introns
  .concat(Bam_Nr_Exons, Psl_Nr_Introns, Psl_Nr_Exons, Bam_Hints_long)
  .set { IntronExon_Gffs }


// merge intron and exon hints
if (params.repeats_gff) {

    process generate_hints_file_w_rep {

      executor = "local"
      cpus = 1
      maxForks = params.threads

      // executor = "${params.executor}"
      // module = []
      // clusterOptions = "-N 1 --ntasks-per-node 16 --account ${params.account} --partition ${params.partition} --qos ${params.qos}"

      publishDir "${WD}/augustus/hints", mode: 'copy', pattern: "merged.hints.gff"

      input:
        file repeat_hints
        file all_gffs from IntronExon_Gffs.collect()

      output:
        file "merged.hints.gff" into hints_file

      script:
        """
        cat ${repeat_hints} ${all_gffs} | \
        tr -s " " "\t" | tr -s "\t" "\t" | \
        sort -k1,1 -k4n,4 -k5nr,5 -T . \
        > merged.hints.gff
        """
    }

} else {

    process generate_hints_file {

      executor = "local"
      cpus = 1
      maxForks = params.threads

      // executor = "${params.executor}"
      // module = []
      // clusterOptions = "-N 1 --ntasks-per-node 16 --account ${params.account} --partition ${params.partition} --qos ${params.qos}"

      publishDir "${WD}/augustus/hints", mode: 'copy', pattern: "merged.hints.gff"

      input:
        file all_gffs from IntronExon_Gffs.collect()

      output:
        file "merged.hints.gff" into hints_file

      script:
        """
        cat ${all_gffs} | \
        tr -s " " "\t" | tr -s "\t" "\t" | \
        sort -k1,1 -k4n,4 -k5nr,5 -T . \
        > merged.hints.gff
        """
    }
}


// split genome file
Channel
  .fromPath("${params.genome_fasta}")
  .set{ Fasta_file }

process split_genome {

  executor = "local"
  cpus = 1
  maxForks = params.threads

  // executor = "${params.executor}"
  // module = []
  // clusterOptions = "-N 1 --ntasks-per-node 16 --account ${params.account} --partition ${params.partition} --qos ${params.qos}"

  input:
    file fasta from Fasta_file

  output:
    file "*.part_[0-9]*.fa" into Genome_ch

  script:
    """
    ${PYTHON3} \
    ${params.source_dir}/split_sequences.py \
    --infile ${fasta} \
    --num-files ${params.threads}
    """
  }

// convert list of files to channel
Genome_ch.flatMap().set{ Aug_genome }

// Run Augustus
process run_augustus {

  executor = "local"
  cpus = 1
  maxForks = params.threads

  // executor = "${params.executor}"
  // module = []
  // clusterOptions = "-N 1 --ntasks-per-node 16 --account ${params.account} --partition ${params.partition} --qos ${params.qos}"

  publishDir "${WD}/process/split_prediction",
  mode: 'copy',
  pattern: '*.gff'

  input:
   file hints_file
   file genome_fasta from Aug_genome

  output:
   file "${params.species_name}.${genome_fasta}.gff" into Augustus_gff

  script:
   """
   ${AUGUSTUS} \
   --species=${params.species_name} \
   --strand=both \
   --genemodel=${params.gene_model} \
   --AUGUSTUS_CONFIG_PATH=${params.augustus_config_path} \
   --hintsfile=${hints_file} \
   --extrinsicCfgFile=${params.extrinsic} \
   --UTR=on \
   --gff3=on \
   --alternatives-from-evidence=true \
   --allow_hinted_splicesites=atac \
   ${genome_fasta} \
   > ${params.species_name}.${genome_fasta}.gff
   """
}


process join_aug_pred {

  executor = "local"
  cpus = 1
  maxForks = params.threads

  // executor = "${params.executor}"
  // module = []
  // clusterOptions = "-N 1 --ntasks-per-node 16 --account ${params.account} --partition ${params.partition} --qos ${params.qos}"

  publishDir "${WD}/augustus/proc.raw_gene_set",
  mode: "copy",
  pattern: "*gff"

  input:
    file all_augustus_gffs from Augustus_gff.collect()

  output:
    file "${params.species_name}.gff" into aug_joined

  script:
  """
  cat ${all_augustus_gffs} | \
  ${PERL} \
  ${params.augustus_scripts}/join_aug_pred.pl \
  > ${params.species_name}.gff
  """
}


process get_anno_fasta {

  executor = "local"
  cpus = 1
  maxForks = params.threads

  // executor = "${params.executor}"
  // module = []
  // clusterOptions = "-N 1 --ntasks-per-node 16 --account ${params.account} --partition ${params.partition} --qos ${params.qos}"

  publishDir "${WD}/augustus/proc.raw_gene_set",
  mode: "copy"

  input:
    file aug_joined

  output:
    file "${params.species_name}.aa" into aug_aa
    file "${params.species_name}.cdsexons" into aug_cdsexons
    file "${params.species_name}.codingseq" into aug_codingseq
    file "${params.species_name}.mrna" into aug_mrna

  script:
  """
  ${PERL} \
  ${params.augustus_scripts}/getAnnoFasta.pl \
  --seqfile=${params.genome_fasta} \
  ${aug_joined}
  """
}


process add_evidence_to_gff3 {

  executor = "local"
  cpus = 1
  maxForks = params.threads

  // executor = "${params.executor}"
  // module = []
  // clusterOptions = "-N 1 --ntasks-per-node 16 --account ${params.account} --partition ${params.partition} --qos ${params.qos}"

  publishDir "${WD}/augustus/proc.raw_gene_set",
  mode: "copy"

  input:
    file aug_joined

  output:
    file "${params.species_name}.evidence.gff" into aug_evidence

  script:
  """
  ${PYTHON3} \
  ${launchDir}/src/add_evidence_to_gff3_file.py \
  --input ${aug_joined} \
  --output ${params.species_name}.evidence.gff
  """
}


process get_ultrashort_peptides {

  executor = "local"
  cpus = 1
  maxForks = params.threads

  // executor = "${params.executor}"
  // module = []
  // clusterOptions = "-N 1 --ntasks-per-node 16 --account ${params.account} --partition ${params.partition} --qos ${params.qos}"

  publishDir "${WD}/augustus/proc.raw_gene_set",
  mode: "copy"

  input:
    file aug_aa

  output:
    file "list.ultrashort_peptides.txt" into short_genes

  script:
    """
    ${BIOAWK} \
    -c fastx \
    '{if (length(\$seq) < ${params.min_peptide_length}) {print \$name}}' \
    ${aug_aa} | \
    sed -e 's/.t[0-9]*//' | \
    sort -k1.2n,1.20 | uniq \
    > list.ultrashort_peptides.txt
    """
}


process genes_over_repeats {

  executor = "local"
  cpus = 1
  maxForks = params.threads

  // executor = "${params.executor}"
  // module = []
  // clusterOptions = "-N 1 --ntasks-per-node 16 --account ${params.account} --partition ${params.partition} --qos ${params.qos}"

  publishDir "${WD}/augustus/proc.raw_gene_set",
  mode: "copy",
  pattern: "*.txt"

  input:
    file aug_joined
    file repeat_hints

  output:
    file "list.rep_overlapping_genes.txt" into rep_genes

  script:
  """
  cat ${aug_joined} | \
  ${AWK} '\$3=="${params.rep_overlapping_feature}"' \
  > features.gff &&
  ${BEDTOOLS} intersect \
  -u \
  -a features.gff \
  -b ${repeat_hints} | \
  cut -f 9 | \
  sed -e 's/.*g/g/' | \
  sed -e 's/.t.*//' | \
  sort | \
  uniq \
  > list.rep_overlapping_genes.txt
  """
}


process filter_genes {

  executor = "local"
  cpus = 1
  maxForks = params.threads

  // executor = "${params.executor}"
  // module = []
  // clusterOptions = "-N 1 --ntasks-per-node 16 --account ${params.account} --partition ${params.partition} --qos ${params.qos}"

  publishDir "${WD}/augustus",
  mode: "copy"

  input:
    file rep_genes
    file short_genes
    file aug_evidence

  output:
    file "${params.species_name}.clean.gff" into aug_clean
    file "removed_genes.txt" into removed_genes

  script:
    """
    cat ${short_genes} ${rep_genes} | \
    sort -k1.2n,1.20 | uniq \
    > removed_genes.txt &&
    ${PYTHON3} \
    ${launchDir}/src/filter_raw_gene_set.py \
    --input ${aug_evidence} \
    --remove-genes removed_genes.txt \
    --evidence ${params.min_evidence} \
    > ${params.species_name}.clean.gff
    """
}


process make_copy_for_gbrowse {

  executor = "local"
  cpus = 1
  maxForks = params.threads

  // executor = "${params.executor}"
  // module = []
  // clusterOptions = "-N 1 --ntasks-per-node 16 --account ${params.account} --partition ${params.partition} --qos ${params.qos}"

  publishDir "${WD}/augustus",
  mode: "copy",
  pattern: "*gff"

  input:
    file aug_clean

  output:
    file "${params.species_name}.clean.gbrowse_formatted.gff" into aug_clean_gbrowse

  script:
    """
    cat ${aug_clean} | \
    ${PERL} \
    ${params.augustus_scripts}/augustus2gbrowse.pl \
    > ${params.species_name}.clean.gbrowse_formatted.gff
    """
}


process get_clean_anno_fasta {

  executor = "local"
  cpus = 1
  maxForks = params.threads

  // executor = "${params.executor}"
  // module = []
  // clusterOptions = "-N 1 --ntasks-per-node 16 --account ${params.account} --partition ${params.partition} --qos ${params.qos}"

  publishDir "${WD}/augustus",
  mode: "copy"

  input:
    file aug_clean
    file aug_aa
    file aug_mrna
    file aug_codingseq

  output:
    file "${params.species_name}.clean.aa" into aug_clean_aa
    file "${params.species_name}.clean.codingseq" into aug_clean_codingseq
    file "${params.species_name}.clean.mrna" into aug_clean_mrna

  script:
  """
  cat ${aug_clean} | \
  ${AWK} '\$3=="transcript"' | \
  cut -f 9 | \
  cut -d ";" -f 1 | \
  cut -d "=" -f 2 \
  > ${params.species_name}.clean.transcript_names.txt &&
  ${PYTHON3} \
  ${launchDir}/src/select-sequences-from-filename.py \
  -in ${aug_aa} \
  -out ${params.species_name}.clean.aa \
  -n ${params.species_name}.clean.transcript_names.txt &&
  ${PYTHON3} \
  ${launchDir}/src/select-sequences-from-filename.py \
  -in ${aug_mrna} \
  -out ${params.species_name}.clean.mrna \
  -n ${params.species_name}.clean.transcript_names.txt &&
  ${PYTHON3} \
  ${launchDir}/src/select-sequences-from-filename.py \
  -in ${aug_codingseq} \
  -out ${params.species_name}.clean.codingseq \
  -n ${params.species_name}.clean.transcript_names.txt
  """
}


process filtering_stats {

  executor = "local"
  cpus = 1
  maxForks = params.threads

  // executor = "${params.executor}"
  // module = []
  // clusterOptions = "-N 1 --ntasks-per-node 16 --account ${params.account} --partition ${params.partition} --qos ${params.qos}"

  publishDir "${WD}/augustus",
  mode: "copy",
  pattern: "*clean.stats.tsv"

  input:
    file aug_joined
    file aug_clean
    file aug_aa
    file aug_mrna
    file aug_codingseq
    file aug_clean_aa
    file aug_clean_mrna
    file aug_clean_codingseq

  output:
    file "${params.species_name}.clean.stats.tsv" into aug_clean_stats

  script:
    """
    cat ${aug_joined} | egrep -v "^#" | ${AWK} '{if (\$3=="gene") {print \$9}}' | cut -d "=" -f 2 | sort | uniq | wc -l > raw.genes.count.txt &&
    cat ${aug_clean} | egrep -v "^#" | ${AWK} '{if (\$3=="gene") {print \$9}}' | cut -d "=" -f 2 | sort | uniq | wc -l > filt.genes.count.txt &&
    ${BIOAWK} -c fastx '{print \$name}' ${aug_aa} | wc -l > raw.aa.count.txt &&
    ${BIOAWK} -c fastx '{print \$name}' ${aug_clean_aa} | wc -l > filt.aa.count.txt &&
    ${BIOAWK} -c fastx '{print \$name}' ${aug_mrna} | wc -l > raw.mrna.count.txt &&
    ${BIOAWK} -c fastx '{print \$name}' ${aug_clean_mrna} | wc -l > filt.mrna.count.txt &&
    ${BIOAWK} -c fastx '{print \$name}' ${aug_codingseq} | wc -l > raw.codingseq.count.txt &&
    ${BIOAWK} -c fastx '{print \$name}' ${aug_clean_codingseq} | wc -l > filt.codingseq.count.txt &&
    echo -e 'Feature\tRaw\tFiltered' > ${params.species_name}.clean.stats.tsv &&
    echo -e "Genes\t\$(cat raw.genes.count.txt)\t\$(cat filt.genes.count.txt)" >> ${params.species_name}.clean.stats.tsv &&
    echo -e "Transcripts\t\$(cat raw.mrna.count.txt)\t\$(cat filt.mrna.count.txt)" >> ${params.species_name}.clean.stats.tsv &&
    echo -e "Proteins\t\$(cat raw.aa.count.txt)\t\$(cat filt.aa.count.txt)" >> ${params.species_name}.clean.stats.tsv
    """
}


process evidence_table {

  executor = "local"
  cpus = 1
  maxForks = params.threads

  // executor = "${params.executor}"
  // module = []
  // clusterOptions = "-N 1 --ntasks-per-node 16 --account ${params.account} --partition ${params.partition} --qos ${params.qos}"

  publishDir "${WD}/augustus",
  mode: "copy",
  pattern: "*clean.stats.tsv"

  input:
    file aug_evidence
    file aug_clean

  output:
    file "${params.species_name}.clean.evidence.tsv" into aug_evidence_stats

  script:
  """
  cat ${aug_evidence} | awk '\$3=="transcript"' | cut -f 9 | cut -d "_" -f 2 | tr -d "%" | awk '\$1 == 100' | wc -l > raw.100.txt
  cat ${aug_evidence} | awk '\$3=="transcript"' | cut -f 9 | cut -d "_" -f 2 | tr -d "%" | awk '\$1 < 100 && \$1 >= 1' | wc -l > raw.1-99.txt
  cat ${aug_evidence} | awk '\$3=="transcript"' | cut -f 9 | cut -d "_" -f 2 | tr -d "%" | awk '\$1 < 1' | wc -l > raw.0.txt
  cat ${aug_clean} | awk '\$3=="transcript"' | cut -f 9 | cut -d "_" -f 2 | tr -d "%" | awk '\$1 == 100' | wc -l > filt.100.txt
  cat ${aug_clean} | awk '\$3=="transcript"' | cut -f 9 | cut -d "_" -f 2 | tr -d "%" | awk '\$1 < 100 && \$1 >= 1' | wc -l > filt.1-99.txt
  cat ${aug_clean} | awk '\$3=="transcript"' | cut -f 9 | cut -d "_" -f 2 | tr -d "%" | awk '\$1 < 1' | wc -l > filt.0.txt
  echo -e 'Evidence\tRaw\tFiltered' > ${params.species_name}.clean.evidence.tsv &&
  echo -e "100%\t\$(cat raw.100.txt)\t\$(cat filt.100.txt)" >> ${params.species_name}.clean.evidence.tsv &&
  echo -e "1-99%\t\$(cat raw.1-99.txt)\t\$(cat filt.1-99.txt)" >> ${params.species_name}.clean.evidence.tsv &&
  echo -e "0%\t\$(cat raw.0.txt)\t\$(cat filt.0.txt)" >> ${params.species_name}.clean.evidence.tsv
  """
}
