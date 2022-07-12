#!/usr/bin/env nextflow

// ##############################################
// HELP SECTION
// ----------------------------------------------

if (params.help) {
  println """

  ------------------------------------------------------------------------------
  Gene prediction from raw illumina/pacbio/nanopore reads

  Written and maintained by Matteo Schiavinato
  Barcelona Supercomputing Center (BSC-CNS)
  2020-2022
  ------------------------------------------------------------------------------

  OPTIONS

  ### input ###

  --help                      This help message
  --threads                   N. of parallel threads for mapping step
  --input_data                TAB-separated file with the input sequencing data
                              (reads can be SE, PE, PB, ONT type)
  --paired_end                --input_data contains some PE samples
  --single_end                --input_data contains some SE samples
  --pacbio                    --input_data contains some PB samples
  --nanopore                  --input_data contains some ONT samples
  --output_dir                Create this directory and put everything in there
  --genome_fasta              Reference FASTA genome
  --adapters                  FASTA file w/ adapter seq. to use in PE/SE trimming
  --extrinsic                 "extrinsic.cfg" file for gene prediction parameters
  --repeats_gff               GFF file containing the annotation of repeats

                              Warning:
                              Repeat annotations must be of type "nonexonpart"
                              "nonexonpart" must be included in --extrinsic file

  ### gene prediction ###

  --gene_model                "complete" or "partial"                           [complete]
  --species_name              See augustus --species=help                       [!]
  --min_peptide_length        Min. peptide length to retain sequence            [10]
  --min_evidence              Min. evidence from RNASeq data to retain gene     [1]

  ### other parameters ###

  Edit directly in the nextflow.config file
  In most cases, defaults are fine

  You may have to generate your extrinsic.cfg file if your species is not in
  --species=help

  ------------------------------------------------------------------------------

  """
  exit 0
}


// -----------------------------------------------------------------------------
// READ INPUT DATA
// -----------------------------------------------------------------------------

// setting output directory
WD="${params.output_dir}"

// reading input data
Channel
  .fromPath("${params.input_data}")
  .splitCsv(header: true, sep: "\t")
  .map{ it -> [it.sample_id, it.read_type, it.reads_for, it.reads_rev] }
  .into{ Input_data_1; Input_data_2; Input_data_3; Input_data_4 }


// paired end reads

if (params.paired_end) {

  // form channel

  Input_data_1
    .filter{ it[1] == "PE" }
    .map{ it -> [it[0], file(it[2]), file(it[3])] }
    .set{ Input_data_PE }

  // trim reads

  process trim_PE_reads {

        executor = "local"
        cpus = 4
        maxForks = 16

        publishDir "${params.output_dir}/reads/fastq/PE", mode: "copy", pattern: "*.{P1,P2,U}.fastq"

        input:
          tuple val(sample_id), file(reads_for), file(reads_rev) from Input_data_PE

        output:
          tuple \
          val(sample_id), file("${sample_id}.P1.fastq"), \
          file("${sample_id}.P2.fastq"), file("${sample_id}.U.fastq") \
          into Trimmed_PE

        script:
          """
          ${TRIMMOMATIC} \
          PE \
          -threads 4 \
          ${reads_for} ${reads_rev} \
          ${sample_id}.P1.fastq ${sample_id}.U1.fastq \
          ${sample_id}.P2.fastq ${sample_id}.U2.fastq \
          ILLUMINACLIP:${params.adapters}:2:30:10 \
          LEADING:${params.leading} \
          TRAILING:${params.trailing} \
          SLIDINGWINDOW:${params.sliding_window} \
          AVGQUAL:${params.avgqual} \
          MINLEN:${params.minlen} &&
          cat ${sample_id}.U1.fastq ${sample_id}.U2.fastq \
          > ${sample_id}.U.fastq
          """

        stub:
          """
          touch ${sample_id}.P1.fastq ${sample_id}.P2.fastq ${sample_id}.U.fastq
          """
  }

  // map reads

  process map_PE_reads {

    executor = "local"
    cpus = params.threads
    maxForks = 1

    publishDir  "${params.output_dir}/mapping/raw/PE",
                mode: "copy", pattern: "*.sam"

    input:
      tuple \
      val(sample_id), file(reads_for), file(reads_rev), file(reads_U) \
      from Trimmed_PE

    output:
      tuple \
      val(sample_id), val("short"), file("${sample_id}.sam") \
      into Sams_PE

    script:
      """
      ${HISAT2_BUILD} -p ${params.threads} ${params.genome_fasta} ref &&
      ${HISAT2} -p ${params.threads} \
      --score-min L,0.0,-0.2 \
      -x ref \
      -1 ${reads_for} -2 ${reads_rev} -U ${reads_U} \
      -S ${sample_id}.sam
      """

    stub:
      """
      touch ${sample_id}.sam
      """
  }
} else {
  Sams_PE = Channel.empty()
}


// single end reads

if (params.single_end) {

  // form channel

  Input_data_2
    .filter{ it[1] == "SE" }
    .map{ it -> [it[0], file(it[2])] }
    .set{ Input_data_SE }

  // trim reads

  process trim_SE_reads {

        executor = "local"
        cpus = 4
        maxForks = 16

        publishDir "${params.output_dir}/reads/fastq/SE", mode: "copy", pattern: "*.se.fastq"

        input:
          tuple val(sample_id), file(reads) from Input_data_SE

        output:
          tuple \
          val(sample_id), file("${sample_id}.se.fastq") \
          into Trimmed_SE

        script:
          """
          ${TRIMMOMATIC} \
          SE \
          -threads 4 \
          ${reads} \
          ${sample_id}.se.fastq \
          ILLUMINACLIP:${params.adapters}:2:30:10 \
          LEADING:${params.leading} \
          TRAILING:${params.trailing} \
          SLIDINGWINDOW:${params.sliding_window} \
          AVGQUAL:${params.avgqual} \
          MINLEN:${params.minlen}
          """

        stub:
          """
          touch ${sample_id}.se.fastq
          """
  }

  // map reads

  process map_SE_reads {

    executor = "local"
    cpus = params.threads
    maxForks = 1

    publishDir  "${params.output_dir}/mapping/raw/SE",
                mode: "copy", pattern: "*.sam"

    input:
      tuple \
      val(sample_id), file(reads) \
      from Trimmed_SE

    output:
      tuple \
      val(sample_id), val("short"), file("${sample_id}.sam") \
      into Sams_SE

    script:
      """
      ${HISAT2_BUILD} -p ${params.threads} ${params.genome_fasta} ref &&
      ${HISAT2} -p ${params.threads} \
      --score-min L,0.0,-0.2 \
      -x ref \
      -U ${reads} \
      -S ${sample_id}.sam
      """

    stub:
      """
      touch ${sample_id}.sam
      """
  }
} else {
  Sams_SE = Channel.empty()
}


// nanopore reads

if (params.nanopore) {

  // form channel

  Input_data_3
    .filter{ it[1] == "ONT" }
    .map{ it -> [it[0], file(it[2])] }
    .set{ Input_data_ONT }

  // convert ONT reads to FASTA

  process ONT_to_FASTA {

    executor = "local"
    cpus = params.threads
    maxForks = 1

    publishDir  "${params.output_dir}/reads/fasta/ONT",
                mode: "copy", pattern: "*.fasta"

    input:
      tuple \
      val(sample_id), file(reads) \
      from Input_data_ONT

    output:
      tuple \
      val(sample_id), file("${sample_id}.fasta") \
      into Fasta_ONT

    script:
      """
      ${BIOAWK} -c fastx '{print \">\"\$name"\\n\"\$seq}' ${reads} \
      > ${sample_id}.fasta
      """

    stub:
      """
      touch ${sample_id}.fasta
      """
  }

  // ONT reads map

  process map_ONT_reads {

    executor = "local"
    cpus = params.threads
    maxForks = 1

    publishDir  "${params.output_dir}/mapping/raw/ONT",
                mode: "copy", pattern: "*.sam"

    input:
      tuple \
      val(sample_id), file(reads) \
      from Fasta_ONT

    output:
      tuple \
      val(sample_id), val("long"), file("${sample_id}.sam") \
      into Sams_ONT

    script:
      """
      ${MINIMAP2} \
      -x map-ont \
      -H \
      -a \
      -t ${params.threads} \
      -o ${sample_id}.sam \
      ${params.genome_fasta} \
      ${reads} \
      """

    stub:
      """
      touch ${sample_id}.sam
      """
  }
} else {
  Sams_ONT = Channel.empty()
}

// pacbio reads

if (params.pacbio) {

  // form channel

  Input_data_4
    .filter{ it[1] == "PB" }
    .map{ it -> [it[0], file(it[2])] }
    .set{ Input_data_PB }

  // convert PB reads to FASTA

  process PB_to_FASTA {

    executor = "local"
    cpus = params.threads
    maxForks = 1

    publishDir  "${params.output_dir}/reads/fasta/PB",
                mode: "copy", pattern: "*.fasta"

    input:
      tuple \
      val(sample_id), file(reads) \
      from Input_data_PB

    output:
      tuple \
      val(sample_id), file("${sample_id}.fasta") \
      into Fasta_PB

    script:
      """
      ${BIOAWK} \
      -c fastx \
      '{print \">\"\$name"\\n\"\$seq}' ${reads} \
      > ${sample_id}.fasta
      """
    stub:
      """
      touch ${sample_id}.fasta
      """
  }

  // PB reads map

  process map_PB_reads {

    executor = "local"
    cpus = params.threads
    maxForks = 1

    publishDir  "${params.output_dir}/mapping/raw/PB",
                mode: "copy", pattern: "*.sam"

    input:
      tuple \
      val(sample_id), file(reads) \
      from Fasta_PB

    output:
      tuple \
      val(sample_id), val("long"), file("${sample_id}.sam") \
      into Sams_PB

    script:
      """
      ${MINIMAP2} \
      -x map-pb \
      -H \
      -a \
      -t ${params.threads} \
      -o ${sample_id}.sam \
      ${params.genome_fasta} \
      ${reads} \
      """

    stub:
      """
      touch ${sample_id}.sam
      """
  }
} else {
  Sams_PB = Channel.empty()
}

// -----------------------------------------------------------------------------
// COMBINE AND PROCESS BAM FILES
// -----------------------------------------------------------------------------

// combine all bam files

Sams_PE
  .concat( Sams_SE, Sams_PB, Sams_ONT )
  .set{ Sams }


// filter bam files

process filter_and_sort {

  executor = "local"
  cpus = params.threads
  maxForks = params.threads

  publishDir  "${params.output_dir}/mapping/filtered",
              mode: "copy", pattern: "*.{fs.bam,fs.bam.bai}"

  input:
    tuple \
    val(sample_id), val(read_type), file(sam) \
    from Sams

  output:
    tuple \
    val(sample_id), val(read_type), file("${sample_id}.fs.bam"), file("${sample_id}.fs.bam.bai") \
    into Filt_bams

  script:
    """
    ${SAMTOOLS} view --threads ${params.threads} -F 0x0100 -F 0x4 -h -b \
    -o ${sample_id}.f.bam ${sam} &&
    ${SAMTOOLS} sort --threads ${params.threads} --output-fmt bam -T ${sample_id} \
    -o ${sample_id}.fs.bam \
    ${sample_id}.f.bam &&
    ${SAMTOOLS} index -@ ${params.threads} -b \
    ${sample_id}.fs.bam
    """

  stub:
    """
    touch ${sample_id}.f.bam ${sample_id}.fs.bam ${sample_id}.fs.bam.bai
    """
}

// resplit channels by long and short

Filt_bams.into{ Filt_bams_1; Filt_bams_2 }

Filt_bams_1
  .filter{ it[1] == "short" }
  .set{ Filt_bams_short }

Filt_bams_2
  .filter{ it[1] == "long" }
  .set{ Filt_bams_long }


// -----------------------------------------------------------------------------
// GENERATE HINTS
// -----------------------------------------------------------------------------

if (Filt_bams_long) {

  process bam2hints_long {

    executor = "local"
    cpus = 1
    maxForks = params.threads

    publishDir "${params.output_dir}/augustus/hints/long", mode: 'copy', pattern: '*.long.hints.gff'

    input:
      tuple val(sample_id), val(read_type), file(bam), file(bai) from Filt_bams_long

    output:
      tuple val(sample_id), val(read_type), file("${sample_id}.${read_type}.hints.gff") into Bam_Hints_long

    script:
      """
      bam2hints \
      --source=PB \
      --maxintronlen=${params.max_intron_size} \
      --priority=${params.long_read_priority} \
      --ep_cutoff=20 \
      --exonhints \
      --in=${bam} \
      --out=${sample_id}.${read_type}.hints.gff
      """

    stub:
      """
      touch ${sample_id}.${read_type}.hints.gff
      """
  }
} else {
  Bam_Hints_long = Channel.empty()
}


// BAM FILES

if (Filt_bams_short) {

  // split streams

  Filt_bams_short.into{ Filt_bams_intron; Filt_bams_exon }

  // get intron hints

  process bam2hints {

    executor = "local"
    cpus = 1
    maxForks = params.threads

    publishDir "${params.output_dir}/augustus/hints/short/bam2hints", mode: 'copy', pattern: '*.introns.gff'

    input:
      tuple val(sample_id), val(read_type), file(bam), file(bai) from Filt_bams_intron

    output:
      tuple val(sample_id), val(read_type), file("${sample_id}.${read_type}.introns.gff") into Bam2Hints

    script:
      """
      bam2hints \
      --intronsonly \
      --source=E \
      --priority=${params.bam_priority} \
      --maxintronlen=${params.max_intron_size} \
      --in=${bam} \
      --out=${sample_id}.${read_type}.introns.gff
      """

    stub:
      """
      touch ${sample_id}.${read_type}.introns.gff
      """
  }

  // get wig profiles

  process bam2wig {

    executor = "local"
    cpus = 1
    maxForks = params.threads

    publishDir "${params.output_dir}/augustus/hints/short/bam2wig", mode: 'copy', pattern: '*.exons.wig'

    input:
      tuple val(sample_id), val(read_type), file(bam), file(bai) from Filt_bams_exon

    output:
      tuple val(sample_id), val(read_type), file("${sample_id}.${read_type}.exons.wig") into Bam2Wig

    script:
      """
      bam2wig \
      ${bam} \
      > ${sample_id}.${read_type}.exons.wig
      """

    stub:
      """
      touch ${sample_id}.${read_type}.exons.wig
      """
  }


  // join channels based on baseName
  // output values in channel are like this:
  // [ baseName , file_gff , file_wig ]
  Bam2Hints
    .join(Bam2Wig, by: [0,1])
    .set { Bam_IntronWigs }


  process bam_noise_reduction {

    executor = "local"
    cpus = 1
    maxForks = params.threads

    publishDir "${params.output_dir}/augustus/hints/short/nr", mode: 'copy', pattern: '*.{introns,exons}.nr.{gff,wig}'

    input:
      tuple val(sample_id), val(read_type), file(intron), file(wig) from Bam_IntronWigs

    output:
      file "${sample_id}.${read_type}.introns.nr.gff" into Bam_Nr_Introns
      tuple val(sample_id), val(read_type), file("${sample_id}.${read_type}.exons.nr.wig") into Bam_Nr_Wigs

    script:
      """
      ${PERL} \
      ${projectDir}/src/RNA-seq-noise-reduction.v2.0.pl \
      ${intron} \
      ${wig} \
      ${sample_id}.${read_type}.introns.nr.gff \
      ${sample_id}.${read_type}.exons.nr.wig
      """

    stub:
      """
      touch ${sample_id}.${read_type}.introns.nr.gff
      touch ${sample_id}.${read_type}.exons.nr.wig
      """
  }

  process bam_wig2hints {

    executor = "local"
    cpus = 1
    maxForks = params.threads

    publishDir "${params.output_dir}/augustus/hints/short/wig2hints", mode: 'copy', pattern: '*.exons.nr.gff'

    input:
      tuple val(sample_id), val(read_type), file(wig) from Bam_Nr_Wigs

    output:
      file "${sample_id}.${read_type}.exons.nr.gff" into Bam_Nr_Exons

    script:
    """
    cat ${wig} | \
    ${PERL} \
    ${params.augustus_scripts}/wig2hints.pl \
    --width=10 --margin=10 --minthresh=2 \
    --minscore=4 --prune=0.1 --src=W --type=ep \
    --radius=4.5 --pri=${params.bam_priority} \
    --strand="." \
    --UCSC=${sample_id}.${read_type}.unstranded.track \
    > ${sample_id}.${read_type}.exons.nr.gff
    """

    stub:
      """
      touch ${sample_id}.${read_type}.exons.nr.gff
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

    publishDir "${params.output_dir}/augustus/hints/repeats", mode: 'copy', pattern: "repeat_hints.gff"

    output:
      file "repeat_hints.gff" into Repeat_hints

    script:
      """
      cat ${params.repeats_gff} | \
      ${AWK} '{print \$1"\t"\$2"\tnonexonpart\t"\$4"\t"\$5"\t"\$6"\t"\$7"\t"\$8"\t"\$9}' | \
      cut -f 1-8 | \
      ${AWK} '{print \$0"\tsrc=RM;pri=${params.rep_priority};"}' \
      > repeat_hints.gff
      """

    stub:
      """
      touch repeat_hints.gff
      """
  }
} else {
  Repeat_hints = Channel.empty()
}

// merge hints

Bam_Nr_Introns
  .concat(Bam_Nr_Exons, Bam_Hints_long, Repeat_hints)
  .set { IntronExon_Gffs }


// // generate hints file

process generate_hints_file {

  executor = "local"
  cpus = 1
  maxForks = params.threads

  publishDir "${params.output_dir}/augustus/hints", mode: 'copy', pattern: "merged.hints.gff"

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

  stub:
    """
    touch merged.hints.gff
    """
}

// split genome file

Channel
  .fromPath("${params.genome_fasta}")
  .set{ Fasta_file }

process split_genome {

  executor = "local"
  cpus = 1
  maxForks = params.threads

  publishDir "${params.output_dir}/augustus/split_genome", mode: 'copy', pattern: "*.[0-9].fa"

  input:
    file fasta from Fasta_file

  output:
    file "${fasta}.[0-9]*.fa" into Aug_genome

  script:
    """
    ${PYTHON3} \
    ${projectDir}/src/split_sequences.py \
    --input-fasta ${fasta} \
    --num-out-files ${params.threads}
    """

  stub:
    """
    ${PYTHON3} \
    ${projectDir}/src/split_sequences.py \
    --input-fasta ${fasta} \
    --num-out-files ${params.threads}
    """
  }


Aug_genome.flatten().set{ Aug_genome }


// Run Augustus
process run_augustus {

  executor = "local"
  cpus = 1
  maxForks = params.threads

  publishDir "${params.output_dir}/augustus/split_prediction", mode: 'copy', pattern: "*.gff"

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
   --UTR=${params.utr} \
   --gff3=on \
   --alternatives-from-evidence=true \
   --allow_hinted_splicesites=atac \
   ${genome_fasta} \
   > ${params.species_name}.${genome_fasta}.gff
   """

  stub:
    """
    touch ${params.species_name}.${genome_fasta}.gff
    """
}


process join_aug_pred {

  executor = "local"
  cpus = 1
  maxForks = params.threads

  publishDir "${params.output_dir}/augustus/proc.raw_gene_set",
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

  publishDir "${params.output_dir}/augustus/proc.raw_gene_set",
  mode: "copy"

  input:
    file aug_joined

  output:
    file "${params.species_name}.aa" into aug_aa
    file "${params.species_name}.cdsexons" into aug_cdsexons
    file "${params.species_name}.codingseq" into aug_codingseq

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

  publishDir "${params.output_dir}/augustus/proc.raw_gene_set",
  mode: "copy"

  input:
    file aug_joined

  output:
    file "${params.species_name}.evidence.gff" into aug_evidence

  script:
  """
  ${PYTHON3} \
  ${projectDir}/src/add_evidence_to_gff3_file.py \
  --input ${aug_joined} \
  --output ${params.species_name}.evidence.gff
  """
}



process get_ultrashort_peptides {

  executor = "local"
  cpus = 1
  maxForks = params.threads

  publishDir "${params.output_dir}/augustus/proc.raw_gene_set",
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

process filter_genes {

  executor = "local"
  cpus = 1
  maxForks = params.threads

  publishDir "${params.output_dir}/augustus",
  mode: "copy"

  input:
    file short_genes
    file aug_evidence

  output:
    file "${params.species_name}.clean.gff" into aug_clean
    file "removed_genes.txt" into removed_genes

  script:
    """
    cat ${short_genes} | \
    sort -k1.2n,1.20 | uniq \
    > removed_genes.txt &&
    ${PYTHON3} \
    ${projectDir}/src/filter_raw_gene_set.py \
    --input ${aug_evidence} \
    --remove-genes removed_genes.txt \
    --min-evidence ${params.min_evidence} \
    > ${params.species_name}.clean.gff
    """
}


process make_copy_for_gbrowse {

  executor = "local"
  cpus = 1
  maxForks = params.threads

  publishDir "${params.output_dir}/augustus",
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

  publishDir "${params.output_dir}/augustus",
  mode: "copy"

  input:
    file aug_clean
    file aug_aa
    file aug_codingseq

  output:
    file "${params.species_name}.clean.aa" into aug_clean_aa
    file "${params.species_name}.clean.codingseq" into aug_clean_codingseq

  script:
  """
  cat ${aug_clean} | \
  ${AWK} '\$3=="transcript"' | \
  cut -f 9 | \
  cut -d ";" -f 1 | \
  cut -d "=" -f 2 \
  > ${params.species_name}.clean.transcript_names.txt &&
  ${PYTHON3} \
  ${projectDir}/src/select-sequences-from-filename.py \
  -in ${aug_aa} \
  -out ${params.species_name}.clean.aa \
  -n ${params.species_name}.clean.transcript_names.txt &&
  ${PYTHON3} \
  ${projectDir}/src/select-sequences-from-filename.py \
  -in ${aug_codingseq} \
  -out ${params.species_name}.clean.codingseq \
  -n ${params.species_name}.clean.transcript_names.txt
  """
}


process filtering_stats {

  executor = "local"
  cpus = 1
  maxForks = params.threads

  publishDir "${params.output_dir}/augustus",
  mode: "copy",
  pattern: "*clean.stats.tsv"

  input:
    file aug_joined
    file aug_clean
    file aug_aa
    file aug_codingseq
    file aug_clean_aa
    file aug_clean_codingseq

  output:
    file "${params.species_name}.clean.stats.tsv" into aug_clean_stats

  script:
    """
    cat ${aug_joined} | egrep -v "^#" | ${AWK} '{if (\$3=="gene") {print \$9}}' | cut -d "=" -f 2 | sort | uniq | wc -l > raw.genes.count.txt &&
    cat ${aug_clean} | egrep -v "^#" | ${AWK} '{if (\$3=="gene") {print \$9}}' | cut -d "=" -f 2 | sort | uniq | wc -l > filt.genes.count.txt &&
    ${BIOAWK} -c fastx '{print \$name}' ${aug_aa} | wc -l > raw.aa.count.txt &&
    ${BIOAWK} -c fastx '{print \$name}' ${aug_clean_aa} | wc -l > filt.aa.count.txt &&
    ${BIOAWK} -c fastx '{print \$name}' ${aug_codingseq} | wc -l > raw.codingseq.count.txt &&
    ${BIOAWK} -c fastx '{print \$name}' ${aug_clean_codingseq} | wc -l > filt.codingseq.count.txt &&
    echo -e 'Feature\tRaw\tFiltered' > ${params.species_name}.clean.stats.tsv &&
    echo -e "Genes\t\$(cat raw.genes.count.txt)\t\$(cat filt.genes.count.txt)" >> ${params.species_name}.clean.stats.tsv &&
    echo -e "Proteins\t\$(cat raw.aa.count.txt)\t\$(cat filt.aa.count.txt)" >> ${params.species_name}.clean.stats.tsv
    """
}


process evidence_table {

  executor = "local"
  cpus = 1
  maxForks = params.threads

  publishDir "${params.output_dir}/augustus",
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
