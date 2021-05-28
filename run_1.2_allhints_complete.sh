#!/usr/bin/env sh

WD="/gpfs/data/fs71579/schmat90/CF/projects/genomics/G00002_Beet_gene_prediction"
cd ${WD}
OUTDIR="RefBeet_1.2_allhints_complete"
if [ ! -d ${OUTDIR} ] ; then mkdir ${OUTDIR} ; fi

cd ${WD}/scripts

echo """\
#!/usr/bin/env sh
#SBATCH -N 1
#SBATCH --ntasks-per-node 10
#SBATCH -J 1.2_ALL
#SBATCH --account p71579
#SBATCH --partition mem_0096
#SBATCH --qos mem_0096
#SBATCH --mail-type FAIL,END
#SBATCH --mail-user matteo.schiavinato@boku.ac.at

cd ${WD}/scripts

nextflow \
main.nf \
-resume \
-work-dir ${WD}/${OUTDIR}/work \
-with-report ${WD}/${OUTDIR}/cmd.sbatch.report.html \
-with-timeline ${WD}/${OUTDIR}/cmd.sbatch.timeline.html \
-with-dag ${WD}/${OUTDIR}/cmd.sbatch.dag.png \
--account p71579 \
--partition mem_0096 \
--qos mem_0096 \
--executor slurm \
--species_name beta_vulgaris \
--output_directory ${WD}/${OUTDIR} \
--threads 10 \
--genome_fasta ${WD}/raw_data/RefBeet_1.2/RefBeet-1.2.masked.fa \
--repeats_gff ${WD}/raw_data/RefBeet_1.2/RefBeet-1.2.masked.gff \
--gene_model complete \
--psl ${WD}/RefBeet_1.2_run/filter_PSL \
--bam ${WD}/RefBeet_1.2_run/filter_BAM \
--bam_long ${WD}/RefBeet_1.2_run/filter_BAM_long \
--min_peptide_length 10 \
--rep_overlapping_feature exon \
--min_evidence 1 \

""" > ${WD}/${OUTDIR}/cmd.sbatch

sbatch \
--output ${WD}/${OUTDIR}/cmd.sbatch.out \
${WD}/${OUTDIR}/cmd.sbatch
