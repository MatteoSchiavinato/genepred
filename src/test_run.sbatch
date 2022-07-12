#!/usr/bin/env sh
#SBATCH -N 1
#SBATCH -c 48
#SBATCH -J genepred

# ------------------------------------------------------------------------------
# EDIT FROM HERE
# ------------------------------------------------------------------------------

OUTDIR="/gpfs/projects/bsc40/project/pipelines/genepred/test_data/test_run"
GENOME_FASTA="/gpfs/projects/bsc40/project/pipelines/genepred/test_data/S_cerevisiae.fa"
SPECIES_NAME="saccharomyces_cerevisiae_S288C"
MAIN_FILE="/gpfs/projects/bsc40/project/pipelines/genepred/main.nf"
CONFIG_FILE="/gpfs/projects/bsc40/project/pipelines/genepred/nextflow.config"
INPUT_DATA="/gpfs/projects/bsc40/project/pipelines/genepred/test_data/input_test_data.tsv"

# ------------------------------------------------------------------------------
# STOP EDITING HERE
# ------------------------------------------------------------------------------

module load impi gcc intel mkl
module load gcc/7.2.0 hdf5/1.10.1 suitesparse/5.2.0 gsl/2.4 lpsolve/5.5.2.11
module load boost/1.68.0 htslib/1.15.1 samtools/1.15 bamtools/2.5.1 augustus/3.4.0
module load trimmomatic/0.39 hisat2/2.2.1 graphviz/2.40.1 bioawk/1.0 minimap2/2.24-r1122

if [ ! -d ${OUTDIR} ]; then mkdir ${OUTDIR}; fi

cd ${OUTDIR}
if [ ! -d nf ]; then mkdir nf; fi

nextflow \
-C ${CONFIG_FILE} \
run \
${MAIN_FILE} \
-resume \
-work-dir ${OUTDIR}/nf/nf-work \
-with-report ${OUTDIR}/nf/cmd.sbatch.report.html \
-with-timeline ${OUTDIR}/nf/cmd.sbatch.timeline.html \
-with-dag ${OUTDIR}/nf/cmd.sbatch.dag.png \
--input_data ${INPUT_DATA} \
--single_end --paired_end --nanopore \
--species_name ${SPECIES_NAME} \
--output_dir ${OUTDIR} \
--threads 40 \
--genome_fasta ${GENOME_FASTA} \
--gene_model complete \
--min_peptide_length 10 \
--rep_overlapping_feature exon \
--min_evidence 1 \
--utr off 
