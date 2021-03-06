#!/usr/bin/env sh
#SBATCH -N 1
#SBATCH -c 48
#SBATCH -J genepred

# ------------------------------------------------------------------------------
# HOW TO RUN
# ------------------------------------------------------------------------------

#  1. Login through the mn0 node (with internet), run module load nextflow/21.04.1
#  2. Run "nextflow run hello" to set up nextflow (only on first usage, requires internet)
#  3. Create an output directory in your "/gpfs/projects/bsc40/current" folder
#  4. Copy this file to the output directory you just created
#  5. Put its path as $OUTDIR here below
#  6. Copy the "src/input_data_example.tsv" file inside $OUTDIR, and edit it with your files
#  7. Put its path as $INPUT_DATA here below
#  8. Edit the $READ_TYPES variable leaving only the read types you used in $INPUT_DATA
#  9. Put the path to the genome fasta you want to use in the prediction as $GENOME_FASTA here below
# 10. Select the species you are working with from "ls config/species"
# 11. Put the species folder name (e.g. MONTSE19P8) as $SPECIES_NAME here below
# 12. Copy the "main.nf" and the "nextflow.config" files to $OUTDIR
# 13. Put their path as $MAIN_FILE and $CONFIG_FILE here below

# ------------------------------------------------------------------------------
# EDIT FROM HERE
# ------------------------------------------------------------------------------

OUTDIR="path to an existing output directory"
INPUT_DATA="copy the src/input_data_example.tsv table to your folders, edit it, and pass the path here"
READ_TYPES="--single_end --paired_end --pacbio --nanopore"
GENOME_FASTA="path to an existing FASTA genome sequence"
SPECIES_NAME="see augustus --species=help"
MAIN_FILE="path to your copied main.nf"
CONFIG_FILE="path to your copied nextflow.config"

# ------------------------------------------------------------------------------
# STOP EDITING HERE
# ------------------------------------------------------------------------------

module purge
module load nextflow/21.04.1
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
${READ_TYPES} \
--species_name ${SPECIES_NAME} \
--output_dir ${OUTDIR} \
--threads 40 \
--genome_fasta ${GENOME_FASTA} \
--gene_model complete \
--min_peptide_length 10 \
--rep_overlapping_feature exon \
--min_evidence 1 \
--utr off
