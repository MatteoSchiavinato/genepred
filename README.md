# genepred
A Nextflow pipeline for gene prediction and functional annotation
Matteo Schiavinato
2020-2021
BOKU - University of Natural Resources and Life Sciences 

SHOW PIPELINE HELP SECTION 

`
nextflow run main.nf --help 
`

BEFORE YOU BEGIN: 
- This pipeline was run on the Vienna Scientific Cluster 3 (VSC3) and 4 (VSC4)
- Make sure you have the dependencies installed (run ./check_dependencies.py contained in the git repo)
- Make sure you change the paths and executables inside "nextflow.config" according to your own environment 
- Average running time on a 1 GB genome size: ~15 hours with 16 parallel threads 

INFORMATION:
- The pipeline accepts two types of mapping files: BAM or PSL (format-wise)
- Besides that, it accepts also long read mapping records in BAM format (--bam_long option) 
- The pipeline works even if you pass only one type of information (say only PSL files)
- If you pass both BAM and PSL information, they will be combined prioritizing results in the more recent BAM format 
- If you pass BAM, PSL and BAM long-read information, the ranking will be long reads > short reads BAM > short reads PSL 
- If you pass repetitive elements annotated in GFF format, they will be used as repeat hints and prioritized over read mapping hints
- You can adjust priority settings from "nextflow.config" if you don't like what you just read above 
- The pipeline doesn't filter the PSL / BAM files, so you must provided files that have been filtered and sorted according to the official guidelines (would you like this to be implemented?) 

FOR MORE INFO, SEE OFFICIAL GUIDELINES:

https://bioinf.uni-greifswald.de/augustus/binaries/readme.rnaseq.html

https://vcru.wisc.edu/simonlab/bioinformatics/programs/augustus/docs/tutorial2015/prediction.html


Example usage: 

`
nextflow run main.nf \
--species_name <...> \
--genome_fasta <your_FASTA_file> \
--repeats_gff <repeat_GFF_annotation> \
--bam_long <dir_containing_long-read_BAM_files> \
--bam <short_read_BAM_files_dir> \
--psl <short_read_PSL_files_dir> \
`

For an overview on how I used it last time, see the `run_****.sh` script contained within repo


### Tutorial 


