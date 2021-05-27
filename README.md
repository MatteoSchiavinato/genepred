# genepred
A Nextflow pipeline for gene prediction and functional annotation
Matteo Schiavinato
2020-2021
BOKU - University of Natural Resources and Life Sciences 

BEFORE YOU BEGIN: 
- This pipeline was run on the Vienna Scientific Cluster 3 (VSC3)
- Make sure you change the executor (currently "slurm") to "local" inside each process of "main.nf" if you're running this pipeline without slurm 
- If you run it in local mode, then change "slurm" to "local" inside the "nextflow.config" file too 
- Make sure you change the --account, --partition and --qos settings in each process (in "main.nf") and inside "nextflow.config" according to your user 
- Make sure you have the dependencies installed (see "nextflow.config", section env{} )
- Make sure you change the paths and executables inside "nextflow.config" according to your own environment 
- The "blast" section has not been tested thoroughly. Use --skip_blast to avoid it if you're in a hurry!
- Average running time on a 1 GB genome size: ~25-30 hours with 16 parallel threads 

INFORMATION:
- The pipeline works even if you pass only one type of information (say only PSL files)
- If you pass both BAM and PSL information, they will be combined prioritizing results in the more recent BAM format 
- If you pass repetitive elements annotated in GFF format, they will be used as repeat hints and prioritized over BAM and PSL hints
- You can adjust priority settings from "nextflow.config" 
- The pipeline doesn't filter the PSL / BAM files, so you must provided files that have been filtered and sorted according to the official guidelines

FOR MORE INFO, SEE OFFICIAL GUIDELINES:

https://bioinf.uni-greifswald.de/augustus/binaries/readme.rnaseq.html

https://vcru.wisc.edu/simonlab/bioinformatics/programs/augustus/docs/tutorial2015/prediction.html


Basic usage (without blast): 

`nextflow run main.nf [OPTIONS] --skip_blast`

For an overview on how I used it last time, see `run.sh` (contained within repo) 
