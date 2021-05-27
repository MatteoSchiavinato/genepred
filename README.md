# genepred

A Nextflow pipeline for gene prediction and functional annotation

Compiled by Matteo Schiavinato

2020-2021

University of Natural Resources and Life Sciences (BOKU), Vienna, Austria


### show pipeline help section

```
nextflow run main.nf --help
```

### example usage

```
nextflow run main.nf \
--species_name <...> \
--genome_fasta <your_FASTA_file> \
--repeats_gff <repeat_GFF_annotation> \
--bam_long <dir_containing_long-read_BAM_files> \
--bam <short_read_BAM_files_dir> \
--psl <short_read_PSL_files_dir> \
```

For an overview on how I used it last time, see the `run_****.sh` script contained within repo


# Tutorial

### Before you begin

To use this pipeline, you will need:
- **BAM** or **PSL** files obtained by mapping RNASeq reads onto the genome sequence where the prediction has to be done
- A **FASTA** genome sequence
- **Augustus parameters** trained for the species that you're doing a prediction for. Parameters are in a file

Consult the official Augustus guidelines to understand what each step is doing:

### official Augustus guidelines

https://bioinf.uni-greifswald.de/augustus/binaries/readme.rnaseq.html

https://vcru.wisc.edu/simonlab/bioinformatics/programs/augustus/docs/tutorial2015/prediction.html



### RNASeq reads mapping

If you don't have BAM or PSL mapping files from your RNASeq reads, then you have to generate them yourself. To do so, you can use any program you like, for example hisat2 (for BAM files) or BLAT (for PSL files). The original Augustus guidelines use BLAT, which is however an old tool (2002). Hence, they progressively moved towards BAM files over time. In my pipeline, you can use both formats at the same time, so there's no need to choose if you want to use both.

### Long-read mapping

Long-read RNA reads are becoming established in the world of research, so I decided to include them in the pipeline since we wanted to use them for a project. After consulting with Katharina Hoff, one of the Augustus authors, I received instructions on how they would handle long reads in the Augustus pipeline (which is how I implemented it in the code).
