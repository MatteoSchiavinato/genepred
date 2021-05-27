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


### Before you begin

To use this pipeline, you will need:
- **BAM** or **PSL** files obtained by mapping RNASeq reads onto the genome sequence where the prediction has to be done. The pipeline can handle both at the same time, so you don't have to choose.
- A **FASTA** genome sequence which is the same where you mapped the reads.
- **Augustus parameters** trained for the species that you're doing a prediction for. Parameters are multiple files contained in a directory called *exactly* like you named your species with `--species_name`. This directory is contained inside `config/species`, which in turn is inside the main directory of your Augustus executable.
- Optional, but recommended: a **GFF** repeat annotation performed with (for example) the RepeatMasker pipeline

Consult the official Augustus guidelines to understand what each step is doing:

### official Augustus guidelines

Rnaseq:
https://bioinf.uni-greifswald.de/augustus/binaries/readme.rnaseq.html

Tutorial:
https://vcru.wisc.edu/simonlab/bioinformatics/programs/augustus/docs/tutorial2015/prediction.html

Manual:
https://math-inf.uni-greifswald.de/storages/uni-greifswald/fakultaet/mnf/mathinf/stanke/augustus_wrp.pdf


# Tutorial

This tutorial brings you from the beginning to the end of a single gene prediction run. The first part of it shows how to generate mapping files the right way, but if you already have them, you can skip it.


### Read mapping

To generate mapping files from your short RNASeq reads you can use any program you like, for example hisat2 (for BAM files) or BLAT (for PSL files). The original Augustus guidelines use BLAT, which is however an old tool (2002). Hence, they progressively moved towards BAM files over time. In my pipeline, you can use both formats at the same time, so there's no need to choose if you want to use both.

Here you can find basic instructions to map your reads with these two programs.

##### HISAT2

Hisat2 requires that you first build an index:

```
hisat2-build \
-p <threads> \
<genome_FASTA> \
<output_prefix>
```

The genome FASTA must be the one you want to predict your genes on. Usually, I use the same name as the genome FASTA also for the output prefix, this way the program will generate index files ending in `*.ht2` right next to the genome sequence. However, this is your choice. All hisat2 needs to map reads is the index, once it's generated.

Once the index is generated, you can then map your reads. Here you can see how I ran it, with the parameters I set. Check parameters meaning and usage on their main manual page.

Here's how to map paired-end reads:

```
hisat2 \
-p <threads> \
-f \
--max-intronlen 500000 \
--mp 6,2 \
--rdg 5,3 \
--rfg 5,3 \
--score-min L,0.0,-0.5 \
-k 5 \
-1 <read_file_1> \
-2 <read_file_2> \
-x <output_prefix> \
-S <output_sam_file>
```

In case of single-end reads, you just have to switch the `-1 / -2` parameters to a single `-U <read_file>`.

Once reads are mapped, the mapping records must be filtered and sorted prior to being used in the pipeline. I do that with samtools using the command below. The `LC_ALL=C` command is important because you probably won't be able to sort the files properly if the normal `locale` is set. To understand what this means, type `locale` in your terminal and enter. Likely, `LC_ALL=` won't be set. This procedure is the same for both single-end and paired-end reads.

```
LC_ALL=C &&
samtools view -h -F 0x4 -F 0x0100 <raw_sam_file> | \
samtools sort -O bam -@ 16 \
> <filt_and_sorted_bam_file>
LC_ALL=''
```

The output files of this command are ready to be used in my pipeline. Try to name them with filenames that represent their sample ID, so that the pipeline can pick up their name properly. Suggestion: `sample_id.pe.fs.bam` (paired-end) or `sample_id.se.fs.bam` (for single end). Should work irregardless of this naming system, but 1) tidy naming is proper naming, and 2) you never know.


##### BLAT

BLAT is a Blast-Like-Alignment-Tool which was released the first time in 2002. It can be painfully slow (i.e. weeks) if you don't set parameters with a grain of salt. After half a PhD spent on gene prediction, I figured these are the settings that work for me. They may not work for you, so understand them well (through their manual) and set them the way you want. A few years later they released **pblat**, which is blat with parallelization.

A characteristic of BLAT is that it cannot map paired-end reads. However, there is a workaround which is to map R1 and R2 separately and then join them afterwards. Hence, for BLAT mapping treat them as two independent single end files. Here, the example uses pblat, which works exactly like blat but has a `-threads` option.

```
pblat \
-noHead \
-t=dna \
-q=rna \
-threads=16 \
-tileSize=10 \
-stepSize=10 \
-minMatch=2 \
-minIdentity=93 \
-maxGap=2 \
-out=psl \
-maxIntron=100000 \
<genome_fasta> \
<read_file> \
<output_psl_file>
```

The `-noHead` option prevents the addition of a header. This may not help the human reader but helps the filtering afterwards so use it. Another option that requires discussion is `-minIdentity`. This option sets the minimum sequence identity between a read and the target genome sequence. This should be adjusted according to your estimation of divergence between reads and reference. If you set it close to 99, the program will be faster. However, to leverage the real power of BLAT (the heuristic approach) you should not keep it too tight. In fact, BLAT can map reads quite well at low sequence identity, due to the fact that it doesn't have an index of k-mers to read from. Do some testing! Finally, the `-maxIntron` option sets how long an intron can be (i.e. a split in a read). Increasing this value makes the program become horribly slow, so handle with care.

Once the mapping is done, you have to filter the files. The filtering command below is for paired-end reads, but for single-end reads it works the same, just substitute `cat <psl_file_R1> <psl_file_R2>` with `cat <psl_file>` because there is only a single file, and remove the `--paired` option from the `filterPSL.pl` step.

```
LC_ALL=C &&
cat <psl_file_R1> <psl_file_R2> | \
sort -k 10,10 -T tmp_sort | \
filterPSL.pl --uniq --paired | \
sort -n -k 16,16 -T tmp_sort | \
sort -s -k 14,14 -T tmp_sort \
> <filt_and_sorted_psl_file>
LC_ALL=''
```

In this command, the `LC_ALL` setting is the same as specified above so I won't explain it again. The `sort -k 10,10` sorts by read name, because the mapping results are usually by genome position. The `filterPSL.pl` script is contained inside the augustus scripts, which you can find in the Augustus main directory. The options `--uniq` tells the script to only keep one record per read, which is of paramount importance. The `--paired` option must be used with paired-end reads only and tells the script to consider two lines at the same time (R1 and R2) when doing the filtering. The following two sort commands re-sort the file by genome coordinate, which makes it ready for gene prediction.

just navigate
##### Long-read mapping

Long-read RNA reads are becoming established in the world of research, so I decided to include them in the pipeline since we wanted to use them for a project. After consulting with Katharina Hoff, one of the Augustus authors, I received instructions on how they would handle long reads in the Augustus pipeline (which is how I implemented it in the code). Here, you can see how I mapped the reads. For the mapping I use minimap2, a very common tool for long read mapping. It's very easy and fast so I would strongly recommend using this. If not, any tool producing a BAM file will do.

```
minimap2 \
-H \
-k 15 \
-d <sample_id>.mm2_index \
-G 500k \
-F 1500 \
-A 2 -B 4 -O 4,24 -E 2,1 \
-u f \
-a \
-o <output_sam_file> \
-t 16 \
-x splice:hq \
<genome_fasta> \
<read_file>
```

The `-k` option sets the kmer size used for the mapping. See the minimap manual for understanding how this parameter works. The file specified with `-d` is the index created by minimap2, so you can use any name. You can read the description of the other options on their manual. What is important is that you specify an output BAM file with `-o`, providing a genome sequence and your read file in FASTA format.


### Running the pipeline

Once you have filtered and sorted mapping records from at least one source (be it hisat2, blat, or minimap2), you can use the pipeline.

##### species parameters

To use it, you will need to have trained Augustus parameters for the species you want to use. You can check the existence of the species you want (or a close relative) among your augustus species library just by running the following on your terminal:

```
augustus --species=help
```

When using the pipeline, there is an option called `--species_name` which you must set with the corresponding identifier of the species you want which you will see when running the command above. If your species is not there, you can generate your parameters by using the augustus guidelines to train parameters. In case you work with beta vulgaris, you can find the parameters in the github repository. In that case, just do the following:

- navigate to the Augustus main directory
- navigate to the `config/species` subdirecory
- create a directory called `beta_vulgaris`
- copy all the parameter files contained in this github repo inside that folder
- run `augustus --species=help` again and see if `beta_vulgaris` is there

For all other species that will require your parameter training, you can use the same strategy but with a different name (once you get the parameters).

##### extrinsic cfg file

Another file that Augustus wants is the `extrinsic.cfg` file. This file contains configuration parameters that Augustus uses. I provided one of these files inside this github repository which was optimized with Beta vulgaris reads. If you open it, you'll see that there are different values associated with different features.

Read the Augustus official manual to understand what they mean (too long to put here).


##### The nextflow.config file

Nextflow needs a config file with associated parameters (which are the same you can pass via command line using the `--`). This file is contained in the same directory where `main.nf` is. You can modify the way you want, but before you run the pipeline, open it in a text editor and check that all the paths that are contained there exist in your system. If not, fix them.


##### run the pipeline

To run the pipeline, here is a basic command:

```
nextflow run main.nf \
--species_name <species_name> \
--genome_fasta <your_FASTA_file> \
--repeats_gff <repeat_GFF_annotation> \
--bam_long <dir_containing_long-read_BAM_files> \
--bam <short_read_BAM_files_dir> \
--psl <short_read_PSL_files_dir> \
```

We discussed how to obtain the species name and the associated parameters, and how to generate the mapping records. Everything else you should have already. To see further command line options, just run `nextflow run main.nf --help`.

An example run is found in a `run**.sh` script that I provided within this github repo. Have a look at that script carefully to see how I set things. Pay attention to the work directory that I set inside there, you may want to read about that directory on the nextflow guidelines. It can become really big if you run multiple runs, so remove it once the run is complete (it's mostly temp files). If the pipeline crashes, you can restart it with `nextflow run main.nf -resume`. However, it will only resume if the `work` directory is there. If you removed it, it will start again from the beginning.
