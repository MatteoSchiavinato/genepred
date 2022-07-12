# hands-on tutorial

This tutorial brings you from the beginning to the end of a single gene prediction run. The first part of it shows how to generate mapping files the right way, but if you already have them, you can skip it.


### Read mapping

To generate mapping files from your short RNASeq reads you can use any program you like, for example [hisat2](http://daehwankimlab.github.io/hisat2/manual "HISAT2 manual") (for BAM files) or [BLAT](https://genome.ucsc.edu/goldenpath/help/blatSpec.html "BLAT manual") (for PSL files). The original Augustus guidelines use [BLAT](https://genome.ucsc.edu/goldenpath/help/blatSpec.html "BLAT manual"), which is however an old tool (2002). Hence, they progressively moved towards BAM files over time. In my pipeline, you can use both formats at the same time, so there's no need to choose if you want to use both.

Here you can find basic instructions to map your reads with these two programs.

##### HISAT2

[hisat2](http://daehwankimlab.github.io/hisat2/manual "HISAT2 manual") requires that you first build an index:

```
hisat2-build \
-p <threads> \
<genome_FASTA> \
<output_prefix>
```

The genome FASTA must be the one you want to predict your genes on. Usually, I use the same name as the genome FASTA also for the output prefix, this way the program will generate index files ending in `*.ht2` right next to the genome sequence. However, this is your choice. All [hisat2](http://daehwankimlab.github.io/hisat2/manual "HISAT2 manual") needs to map reads is the index, once it's generated.

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

[BLAT](https://genome.ucsc.edu/goldenpath/help/blatSpec.html "BLAT manual") is a Blast-Like-Alignment-Tool which was released the first time in 2002. It can be painfully slow (i.e. weeks) if you don't set parameters with a grain of salt. After half a PhD spent on gene prediction, I figured these are the settings that work for me. They may not work for you, so understand them well (through their manual) and set them the way you want. A few years later they released **pblat**, which is [BLAT](https://genome.ucsc.edu/goldenpath/help/blatSpec.html "BLAT manual") with parallelization.

A characteristic of [BLAT](https://genome.ucsc.edu/goldenpath/help/blatSpec.html "BLAT manual") is that it cannot map paired-end reads. However, there is a workaround which is to map R1 and R2 separately and then join them afterwards. Hence, for [BLAT](https://genome.ucsc.edu/goldenpath/help/blatSpec.html "BLAT manual") mapping treat them as two independent single end files. Here, the example uses pblat, which works exactly like [BLAT](https://genome.ucsc.edu/goldenpath/help/blatSpec.html "BLAT manual") but has a `-threads` option.

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

The `-noHead` option prevents the addition of a header. This may not help the human reader but helps the filtering afterwards so use it. Another option that requires discussion is `-minIdentity`. This option sets the minimum sequence identity between a read and the target genome sequence. This should be adjusted according to your estimation of divergence between reads and reference. If you set it close to 99, the program will be faster. However, to leverage the real power of [BLAT](https://genome.ucsc.edu/goldenpath/help/blatSpec.html "BLAT manual") (the heuristic approach) you should not keep it too tight. In fact, [BLAT](https://genome.ucsc.edu/goldenpath/help/blatSpec.html "BLAT manual") can map reads quite well at low sequence identity, due to the fact that it doesn't have an index of k-mers to read from. Do some testing! Finally, the `-maxIntron` option sets how long an intron can be (i.e. a split in a read). Increasing this value makes the program become horribly slow, so handle with care.

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


##### Long-read mapping

Long-read RNA reads are becoming established in the world of research, so I decided to include them in the pipeline since we wanted to use them for a project. After consulting with Katharina Hoff, one of the Augustus authors, I received instructions on how they would handle long reads in the Augustus pipeline (which is how I implemented it in the code). Here, you can see how I mapped the reads. For the mapping I use [minimap2](https://lh3.github.io/minimap2/minimap2.html "minimap2 manual"), a very common tool for long read mapping. It's very easy and fast so I would strongly recommend using this. If not, any tool producing a BAM file will do.

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

The `-k` option sets the kmer size used for the mapping. See the minimap manual for understanding how this parameter works. The file specified with `-d` is the index created by [minimap2](https://lh3.github.io/minimap2/minimap2.html "minimap2 manual"), so you can use any name. You can read the description of the other options on their manual. What is important is that you specify an output BAM file with `-o`, providing a genome sequence and your read file in FASTA format.
