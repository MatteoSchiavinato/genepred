# genepred

A Nextflow pipeline for gene prediction and functional annotation

Written and maintained by Matteo Schiavinato
Barcelona Supercomputing Center (BSC-CNS)
2020-2022

### Table of contents

- [genepred](#genepred)
    + [Table of contents](#table-of-contents)
    + [What is this for?](#what-is-this-for-)
    + [Quick setup and run](#quick-setup-and-run)
      - [Installation](#installation)
        * [Something didn't work?](#something-didn-t-work-)
      - [Data preparation](#data-preparation)
        * [official Augustus guidelines](#official-augustus-guidelines)
      - [Run the pipeline](#run-the-pipeline)
        * [species parameters](#species-parameters)
        * [extrinsic cfg file](#extrinsic-cfg-file)
      - [Resume / manage the pipeline](#resume---manage-the-pipeline)


### What is this for?

This pipeline takes the following input:
- Raw sequencing reads (pacbio, nanopore, illumina) in a table (`--input_data`)
- A genome sequence
- Augustus parameters
- Optional: modelled repetitive elements in GFF format

It automates the following steps:
- read preprocessing
- read mapping
- exon hints generation
- intron hints generation
- hints integration from multiple sources (PSL, BAM from short reads, BAM from long reads)
- hints prioritization
- parallelized gene prediction with Augustus
- gene set filtering based on evidence (0%, 1-99%, 100%)
- mRNA and protein sequence extraction
- removal of genes overlapping repetitive elements
- removal of ultrashort peptides
- production of a GBrowse-ready GFF file
- evaluation of run statistics and metrics in a file

It produces the following output:
- Raw and high quality gene sets predicted on the genome sequence
- Protein and transcript sequences of the predicted genes
- Run statistics and metrics

### Quick setup and run

#### Installation

First, clone the repository:

```
git clone https://github.com/MatteoSchiavinato/genepred.git
```

Then, check the missing dependencies, installing any missing dependency:

```
./check_dependencies.py
```

If missing, install nextflow from a terminal with access to an internet connection:

```
java -version     # must be 8 or higher
cd /your/bin/directory
curl -s https://get.nextflow.io | bash
./nextflow run hello
```

Show the pipeline's help section, to assess that everything is in place:

```
nextflow run main.nf --help
```

Edit the `nextflow.config` file, inserting your own executables and paths.

```
nano nextflow.config
```


##### Something didn't work?

Check out the suggestions on fixing issues in [this document](docs/ISSUES.md).


#### Data preparation

Train your augustus parameters if your species is not one of those contained inside `config/species` folder of the Augustus installation folder:

https://vcru.wisc.edu/simonlab/bioinformatics/programs/augustus/docs/tutorial2015/training.html

You can check the existence of the species you want (or a close relative) among your augustus species library just by running the following on your terminal:

```
augustus --species=help
```

Then, map your reads using any of the strategies described in the [tutorial document](docs/TUTORIAL.md).

##### official Augustus guidelines

Rnaseq:
https://bioinf.uni-greifswald.de/augustus/binaries/readme.rnaseq.html

Tutorial:
https://vcru.wisc.edu/simonlab/bioinformatics/programs/augustus/docs/tutorial2015/prediction.html

Manual:
https://math-inf.uni-greifswald.de/storages/uni-greifswald/fakultaet/mnf/mathinf/stanke/augustus_wrp.pdf

#### Run the pipeline

To run the pipeline, here is a basic command:

```
nextflow run main.nf \
--species_name <species_name> \
--genome_fasta <your_FASTA_file> \
--repeats_gff <repeat_GFF_annotation> \
--input_data <input_data_table>
```

An example run is found in a `run**.sh` script that I provided within this github repo. Have a look at that script carefully to see how I set things.


##### species parameters

When using the pipeline, there is an option called `--species_name` which you must set with the corresponding identifier of the species you want which you will see when running the command above. This is either a pre-existing species for which parameters are available, or one that you generated yourself when you trained your augustus parameters (see above).

##### extrinsic cfg file

Another file that Augustus wants is the `extrinsic.cfg` file. This file contains configuration parameters that Augustus uses. I provided one of these files inside this github repository. If you open it, you'll see that there are different values associated with different features.

Read the Augustus official manual to understand what they mean:

https://math-inf.uni-greifswald.de/storages/uni-greifswald/fakultaet/mnf/mathinf/stanke/augustus_wrp.pdf


#### Resume / manage the pipeline

You can find all the information on this in the [dedicated document](docs/NEXTFLOW.md).
