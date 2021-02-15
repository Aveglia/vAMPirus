![vAMPirus logo](https://sync.palmuc.org/index.php/s/gscg4m4PtdN5ZwY/preview)

                            An automated virus amplicon sequencing analysis pipeline

# Table of contents

* [Quick intro](#Quick-intro)
  * [Contact/support](#Contact/support)  
* [Getting started](#Getting-started)
  * [Order of operations](##Order-of-operations)
    * [Dependencies](###Dependencies-(See-How-to-cite))
* [Installing vAMPirus](#Installing-vAMPirus)
  * [Setting up vAMPirus dependencies](##Setting-up-vAMPirus-dependencies)
* [Testing vAMPirus install](#Testing-vAMPirus-installation)
* [Running vAMPirus](#Running-vAMPirus)
* [Who to cite](#Who-to-cite)


# Quick intro

Here we present vAMPirus (https://github.com/Aveglia/vAMPirus.git), an automated virus amplicon sequencing analysis program.

The vAMPirus program contains two different pipelines:

1. DataCheck pipeline: provides the user an interactive html report file containing information regarding sequencing success per sample as well as a preliminary look into the clustering behavior of the data which can be leveraged by the user to inform future analyses

2. Analyze pipeline: a comprehensive analysis of the provided data producing a wide range of results and outputs which includes an interactive report with figures and statistics.


NOTE => This is a more brief overview of how to install and set up vAMPirus, for more detail see the wiki.


## Contact/support:

If you have a feature request or any feedback/questions, feel free to email vAMPirusHelp@gmail.com or you can open an Issue on GitHub.


# Getting started

## Quick order of operations

1. Clone vAMPirus from github
2. Execute the vampirus_startup.sh script to install dependencies and any databases specified
3. Test installation with supplied test dataset
5. Launch the DataCheck pipeline with your dataset and adjust parameters if necessary
6. Launch the Analyze pipeline with your dataset

### Dependencies (see Who to cite section)    

1. Python versions 3.6/2.7
2. Diamond version 0.9.30
3. FastQC version 0.11.9
4. fastp version 0.20.1
5. Clustal Omega version 1.2.4
6. IQ-TREE version 2.0.3
7. ModelTest-NG version 0.1.6
8. MAFFT version 7.464
9. vsearch version 2.14.2
10. BBMap version 38.79
11. trimAl version 1.4.1
12. CD-HIT version 4.8.1
13. EMBOSS version 6.5.7.0
14. seqtk version 1.3

If you plan on using Conda to run vAMPirus, all dependencies will be installed as a Conda environment automatically with the vampirus_startup.sh script.

If you plan to use a container engine like Singularity, the dependencies of vAMPirus have been built as a Docker container (aveglia/vAMPirus) that's stored in the vAMPirus directory and will be set up by Nextflow upon the initial launch of a vAMPirus pipeline.


# Installing vAMPirus

Clone the most recent version of vAMPirus from GitHub using:

  git clone https://github.com/Aveglia/vAMPirus.git


## Setting up vAMPirus dependencies

To launch vAMPirus, you will need to have Nextflow (https://www.nextflow.io/) installed on your system. If you do not already, the vampirus_startup.sh script will install it for you.

You will also need to decide whether you would like to set up and run vAMPirus with Anaconda/Miniconda OR with a container engine like Singularity.

If you plan on using Anaconda, the startup script will install conda (if necessary) and build a vAMPirus environment with all necessary dependencies.

If you plan to use Singularity, it will need to be installed separately prior to launching the vAMPirus.nf script.

To set up and install vAMPirus dependencies, simply move to the vAMPirus directory and run the vampirus_startup.sh script.

    cd ./vAMPirus; ./vampirus_startup.sh -h

The start up script will check your system for Nextflow and Anaconda/Miniconda (can be skipped) and if they are not present, the script will ask if you would like to install these programs. If you answer with 'y', the script will install the missing programs and will build the vAMPirus conda environment and the installation is complete.

You can also use the startup script to install different databases to use for vAMPirus analyses, these include:

1. NCBIs Viral protein RefSeq database
2. The proteic version of the Reference Viral DataBase (RVDB) (See https://f1000research.com/articles/8-530)
3. The complete NCBI NR protein database

To use the vampirus_startup.sh script to download any or all of these databases listed above you just need to use the "-d" option.
If we look at the script usage:

    General execution:

    vampirus_startup.sh -h [-d 1|2|3|4] [-s]

        Command line options:

            [ -h ]                       	Print help information

            [ -d 1|2|3|4 ]                Set this option to create a database directiory within the current working directory and download the following databases for taxonomy assignment:

                                                        1 - Download the proteic version of the Reference Viral DataBase (See the paper for more information on this database: https://f1000research.com/articles/8-530)
                                                        2 - Download only NCBIs Viral protein RefSeq database
                                                        3 - Download only the complete NCBI NR protein database
                                                        4 - Download all three databases

            [ -s ]                       Set this option to skip conda installation and environment set up (you can use if you plan to run with Singularity and the vAMPirus Docker container)

For example, if you would like to install Nextflow, download NCBIs Viral protein RefSeq database, run:

    `./vampirus_startup.sh -d 1`

and if we wanted to do the same thing as above but skip the Conda check/installation, run:

    `./vampirus_startup.sh -d 1 -s`

### For MacOS users

To be safe, it is recommended that MacOS users set up and install Singularity to use the vAMPirus Docker image. We have noticed errors caused by m

# Testing vAMPirus installation

The startup script will generate a text file (STARTUP_HELP.txt) that has instructions and example commands to test the installation.

NOTE => If using Singularity, when you run the test command calling for singularity (-profile singularity) Nextflow will set up the dependencies.

Launch commands for testing (you do not need to edit anything in the config files for test commands):

### DataCheck test =>

      `/path/to/nextflow run /path/to/vAMPirus.nf -c /path/to/vampirus.config -profile conda,test --DataCheck`

OR

      `nextflow run vAMPirus.nf -c vampirus.config -profile singularity,test --DataCheck`

### Analyze test =>

      `/path/to/nextflow run /path/to/vAMPirus.nf -c /path/to/vampirus.config -profile conda,test --Analyze --ncASV --pcASV --stats run`

OR

      `nextflow run vAMPirus.nf -c vampirus.config -profile singularity,test --Analyze --ncASV --pcASV --stats run`


# Running vAMPirus

If you done the setup and confirmed installation success with the test commands, you are good to get going with your own data. Before getting started edit the configuration file with the parameters and other options you plan to use.

Here are some example vAMPirus launch commands:
### DataCheck pipeline =>

Example 1. Launching the vAMPirus DataCheck pipeline using conda

      `nextflow run vAMPirus.nf -c vampirus.config -profile conda --DataCheck`

Example 2. Launching the vAMPirus DataCheck pipeline using Singularity and multiple primer removal with the path to the fasta file with the primer sequences set in the launch command

      `nextflow run vAMPirus.nf -c vampirus.config -profile singularity --DataCheck --multi --primers /PATH/TO/PRIMERs.fa`

Example 3. Launching the vAMPirus DataCheck pipeline with primer removal by global trimming of 20 bp from forward reads and 26 bp from reverse reads

      `nextflow run vAMPirus.nf -c vampirus.config -profile conda --DataCheck --GlobTrim 20,26`


### Analyze pipeline =>

Example 4. Launching the vAMPirus Analyze pipeline with singularity with ASV and AminoType generation with all accesory analyses (taxonomy assignment, EMBOSS, IQTREE, statistics)

      `nextflow run vAMPirus.nf -c vampirus.config -profile singularity --Analyze --stats run`

Example 5. Launching the vAMPirus Analyze pipeline with conda to perform multiple primer removal and protein-based clustering of ASVs, but skip most of the extra analyses

      `nextflow run vAMPirus.nf -c vampirus.config -profile conda --Analyze --pcASV --skipPhylogeny --skipEMBOSS --skipTaxonomy --skipReport`

Example 6. Launching vAMPirus Analyze pipeline with conda to produce only ASV-related results

      `nextflow run vAMPirus.nf -c vampirus.config -profile conda --Analyze --skipAminoTyping --stats run`


## Resuming analyses =>

If an analysis is interrupted, you can use Nextflows "-resume" option that will start from the last cached "check point".

For example if the analysis launched with the command from Example 6 above was interrupted, all you would need to do is add the "-resume" to the end of the command like so:

      `nextflow run vAMPirus.nf -c vampirus.config -profile conda --Analyze --skipAminoTyping --stats run -resume`


# Who to cite:

If you do use vAMPirus for your analyses, please cite the following ->

1. vAMPirus - Veglia, A.J., Rivera Vicens, R., Grupstra, C., Howe-Kerr, L., and Correa A.M.S. (2020) vAMPirus: An automated virus amplicon sequencing analysis pipeline. Zenodo. DOI:

2. Diamond - Buchfink B, Xie C, Huson DH. (2015) Fast and sensitive protein alignment using DIAMOND. Nat Methods. 12(1):59-60. doi:10.1038/nmeth.3176

3. FastQC - Andrews, S. (2010). FastQC:  A Quality Control Tool for High Throughput Sequence Data [Online]. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

4. fastp - Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics, 34(17), i884-i890.

5. Clustal Omega - Sievers, F., Wilm, A., Dineen, D., Gibson, T.J., Karplus, K., Li, W., Lopez, R., McWilliam, H., Remmert, M., Söding, J. and Thompson, J.D., 2011. Fast, scalable generation of high‐quality protein multiple  sequence alignments using Clustal Omega. Molecular systems biology, 7(1), p.539.

6. IQ-TREE - Minh, B. Q., Schmidt, H. A., Chernomor, O., Schrempf, D., Woodhams, M. D., Von Haeseler, A., & Lanfear, R. (2020). IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era. Molecular Biology and Evolution, 37(5), 1530-1534.

7. ModelTest-NG - Darriba, D., Posada, D., Kozlov, A. M., Stamatakis, A., Morel, B., & Flouri, T. (2020). ModelTest-NG: a new and scalable tool for the selection of DNA and protein evolutionary models. Molecular biology and evolution, 37(1), 291-294.

8. MAFFT - Katoh, K., & Standley, D. M. (2013). MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Molecular biology and evolution, 30(4), 772-780.

9. vsearch - Rognes, T., Flouri, T., Nichols, B., Quince, C., & Mahé, F. (2016). VSEARCH: a versatile open source tool for metagenomics. PeerJ, 4, e2584.

10. BBMap - Bushnell, B. (2014). BBTools software package. URL http://sourceforge. net/projects/bbmap.

11. trimAl - Capella-Gutiérrez, S., Silla-Martínez, J. M., & Gabaldón, T. (2009). trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. Bioinformatics, 25(15), 1972-1973.

12. CD-HIT - Fu, L., Niu, B., Zhu, Z., Wu, S., & Li, W. (2012). CD-HIT: accelerated for clustering the next-generation sequencing data. Bioinformatics, 28(23), 3150-3152.

13. EMBOSS - Rice, P., Longden, I., & Bleasby, A. (2000). EMBOSS: the European molecular biology open software suite.

14. seqtk - Li, H. (2012). seqtk Toolkit for processing sequences in FASTA/Q formats. GitHub, 767, 69.

15. UNOISE algorithm - R.C. Edgar (2016). UNOISE2: improved error-correction for Illumina 16S and ITS amplicon sequencing, https://doi.org/10.1101/081257
