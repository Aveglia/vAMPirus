![vAMPirus logo](https://sync.palmuc.org/index.php/s/gscg4m4PtdN5ZwY/preview)

                            An automated virus amplicon sequencing analysis pipeline

# Table of contents

* [Quick intro](#Quick-intro)
  * [Contact/support](#Contact/support)  
                                
# Quick intro

Here we present vAMPirus (https://github.com/Aveglia/vAMPirus.git), an automated and easy-to-use virus amplicon sequencing analysis program.

The vAMPirus program contains two different pipelines:

1. DataCheck pipeline: provides the user an interactive html report file containing information regarding sequencing success per sample as well as a preliminary look into the clustering behavior of the data which can be leveraged by the user to inform future analyses

2. Analyze pipeline: a comprehensive analysis of the provided data producing a wide range of results and outputs which includes an interactive report with figures and statistics.


NOTE => This is a more brief overview of how to install and set up vAMPirus, for more detail see the wiki.

## Contact/support:

If you have a feature request or any feedback/questions, feel free to email vAMPirusHelp@gmail.com or you can open an Issue on GitHub.

# Getting started

## Order of operations

    1. Clone vAMPirus from github  -  `git clone https://github.com/Aveglia/vAMPirus.git`

    2. Run the vAMPirus start up script to install dependencies and any databases specified

    3. Test vAMPirus pipelines with the supplied test dataset

    4. Edit parameters in the vAMPirus configuration file

    5. Run DataCheck mode with dataset

    6. Run Analyze mode with desired clustering technique and %ID

**Dependencies (will be built as a conda environment) ->**    

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

## Contact/support:

Please contact Alex Veglia at ajv5@rice.edu with any feedback or questions. Any kind of input from the community is welcomed and encouraged!

# Getting started

## Installing vAMPirus

Clone the most recent version of vAMPirus from github using:

    `git clone https://github.com/Aveglia/vAMPirus.git`

## Setting up vAMPirus dependencies and checking installation

To deploy vAMPirus, you need to have Nextflow and Anaconda/Miniconda installed on your system. To set up and install vAMPirus
dependencies, simply move to the vAMPirus directory and run the vampirus_startup.sh script.

    `cd ./vAMPirus`; `./vampirus_startup.sh`

The start up script will check your system for Nextflow and Anaconda/Miniconda and if they are not present, the script will ask
if you would like to install these programs. If you answer with 'y', the script will install the missing programs and will build
the vAMPirus conda environment and the installation is complete.

You can also use the startup script to install different databases to use for vAMPirus analyses, these include:

    1. NCBIs Viral protein RefSeq database
    2. The proteic version of the Reference Viral DataBase (RVDB) (See https://f1000research.com/articles/8-530)
    3. The complete NCBI NR protein database

To use the vampirus_startup.sh script to download any or all of these databases listed above you just need to use the "-d" option.
If we look at the script usage:


    General execution:

    vampirus_startup.sh -h -d [1|2|3|4]


    Command line options:

            [ -h ]                       	Print help information

            [ -d 1|2|3|4 ]                  Set this option to create a database directiory within the current working directory and download the following databases for taxonomy assignment:

                                                        1 - Download only NCBIs Viral protein RefSeq database
                                                        2 - Download the proteic version of the Reference Viral DataBase (See the paper for more information on this database: https://f1000research.com/articles/8-530)
                                                        3 - Download only the complete NCBI NR protein database
                                                        4 - Download all three databases    

So, if we wanted to download NCBIs Viral protein RefSeq database, we would just need to run:

    `./vampirus_startup.sh -d 1`

[This will check for Nextflow and Anaconda/Miniconda installations, then download specified database(s), which in this case is NCBIs RefSeq database.]  

## Testing vAMPirus installation

    A test dataset is provided in the vAMPirus/example_data. To ensure that vAMPirus is set up properly before running with your own data, you can run:

    `nextflow run vAMPirusv0.1.0.nf -c ./example_data/vampirus_test.config -with-conda /PATH/TO/miniconda3/env/vAMPirus --Analyze --testing`

## Example vAMPirus launch commands

DataCheck pipeline =>

Example 1. Launching the vAMPirus DataCheck pipeline using conda

    `nextflow run vAMPirus.nf -c vampirus.config -profile conda --DataCheck`

Example 2. Launching the vAMPirus DataCheck pipeline using Singularity and multiple primer removal with the path to the fasta file with the primer sequences set in the launch command

    `nextflow run vAMPirus.nf -c vampirus.config -profile singularity --DataCheck --multi --primers /PATH/TO/PRIMERs.fa`

Example 3. Launching the vAMPirus DataCheck pipeline with primer removal by global trimming of 20 bp from forward reads and 26 bp from reverse reads

    `nextflow run vAMPirus.nf -c vampirus.config -profile conda --DataCheck --GlobTrim 20,26`


Analyze pipeline =>

Example 4. Launching the vAMPirus Analyze pipeline with singularity with ASV and AminoType generation with all accesory analyses (taxonomy assignment, EMBOSS, IQTREE, statistics)

    `nextflow run vAMPirus.nf -c vampirus.config -profile singularity --Analyze --stats run`

Example 5. Launching the vAMPirus Analyze pipeline with conda to perform multiple primer removal and protein-based clustering of ASVs, but skip most of the extra analyses

    `nextflow run vAMPirus.nf -c vampirus.config -profile conda --Analyze --pcASV --skipPhylogeny --skipEMBOSS --skipTaxonomy --skipReport`

Example 6. Launching vAMPirus Analyze pipeline with conda to produce only ASV-related results

    `nextflow run vAMPirus.nf -c vampirus.config -profile conda --Analyze --skipAminoTyping --stats run`


Resuming analyses =>

If an analysis is interupted, you can use Nextflows "-resume" option that will start from the last cached "check point".

For example if the analysis launched with the command from Example 6 above was interupted, all you would need to do is add the "-resume" to the end of the command like so:

    `nextflow run vAMPirus.nf -c vampirus.config -profile conda --Analyze --skipAminoTyping --stats run -resume`


## How to cite:

    If you do use vAMPirus for your analyses, please cite the following ->

        1. vAMPirus - Veglia, A.J., Rivera Vicens, R., Grupstra, C., Howe-Kerr, L., and Correa A.M.S. (2020) vAMPirus: An automated virus amplicon sequencing analysis pipeline. Zenodo. DOI:

        2. Diamond - Buchfink B, Xie C, Huson DH. (2015) Fast and sensitive protein alignment using DIAMOND. Nat Methods. 12(1):59-60. doi:10.1038/nmeth.3176

        3. FastQC - Andrews, S. (2010). FastQC:  A Quality Control Tool for High Throughput Sequence Data [Online]. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

        4. fastp - Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics, 34(17), i884-i890.

        5. Clustal Omega - Sievers, F., Wilm, A., Dineen, D., Gibson, T.J., Karplus, K., Li, W., Lopez, R., McWilliam, H., Remmert, M., Söding, J. and Thompson, J.D., 2011. Fast, scalable generation of high‐quality protein multiple
                           sequence alignments using Clustal Omega. Molecular systems biology, 7(1), p.539.

        6. IQ-TREE - Minh, B. Q., Schmidt, H. A., Chernomor, O., Schrempf, D., Woodhams, M. D., Von Haeseler, A., & Lanfear, R. (2020). IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era.
                     Molecular Biology and Evolution, 37(5), 1530-1534.

        7. ModelTest-NG - Darriba, D., Posada, D., Kozlov, A. M., Stamatakis, A., Morel, B., & Flouri, T. (2020). ModelTest-NG: a new and scalable tool for the selection of DNA and protein evolutionary models. Molecular biology and
                          evolution, 37(1), 291-294.

        8. MAFFT - Katoh, K., & Standley, D. M. (2013). MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Molecular biology and evolution, 30(4), 772-780.

        9. vsearch - Rognes, T., Flouri, T., Nichols, B., Quince, C., & Mahé, F. (2016). VSEARCH: a versatile open source tool for metagenomics. PeerJ, 4, e2584.

        10. BBMap - Bushnell, B. (2014). BBTools software package. URL http://sourceforge. net/projects/bbmap.

        11. trimAl - Capella-Gutiérrez, S., Silla-Martínez, J. M., & Gabaldón, T. (2009). trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. Bioinformatics, 25(15), 1972-1973.

        12. CD-HIT - Fu, L., Niu, B., Zhu, Z., Wu, S., & Li, W. (2012). CD-HIT: accelerated for clustering the next-generation sequencing data. Bioinformatics, 28(23), 3150-3152.

        13. EMBOSS - Rice, P., Longden, I., & Bleasby, A. (2000). EMBOSS: the European molecular biology open software suite.

        14. seqtk - Li, H. (2012). seqtk Toolkit for processing sequences in FASTA/Q formats. GitHub, 767, 69.

        15. UNOISE algorithm - R.C. Edgar (2016). UNOISE2: improved error-correction for Illumina 16S and ITS amplicon sequencing, https://doi.org/10.1101/081257
