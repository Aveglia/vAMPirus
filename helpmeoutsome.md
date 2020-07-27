                                                ,---.,-.-.,---.o
                                          .    ,|---|| | ||---'.,---..   .,---.
                                           \  / |   || | ||    ||    |   |`---.
                                            `'  `   `` ` ``    ``    `---``---`
                                An automated virus amplicon sequencing analysis pipeline

## Overview


## Getting started

### Installing vAMPirus

Clone the most recent version of vAMPirus from github using:

`git clone https://github.com/Aveglia/vAMPirus.git`

### Setting up vAMPirus dependencies and checking installation

To deploy vAMPirus, you need to have Nextflow and Anaconda/Miniconda installed on your system. To set up and install vAMPirus
dependencies, simply move to the vAMPirus directory and run the vampirus_startup.sh script.

`cd ./vAMPirus`
`./vampirus_startup.sh`

The start up script will check your system for Nextflow and Anaconda/Miniconda and if they are not present, the script will ask
if you would like to install these programs. If you answer with 'y', the script will install the missing programs and will build
the vAMPirus Conda environment.

You can also use the script to install different databases which include:

1. NCBIs Viral protein RefSeq database
2. The proteic version of the Reference Viral DataBase (RVDB) (See https://f1000research.com/articles/8-530)
3. The complete NCBI NR protein database

To use the vampirus_startup.sh script to download any or all of these databases listed above you just need to use the "-d" option.
If we look at the script usage guidelines:


*General execution:

vampirus_startup.sh -h -d [1|2|3|4]


    Command line options:

        [ -h ]                       	Print help information

        [ -d 1|2|3|4 ]                  Set this option to create a database directiory within the current working directory and download the following databases for taxonomy assignment:

                                                    1 - Download only NCBIs Viral protein RefSeq database
                                                    2 - Download the proteic version of the Reference Viral DataBase (See the paper for more information on this database: https://f1000research.com/articles/8-530)
                                                    3 - Download only the complete NCBI NR protein database
                                                    4 - Download all three databases

**

So, if we wanted to download NCBIs Viral protein RefSeq database, we would just need to run:

`./vampirus_startup.sh -d 1`
[^1] This will check for Nextflow and Conda installations, then download specified database(s) which in this case is NCBIs RefSeq database.

It should be noted, that any database can be used, but it needs to be in fasta format and the headers for reference sequences need to match
one of two patterns:

RVDB format (default) -> ">acc|GENBANK|AYD68780.1|GENBANK|MH171300|structural polyprotein [Marine RNA virus BC-4]"

NCBI NR/RefSeq format -> ">KJX92028.1 hypothetical protein TI39_contig5958g00003 [Zymoseptoria brevis]"

During Taxonomy Assignment, vAMPirus extacts results using the information stored in the headers of reference sequence that are matched. If the
database headers do not match these patterns, you are bound to see errors in the naming of files created during the Taxonomy Assignment phase of vAMPirus.

### Testing vAMPirus installation

A test dataset is provided in the vAMPirus/example_data. To ensure that vAMPirus is set up properly before running with your own data, you can run:

`nextflow run vAMPirusv0.1.0.nf -c ./example_data/vampirus_test.config -with-conda /PATH/TO/miniconda3/env/vAMPirus`

## Running vAMPirus

vAMPirus is deployed using Nextflow which "enables scalable and reproducible scientific workflows using software containers. It allows the adaptation of
pipelines written in the most common scripting languages. Its fluent DSL simplifies the implementation and the deployment of complex parallel and reactive
workflows on clouds and clusters."

To learn more about Nextflow, visit Nextflow enables scalable an

### Understanding the vAMPirus config file and

The benefit of using Nextflow
