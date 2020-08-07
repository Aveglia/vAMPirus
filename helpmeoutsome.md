****************************************************************************************************
                                                            ,---.,-.-.,---.o
                                                      .    ,|---|| | ||---'.,---..   .,---.
                                                       \  / |   || | ||    ||    |   |`---.
                                                        `'  `   `` ` ``    ``    `---``---`
                                             An automated virus amplicon sequencing analysis pipeline
*****************************************************************************************************
# Introduction to vAMPirus

The main motive behind vAMPirus is to provide a robust and easy-to-use bioinformatics workflow for virus amplicon sequencing analysis. The vAMPirus workflow
allows easy reproducibility of project-specific analyses and is flexible enough to tailor your analysis to your own data.

## Contact/support:

Please contact Alex Veglia at ajv5@rice.edu with any feedback or questions. Any kind of input from the community is welcomed and encouraged!

## How to cite:

If you do use vAMPirus for your analyses, please cite using the following ->

Veglia, A.J. *et.al.,* 2020



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

[This will check for Nextflow and Anaconda/Miniconda installations, then download specified database(s), which in this case is NCBIs RefSeq database.]

It should be noted, that any database can be used, but it needs to be in fasta format and the headers for reference sequences need to match
one of two patterns:

RVDB format (default) -> ">acc|GENBANK|AYD68780.1|GENBANK|MH171300|structural polyprotein [Marine RNA virus BC-4]"

NCBI NR/RefSeq format -> ">KJX92028.1 hypothetical protein TI39_contig5958g00003 [Zymoseptoria brevis]"

During Taxonomy Assignment, vAMPirus infers results by extracting the information stored in the reference sequence headers. If the database sequence headers do not match these
patterns, you are bound to see errors in the naming of files created during the Taxonomy Assignment phase of vAMPirus.

The default is that vAMPirus assumes that the database headers are in RVDB format, to change this assumption, you would need to edit the configuration file at line 78 where "refseq=F". Change the "F" to "T" and
you are good to go! You could also change this within the launch command with adding "--refseq T", but setting parameters will be discussed further in a section later.


## Testing vAMPirus installation

A test dataset is provided in the vAMPirus/example_data. To ensure that vAMPirus is set up properly before running with your own data, you can run:

`nextflow run vAMPirusv0.1.0.nf -c ./example_data/vampirus_test.config -with-conda /PATH/TO/miniconda3/env/vAMPirus`

# Quick Notes Before Running vAMPirus

## The Nextflow workflow manager and the "launch command"

vAMPirus is deployed using the Nextflow workflow manager which "enables scalable and reproducible scientific workflows using software containers. It allows the
adaptation of pipelines written in the most common scripting languages. Its fluent DSL simplifies the implementation and the deployment of complex parallel and reactive
workflows on clouds and clusters." To learn more about Nextflow and to learn more how to monitor your submitted jobs from a web portal with Nextflow Tower, visit nextflow.io.

The great thing about vAMPirus being integrated into Nextflow is that it is just as easy to run vAMPirus on a HPC ad it is to run locally on your personal machine. This is
where the configuration file comes into play and we will discuss this in more detail in the next section. First, lets understand the basics of how to launch vAMPirus using
Nextflow.

Here is a basic "launch command" to deploy the vAMPirus pipeline (we will talk more about the mandatory/optional arguments of vAMPirus later):

`nextflow run vAMPirusv0.1.0.nf -c vampirus.config -with-conda /PATH/TO/miniconda3/env/vAMPirus --Analyze`

In the command above, there are five necessary pieces of information needed to successfully launch the vAMPirus workflow:

    1. The first is the "nextflow" executable in the beginning. Nextflow is responsible for launching the vAMPirus processes and this is why nextflow must be either added
       to your $PATH variable, specify the path to your nextflow executable in the command itself, or copy the executable to your working directory.

    2. Second, you must tell Nextflow to "run" the vAMPirus program which is described in the "vAMPirusv0.1.0.nf" file. Depending on where you plan to submit this command,
       you may have to specify the path to the vAMPirusv0.1.0.nf file or you can copy the file to your working directory.

    3. Next, we need to tell Nextflow what configuration file we would like to use for the vAMPirus run which is illustrated by the "-c vampirus.config" segment of the command.
       We will talk more about the configuration file in the next section.

    4. The next piece of information Nextflow needs is the environment to use for the vAMPirus workflow. The dependencies of vAMPirus are stored as a conda environment so we need
       to tell Nextflow that we would like to run vAMPirus with conda and specifically the vAMPirus environment that was built when running the vampirus_startup.sh script during
       installation. To find the path to the vAMPirus conda environment, you can run "conda info --envs" and this will give you the path to your conda environments.

    5. And finally, there are a few different mandatory arguments that one can set when running vAMPirus to specify the type of analysis the user would like to run (we will go over these
       later in this documentation). In this example, we are using "--Analyze" which includes ASV and AminoType generation and all subsequent analyses (Taxonomy Assignment, Phylogeny, EMBOSS, etc.).

Now that we have an understanding on how to deploy vAMPirus with Nextflow, lets look at how to set both analysis- and resource-related parameters for your vAMPirus runs.

## Understanding the vAMPirus config file and setting parameters

### The configuration file (vampirus.config)

Nextflow deployment of vAMPirus relies on the use of the configuration file (vampirus.config) that is found in the vAMPirus program directory. The configuration file is a great way to store parameters/options
used in your analyses. It also makes it pretty easy to set and keep track of multiple parameters as well as storing custom default values that you feel work best for your data. You can also have multiple copies
of vampirus configuration files with different parameters, you would just have to specify the correct file with the "-c" argument shown in the section before.

Furthermore, the configuration file contains analysis-specific parameters AND resource-specific Nextflow launching parameters. A benefit of Nextflow integration, is that you can run the vAMPirus workflow on a large
HPC just as easily as you could on your local machine. If you look at line 151 and greater in the vampirus.config file, you will see resource-specific parameters that you can alter before any run. Nexflow is capable
of submitting jobs automatically using slurm and PBS, check out the Nextflow docs to learn more (https://www.nextflow.io/docs/latest/index.html)!

### Setting parameter values

There are two ways to set parameters with Nextflow and vAMPirus:

    1. Just edit the config file!

        Here we have a block from the vampirus.config file that stores information related to your run:

            // Project/analyses- specific information
                // Project name - Name that will be used as a prefix for namng files by vAMPirus
                     projtag="vAMPrun"
                // Path to metadata spreadsheet file to be used for plot
                     metadata="/data/alex/PVID_dinorna/AMPS/testvamp/vAMPirus/fisces_metdata.csv"
                // Minimum number of hit counts for a sample to have to be included in the downstream analyses and report generation
        	         minimumCounts="1000"
        	    // PATH to current working directory that contains (i) the vAMPirus.nf script, (ii)
                     mypwd="/data/alex/PVID_dinorna/AMPS/testvamp/vAMPirus"
                     email="your_email@web.com"
                // reads directory
                     reads="/data/alex/PVID_dinorna/AMPS/testvamp/vAMPirus/reads/R{1,2}_001.fastq.gz"
                // Directory to store output of vAMPirus analyses
                     outdir="results"

        The first one in the block is the project tag or "projtag" which by default, if unchanged, will use the prefix "vAMPrun". To change this value, and any other parameter value, just edit right in the configuration
        file so if you wanted to call the run "VirusRun1" you would edit the line to:

        // Project/analyses- specific information
            // Project name - Name that will be used as a prefix for naming files by vAMPirus
                 projtag="VirusRun1"


    2. Set the value within the launch command itself!

        Instead of editing the configuration file directly, you could set parmeters within the launching command itself. So, for example, if we wanted to run the analysis with nucletide-based clustering of ASVs at 95%
        similarity, you would do so like this:

                `nextflow run vAMPirusv0.1.0.nf -c vampirus.config -with-conda /PATH/TO/miniconda3/env/vAMPirus --Analyze --nOTU --clusterNuclID .95`

        Here we use the "--Analyze" option that tells vAMPirus that we are ready to analyze soem data. Then the "--nOTU" argument with the "--clisterNuclID .95" tells vAMPirus we would like to cluster our ASVs based on
        95% nucleotide similarity. The default ID value is stored at line 51 in the vampirus.config file (currently 85%), but as soon as you specify and provide a value in the command, the default value is overwritten.

### vAMPirus skip options

To specify certain parts of the vAMPirus workflow to perform in a given run, you can use skip options to have vAMPirus ignore certain processes. Here are the current skip options you can specify within the launch command:

    *// Skip options
        // Skip all Read Processing
            skipReadProcessing=false
        // Skip quality control processes only
            skipFastQC = false
        // Skip adapter removal process only
            skipAdapterRemoval=false
        // Skip primer removal process only
            skipPrimerRemoval=false
        // Skip AminoTyping
            skipAminoTyping=false
        // Skip Taxonomy
            skipTaxonomy=false
        // Skip phylogeny
            skipPhylogeny = false
        // Skip EMBOSS analyses
            skipEMBOSS = false*

To utilize these skip options is pretty simple where you would just add it to the launch command like so:

    `nextflow run vAMPirusv0.1.0.nf -c vampirus.config -with-conda /PATH/TO/miniconda3/env/vAMPirus --Analyze --nOTU --clusterNuclID .95 --skipPhylogeny --skipTaxonomy`

With this launch command, vAMPirus will perform ASV generation and nucleotide-based clustering to produce nOTUs, then will generate counts tables, matrices and the final report for you.


# Running the vAMPirus workflow

## For the impatient

Once you have everything set up and you have edited the parameters of interest in your configuration file you can run the following launch command for a full analysis:

    `nextflow run vAMPirusv0.1.0.nf -c vampirus.config -with-conda /PATH/TO/miniconda3/env/vAMPirus --Analyze --nOTU --pOTU`

This launch command will run all aspects of the vAMPirus workflow on your data and spit out final reports for each clustering technique.

## Necessary input

### Sequencing reads

Input can be raw or processed gzipped fastq files with names containing *"_R1"* or *"_R2"*. You can specify the directory containing your reads in line 25 of the vampirus.config file.

NOTE: Sample names are extracted from read library names by using the string to the left of the *"_R"* in the filename automatically.

### Metadata file

For every analysis, vAMPirus generates a final report and uses a user supplied metadata file with sample names and treatment. Treatment is how vAMPirus groups samples in downstream statistical
analyses performed to generate for the final report. For example, if comparing samples from different species of corals, you would set up a metadata file like so:

        sample,treatment
        Coral1,Ofaveolata
        Coral2,Ofaveolata
        Coral3,Ofaveolata
        Coral4,Mcavernosa
        Coral5,Mcavernosa
        Coral6,Mcavernosa

The metadata file needs to be comma separated with the first column being "sample" and the second column must be "treatment". These species names could easily be replaced with "Heat"/"Control".

## The mandatory arguments

To run the vAMPirus workflow, you must specify one or two mandatory arguments:

    1. "--DataCheck"

        Usage example:

            `nextflow run vAMPirusv0.1.0.nf -c vampirus.config -with-conda /PATH/TO/miniconda3/env/vAMPirus --DataCheck`

        The DataCheck feature of vAMPirus is meant to give the user some information about their data so they can tailor their final analysis appropriately. In DataCheck mode, vAMPirus
        performs all read processing operations then generates ASVS and performs nucleotide- and p rotein-based clustering at 24 different clustering percentages ranging from 55-99% ID.
        vAMPirus then generates an html report that displays and visualizes read processing and clustering stats. It is recommended that before running any dataset through vAMPirus, you
        run the data through the DataCheck to allow informed decisions on clustering percentage, max read length, etc., etc.

        Here is how Nextflow will display the list of processes vAMPirus will execute during DataCheck (executed with the launch command above):

                executor >  local (57)
                [8a/75e048] process > Build_database                                  [100%] 1 of 1 ✔
                [b4/c08216] process > QualityCheck_1DC                                [100%] 9 of 9 ✔
                [5c/ed7618] process > Adapter_Removal_DC                              [100%] 9 of 9 ✔
                [25/56a27d] process > Primer_Removal_DC                               [100%] 9 of 9 ✔
                [fa/736d66] process > QualityCheck_2_DC                               [100%] 9 of 9 ✔
                [d8/4b0c4e] process > Read_Merging_DC                                 [100%] 9 of 9 ✔
                [93/f5d315] process > Compile_Reads_DC                                [100%] 1 of 1 ✔
                [4d/8d83dd] process > Compile_Names_DC                                [100%] 1 of 1 ✔
                [0a/fdd8e8] process > Length_Filtering_DC                             [100%] 1 of 1 ✔
                [b6/097dd8] process > Extract_Uniques_DC                              [100%] 1 of 1 ✔
                [1b/1c4476] process > Identify_ASVs_DC                                [100%] 1 of 1 ✔
                [2e/9101c3] process > Chimera_Check_DC                                [100%] 1 of 1 ✔
                [14/365745] process > NucleotideBased_ASV_clustering_DC               [100%] 1 of 1 ✔
                [99/938f26] process > Translation_For_ProteinBased_Clustering_DC      [100%] 1 of 1 ✔
                [34/9eb77a] process > Protein_clustering_DC                           [100%] 1 of 1 ✔
                [26/1143ba] process > combine_csv_DC                                  [100%] 1 of 1 ✔
                [2e/e5fea3] process > Report_DataCheck                                [100%] 1 of 1 ✔

        This screen is displayed uponEach column contains different information related to progress and execution of the analysis.




    2. "--Analyze" -- We've seen this one before

        Usage example:

        *nextflow run vAMPirusv0.1.0.nf -c vampirus.config -with-conda /PATH/TO/miniconda3/env/vAMPirus --Analyze*

        The Analyze option allows vAMPirus to know that you plan to analyze your data with the given parameters either within the launch command or sourced from the configuration file. On
        its own, "--Analyze" will run all read processing operations, generate ASVs, ASV counts files/matrices, ASV phylogeny, ASV taxonomy assignment, generate Aminotypes, Aminotype counts/matrices,
        Aminotype phylogeny, Aminotype taxonomy assignment and EMBOSS analyses. vAMPirus will also produce final reports for ASV and Aminotype analyses.

        To generate nOTUs (nucleotide-based OTUs) or pOTUs (protein-based OTUs) and run all subsequent analyses with them (phylogeny, taxonomy assignment), you would just need to add
