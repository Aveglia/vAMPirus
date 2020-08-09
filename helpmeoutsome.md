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

## Order of operations

    1. Clone vAMPirus from github
    2. Run the vAMPirus start up script to download Nextflow and create conda environment
    3. Edit vAMPirus configuration file
    4. Run DataCheck mode with dataset
    5. Run Analyze mode with desired clustering technique and %ID

## Contact/support:

Please contact Alex Veglia at ajv5@rice.edu with any feedback or questions. Any kind of input from the community is welcomed and encouraged!

## How to cite:

If you do use vAMPirus for your analyses, please cite the following ->

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

`nextflow run vAMPirusv0.1.0.nf -c ./example_data/vampirus_test.config -with-conda /PATH/TO/miniconda3/env/vAMPirus --Analyze --testing`

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

## The Nextflow output

When submitting a launch command to start your vAMPirus run, Nextflow will spit out something that looks like this (output for DataCheck mode):

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

Nextflow allows for interactive monitoring of submitted workflows, so in this example, we see the left column containing working directories for each process being executed, next to that we see the process name,
and the final column on the right contains the status and success of each process. In this example each process has been executed successfully and has been cached. The amazing thing about Nextflow is that say you received
an error or decide you would like to change a parameter/add a type of clustering you can use the Nextflow "-resume" option that means you don't re-run already completed processes.

## Understanding the vAMPirus config file and setting parameters

### The configuration file (vampirus.config)

Nextflow deployment of vAMPirus relies on the use of the configuration file (vampirus.config) that is found in the vAMPirus program directory. The configuration file is a great way to store parameters/options
used in your analyses. It also makes it pretty easy to set and keep track of multiple parameters as well as storing custom default values that you feel work best for your data. You can also have multiple copies
of vAMPirus configuration files with different parameters, you would just have to specify the correct file with the "-c" argument shown in the section before.

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

NOTE: Nextflow also has options in the launch command. To tell them apart, Nextflow options uses a single dash (e.g. -with-conda) while vAMPirus options are always with a double dash (e.g. --Analyze)


### vAMPirus skip options

To specify certain parts of the vAMPirus workflow to perform in a given run, you can use skip options to have vAMPirus ignore certain processes. Here are the current skip options you can specify within the launch command:

            // Skip options
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
                    skipEMBOSS = false
                // Skip Reports
                    skipReports = false

To utilize these skip options is pretty simple where you would just add it to the launch command like so:

    `nextflow run vAMPirusv0.1.0.nf -c vampirus.config -with-conda /PATH/TO/miniconda3/env/vAMPirus --Analyze --nOTU --clusterNuclID .95 --skipPhylogeny --skipTaxonomy`

With this launch command, vAMPirus will perform ASV generation and nucleotide-based clustering to produce nOTUs, then will generate counts tables, matrices and the final report for you.


# Running the vAMPirus workflow

## Recommended order of operations

    1. Clone vAMPirus from github
    2. Run the vAMPirus start up script to download Nextflow and create conda environment
    3. Edit vAMPirus configuration file
    4. Run DataCheck mode with dataset
    5. Run Analyze mode with desired clustering technique and %ID

## For the impatient

Once you have everything set up and you have edited the parameters of interest in your configuration file you can run the following launch command for a full analysis:

    `nextflow run vAMPirusv0.1.0.nf -c vampirus.config -with-conda /PATH/TO/miniconda3/env/vAMPirus --Analyze --nOTU --pOTU`

This launch command will run all aspects of the vAMPirus workflow on your data and spit out final reports for each clustering %ID and technique.

## Necessary input

### Sequencing reads

Input can be raw or processed gzipped fastq files with names containing "\_R1" or "\_R2". You can specify the directory containing your reads in line 25 of the vampirus.config file.

NOTE: Sample names are extracted from read library names by using the string to the left of the "\_R" in the filename automatically.

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

The metadata file needs to be comma separated with the first column being "sample" and the second column must be "treatment". These species names could easily be replaced with "Heat"/"Control"
or any other way you would like to categorize the samples.

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

        Every time you launch vAMPirus with Nextflow, you will see this kind of output that refreshes with the status of the different processes during the run.


    2. "--Analyze" -- We've seen this one before

        Usage example:

        `nextflow run vAMPirusv0.1.0.nf -c vampirus.config -with-conda /PATH/TO/miniconda3/env/vAMPirus --Analyze`

        Example Nextflow output for this launch command:

                executor >  local (8)
                [8a/75e048] process > Build_database                      [100%] 1 of 1 ✔
                [8f/5bf47f] process > QualityCheck_1                      [100%] 9 of 9 ✔
                [1e/d40d5f] process > Adapter_Removal                     [100%] 9 of 9 ✔
                [-        ] process > Primer_Removal                      [100%] 9 of 9 ✔
                [-        ] process > QualityCheck_2                      [100%] 9 of 9 ✔
                [-        ] process > Read_Merging                        [100%] 9 of 9 ✔
                [-        ] process > Compile_Reads                       [100%] 1 of 1 ✔
                [-        ] process > Compile_Names                       [100%] 1 of 1 ✔
                [-        ] process > Length_Filtering                    [100%] 1 of 1 ✔
                [-        ] process > Extract_Uniques                     [100%] 1 of 1 ✔
                [-        ] process > Identify_ASVs                       [100%] 1 of 1 ✔
                [-        ] process > Chimera_Check                       [100%] 1 of 1 ✔
                [-        ] process > ASV_Taxonomy_Assignment             [100%] 1 of 1 ✔
                [-        ] process > Generate_ASV_Counts_Tables          [100%] 1 of 1 ✔
                [-        ] process > Generate_ASV_Matrix                 [100%] 1 of 1 ✔
                [-        ] process > ASV_Phylogeny                       [100%] 1 of 1 ✔
                [-        ] process > Translate_For_AminoTyping           [100%] 1 of 1 ✔
                [-        ] process > Generate_AminoTypes                 [100%] 1 of 1 ✔
                [-        ] process > Generate_AminoType_Matrix           [100%] 1 of 1 ✔
                [-        ] process > AminoType_EMBOSS_Analyses           [100%] 1 of 1 ✔
                [-        ] process > Taxonomy_Assignment_AminoTypes      [100%] 1 of 1 ✔
                [-        ] process > AminoType_Phylogeny                 [100%] 1 of 1 ✔
                [-        ] process > Generate_AminoTypes_Counts_Table    [100%] 1 of 1 ✔
                [-        ] process > combine_csv                         [100%] 1 of 1 ✔
                [-        ] process > Report_ASVs                         [100%] 1 of 1 ✔
                [-        ] process > Report_AmynoType                    [100%] 1 of 1 ✔


        The Analyze option allows vAMPirus to know that you plan to analyze your data with the given parameters either within the launch command or sourced from the configuration file. On
        its own, "--Analyze" will run all read processing operations, generate ASVs, ASV counts files/matrices, ASV phylogeny, ASV taxonomy assignment, generate Aminotypes, Aminotype counts/matrices,
        Aminotype phylogeny, Aminotype taxonomy assignment and EMBOSS analyses. vAMPirus will also produce final reports for ASV and Aminotype analyses.

        To generate nOTUs (nucleotide-based OTUs) or pOTUs (protein-based OTUs) and run all subsequent analyses with them (phylogeny, taxonomy assignment), you would just need to add the "--nOTU" and "--pOTU"
        options like so:

        `nextflow run vAMPirusv0.1.0.nf -c vampirus.config -with-conda /PATH/TO/miniconda3/env/vAMPirus --Analyze --nOTU --pOTU`

        Here is what the Nextflow output would look like for this launch command:

                executor >  local (8)
                [-        ] process > Build_database                      [  0%] 0 of 1
                [-        ] process > QualityCheck_1                      -
                [-        ] process > Adapter_Removal                     -
                [-        ] process > Primer_Removal                      -
                [-        ] process > QualityCheck_2                      -
                [-        ] process > Read_Merging                        -
                [-        ] process > Compile_Reads                       -
                [-        ] process > Compile_Names                       -
                [-        ] process > Length_Filtering                    -
                [-        ] process > Extract_Uniques                     -
                [-        ] process > Identify_ASVs                       -
                [-        ] process > Chimera_Check                       -
                [-        ] process > NucleotideBased_ASV_clustering      -
                [-        ] process > Nucleotide_Taxonomy_Assignment      -
                [-        ] process > Generate_Counts_Tables_Nucleotide   -
                [-        ] process > Generate_Nucleotide_Matrix          -
                [-        ] process > Nucleotide_Phylogeny                -
                [-        ] process > Translate_For_AminoTyping           -
                [-        ] process > Generate_AminoTypes                 -
                [-        ] process > Generate_AminoType_Matrix           -
                [-        ] process > AminoType_EMBOSS_Analyses           -
                [-        ] process > Taxonomy_Assignment_AminoTypes      -
                [-        ] process > AminoType_Phylogeny                 -
                [-        ] process > Generate_AminoTypes_Counts_Table    -
                [-        ] process > Translation_For_pOTU_Generation     -
                [-        ] process > Generate_pOTUs                      -
                [-        ] process > pOTU_Nucleotide_Taxonomy_Assignment -
                [-        ] process > Generate_Nucleotide_pOTU_Counts     -
                [-        ] process > Generate_pOTU_Nucleotide_Matrix     -
                [-        ] process > pOTU_Nucleotide_Phylogeny           -
                [-        ] process > pOTU_AminoAcid_Matrix               -
                [-        ] process > pOTU_EMBOSS_Analyses                -
                [-        ] process > pOTU_AminoAcid_Taxonomy_Assignment  -
                [-        ] process > pOTU_Protein_Phylogeny              -
                [-        ] process > Generate_pOTU_Protein_Counts        -
                [-        ] process > combine_csv                         -
                [-        ] process > Report_ASV                          -
                [-        ] process > Report_nOTU                         -
                [-        ] process > Report_AmynoTypes                   -
                [-        ] process > Report_pOTU_AminoAcid               -
                [-        ] process > Report_pOTU_Nucleotide              -

        You can notice that there are a few more processes now compared to the output of the previous launch command which is what we expect since we are asking vAMPirus to do a little bit more work for us :).

# Breaking it down: The vAMPirus workflow

## Read processing

The read processing segment of vAMPirus includes FastQC report generation, adapter removal with fastp, primer removal with bbduk.sh, read merging with vsearch, and a final length filtering/global trimming with fastp and bbduk.sh.

### Adapter removal with fastp

Adapter contamination in reads are automatically detected with fastp. Overrepresentation analysis and quality filtering is also ran during this process prior to primer removal. Adapter contamination removal can be skipped using the
"--skipAdapterRemoval" option, but is not recommended as even if you already performed adapter removal, it doesn't hurt to check again.

### Primer removal with bbduk.sh

There are two ways that vAMPirus is able to remove primer sequences from reads with bbduk.sh:

    1. Primer removal by chopping off specified number of bases from each read -

            This is the default action of vAMPirus if no primer sequences are specified. To use this method of primer removal, you must specify "--GlobTrim" and either specify the number of bases in the launch command like so:

                `nextflow run vAMPirusv0.1.0.nf -c vampirus.config -with-conda /PATH/TO/miniconda3/env/vAMPirus --Analyze --nOTU --pOTU --GlobTrim 23,26`

            In this situation, we tell vAMPirus we would like to remove primers by trimming 23 bases from the forward reads and 26 bases from the reverse reads. The other way to initiate this method of primer removal is to add
             the same information at line 39 in the configuration file:

             // Primer Removal parameters
                 // If not specifying primer sequences, forward and reverse reads will be trimmed by number of bases specified using --GlobTrim #basesfromforward,#basesfromreverse
                     GlobTrim="23,26"

            By adding the information to line 39, vAMPirus will automatically use this method and these parameters for primer removal. If you want to change the number of bases without editing the configuration file, all you would
            need to do is then specify in the launch command with "--GlobTrim 20,27" and vAMPirus will ignore the "23,26" in the configuration file.

            NOTE: Specifying global trimming by editing line 39 in the config file or using "--GlobTrim" in the launch command will override the use of primer sequences for removal if both are specified  

    2. Primer removal by specifying primer sequences

            To signal to vAMPirus to utilize this method of primer removal, similar to the previous method, you could edit the configuration file or specify within the launch command. To specify in launch command, we would need to used the "--fwd" and "--rev" options:

                `nextflow run vAMPirusv0.1.0.nf -c vampirus.config -with-conda /PATH/TO/miniconda3/env/vAMPirus --Analyze --nOTU --pOTU --fwd FWDPRIMER --rev REVPRIMER`

            vAMPirus will then provide these sequences to bbduk.sh for them to be removed. The primer sequences could also be stored in the configuration file in lines 40-43:

                // Specific primer sequence on forward reads to be removed
                    fwd="FWDPRIMER"
                // Reverse primer sequence
                    rev="REVPRIMER"


### Read merging and length filtering

Read merging in the vAMPirus workflow is performed by vsearch and afterwards, reads are trimmed to the expected amplicon length (--maxLen) and any reads with lengths below the user specified minimum read length (--minLen). There are three parameters that you can edit to influence this segment of vAMPirus. If we look at lines 29-35:

    // Merged read length filtering parameters
        // Minimum merged read length - reads below the specified maximum read length will be used for counts only
            minLen="400"
        // Maximum merged read length - reads with length equal to the specified max read length will be used to generate uniques and ASVs
            maxLen="422"
        // Maximum expected error for vsearch merge command
            maxEE="1"

The user can edit the minimum length (--minLen) for reads to be used for counts table generation, maximum length (--maxLen) for reads used to generate uniques and subsequent ASVs, and the expected error rate (--maxEE) for overlapping region of reads during read merging with vsearch. The values above are default and should be edited before running your data with --Analyze.

This is where the DataCheck report is very useful, you can review the report and see the number of reads that merge per library and you can edit the expected error value to be less stringent if needed. The DataCheck report also
contains a read length distribution that you can use to select an ideal maximum/minimum read length.


## Amplicon Sequence Variants, AminoTypes and Operational Taxonomic Units

The goal of vAMPirus was to make it easy for the user to analyze their data is many different ways to potentially reveal patterns that would have been missed if limited to one method/pipeline. A major and sometimes difficult step in analyzing virus amplicon sequence data is deciding the method to use for identifying or defining different viral "species" in the data. To aid this process, vAMPirus has the DataCheck mode discussed above and has several different options for sequence clustering/analysis for the user to decide between.  

 vAMPirus relies on vsearch using the UNOISE3 algorithm to generate Amplicon Sequencing Variants (ASVs) from dereplicated amplicon reads. ASVs are always generated by default and there are two parameters that the user can specify either in the launch command or by editing the configuration file at lines 45-49:

     // ASV generation and clustering parameters
         // Alpha value for denoising - the higher the alpha the higher the chance of false positives in ASV generation (1 or 2)
             alpha="1"
         // Minimum size or representation for sequence to be considered in ASV generation
             minSize="8"

The smaller the alpha value (--alpha) the more stringent vsearch is ASV generation while minimum size is the minimum representation of a unique sequence to have to be considered in the ASV generation.

Now, onto clustering ASVs into Operation Taxonomic Units (OTUs). vAMPirus is able to use two different techniques for generating OTUs for the user:


    1. AminoTyping -

        vAMPirus by default, unless the --skipAminoTyping option is set, will generate unique amino acid sequences or  "AminoTypes" from generated ASVs. These AminoTypes, barring any skip options set, will run through all the
        same analyses as ASVs.

        vAMPirus will translate the ASVs with Virtual Ribosome and relies on the user to specify the expected or minimum amino acid sequence length (--minAA) to be used for AminoTyping and pOTU generation (discussed below). For example, if you amplicon size is ~422 bp long, you would expect the amino acid translations to be ~140. Thus, you would either edit the --minAA value to be 140 in the configuration file (line 59) or in the launch command. You can make it shorter if you would like, but based on personal observation, a shorter translation is usually the result of stop codons present which would usually be removed from subsequent analysis. If there are any sequences below the minimum amino acid sequence length the problematic sequence(s) and its translation will be stored in a directory for you to review.     

        Example launch command:

            `nextflow run vAMPirusv0.1.0.nf -c vampirus.config -with-conda /PATH/TO/miniconda3/env/vAMPirus --Analyze --minAA 140`


    2. Nucleotide-based clustering into OTUS (nOTUs) -

        In this technique, as the name infers, the clustering of ASVs into OTUs is based on nucleotide identity. To signal vAMPirus to generate nOTUs, we just add the "--nOTU" option to our launch command:

                `nextflow run vAMPirusv0.1.0.nf -c vampirus.config -with-conda /PATH/TO/miniconda3/env/vAMPirus --Analyze --nOTU`

        With nOTU analysis specified, vAMPirus will still generate ASVs and produce ASV-related results with a report, so you are not losing any results when you the clustering and this way, you can compare the reports afterwards.

        vAMPirus allows you to specify either one or multiple clustering percentages that can be specified in the configuration file or the launch command. For each ID% specified, nOTUs will be generated and will be ran through
        all subsequent processes (e.g. counts, phylogeny, taxonomy, report generation, etc.)

            // Percent similarity to cluster nucleotide ASV sequences
                clusterNuclID=".85"
            // List of percent similarities to cluster nucleotide ASV sequences - must be separated by a comma (ex. ".95,.96")
                clusterNuclIDlist=""

        Above, the specified clustering percentage (--clusterNuclID) for nOTU generation is 85%. The value must always be in decimal format and can be put into the launch command like so:

            `nextflow run vAMPirusv0.1.0.nf -c vampirus.config -with-conda /PATH/TO/miniconda3/env/vAMPirus --Analyze --nOTU --clusterNuclID .85`

        You could also have a list of percentages to cluster by:

            `nextflow run vAMPirusv0.1.0.nf -c vampirus.config -with-conda /PATH/TO/miniconda3/env/vAMPirus --Analyze --nOTU --clusterNuclIDlist .85,.90,.96`

        Using "--clusterNuclIDlist" will override the single percentage clustering and vAMPirus will generate separate nOTUs fastas based on 85%, 90% and 96% nucleotide similarity, you could theoretically cluster by any amount
        of percentages between 1-100 that your data requires or your heart desires.


    3. Protein-based clustering into OTUs (pOTUs)

        For this clustering method, ASVs are translated into amino acid sequences and clustered based on shared amino acid similarity, %ID is user specified in a similar manner to nOTU cluster options:

            Lines 54-57 in vampirus.config ->

            // Default percent similarity to cluster aminoacid sequences
                clusterAAID=".97"
            // List of percent similarities to cluster aminoacid sequences - must be separated by ".95,.96"
                clusterAAIDlist=""

        When clustering for pOTUs, vAMPirus will create nucleotide and protein versions and both will go through the rest of the analyses within the pipeline (phylogeny, taxonomy, etc.). It should be noted that the --minAA option
        described above in the AminoType section applies also to pOTU generation. If there are any sequences below the minimum amino acid sequence length, the problematic sequence(s) and its translation will be stored in a file for
        you to review.

        Example launch command:

            `nextflow run vAMPirusv0.1.0.nf -c vampirus.config -with-conda /PATH/TO/miniconda3/env/vAMPirus --Analyze --pOTU --clusterAAIDlist .85,.90,.96`

## Counts tables and percent ID matrices

vAMPirus generate nucleotide-based counts tables using vsearch and protein-based counts tables using Diamond. Counts tables and percent ID matrices are always produced for each ASV, AminoType and all OTU fasta files produced. Here
are the parameters you can edit at lines 61-70:

    // Counts table generation parameters
        // Similarity ID to use for ASV counts table
            asvcountID=".97"
        // Protein counts table parameters
            // Minimum Bitscore for counts
                ProtCountsBit="50"
            // Minimum aminoacid sequence similarity for hit to count
                ProtCountID="85"
            // Minimum alignment length for hit to count
                ProtCountsLength="50"

The "--asvcountID" is the percent ID during global alignment that vsearch would use to classify a hit to a given ASV. For all OTUs, the specified clustering percentage is used for counts. Protein-based counts file generation has a
few more parameters the user can alter: "--ProtCountsBit" is the minimum bitscore for an alignment to be recorded, "--ProtCountID" is the minimum percent amino acid similarity an alignment needs to have to be recorded, and
"--ProtCountsLength" is the minimum alignment length for a hit to be recorded.


## Phylogenetic analysis and model testing

Phylogenetic trees are produced automatically for ASVs (unless --nOTU specified), nOTUs, pOTUs and AminoTypes using IQ-TREE. All produced sequence fastas are aligned using the MAFFT algorithm then alignments are trimmed
automatically using TrimAl. Post alignment and trimming, there is some flexibility in this process where you can specify a few different aspects of the analysis:

### Substitution model testing

ModelTest-NG is always ran to determine the substitution model and generate a starting tree to be used in the IQ-TREE command for all fasta files generated in the pipeline. The user decides whether they would like to use the model
determined by ModelTest-NG in the IQ-TREE command or they can have ModelFinder Plus determine the model instead. By default, the ModelFinder Plus determined model is used, however, to trigger the use of the model determined by
ModelTest-NG for both nucleotide and amino acid model testing, you would use "--ModelTnt" and "--ModelTaa". So:

        Example launch command using ModelTest-NG ->

        `nextflow run vAMPirusv0.1.0.nf -c vampirus.config -with-conda /PATH/TO/miniconda3/env/vAMPirus --Analyze --pOTU --clusterAAIDlist .85,.90,.96 --ModelTnt --ModelTaa`

If you didn't want to use it for both, you would just use the "--ModelTnt" or "--ModelTaa" alone in the launch command.

### Bootstrapping

IQ-TREE is capable of performing parametric or non-parametric bootstrapping. You can specify which one using "--parametric" or "--nonparametric" and to set how many boostraps to perform, you would use "--boots #ofbootstraps" or edit
 line 92 in the vampirus.config file. Here is an example for creating a tree using the model determined by ModelTest-NG, non-parametric boostrapping and 500 bootstraps:

        `nextflow run vAMPirusv0.1.0.nf -c vampirus.config -with-conda /PATH/TO/miniconda3/env/vAMPirus --Analyze --pOTU --clusterAAIDlist .85,.90,.96 --ModelTnt --ModelTaa --nonparametric --boots 500`

### Custom IQ-TREE command

The default IQ-TREE command looks like this:

        `iqtree -s {alignment} --prefix {samplename} -m {model} --redo -t {startingtree} -nt auto -b {#bootstraps}`

If you were to input your custom IQ-TREE command, you would add on to:

        `iqtree -s {alignment} --prefix {samplename} --redo -t {startingtree} -T auto {YOUR_CUSTOM_COMMAND_HERE}`

To add your custom part of the IQ-TREE command, you would edit lines 82 and 83 for nucleotide and amino acid tree command, respectively:

        // Phylogeny analysis parameters
            // Customs options for IQ-TREE (Example: "-option1 A -option2 B -option3 C -option4 D")
                iqCustomnt="-your Custom -IQ tree -command here"
                iqCustomaa="-your Custom -IQ tree -command here"


## Taxonomy assignment

vAMPirus uses Diamond blastx/blastp and the provided protein database to assign taxonomy to ASVs/OTUs/AminoTypes. There are summary files generated, one as a phyloseq object and the other as a .tsv with information in a different
 arrangement. Results are also visualized as a donut graph in the final reports. You can adjust the following parameters:

    // Taxonomy assignment parameters
        // Specify name of database to use for analysis
            dbname="DATABASENAME.FASTA"
        // Path to Directory where database is being stored
            dbdir="/PATH/TO/DATABASE/DIRECTORY"
        // Toggle use of RefSeq header format; default is Reverence Viral DataBase (RVDB)
            refseq="F"

This was mentioned in an earlier section of the docs, but before running vAMPirus without the "--skipTaxonomy" option set, you should add the name and path of the database you would like to use at lines 74 and 76 in the config file,
 respectively. The database, currently, needs to be protein sequences and be in fasta format. The database can be a custom database but for proper reporting of the results, the headers need to follow either RVDB or RefSeq
 header formats:


    1. RVDB format (default) -> ">acc|GENBANK|AYD68780.1|GENBANK|MH171300|structural polyprotein [Marine RNA virus BC-4]"

    2. NCBI NR/RefSeq format -> ">KJX92028.1 hypothetical protein TI39_contig5958g00003 [Zymoseptoria brevis]"

NOTE: By default, vAMPirus assumes the headers are in RVDB format, to trigger the use of NCBI RefSeq format, edit the "F" to "T" at line 78 in the config file or add "--refseq T" to the launch command.

## EMBOSS Analyses

For pOTU protein files and AminoTypes, vAMPirus will run different protein property analyses using scripts within EMBOSS. To skip this process, just add "--skipEMBOSS" to the launch command.

# vAMPirus output

Output files produced throughout the vAMPirus pipeline is organized in directories within the specified results/output directory (e.g. )
