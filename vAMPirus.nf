#!/usr/bin/env nextflow

/*
========================================================================================
                                vAMPirus
========================================================================================
                       Virus Amplicon Sequencing Analysis Pipeline
                       Author: Alex J. Veglia and Ramón Rivera-Vicéns
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info """
    ==============================================================================================================================================================================================
                                                                         Quick help, use --fullHelp for usage examples for vAMPirus v${workflow.manifest.version}
    ==============================================================================================================================================================================================

            Steps:
                1- Before launching the vAMPirus.nf, be sure to run the vampirus_startup.sh script to install dependencies and/or databases

                2- Test the vAMPirus installation with the provided test dataset (if you have ran the start up script, you can see STARTUP_HELP.txt for test commands and other examples)

                3. Edit parameters in vampirus.config file

                4. Launch the DataCheck pipeline to get summary information about your dataset

                5. Change any parameters in vampirus.config file that might aid your analysis (e.g. clustering ID, maximum merged read length)

                6. Launch the Analyze pipeline to perform a comprehensive analysis with your dataset

                7. Explore results directories and produced final reports


            Usage:

                    nextflow run vAMPirus.nf -c vampirus.config -profile [conda|singularity] --[Analyze|DataCheck] [--ncASV] [--pcASV] [--asvMED] [--aminoMED] [--stats]


            --Help options--

                    --help                          Print help information

                    --fullHelp                      Print even more help information


            --Mandatory arguments (choose one)--

                    --Analyze                       Run absolutely everything

                    --DataCheck                     Assess how data performs with during processing and clustering


            --ASV clustering arguments--

                    --ncASV                          Set this option to have vAMPirus cluster nucleotide amplicon sequence variants (ASVs) into nucleotide-based operational taxonomic units (ncASVs) - See options below to define a single percent similarity or a list

                    --pcASV                          Set this option to have vAMPirus cluster nucleotide and translated ASVs into protein-based operational taxonomic units (pcASVs) - See options below to define a single percent similarity or a list


            --Phylogeny-based clustering--

                    --asvTClust                        Set this option to perform phylogeny-based clustering of ASV sequences, see manual for more information.

                    --aminoTClust                      Set this option to perform phylogeny-based clustering of AminoType sequences, see manual for more information.


            --Minimum Entropy Decomposition arguments--

                    --asvMED                        Set this option to perform Minimum Entropy Decomposition on ASV sequences, see manual for more information. You will need to set a value for --asvC to perform this analysis

                    --aminoMED                     Set this option to perform Minimum Entropy Decomposition on AminoType sequences, see manual for more information. You will need to set a value for --aminoC to perform this analysis

            --Skip options--

                    --skipReadProcessing            Set this option to skip all read processing steps in the pipeline

                    --skipFastQC                    Set this option to skiip FastQC steps in the pipeline

                    --skipAdapterRemoval            Set this option to skip adapter removal in the pipeline

                    --skipPrimerRemoval             Set this option to skip primer removal process

                    --skipMerging                   Set this option to skip read merging

                    --skipAminoTyping               Set this option to skip AminoTyping processes

                    --skipTaxonomy                  Set this option to skip taxonomy assignment processes

                    --skipPhylogeny                 Set this option to skip phylogeny processes

                    --skipEMBOSS                    Set this option to skip EMBOSS processes

                    --skipReport                    Set this option to skip html report generation

      **NOTE** Most opitons below can be set using the configuration file (vampirus.config) to avoid a lengthy launch command.

           --Project/analysis information--

                    --projtag                       Set project name to be used as a prefix for output files

                    --metadata                      Set path to metadata spreadsheet file to be used for report generation (must be defined if generating report)

                    --reads                         Path to directory containing read libraries, must have *R{1,2}* in the library names

                    --workingdir                    Path to working directory where Nextflow will put all Nextflow and vAMPirus generated output files

                    --outdir                        Name of results directory containing all output from the chosen pipeline (will be made within the working directory)


            --Merged read length filtering--

                    --minLen                        Minimum merged read length - reads below the specified maximum read length will be used for counts only

                    --maxLen                        Maximum merged read length - reads with length equal to the specified max read length will be used to identifying unique sequences and  subsequent Amplicon Sequence Variant (ASV) analysis

                    --maxEE                         Use this option to set the maximum expected error rate for vsearch merging. Default is 1.

                    --diffs                         Maximum number of non-matching nucleotides allowed in overlap region.

                    --maxn                          Maximum number of "N"'s in a sequence - if above the specified value, sequence will be discarded.

                    --minoverlap                    Minimum length of overlap for sequence merging to occur for a pair.


            --Primer removal--

                General primer removal parameters

                    --primerLength                  Use this option to set the max primer length to restrict bbduk.sh primer trimming to the first x number of bases

                    --maxkmer                       Maximum kmer length for bbduk.sh to use for primer detection and removal (must be shorter than your primer length; default = 13)

                    --minkmer                       Minimum kmer length for primer removal (default = 3)

                    --minilen                       Minimum non-merged read length after adapter and primer removal (default = 100)

                Single primer set removal-

                    --gtrim                      Set this option to perform global trimming to reads to remove primer sequences. Example usage "--gtrim #basesfromforward,#basesfromreverse"

                    --fwd                           Forward primer sequence for reads to be detected and removed from reads (must specify reverse sequence if providing forward)

                    --rev                           Reverse primer sequence for reads to be detected and removed from reads (must specify forward sequence if providing reverse)

                Multiple primer set removal-

                    --multi                         Use this option to signal multiple primer sequence removal within the specified pipeline

                    --primers                       Use this option to set the path to a fasta file with all of the primer sequences to be detected and removed from reads


            --Amplicon Sequence Variant (ASV) genration and clustering--

                    --alpha                         Alpha value for denoising - the higher the alpha the higher the chance of false positives in ASV generation (1 or 2)

                    --minSize                       Minimum size or representation in the dataset for sequence to be considered in ASV generation

                    --clusterNuclID                 With --ncASV set, use this option to set a single percent similarity to cluster nucleotide ASV sequences into ncASVs by [ Example: --clusterNuclID .97 ]

                    --clusterNuclIDlist             With --ncASV set, use this option to perform nucleotide clustering with a comma separated list of percent similarities [ Example: --clusterNuclIDlist .95,.96,.97,.98 ]

                    --clusterAAID                   With --pcASV set, use this option to set a single percent similarity for protein-based ASV clustering to generation pcASVs  [ Example: --clusterAAID .97 ]

                    --clusterAAIDlist               With --pcASV set, use this option to perform protein-based ASV clustering to generate pcASVs with a comma separated list of percent similarities [ Example: --clusterAAIDlist .95,.96,.97,.98 ]

                    --minAA                         With --pcASV set, use this option to set the expected or minimum amino acid sequence length of open reading frames within your amplicon sequences

           --Minimum Entropy Decomposition--

                    --asvC                          Number of high entropy positions to use for ASV MED analysis and generate "Groups"

                    --aminoC                        Number of high entropy positions to use for AminoType MED analysis and generate "Groups"

           --Counts table generation--

                    --asvcountID                    Similarity ID to use for ASV counts

                    --ProtCountID                   Minimum amino acid sequence similarity for hit to count

                    --ProtCountsLength              Minimum alignment length for hit to count

                    --ProtCountsBit                 Minimum bitscore for hit to be counted


            --Taxonomy inference parameters--

                    --dbname                       Specify name of database to use for analysis

                    --dbdir                        Path to Directory where database is being stored

                    --headers                      Set taxonomy database header format -> headers= "NCBI" to toggle use of NCBI header format; set to "RVDB" to signal the use of Reverence Viral DataBase (RVDB) headers

                    --dbanno                       Path to directory hmm annotation .txt file - see manual for information on this. Leave as is if not planning on using.

                    --lca                          Set --lca T if you would like to add taxonomic classification to taxonomy results - example: "ASV1, Viruses::Duplodnaviria::Heunggongvirae::Peploviricota::Herviviricetes::Herpesvirales::Herpesviridae::Gammaherpesvirinae::Macavirus"

                    --bitscore                     Set minimum bitscore to allow for best hit in taxonomy assignment

                    --minID                        Set minimum percent amino acid similarity for best hit to be counted in taxonomy assignment

                    --minaln                       Set minimum amino acid alignment length for best hit to be counted in taxonomy assignment


            --Phylogeny analysis parameters--

              Setting customs options for IQ-TREE (Example: "-option1 A -option2 B -option3 C -option4 D") - might be easier to set in the vampirus.config file at lines 108/109

                    --iqCustomnt                   Use option to set custom options to use in all IQTREE analyses with nuceoltide sequences

                    --iqCustomaa                   Use option to set custom options to use in all IQTREE analyses with amino acid sequences

              These options below you can set at the command, for example, to set to use model from ModelTest-NG with parametric bootstrapping --ModelTnt --ModelTaa --parametric

                    --ModelTnt=false               Signal for IQ-TREE to use model determined by ModelTest-NG for all IQTREE analyses with nuceoltide sequences (Default is IQ-TREE will do automatic model testing with ModelFinder Plus)

                    --ModelTaa=false               Signal for IQ-TREE to use model determined by ModelTest-NG for all IQTREE analyses with amino acid sequences

                    --parametric                   Set to use parametric bootstrapping in IQTREE analyses

                    --nonparametric                Set to use parametric bootstrapping in IQTREE analyses

                    --boots                        Number of bootstraps (recommended 1000 for parametric and 100 for non-parametric)


              --Statistics options--

                    --stats                        Set "--stats run" to signal statstical tests to be performed and included in the final report

                    --minimumCounts                Minimum number of hit counts for a sample to have to be included in the downstream statistical analyses and report generation

                    --trymax                       Maximum number of iterations performed by metaMDS

        |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    """.stripIndent()
}
def fullHelpMessage() {
    log.info """
    ==============================================================================================================================================================================================
                                                        THIS IS A LONGER HELP WITH USAGE EXAMPLES vAMPirus v${workflow.manifest.version}
    ==============================================================================================================================================================================================

              Steps:
                  1- Before launching the vAMPirus.nf, be sure to run the vampirus_startup.sh script to install dependencies and/or databases

                  2- Test the vAMPirus installation with the provided test dataset (if you have ran the start up script, you can see STARTUP_HELP.txt for test commands and other examples)

                  3. Edit parameters in vampirus.config file

                  4. Launch the DataCheck pipeline to get summary information about your dataset

                  5. Change any parameters in vampirus.config file that might aid your analysis (e.g. clustering ID, maximum merged read length)

                  6. Launch the Analyze pipeline to perform a comprehensive analysis with your dataset

                  7. Explore results directories and produced final reports


              Usage:

                      nextflow run vAMPirus.nf -c vampirus.config -profile [conda|singularity] --[Analyze|DataCheck] [--ncASV] [--pcASV]


                      --Help options--

                              --help                          Print help information

                              --fullHelp                      Print even more help information


                      --Mandatory arguments (choose one)--

                              --Analyze                       Run absolutely everything

                              --DataCheck                     Assess how data performs with during processing and clustering


                      --ASV clustering arguments--

                              --ncASV                          Set this option to have vAMPirus cluster nucleotide amplicon sequence variants (ASVs) into nucleotide-based operational taxonomic units (ncASVs) - See options below to define a single percent similarity or a list

                              --pcASV                          Set this option to have vAMPirus cluster nucleotide and translated ASVs into protein-based operational taxonomic units (pcASVs) - See options below to define a single percent similarity or a list


                      --Phylogeny-based clustering--

                              --asvTClust                        Set this option to perform phylogeny-based clustering of ASV sequences, see manual for more information.

                              --aminoTClust                      Set this option to perform phylogeny-based clustering of AminoType sequences, see manual for more information.


                      --Minimum Entropy Decomposition arguments--

                              --asvMED                        Set this option to perform Minimum Entropy Decomposition on ASV sequences, see manual for more information. You will need to set a value for --asvC to perform this analysis

                              --aminoMED                     Set this option to perform Minimum Entropy Decomposition on AminoType sequences, see manual for more information. You will need to set a value for --aminoC to perform this analysis


                      --Skip options--

                              --skipReadProcessing            Set this option to skip all read processing steps in the pipeline

                              --skipFastQC                    Set this option to skiip FastQC steps in the pipeline

                              --skipAdapterRemoval            Set this option to skip adapter removal in the pipeline

                              --skipPrimerRemoval             Set this option to skip primer removal process

                              --skipMerging                   Set this option to skip read merging

                              --skipAminoTyping               Set this option to skip AminoTyping processes

                              --skipTaxonomy                  Set this option to skip taxonomy assignment processes

                              --skipPhylogeny                 Set this option to skip phylogeny processes

                              --skipEMBOSS                    Set this option to skip EMBOSS processes

                              --skipReport                    Set this option to skip html report generation


                **NOTE** Most opitons below can be set using the configuration file (vampirus.config) to avoid a lengthy launch command.

                     --Project/analysis information--

                              --projtag                       Set project name to be used as a prefix for output files

                              --metadata                      Set path to metadata spreadsheet file to be used for report generation (must be defined if generating report)

                              --reads                         Path to directory containing read libraries, must have *R{1,2}* in the library names

                              --workingdir                    Path to working directory where Nextflow will put all Nextflow and vAMPirus generated output files

                              --outdir                        Name of results directory containing all output from the chosen pipeline (will be made within the working directory)


                      --Merged read length filtering--

                              --minLen                        Minimum merged read length - reads below the specified maximum read length will be used for counts only

                              --maxLen                        Maximum merged read length - reads with length equal to the specified max read length will be used to identifying unique sequences and  subsequent Amplicon Sequence Variant (ASV) analysis

                              --maxEE                         Use this option to set the maximum expected error rate for vsearch merging. Default is 1.

                              --diffs                         Maximum number of non-matching nucleotides allowed in overlap region.

                              --maxn                          Maximum number of "N"'s in a sequence - if above the specified value, sequence will be discarded.

                              --minoverlap                    Minimum length of overlap for sequence merging to occur for a pair.


                      --Primer removal--

                          General primer removal parameters

                              --primerLength                  Use this option to set the max primer length to restrict bbduk.sh primer trimming to the first x number of bases

                              --maxkmer                       Maximum kmer length for bbduk.sh to use for primer detection and removal (must be shorter than your primer length; default = 13)

                              --minkmer                       Minimum kmer length for primer removal (default = 3)

                              --minilen                       Minimum non-merged read length after adapter and primer removal (default = 100)

                          Single primer set removal-

                              --gtrim                      Set this option to perform global trimming to reads to remove primer sequences. Example usage "--gtrim #basesfromforward,#basesfromreverse"

                              --fwd                           Forward primer sequence for reads to be detected and removed from reads (must specify reverse sequence if providing forward)

                              --rev                           Reverse primer sequence for reads to be detected and removed from reads (must specify forward sequence if providing reverse)

                          Multiple primer set removal-

                              --multi                         Use this option to signal multiple primer sequence removal within the specified pipeline

                              --primers                       Use this option to set the path to a fasta file with all of the primer sequences to be detected and removed from reads


                      --Amplicon Sequence Variant (ASV) genration and clustering--

                              --alpha                         Alpha value for denoising - the higher the alpha the higher the chance of false positives in ASV generation (1 or 2)

                              --minSize                       Minimum size or representation in the dataset for sequence to be considered in ASV generation

                              --clusterNuclID                 With --ncASV set, use this option to set a single percent similarity to cluster nucleotide ASV sequences into ncASVs by [ Example: --clusterNuclID .97 ]

                              --clusterNuclIDlist             With --ncASV set, use this option to perform nucleotide clustering with a comma separated list of percent similarities [ Example: --clusterNuclIDlist .95,.96,.97,.98 ]

                              --clusterAAID                   With --pcASV set, use this option to set a single percent similarity for protein-based ASV clustering to generation pcASVs  [ Example: --clusterAAID .97 ]

                              --clusterAAIDlist               With --pcASV set, use this option to perform protein-based ASV clustering to generate pcASVs with a comma separated list of percent similarities [ Example: --clusterAAIDlist .95,.96,.97,.98 ]

                              --minAA                         With --pcASV set, use this option to set the expected or minimum amino acid sequence length of open reading frames within your amplicon sequences

                     --Minimum Entropy Decomposition--

                              --asvC                          Number of high entropy positions to use for ASV MED analysis and generate "Groups"

                              --aminoC                        Number of high entropy positions to use for AminoType MED analysis and generate "Groups"

                     --Counts table generation--

                              --asvcountID                    Similarity ID to use for ASV counts

                              --ProtCountID                   Minimum amino acid sequence similarity for hit to count

                              --ProtCountsLength              Minimum alignment length for hit to count

                              --ProtCountsBit                 Minimum bitscore for hit to be counted


                      --Taxonomy inference parameters--

                              --dbname                       Specify name of database to use for analysis

                              --dbdir                        Path to Directory where database is being stored

                              --headers                      Set taxonomy database header format -> headers= "NCBI" to toggle use of NCBI header format; set to "RVDB" to signal the use of Reverence Viral DataBase (RVDB) headers

                              --dbanno                       Path to directory hmm annotation .txt file - see manual for information on this. Leave as is if not planning on using.

                              --lca                          Set --lca T if you would like to add taxonomic classification to taxonomy results - example: "ASV1, Viruses::Duplodnaviria::Heunggongvirae::Peploviricota::Herviviricetes::Herpesvirales::Herpesviridae::Gammaherpesvirinae::Macavirus"

                              --bitscore                     Set minimum bitscore to allow for best hit in taxonomy assignment

                              --minID                        Set minimum percent amino acid similarity for best hit to be counted in taxonomy assignment

                              --minaln                       Set minimum amino acid alignment length for best hit to be counted in taxonomy assignment


                      --Phylogeny analysis parameters--

                        Setting customs options for IQ-TREE (Example: "-option1 A -option2 B -option3 C -option4 D") - might be easier to set in the vampirus.config file at lines 108/109

                              --iqCustomnt                   Use option to set custom options to use in all IQTREE analyses with nuceoltide sequences

                              --iqCustomaa                   Use option to set custom options to use in all IQTREE analyses with amino acid sequences

                        These options below you can set at the command, for example, to set to use model from ModelTest-NG with parametric bootstrapping --ModelTnt --ModelTaa --parametric

                              --ModelTnt=false               Signal for IQ-TREE to use model determined by ModelTest-NG for all IQTREE analyses with nuceoltide sequences (Default is IQ-TREE will do automatic model testing with ModelFinder Plus)

                              --ModelTaa=false               Signal for IQ-TREE to use model determined by ModelTest-NG for all IQTREE analyses with amino acid sequences

                              --parametric                   Set to use parametric bootstrapping in IQTREE analyses

                              --nonparametric                Set to use parametric bootstrapping in IQTREE analyses

                              --boots                        Number of bootstraps (recommended 1000 for parametric and 100 for non-parametric)


                        --Statistics options--

                              --stats                        Set "--stats run" to signal statstical tests to be performed and included in the final report

                              --minimumCounts                Minimum number of hit counts for a sample to have to be included in the downstream statistical analyses and report generation

                              --trymax                       Maximum number of iterations performed by metaMDS


        |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        #################################################################################################

                            Various example launch commands for deploying vAMPirus

        #################################################################################################


        DataCheck pipeline =>

        Example 1. Launching the vAMPirus DataCheck pipeline with MED analyses using conda

            nextflow run vAMPirus.nf -c vampirus.config -profile conda --DataCheck --asvMED --aminoMED

        Example 2. Launching the vAMPirus DataCheck pipeline using Singularity and multiple primer removal with the path to the fasta file with the primer sequences set in the launch command

            nextflow run vAMPirus.nf -c vampirus.config -profile singularity --DataCheck --multi --primers /PATH/TO/PRIMERs.fa

        Example 3. Launching the vAMPirus DataCheck pipeline with primer removal by global trimming of 20 bp from forward reads and 26 bp from reverse reads

            nextflow run vAMPirus.nf -c vampirus.config -profile conda --DataCheck --gtrim 20,26


        Analyze pipeline =>

        Example 4. Launching the vAMPirus Analyze pipeline with singularity with ASV and AminoType generation with all accesory analyses (taxonomy assignment, EMBOSS, IQTREE, statistics)

            nextflow run vAMPirus.nf -c vampirus.config -profile singularity --Analyze --stats

        Example 5. Launching the vAMPirus Analyze pipeline with conda to perform multiple primer removal and protein-based clustering of ASVs, but skip most of the extra analyses

            nextflow run vAMPirus.nf -c vampirus.config -profile conda --Analyze --pcASV --skipPhylogeny --skipEMBOSS --skipTaxonomy --skipReport

        Example 6. Launching vAMPirus Analyze pipeline with conda to produce only ASV-related results

            nextflow run vAMPirus.nf -c vampirus.config -profile conda --Analyze --skipAminoTyping --stats

        Example 7. Launching vAMPirus Analyze pipeline with conda to perform ASV analyses with Minimum Entropy Decomposition to form "Groups"

            nextflow run vAMPirus.nf -c vampirus.config -profile conda --Analyze --skipAminoTyping --stats --asvMED --asvC 24


        Resuming analyses =>

        If an analysis is interupted, you can use Nextflows "-resume" option that will start from the last cached "check point".

        For example if the analysis launched with the command from Example 6 above was interupted, all you would need to do is add the "-resume" to the end of the command like so:

            nextflow run vAMPirus.nf -c vampirus.config -profile conda --Analyze --skipAminoTyping --stats run -resume


    """.stripIndent()
}
// Show help message
if (params.help) {
    helpMessage()
    exit 0
} else if (params.fullHelp) {
    fullHelpMessage()
    exit 0
}

log.info """\
================================================================================================================================================
                              vAMPirus v${workflow.manifest.version} - Virus Amplicon Sequencing Analysis Pipeline
================================================================================================================================================

               -------------------------------------------------Project details---------------------------------------------
                                                     Project name:          ${params.projtag}
                                                Working directory:          ${params.workingdir}
                                                Results directory:          ${params.outdir}
                                                    Metadata file:          ${params.metadata}

               ---------------------------------------------------Run details------------------------------------------------
                                                      Single input:         ${params.single}
                                        Minimum merged read length:         ${params.maxLen}
                                                     ASV filtering:         ${params.filter}
                                                Database directory:         ${params.dbdir}
                                                     Database name:         ${params.dbname}
                                                     Database type:         ${params.dbtype}
                                                             ncASV:         ${params.ncASV}
                                                             pcASV:         ${params.pcASV}
                                                           ASV MED:         ${params.asvMED}
                                                     AminoType MED:         ${params.aminoMED}
                                    Phylogeny-based ASV clustering:         ${params.asvTClust}
                              Phylogeny-based AminoType clustering:         ${params.aminoTClust}
                                                       Skip FastQC:         ${params.skipFastQC}
                                              Skip read processing:         ${params.skipReadProcessing}
                                              Skip adapter removal:         ${params.skipAdapterRemoval}
                                               Skip primer removal:         ${params.skipPrimerRemoval}
                                                 Skip read merging:         ${params.skipMerging}
                                                  Skip AminoTyping:         ${params.skipAminoTyping}
                                                     Skip Taxonomy:         ${params.skipTaxonomy}
                                                    Skip phylogeny:         ${params.skipPhylogeny}
                                                       Skip EMBOSS:         ${params.skipEMBOSS}
                                                      Skip Reports:         ${params.skipReport}
        """.stripIndent()

if (params.readsTest) {
    println("\n\tRunning vAMPirus with TEST dataset\n")
    Channel
        .fromFilePairs(params.readsTest)
        .ifEmpty{ exit 1, "params.readTest was empty - no input files supplied" }
        .into{ reads_ch; reads_qc_ch; reads_processing }
} else if (params.single) {
    println("\n\tLocating single read files...\n")
    Channel
        .fromFilePairs("${params.reads}", size: -1, checkIfExists: true)
        .ifEmpty{ exit 1, "params.reads was empty - no input files supplied" }
        .into{ reads_ch; reads_qc_ch; reads_processing }
} else {
    println("\n\tLocating paired read files...\n")
    Channel
        .fromFilePairs("${params.reads}", checkIfExists: true)
        .ifEmpty{ exit 1, "params.reads was empty - no input files supplied" }
        .into{ reads_ch; reads_qc_ch; reads_processing }
}
if (params.clusterNuclIDlist == "") {
    a=params.clusterNuclID
    slist=[a]
    nnuc=slist.size()
} else {
    msize=params.clusterNuclIDlist
    slist=msize.split(',').collect{it as int}
    nnuc=slist.size()
}
if (params.clusterAAIDlist == "") {
    b=params.clusterAAID
    slist2=[b]
    naa=slist2.size()
} else {
    msize2=params.clusterAAIDlist
    slist2=msize2.split(',').collect{it as int}
    naa=slist2.size()
}

if (params.DataCheck || params.Analyze) {

    println("\n\tRunning vAMPirus - This might take a while depending on the mode and dataset size, check out Nextflow tower (tower.nf) to remotely monitor the run.\n")

    if (!params.skipTaxonomy && params.Analyze) {

        process Database_Check {

            conda (params.condaActivate ? "-c conda-forge bioconda::diamond=2.0.15=hb97b32f_1" : null)

            container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/diamond:2.0.15--hb97b32f_1" : "quay.io/biocontainers/diamond:2.0.15--hb97b32f_1")

            script:
                """
                cd ${params.workingdir}
                echo -e "-- Checking if specified database directory exists --\\n"
                if [ ! -d ${params.dbdir} ];then
                    echo "-- Directory does not exist where you say it does..."
                    echo "Maybe it was a mistake, looking for the Databases directory in vAMPdir"
                    if [ ! -d ${params.vampdir}/Databases ];then
                        echo "-- Databases directory is not present, check the configuration file and help documentation to set path to databases"
                        exit 1
                    else
                        echo "Ok, you have the database directory present in ${params.vampdir}"
                        echo "Checking now for any databases that matched the database name specified.."
                        if [ ! -e ${params.vampdir}/Databases/${params.dbname} ];then
                            echo "Nope, no database by that name here. Check the configuration file and help documentation to set path to the database."
                            exit 1
                        else
                            echo "Ok, found the database specified here. Lets move it to the directory you wanted."
                            mkdir ${params.dbdir}
                            cp ${params.vampdir}/Databases/${params.dbname}* ${params.dbdir}/
                            if [ ! -e ${params.dbdir}/${params.dbname}.dmnd ];then
                                echo "It needs to be built upp, doing it now"
                                if [[ ${params.ncbitax} == "true" && ${params.dbtype} == "NCBI" ]]
                                then    diamond makedb --in ${params.dbdir}/${params.dbname} -d ${params.dbdir}/${params.dbname} --taxonmap ${params.dbdir}/NCBItaxonomy/prot.accession2taxid.FULL --taxonnodes ${params.dbdir}/NCBItaxonomy/nodes.dmp --taxonnames ${params.dbdir}/NCBItaxonomy/names.dmp
                                else    diamond makedb --in ${params.dbdir}/${params.dbname} -d ${params.dbdir}/${params.dbname}
                                fi
                                export virdb=${params.dbdir}/${params.dbname}
                            else
                                echo "Database looks to be present and built."
                            fi
                        fi
                    fi
                    cd  ${params.workingdir}
                elif [ -d ${params.dbdir} ];then
                    echo -e "-- Directory exists. Checking if specified database is present now.. --\\n"
                    if [ ! -e ${params.dbdir}/${params.dbname} ];then
                        echo "Specified database not present, edit the configuraton file with the database name plz."
                        exit 1
                    else
                        echo "Database present, checking if built now.."
                        if [ ! -e ${params.dbdir}/${params.dbname}.dmnd  ];then
                            echo "Database not built, building now..."
                            if [ ! -e ${params.dbdir}/${params.dbname}.dmnd ];then
                                echo "It needs to be built upp, doing it now"
                                if [[ ${params.ncbitax} == "true" && ${params.dbtype} == "NCBI" ]]
                                then    diamond makedb --in ${params.dbdir}/${params.dbname} -d ${params.dbdir}/${params.dbname} --taxonmap ${params.dbdir}/NCBItaxonomy/prot.accession2taxid.FULL --taxonnodes ${params.dbdir}/NCBItaxonomy/nodes.dmp --taxonnames ${params.dbdir}/NCBItaxonomy/names.dmp
                                else    diamond makedb --in ${params.dbdir}/${params.dbname} -d ${params.dbdir}/${params.dbname}
                                fi
                                export virdb=${params.dbdir}/${params.dbname}
                            else
                                echo "Database looks to be present and built."
                            fi
                        else
                            echo "-- Database is ready to go!"
                            export virdb=${params.dbdir}/${params.dbname}
                        fi
                    fi
                    cd  ${params.workingdir}
                fi
                """
        }
    }

    if (!params.skipReadProcessing || !params.skipMerging ) {

        if (!params.skipFastQC) {

                if (!params.single) {

                    process QualityCheck_1 {

                        label 'low_cpus'

                        tag "${sample_id}"

                        publishDir "${params.workingdir}/${params.outdir}/ReadProcessing/FastQC/PreClean", mode: "copy", overwrite: true

                        conda (params.condaActivate ? "-c conda-forge bioconda::fastqc=0.11.9=0" : null)

                        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0" : "quay.io/biocontainers/fastqc:0.11.9--0")

                        input:
                            tuple sample_id, file(reads) from reads_qc_ch

                        output:
                            tuple sample_id, file("*_fastqc.{zip,html}") into fastqc_results

                        script:
                            """
                            fastqc --quiet --threads ${task.cpus} ${reads}
                            """
                    }
                } else if (params.single) {

                        process QualityCheck_se1 {

                        label 'low_cpus'

                        tag "${sample_id}"

                        publishDir "${params.workingdir}/${params.outdir}/ReadProcessing/FastQC/PreClean", mode: "copy", overwrite: true

                        conda (params.condaActivate ? "-c conda-forge bioconda::fastqc=0.11.9=0" : null)

                        conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "bioconda::fastp=0.20.1=h8b12597_0" : null)

                        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0" : "quay.io/biocontainers/fastqc:0.11.9--0")

                        input:
                            tuple sample_id, file(reads) from reads_qc_ch

                        output:
                            tuple sample_id, file("*_fastqc.{zip,html}") into fastqc_results

                        script:
                            """
                            fastqc --quiet --threads ${task.cpus} ${reads}
                            """
                }
            }
        }

        if (!params.skipAdapterRemoval ) {

            if (!params.single) {

                    process Adapter_Removal {

                        label 'norm_cpus'

                        tag "${sample_id}"

                        publishDir "${params.workingdir}/${params.outdir}/ReadProcessing/AdapterRemoval", mode: "copy", overwrite: true, pattern: "*.filter.fq"
                        publishDir "${params.workingdir}/${params.outdir}/ReadProcessing/AdapterRemoval/fastpOut", mode: "copy", overwrite: true, pattern: "*.fastp.{json,html}"

                        conda (params.condaActivate ? "bioconda::fastp=0.20.1=h8b12597_0" : null)

                        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/fastp:0.23.2--hb7a2d85_2" : "quay.io/biocontainers/fastp:0.23.2--hb7a2d85_2")

                        input:
                            tuple sample_id, file(reads) from reads_ch

                        output:
                            tuple sample_id, file("*.fastp.{json,html}") into fastp_results
                            tuple sample_id, file("*.fastp.json") into fastp_json
                            tuple sample_id, file("*.filter.fq") into reads_fastp_ch


                        script:
                            """
                            echo ${sample_id}

                            fastp -i ${reads[0]} -I ${reads[1]} -o left-${sample_id}.filter.fq -O right-${sample_id}.filter.fq --detect_adapter_for_pe \
                            --average_qual ${params.avQ} -n ${params.mN} -c --overrepresentation_analysis --html ${sample_id}.fastp.html --json ${sample_id}.fastp.json --thread ${task.cpus} \
                            --report_title ${sample_id}

                            """
                        }

            } else if (params.single) {

                process Adapter_Removal_se {

                    label 'norm_cpus'

                    tag "${sample_id}"

                    publishDir "${params.workingdir}/${params.outdir}/ReadProcessing/AdapterRemoval", mode: "copy", overwrite: true, pattern: "*.filter.fq"
                    publishDir "${params.workingdir}/${params.outdir}/ReadProcessing/AdapterRemoval/fastpOut", mode: "copy", overwrite: true, pattern: "*.fastp.{json,html}"

                    conda (params.condaActivate ? "bioconda::fastp=0.20.1=h8b12597_0" : null)

                    container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/fastp:0.23.2--hb7a2d85_2" : "quay.io/biocontainers/fastp:0.23.2--hb7a2d85_2")

                    input:
                        tuple sample_id, file(reads) from reads_ch

                    output:
                        tuple sample_id, file("*.fastp.{json,html}") into fastp_results
                        tuple sample_id, file("*.fastp.json") into fastp_json
                        tuple sample_id, file("*.filter.fq") into reads_fastp_ch

                    script:
                        """
                        echo ${sample_id}

                        fastp -i ${reads} -o ${sample_id}.filter.fq --average_qual ${params.avQ} -n ${params.mN} --overrepresentation_analysis --html ${sample_id}.fastp.html --json ${sample_id}.fastp.json --thread ${task.cpus} --report_title ${sample_id}
                        """
                }
            }

            process Read_stats_processing {

                label 'norm_cpus'

                tag "${sample_id}"

                publishDir "${params.workingdir}/${params.outdir}/ReadProcessing/AdapterRemoval", mode: "copy", overwrite: true, pattern: "*.filter.fq"
                publishDir "${params.workingdir}/${params.outdir}/ReadProcessing/AdapterRemoval/fastpOut", mode: "copy", overwrite: true, pattern: "*.fastp.{json,html}"

                conda (params.condaActivate ? "conda-forge::jq=1.6" : null)

                container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/jq:1.6" : "quay.io/biocontainers/jq:1.6")

                input:
                    tuple sample_id, file(json) from fastp_json

                output:
                    file("*.csv") into ( fastp_csv_in1, fastp_csv_in2 )

                script:
                    """
                    echo ${sample_id}

                    bash ${params.vampdir}/bin/get_readstats.sh ${json}  ####check with ramon why tf this isnt automatic with tools being exprted correctly
                    """
                }

        } else {
                reads_ch
                    .set{ reads_fastp_ch }
                fastp_results = Channel.empty()
        }

        if (!params.skipPrimerRemoval) {

            if (!params.single) {

                process Primer_Removal {

                    label 'norm_cpus'

                    tag "${sample_id}"

                    publishDir "${params.workingdir}/${params.outdir}/ReadProcessing/PrimerRemoval", mode: "copy", overwrite: true

                    conda (params.condaActivate ? "-c bioconda -c conda-forge  bbmap=39.01" : null)

                    container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/bbamap:38.98--h5c4e2a8_0" : "quay.io/biocontainers/bbmap:38.98--h5c4e2a8_0")

                    input:
                        tuple sample_id, file(reads) from reads_fastp_ch

                    output:
                        tuple sample_id, file("*bbduk*.fastq.gz") into ( reads_bbduk_ch, readsforqc2 )

                    script:
                        // check if we need to check this outside processes
                        if ( params.fwd == "" && params.rev == "" && !params.multi && params.gtrim == "" ) {
                            """
                            bbduk.sh in1=${reads[0]} out=${sample_id}_bb_R1.fastq.gz ftl=${params.defaultFwdTrim} t=${task.cpus}
                            bbduk.sh in=${reads[1]} out=${sample_id}_bb_R2.fastq.gz ftl=${params.defaultRevTrim} t=${task.cpus}
            		            repair.sh in1=${sample_id}_bb_R1.fastq.gz in2=${sample_id}_bb_R2.fastq.gz out1=${sample_id}_bbduk_R1.fastq.gz out2=${sample_id}_bbduk_R2.fastq.gz outs=sing.fq repair
                            """
                        } else if ( params.gtrim != "" && !params.multi ) {
                            """
                            FTRIM=\$( echo ${params.gtrim} | cut -f 1 -d "," )
                            RTRIM=\$( echo ${params.gtrim} | cut -f 2 -d "," )
                            bbduk.sh in=${reads[0]} out=${sample_id}_bb_R1.fastq.gz ftl=\${FTRIM} t=${task.cpus}
                            bbduk.sh in=${reads[1]} out=${sample_id}_bb_R2.fastq.gz ftl=\${RTRIM} t=${task.cpus}
            		        repair.sh in1=${sample_id}_bb_R1.fastq.gz in2=${sample_id}_bb_R2.fastq.gz out1=${sample_id}_bbduk_R1.fastq.gz out2=${sample_id}_bbduk_R2.fastq.gz outs=sing.fq repair
                            """
                        } else if ( params.multi && params.primers ) {
                            """
                            bbduk.sh in=${reads[0]} in2=${reads[1]} out=${sample_id}_bbduk_R1.fastq.gz out2=${sample_id}_bbduk_R2.fastq.gz ref=${params.primers} copyundefined=t t=${task.cpus} restrictleft=${params.primerLength} k=${params.maxkmer} ordered=t mink=${params.minkmer} ktrim=l ecco=t rcomp=t minlength=${params.minilen} tbo tpe
                            """
                        } else {
                            """
                            bbduk.sh in=${reads[0]} in2=${reads[1]} out=${sample_id}_bbduk_R1.fastq.gz out2=${sample_id}_bbduk_R2.fastq.gz literal=${params.fwd},${params.rev} copyundefined=t t=${task.cpus} restrictleft=${params.primerLength} k=${params.maxkmer} ordered=t mink=${params.minkmer} ktrim=l ecco=t rcomp=t minlength=${params.minilen} tbo tpe
                            """
                        }
                }

            } else if (params.single) {

                process Primer_Removal_se {

                    label 'norm_cpus'

                    tag "${sample_id}"

                    publishDir "${params.workingdir}/${params.outdir}/ReadProcessing/PrimerRemoval", mode: "copy", overwrite: true

                    conda (params.condaActivate ? "-c bioconda -c conda-forge  bbmap=39.01" : null)

                    container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/bbamap:38.98--h5c4e2a8_0" : "quay.io/biocontainers/bbmap:38.98--h5c4e2a8_0")

                    input:
                        tuple sample_id, file(reads) from reads_fastp_ch

                    output:
                        tuple sample_id, file("*bbduk*.fastq.gz") into ( reads_bbduk_ch, readsforqc2 )

                    script:
                        // check if we need to check this outside processes
                        if ( params.fwd == "" && params.rev == "" && !params.multi && params.trim == "") {
                            """
                            bbduk.sh in1=${reads} out=${sample_id}_bbduk_SE.fastq.gz ftl=${params.defaultFwdTrim} t=${task.cpus}
                            """
                        } else if ( params.trim != "" && !params.multi ) {
                            """
                            bbduk.sh in=${reads} out=${sample_id}_bbduk_SE.fastq.gz ftl=${trim} t=${task.cpus}

                            """
                        } else if ( params.multi && params.primers ) {
                            """
                            bbduk.sh in=${reads} out=${sample_id}_bbduk_SE.fastq.gz ref=${params.primers} copyundefined=t t=${task.cpus} restrictleft=${params.primerLength} k=${params.maxkmer} ordered=t mink=${params.minkmer} ktrim=l ecco=t rcomp=t minlength=${params.minilen} tbo tpe
                            """
                        } else {
                            """
                            bbduk.sh in=${reads} out=${sample_id}_bbduk_SE.fastq.gz literal=${params.fwd} copyundefined=t t=${task.cpus} restrictleft=${params.primerLength} k=${params.maxkmer} ordered=t mink=${params.minkmer} ktrim=l ecco=t rcomp=t minlength=${params.minilen} tbo tpe
                            """
                        }
                }
            }

        } else {
            reads_fastp_ch
                .set{ reads_bbduk_ch }
        }

        if (!params.skipFastQC && !params.skipPrimerRemoval) {

            if (!params.single) {

                    process QualityCheck_2 {

                        label 'low_cpus'

                        tag "${sample_id}"

                        publishDir "${params.workingdir}/${params.outdir}/ReadProcessing/FastQC/PostClean", mode: "copy", overwrite: true

                        conda (params.condaActivate ? "-c conda-forge bioconda::fastqc=0.11.9=0" : null)

                        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0" : "quay.io/biocontainers/fastqc:0.11.9--0")

                        input:
                            tuple sample_id, file(reads) from readsforqc2

                        output:
                            tuple sample_id, file("*_fastqc.{zip,html}") into fastqc2_results

                        script:
                            """
                            fastqc --quiet --threads ${task.cpus} ${reads}
                            """
                    }

            } else if (params.single) {

                    process QualityCheck_se2 {

                        label 'low_cpus'

                        tag "${sample_id}"

                        publishDir "${params.workingdir}/${params.outdir}/ReadProcessing/FastQC/PostClean", mode: "copy", overwrite: true

                        conda (params.condaActivate ? "-c conda-forge bioconda::fastqc=0.11.9=0" : null)

                        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0" : "quay.io/biocontainers/fastqc:0.11.9--0")

                        input:
                            tuple sample_id, file(reads) from readsforqc2

                        output:
                            file("*_fastqc.{zip,html}") into fastqc2_results

                        script:
                            """
                            fastqc --quiet --threads ${task.cpus} ${reads}
                            """
                }
            }
        }

    } else {
        reads_ch
            .set{ reads_bbduk_ch }
    }

    if (!params.skipMerging && !params.single) {

        process Read_Merging {

            label 'norm_cpus'

            tag "${sample_id}"

            publishDir "${params.workingdir}/${params.outdir}/ReadProcessing/ReadMerging/Individual", mode: "copy", overwrite: true, pattern: "*mergedclean.fastq"
            publishDir "${params.workingdir}/${params.outdir}/ReadProcessing/ReadMerging/Individual/notmerged", mode: "copy", overwrite: true, pattern: "*notmerged*.fastq"

            conda (params.condaActivate ? "-c bioconda -c conda-forge vsearch=2.21.1=hf1761c0_1" : null)

            container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/vsearch:2.21.1--hf1761c0_1" : "quay.io/biocontainers/vsearch:2.21.1--hf1761c0_1")

            input:
                tuple sample_id, file(reads) from reads_bbduk_ch

            output:
                file("*_mergedclean.fastq") into reads_vsearch1_ch
                file("*.name") into names
                file("*notmerged*.fastq") into notmerged

            script:
                """
                vsearch --fastq_mergepairs ${reads[0]} --reverse ${reads[1]} --threads ${task.cpus} --fastqout ${sample_id}_mergedclean.fastq --fastqout_notmerged_fwd ${sample_id}_notmerged_fwd.fastq --fastqout_notmerged_rev ${sample_id}_notmerged_rev.fastq  --fastq_maxdiffs ${params.diffs} --fastq_maxns ${params.maxn} --fastq_allowmergestagger --fastq_maxee ${params.maxEE} --fastq_minovlen ${params.minoverlap} --relabel ${sample_id}.
                echo ${sample_id} > ${sample_id}.name
                """

        }

    } else if (params.single) {

        process Sequence_renaming {

            label 'norm_cpus'

            tag "${sample_id}"

            publishDir "${params.workingdir}/${params.outdir}/ReadProcessing/ReadMerging/Individual", mode: "copy", overwrite: true, pattern: "*mergedclean.fastq"
            publishDir "${params.workingdir}/${params.outdir}/ReadProcessing/ReadMerging/", mode: "copy", overwrite: true, pattern: "*.txt"

            conda (params.condaActivate ? "-c bioconda -c conda-forge vsearch=2.21.1=hf1761c0_1" : null)

            container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/vsearch:2.21.1--hf1761c0_1" : "quay.io/biocontainers/vsearch:2.21.1--hf1761c0_1")

            input:
                tuple sample_id, file(reads) from reads_bbduk_ch

            output:
                file("*_mergedclean.fastq") into reads_vsearch1_ch
                file("*.name") into names
                file("README_IF_SET_SINGLE-END.txt") into lol

            script:
                """
                vsearch --fastq_convert ${reads} --threads ${task.cpus} --fastqout ${sample_id}_mergedclean.fastq --relabel ${sample_id}.
                echo ${sample_id} > ${sample_id}.name
                echo "These are renamed sequences, not merged. Single-end data mode set." >>README_IF_SET_SINGLE-END.txt
                """

        }
    } else {
        reads_bbduk_ch
          .set{ reads_vsearch1_ch }
    }

    process Filtering_Prep1 {

        label 'low_cpus'

        publishDir "${params.workingdir}/${params.outdir}/ReadProcessing/ReadMerging/LengthFiltering", mode: "copy", overwrite: true

        conda (params.condaActivate ? null : null)

        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? null : null)

        input:
            file(reads) from reads_vsearch1_ch
                .collect()

        output:
            file("*_all_merged_preFilt_preClean.fastq") into collect_samples_ch

        script:
            """
            cat ${reads} >>${params.projtag}_all_merged_preFilt_preClean.fastq
	        """
    }

    process Filtering_Prep2 {

        label 'low_cpus'

        publishDir "${params.workingdir}/${params.outdir}/ReadProcessing/ReadMerging", mode: "copy", overwrite: true

        conda (params.condaActivate ? null : null)

        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? null : null)

        input:
            file(names) from names
                .collect()

        output:
            file("*sample_ids.list") into ( samplelist, samplistpotu )

        script:
            """
            cat ${names} >>${params.projtag}_sample_ids.list
            """
    }

    process Length_Filtering {

        label 'norm_cpus'

        publishDir "${params.workingdir}/${params.outdir}/ReadProcessing/ReadMerging/LengthFiltering", mode: "copy", overwrite: true, pattern: "*_merged_preFilt*.fasta"
        publishDir "${params.workingdir}/${params.outdir}/ReadProcessing/ReadMerging", mode: "copy", overwrite: true, pattern: "*Lengthfiltered.fastq"
        publishDir "${params.workingdir}/${params.outdir}/ReadProcessing/ReadMerging/Histograms/pre_length_filtering", mode: "copy", overwrite: true, pattern: "*preFilt_*st.txt"
        publishDir "${params.workingdir}/${params.outdir}/ReadProcessing/ReadMerging/Histograms/post_length_filtering", mode: "copy", overwrite: true, pattern: "*postFilt_*st.txt"

        conda (params.condaActivate ? "-c bioconda -c conda-forge fastp=0.23.2 bbmap=39.01" : null)

        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/mulled-v2-22ce7eae944f2086a1f83e59d4735573352eff58:f8832d34305db510dcb75c90775f8b0ff5aab759-0" : "quay.io/biocontainers/mulled-v2-22ce7eae944f2086a1f83e59d4735573352eff58:f8832d34305db510dcb75c90775f8b0ff5aab759-0")

        input:
            file(reads) from collect_samples_ch

        output:
            file("*_merged_preFilt_clean.fastq") into (  mergeforprotcounts, mergeforpcASVaacounts )
            file("*_merged_preFilt_clean.fasta") into ( nuclCounts_mergedreads_asv_ch, nuclCounts_mergedreads_ncasv_ch, pcASV_mergedreads_ch )
            file("*_merged_clean_Lengthfiltered.fastq") into reads_vsearch2_ch
            file("*preFilt_preClean_baseFrequency_hist.csv") into prefilt_basefreq
            file("*preFilt_preClean_qualityScore_hist.csv") into prefilt_qualityscore
            file("*preFilt_preClean_gcContent_hist.csv") into prefilt_gccontent
            file("*preFilt_preClean_averageQuality_hist.csv") into prefilt_averagequality
            file("*preFilt_preClean_length_hist.csv") into prefilt_length

            file("*postFilt_baseFrequency_hist.csv") into postFilt_basefreq
            file("*postFilt_qualityScore_hist.csv") into postFilt_qualityscore
            file("*postFilt_gcContent_hist.csv") into postFilt_gccontent
            file("*postFilt_averageQuaulity_hist.csv") into postFilt_averagequality
            file("*postFilt_length_hist.csv") into postFilt_length
            file("reads_per_sample_preFilt_preClean.csv") into reads_per_sample_preFilt
            file("read_per_sample_postFilt_postClean.csv") into reads_per_sample_postFilt

        script:
            """
            # from DC
            bbduk.sh in=${reads} bhist=${params.projtag}_all_merged_preFilt_preClean_baseFrequency_hist.txt qhist=${params.projtag}_all_merged_preFilt_preClean_qualityScore_hist.txt gchist=${params.projtag}_all_merged_preFilt_preClean_gcContent_hist.txt aqhist=${params.projtag}_all_merged_preFilt_preClean_averageQuality_hist.txt lhist=${params.projtag}_all_merged_preFilt_preClean_length_hist.txt gcbins=auto
            for x in *preFilt*hist.txt;do
                pre=\$(echo \$x | awk -F ".txt" '{print \$1}')
                cat \$x | tr "\t" "," > \${pre}.csv
                rm \$x
            done
            bbduk.sh in=${reads} out=${params.projtag}_qtrimmed.fastq t=${task.cpus} qtrim=rl trimq=${params.trimq}
            reformat.sh in=${params.projtag}_qtrimmed.fastq out=${params.projtag}_preFilt_preclean.fasta t=${task.cpus}
            echo "sample,reads" >> reads_per_sample_preFilt_preClean.csv
            grep ">" ${params.projtag}_preFilt_preclean.fasta | awk -F ">" '{print \$2}' | awk -F "." '{print \$1}' | sort --parallel=${task.cpus} | uniq -c | sort -brg --parallel=${task.cpus} | awk '{print \$2","\$1}' >> reads_per_sample_preFilt_preClean.csv
            rm ${params.projtag}_preFilt_preclean.fasta
            fastp -i ${reads} -o ${params.projtag}_merged_preFilt_clean.fastq -b ${params.maxLen} -l ${params.minLen} --thread ${task.cpus} -n ${params.maxn}
            reformat.sh in=${params.projtag}_merged_preFilt_clean.fastq out=${params.projtag}_merged_preFilt_clean.fasta t=${task.cpus}
            bbduk.sh in=${params.projtag}_merged_preFilt_clean.fastq out=${params.projtag}_merged_clean_Lengthfiltered.fastq minlength=${params.maxLen} maxlength=${params.maxLen} t=${task.cpus}
            bbduk.sh in=${params.projtag}_merged_clean_Lengthfiltered.fastq bhist=${params.projtag}_all_merged_postFilt_baseFrequency_hist.txt qhist=${params.projtag}_all_merged_postFilt_qualityScore_hist.txt gchist=${params.projtag}_all_merged_postFilt_gcContent_hist.txt aqhist=${params.projtag}_all_merged_postFilt_averageQuaulity_hist.txt lhist=${params.projtag}_all_merged_postFilt_length_hist.txt gcbins=auto
            for x in *postFilt*hist.txt;do
                pre=\$(echo \$x | awk -F ".txt" '{print \$1}')
                cat \$x | tr "\t" "," > \${pre}.csv
                rm \$x
            done
            reformat.sh in=${params.projtag}_merged_clean_Lengthfiltered.fastq out=${params.projtag}_merged_clean_Lengthfiltered.fasta t=${task.cpus}
            echo "sample,reads" >> read_per_sample_postFilt_postClean.csv
            grep ">" ${params.projtag}_merged_clean_Lengthfiltered.fasta | awk -F ">" '{print \$2}' | awk -F "." '{print \$1}' | sort --parallel=${task.cpus} | uniq -c | sort -brg --parallel=${task.cpus} | awk '{print \$2","\$1}' >> read_per_sample_postFilt_postClean.csv
            """
    }

    process Extracting_Uniques {

        label 'low_cpus'

        publishDir "${params.workingdir}/${params.outdir}/ReadProcessing/ReadMerging/Uniques", mode: "copy", overwrite: true

        conda (params.condaActivate ? "-c bioconda -c conda-forge vsearch=2.21.1=hf1761c0_1" : null)

        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/vsearch:2.21.1--hf1761c0_1" : "quay.io/biocontainers/vsearch:2.21.1--hf1761c0_1")

        input:
            file(reads) from reads_vsearch2_ch

        output:
            file("*unique_sequences.fasta") into reads_vsearch3_ch
            file("*unique_sequences.fastq") into reads_fastquniq

        script:
            """
            vsearch --fastx_uniques ${reads} --sizeout --relabel_keep --fastaout ${params.projtag}_unique_sequences.fasta --fastqout ${params.projtag}_unique_sequences.fastq
            """
    }

    process Identify_ASVs {

        label 'norm_cpus'

        publishDir "${params.workingdir}/${params.outdir}/ReadProcessing/ASVs/ChimeraCheck", mode: "copy", overwrite: true

        conda (params.condaActivate ? "-c bioconda -c conda-forge vsearch=2.21.1=hf1761c0_1" : null)

        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/vsearch:2.21.1--hf1761c0_1" : "quay.io/biocontainers/vsearch:2.21.1--hf1761c0_1")

        input:
            file(reads) from reads_fastquniq

        output:
            file("*notChecked.fasta") into reads_vsearch4_ch

        script:
            """
            vsearch --cluster_unoise ${reads} --unoise_alpha ${params.alpha} --relabel ASV --centroids ${params.projtag}_notChecked.fasta --minsize ${params.minSize} --sizeout
            """
    }

    process Chimera_Check {

        label 'low_cpus'

        publishDir "${params.workingdir}/${params.outdir}/ReadProcessing/ASVs", mode: "copy", overwrite: true

        conda (params.condaActivate ? "-c bioconda -c conda-forge vsearch=2.21.1=hf1761c0_1" : null)

        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/vsearch:2.21.1--hf1761c0_1" : "quay.io/biocontainers/vsearch:2.21.1--hf1761c0_1")

        input:
            file(fasta) from reads_vsearch4_ch

        output:
            file("*ASVs.fasta") into asvforfilt

        script:
            """
	        vsearch --uchime3_denovo ${fasta} --relabel ASV --nonchimeras ${params.projtag}_ASVs.fasta
            """
    }

    // UNTIL HERE DEFAULT
    ///multi envs _ diamond,seqtk !!!!!!!!!!!!!!!!!
    if (params.filter) {

        process ASV_Filtering { //CHECK rename_seq.py in container !!!!!!!!!!!

            label 'norm_cpus'

            publishDir "${params.workingdir}/${params.outdir}/ReadProcessing/ASVFiltering", mode: "copy", overwrite: true
            publishDir "${params.workingdir}/${params.outdir}/ReadProcessing/ASVs", mode: "copy", overwrite: true, pattern: "*ASV.fasta"

            conda (params.condaActivate ? "bioconda::diamond=2.0.15=hb97b32f_1 bioconda::seqtk=1.3=h7132678_4" : null)

            container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/mulled-v2-b650c9790d7ca81a2448038cd6bbb974b925ecf2:00a620a1a2f7f99fcc9bb6e73d2cf216ac34bda1-0" : "quay.io/biocontainers/mulled-v2-b650c9790d7ca81a2448038cd6bbb974b925ecf2:00a620a1a2f7f99fcc9bb6e73d2cf216ac34bda1-0")

            input:
                file(asv) from asvforfilt

            output:
                file("*ASV.fasta") into ( reads_vsearch5_ch, reads_vsearch6_ch, asv_med, asv_for_med, nucl2aa, asvsforAminotyping, asvfastaforcounts, asvaminocheck )
                file("*.csv") into ( nothing )
                file("*diamondfilter.out") into ( noth )

            script:
                """
                cp ${params.vampdir}/bin/rename_seq.py .

                #create and rename  filter database
                if [[ ${params.filtDB} != "" ]]
                then    grep ">" ${params.filtDB} | sed 's/ //g' | awk -F ">" '{print \$2}' > filt.head
                        j=1
                        for y in \$( cat filt.head );do
                            echo ">Filt"\$j"" >> filt.headers
                            j=\$(( \${j}+1 ))
                        done
                        ./rename_seq.py ${params.filtDB} filt.headers filterdatabaserenamed.fasta
                        cat filterdatabaserenamed.fasta >> combodatabase.fasta
                        paste -d',' filt.head filt.headers > filtername_map.csv
                        rm filterdatabaserenamed.fasta filt.head
                fi
                #create and rename keep database if available
                if [[ ${params.keepDB} != "" ]]
                then    grep ">" ${params.keepDB} | sed 's/ //g' | awk -F ">" '{print \$2}' > keep.head
                        d=1
                        for y in \$( cat keep.head );do
                        echo ">keep"\$d"" >> keep.headers
                        d=\$(( \${d}+1 ))
                        done
                        ./rename_seq.py ${params.keepDB} keep.headers keepdatabaserenamed.fasta
                        cat keepdatabaserenamed.fasta >> combodatabase.fasta
                        paste -d',' keep.head keep.headers > keepername_map.csv
                        rm keepdatabaserenamed.fasta keep.head
                fi
                #index database
                diamond makedb --in combodatabase.fasta --db combodatabase.fasta
                #run diamond_db
                diamond blastx -q ${asv} -d combodatabase.fasta -p ${task.cpus} --id ${params.filtminID} -l ${params.filtminaln} -e ${params.filtevalue} --${params.filtsensitivity} -o ${params.projtag}_diamondfilter.out -f 6 qseqid qlen sseqid qstart qend qseq sseq length qframe evalue bitscore pident --max-target-seqs 1 --max-hsps 1
                #get asvs
                grep ">" ${asv} | awk -F ">" '{print \$2}' > asv.list

                for x in \$(cat asv.list)
                do  #check for a hit
                    if [[ \$(grep -cw "\$x" ${params.projtag}_diamondfilter.out) -eq 1 ]]
                    then    hit=\$(grep -w "\$x" ${params.projtag}_diamondfilter.out | awk '{print \$3}')
                            #check if hit is to filter
                            if [[ ${params.filtDB} != "" ]]
                            then    if [[ \$(grep -wc "\$hit" filt.headers) -eq 1 ]]
                                    then    echo "\$x,\$hit" >> filtered_asvs_summary.csv
                                    fi
                            fi
                            if [[ ${params.keepDB} != "" ]]
                            then    #check if hit is to keep
                                    if [[ \$(grep -wc "\$hit" keep.headers) -eq 1 ]]
                                    then    echo "\$x" >> kep.list
                                    fi
                            fi
                    else    echo \$x >> nohit.list
                    fi
                done
                #Add ASVs with no hit
                if [[ ${params.keepnohit} == "true" ]] || [[ ${params.keepDB} = "" ]]
                then    cat nohit.list >> kep.list
                        cat kep.list | sort >> keep.list
                        seqtk subseq ${asv} keep.list > kept.fasta
                        u=1
                        for y in \$( cat keep.list );do
                            echo ">ASV\${u}" >> asvrename.list
                            u=\$(( \${u}+1 ))
                        done
                        ./rename_seq.py kept.fasta asvrename.list ${params.projtag}_ASV.fasta
                else
                        cat kep.list | sort > keep.list
                        seqtk subseq ${asv} keep.list > kept.fasta
                        u=1
                        for y in \$( cat keep.list );do
                            echo ">ASV\${u}" >> asvrename.list
                            u=\$(( \${u}+1 ))
                        done
                        ./rename_seq.py kept.fasta asvrename.list ${params.projtag}_ASV.fasta
                fi
                paste -d',' keep.list asvrename.list > ASV_rename_map.csv
                """
        }

    } else {
        asvforfilt
            .into{ reads_vsearch5_ch; reads_vsearch6_ch; asv_med; asv_for_med; nucl2aa; asvsforAminotyping; asvfastaforcounts; asvaminocheck }
    }

    if (params.DataCheck) {

        process NucleotideBased_ASV_clustering_DC {

            label 'norm_cpus'

            publishDir "${params.workingdir}/${params.outdir}/DataCheck/ClusteringTest/Nucleotide", mode: "copy", overwrite: true, pattern: '*{.csv}'

            conda (params.condaActivate ? "-c bioconda -c conda-forge vsearch=2.21.1=hf1761c0_1" : null)

            container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/vsearch:2.21.1--hf1761c0_1" : "quay.io/biocontainers/vsearch:2.21.1--hf1761c0_1")

            input:
                file(fasta) from reads_vsearch5_ch

            output:
                file("number_per_percentage_nucl.csv") into number_per_percent_nucl_plot

            script:
                if (params.datacheckntIDlist) {
                """
                for id in `echo ${params.datacheckntIDlist} | tr "," "\\n"`;do
                    vsearch --cluster_fast ${fasta} --centroids ${params.projtag}_ncASV\${id}.fasta --threads ${task.cpus} --relabel OTU --id \${id}
                done
                for x in *ncASV*.fasta;do
                    id=\$( echo \$x | awk -F "_ncASV" '{print \$2}' | awk -F ".fasta" '{print \$1}')
                    numb=\$( grep -c ">" \$x )
                    echo "\${id},\${numb}" >> number_per_percentage_nucl.csv
                done
                yo=\$(grep -c ">" ${fasta})
    	        echo "1.0,\${yo}" >> number_per_percentage_nucl.csv
                """
                }
        }

        process ASV_Pairwise_Similarity_Check {

                label 'norm_cpus'

                publishDir "${params.workingdir}/${params.outdir}/DataCheck/ClusteringTest/Nucleotide", mode: "copy", overwrite: true, pattern: '*{.matrix}'

                conda (params.condaActivate ? "-c bioconda -c conda-forge clustalo=1.2.4=h87f3376_5" : null)

                container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/clustalo:1.2.4--h87f3376_5" : "quay.io/biocontainers/clustalo:1.2.4--h87f3376_5")

                input:
                    file(fasta) from reads_vsearch6_ch

                output:
                    file("${params.projtag}_ASV_PairwisePercentID.matrix") into asvpdm

                script:
                    """
                    clustalo -i ${fasta} --distmat-out=${params.projtag}_PercentIDq.matrix --percent-id --full --force --threads=${task.cpus}
                    cat ${params.projtag}_PercentIDq.matrix | tr " " "," | sed 's/,,/,/g' | grep "," >${params.projtag}_ASV_PairwisePercentID.matrix
                    rm ${params.projtag}_PercentIDq.matrix
                    """
       }

       process Translation_For_ProteinBased_Clustering_DC {

           label 'norm_cpus'

           publishDir "${params.workingdir}/${params.outdir}/DataCheck/ClusteringTest/Aminoacid/translation", mode: "copy", overwrite: true

           conda (params.condaActivate ? null : null)

           container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? null : null)

           input:
                file(fasta) from nucl2aa

            output:
                file("*ASVtranslations.fasta") into clustering_aa
                file("*_translation_report") into reportaa_VR
                file("*_ASV_all.fasta") into asvfastaforaaclust

            script:
                """
                ${tools}/virtualribosomev2/dna2pep.py ${fasta} -r all -x -o none --fasta ${params.projtag}_ASVprot.fasta --report ${params.projtag}_translation_report
                awk '/^>/ { print (NR==1 ? "" : RS) \$0; next } { printf "%s", \$0 } END { printf RS }' ${params.projtag}_ASVprot.fasta > ${params.projtag}_ASVtranslations.fasta
                cp ${fasta} ${params.projtag}_ASV_all.fasta
                """
        }

        process Protein_clustering_DC { //CHECK rename_seq.py in container !!!!!!!!!!!

            label 'norm_cpus'

            publishDir "${params.workingdir}/${params.outdir}/DataCheck/ClusteringTest/Aminoacid", mode: "copy", overwrite: true, pattern: '*{.csv}'

            conda (params.condaActivate ? "-c bioconda -c conda-forge vsearch=2.21.1 cd-hit=4.8.1 seqtk=1.3 bbmap=39.01" : null)

            container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/mulled-v2-089fb9f3537921c3d6dbcc7521fbc33d82301df5:1e1ccff83e5d9864e7f3c008bd4ece458ffbdb8d-0" : "quay.io/biocontainers/mulled-v2-089fb9f3537921c3d6dbcc7521fbc33d82301df5:1e1ccff83e5d9864e7f3c008bd4ece458ffbdb8d-0")

            input:
                file(fasta) from clustering_aa
                file(asvs) from asvfastaforaaclust

            output:
                file("number_per_percentage_prot.csv") into number_per_percent_prot_plot
                file("*aminoacid_pcASV1.0_noTaxonomy.fasta") into amino_med, aminomatrix

            script:
            // add awk script to count seqs
                """
                set +e
                cp ${params.vampdir}/bin/rename_seq.py .
                for id in `echo ${params.datacheckaaIDlist} | tr "," "\\n"`;do
                    if [ \${id} == ".55" ];then
                        word=3
                    elif [ \${id} == ".65" ];then
                        word=4
                    else
                        word=5
                    fi
                    ####awk 'BEGIN{RS=">";ORS=""}length(\$2)>=${params.minAA}{print ">"\$0}' ${fasta} > ${params.projtag}_filtered_proteins.fasta
                    awk -v RS='>[^\n]+\n' 'length() >= ${params.minAA} {printf "%s", prt \$0} {prt = RT}' ${fasta} > ${params.projtag}_filtered_proteins.fasta
                    cd-hit -i ${params.projtag}_filtered_proteins.fasta -n \${word} -c \${id} -o ${params.projtag}_pcASV\${id}.fasta
                    sed 's/>Cluster />Cluster_/g' ${params.projtag}_pcASV\${id}.fasta.clstr >${params.projtag}_pcASV\${id}.clstr
                    grep ">Cluster_" ${params.projtag}_pcASV\${id}.clstr >temporaryclusters.list
                    y=\$(grep -c ">Cluster_" ${params.projtag}_pcASV\${id}.clstr)
                    echo ">Cluster_"\${y}"" >> ${params.projtag}_pcASV\${id}.clstr
                    t=1
                    b=1
                    for x in \$(cat temporaryclusters.list);do
                        echo "Extracting \$x"
                        name="\$( echo \$x | awk -F ">" '{print \$2}')"
                        clust="pcASV"\${t}""
                        echo "\${name}"
                        awk '/^>'\${name}'\$/,/^>Cluster_'\${b}'\$/' ${params.projtag}_pcASV\${id}.clstr > "\${name}"_"\${clust}"_tmp.list
                        t=\$(( \${t}+1 ))
                        b=\$(( \${b}+1 ))
                    done
                    ls *_tmp.list
                    u=1
                    for x in *_tmp.list;do
                        name="\$(echo \$x | awk -F "_p" '{print \$1}')"
                        echo "\${name}"
                        cluster="\$(echo \$x | awk -F "_" '{print \$3}')"
                        echo "\${cluster}"
                        grep "ASV" \$x | awk -F ", " '{print \$2}' | awk -F "_" '{print \$1}' | awk -F ">" '{print \$2}' > \${name}_\${cluster}_seqs_tmps.list
                        seqtk subseq ${asvs} \${name}_\${cluster}_seqs_tmps.list > \${name}_\${cluster}_nucleotide_sequences.fasta
                        vsearch --cluster_fast \${name}_\${cluster}_nucleotide_sequences.fasta --id 0.2 --centroids \${name}_\${cluster}_centroids.fasta
                        grep ">" \${name}_\${cluster}_centroids.fasta >> \${name}_\${cluster}_tmp_centroids.list
                        for y in \$( cat \${name}_\${cluster}_tmp_centroids.list );do
                            echo ">\${cluster}_type"\$u"" >> \${name}_\${cluster}_tmp_centroid.newheaders
                            u=\$(( \${u}+1 ))
                        done
                        u=1
                        ./rename_seq.py \${name}_\${cluster}_centroids.fasta \${name}_\${cluster}_tmp_centroid.newheaders \${cluster}_types_labeled.fasta
                    done
                    cat *_types_labeled.fasta >> ${params.projtag}_nucleotide_pcASV\${id}_noTaxonomy.fasta
                    grep -w "*" ${params.projtag}_pcASV\${id}.clstr | awk '{print \$3}' | awk -F "." '{print \$1}' >tmphead.list
                    grep -w "*" ${params.projtag}_pcASV\${id}.clstr | awk '{print \$2}' | awk -F "," '{print \$1}' >tmplen.list
                    paste -d"," temporaryclusters.list tmphead.list >tmp.info.csv
                    grep ">" ${params.projtag}_pcASV\${id}.fasta >lala.list
                    j=1
                    for x in \$(cat lala.list);do
                        echo ">${params.projtag}_pcASV\${j}" >>${params.projtag}_aminoheaders.list
                        echo "\${x},>${params.projtag}_pcASV\${j}" >>tmpaminotype.info.csv
                        j=\$(( \${j}+1 ))
                    done
                    rm lala.list
                    awk -F "," '{print \$2}' tmp.info.csv >>tmporder.list
                    for x in \$(cat tmporder.list);do
                        grep -w "\$x" tmpaminotype.info.csv | awk -F "," '{print \$2}' >>tmpder.list
                    done
                    paste -d "," temporaryclusters.list tmplen.list tmphead.list tmpder.list >${params.projtag}_pcASVCluster\${id}_summary.csv
                    ./rename_seq.py ${params.projtag}_pcASV\${id}.fasta ${params.projtag}_aminoheaders.list ${params.projtag}_aminoacid_pcASV\${id}_noTaxonomy.fasta
                    stats.sh in=${params.projtag}_aminoacid_pcASV\${id}_noTaxonomy.fasta gc=${params.projtag}_pcASV\${id}_aminoacid_clustered.gc gcformat=4 overwrite=true
                    stats.sh in=${params.projtag}_nucleotide_pcASV\${id}_noTaxonomy.fasta gc=${params.projtag}_pcASV\${id}_nucleotide_clustered.gc gcformat=4 overwrite=true
                    ####awk 'BEGIN{RS=">";ORS=""}length(\$2)<${params.minAA}{print ">"\$0}' ${fasta} >${params.projtag}_pcASV\${id}_problematic_translations.fasta
                    awk -v RS='>[^\n]+\n' 'length() < ${params.minAA} {printf "%s", prt \$0} {prt = RT}' ${fasta} >${params.projtag}_pcASV\${id}_problematic_translations.fasta
                    if [ `wc -l ${params.projtag}_pcASV\${id}_problematic_translations.fasta | awk '{print \$1}'` -gt 1 ];then
                        grep ">" ${params.projtag}_pcASV\${id}_problematic_translations.fasta | awk -F ">" '{print \$2}' > problem_tmp.list
                        seqtk subseq ${asvs} problem_tmp.list > ${params.projtag}_pcASV\${id}_problematic_nucleotides.fasta
                    else
                       rm ${params.projtag}_pcASV\${id}_problematic_translations.fasta
                    fi
                    rm *.list
                    rm Cluster*
                    rm *types*
                    rm *tmp*
                    rm ${params.projtag}_pcASV\${id}.fast*
                done
                for x in *aminoacid*noTaxonomy.fasta;do
                    id=\$( echo \$x | awk -F "_noTax" '{print \$1}' | awk -F "pcASV" '{print \$2}')
                    numb=\$( grep -c ">" \$x)
                    echo "\${id},\${numb}" >> number_per_percentage_protz.csv
                done
                yesirr=\$( wc -l number_per_percentage_protz.csv | awk '{print \$1}')
                tail -\$(( \${yesirr}-1 )) number_per_percentage_protz.csv > number_per_percentage_prot.csv
                head -1 number_per_percentage_protz.csv >> number_per_percentage_prot.csv
                rm number_per_percentage_protz.csv
                """
        }

        process AminoType_Pairwise_Similarity_Check {

            label 'norm_cpus'

            publishDir "${params.workingdir}/${params.outdir}/DataCheck/ClusteringTest/Aminoacid", mode: "copy", overwrite: true, pattern: '*{.csv}'

            conda (params.condaActivate ? "-c bioconda -c conda-forge clustalo=1.2.4=h87f3376_5" : null)

            container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/clustalo:1.2.4--h87f3376_5" : "quay.io/biocontainers/clustalo:1.2.4--h87f3376_5")

                input:
                    file(fasta) from aminomatrix

                output:
                    file("${params.projtag}_AminoType_PairwisePercentID.matrix") into aminopdm

                script:
                    """
                    clustalo -i ${fasta} --distmat-out=${params.projtag}_PercentIDq.matrix --percent-id --full --force --threads=${task.cpus}
                    cat ${params.projtag}_PercentIDq.matrix | tr " " "," | sed 's/,,/,/g' | grep "," >${params.projtag}_AminoType_PairwisePercentID.matrix
                    rm ${params.projtag}_PercentIDq.matrix
                    """
        }

        if (params.asvMED) {

            process ASV_Shannon_Entropy_Analysis_step1 {

                label 'norm_cpus'

                //publishDir "${params.workingdir}/${params.outdir}/DataCheck/ClusteringTest/Nucleotide/ShannonEntropy", mode: "copy", overwrite: true

                conda (params.condaActivate ? "-c bioconda -c conda-forge muscle=5.1" : null)

                container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/muscle:5.1--h9f5acd7" : "quay.io/biocontainers/muscle:5.1--h9f5acd7")

                input:
                    file(asvs) from asv_med

                output:
                    file("_ASVs_muscleAlign.fasta") into shannon_trim
                    file("*.efa")

                script:
                    """
                    pre=\$(echo ${asvs} | awk -F ".fasta" '{print \$1}' )

                    if [[ ${params.srep} == "true" && ${params.ensemble} == "false" ]];
                    then  if [[ \$( grep -c ">" ${asvs}) -lt 300 ]]
                          then    comm="align"
                          else    comm="super5"
                          fi
                          muscle -"\$comm" ${asvs} -perm ${params.perm} -perturb ${params.pert} -output ${params.projtag}_ASVs_muscleAlign.fasta -threads ${task.cpus} -quiet
                          echo "single replicate alignment chosen; look over muscle5 documentation to learn about ensemble alignment approach" >note.efa
                    elif [[ ${params.srep} == "false" && ${params.ensemble} == "true" ]];
                    then  muscle -align ${asvs} -${params.fied} -output \${pre}_muscle.efa -threads ${task.cpus} -quiet
                          muscle -maxcc \${pre}_muscle.efa -output \${pre}_muscle_raw_ALN.fasta
                    else  if [[ \$( grep -c ">" ${asvs}) -lt 300 ]]
                            then    comm="align"
                            else    comm="super5"
                            fi
                            muscle -"\$comm" ${asvs} -output \${pre}_muscle_raw_ALN.fasta -threads ${task.cpus} -quiet
                            echo "single replicate alignment chosen; look over muscle5 documentation to learn about ensemble alignment approach" >note.efa
                    fi
                    """
            }

            process ASV_Shannon_Entropy_Analysis_step2 {

              label 'norm_cpus'

              //publishDir "${params.workingdir}/${params.outdir}/DataCheck/ClusteringTest/Nucleotide/ShannonEntropy", mode: "copy", overwrite: true

              conda (params.condaActivate ? "bioconda::trimal=1.4.1=h9f5acd7_6" : null)

              container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/trimal:1.4.1--h9f5acd7_6" : "quay.io/biocontainers/trimal:1.4.1--h9f5acd7_6")

              input:
                  file(asvs) from trim

              output:
                  file("*_ASVs_muscleAligned.fasta") into shannon_trim2

              script:
                """
                #set +e
                #trimming
                trimal -in ${asvs} -out ${params.projtag}_ASVs_muscleAligned.fasta  -keepheader -fasta -automated1
                rm ${params.projtag}_ASVs_muscleAlign.fasta
                """
            }

            process ASV_Shannon_Entropy_Analysis_step3 {

              label 'norm_cpus'

              publishDir "${params.workingdir}/${params.outdir}/DataCheck/ClusteringTest/Nucleotide/ShannonEntropy", mode: "copy", overwrite: true

              conda (params.condaActivate ? "${params.vampdir}/bin/yamls/oligotyping.yml" : null)

              container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/oligotyping:2.1--py27_0" : "quay.io/biocontainers/oligotyping:2.1--py27_0")

              input:
                  file(asvs) from shannon_trim2

              output:
                  file("*_ASV_entropy_breakdown.csv") into asv_entro_csv
                  file("*Aligned_informativeonly.fasta-ENTROPY") into asv_entropy
                  file("*ASV*") into entrop

              script:
                """
                #set +e
                o-trim-uninformative-columns-from-alignment ${asvs}  #CHECK
                mv ${params.projtag}_ASVs_muscleAligned.fasta-TRIMMED ./${params.projtag}_ASVs_Aligned_informativeonly.fasta
                #entopy analysis
                entropy-analysis ${params.projtag}_ASVs_Aligned_informativeonly.fasta
                #summarize entropy peaks
                awk '{print \$2}' ${params.projtag}_ASVs_Aligned_informativeonly.fasta-ENTROPY >> tmp_value.list
                for x in \$(cat tmp_value.list)
                do      echo "\$x"
                        if [[ \$(echo "\$x > 0.0"|bc -l) -eq 1 ]];
                        then    echo dope >> above-0.0-.list
                        fi
                        if [[ \$(echo "\$x > 0.1"|bc -l) -eq 1 ]];
                        then    echo dope >> above-0.1-.list
                        fi
                        if [[ \$(echo "\$x > 0.2"|bc -l) -eq 1 ]];
                        then    echo dope >> above-0.2-.list
                        fi
                        if [[ \$(echo "\$x > 0.3"|bc -l) -eq 1 ]];
                        then    echo dope >> above-0.3-.list
                        fi
                        if [[ \$(echo "\$x > 0.4"|bc -l) -eq 1 ]];
                        then    echo dope >> above-0.4-.list
                        fi
                        if [[ \$(echo "\$x > 0.5"|bc -l) -eq 1 ]];
                        then    echo dope >> above-0.5-.list
                        fi
                        if [[ \$(echo "\$x > 0.6"|bc -l) -eq 1 ]];
                        then    echo dope >> above-0.6-.list
                        fi
                        if [[ \$(echo "\$x > 0.7"|bc -l) -eq 1 ]];
                        then    echo dope >> above-0.7-.list
                        fi
                        if [[ \$(echo "\$x > 0.8"|bc -l) -eq 1 ]];
                        then    echo dope >> above-0.8-.list
                        fi
                        if [[ \$(echo "\$x > 0.9"|bc -l) -eq 1 ]];
                        then    echo dope >> above-0.9-.list
                        fi
                        if [[ \$(echo "\$x > 1.0"|bc -l) -eq 1 ]];
                        then    echo dope >> above-1.0-.list
                        fi
                        if [[ \$(echo "\$x > 1.1"|bc -l) -eq 1 ]];
                        then    echo dope >> above-1.1-.list
                        fi
                        if [[ \$(echo "\$x > 1.2"|bc -l) -eq 1 ]];
                        then    echo dope >> above-1.2-.list
                        fi
                        if [[ \$(echo "\$x > 1.3"|bc -l) -eq 1 ]];
                        then    echo dope >> above-1.3-.list
                        fi
                        if [[ \$(echo "\$x > 1.4"|bc -l) -eq 1 ]];
                        then    echo dope >> above-1.4-.list
                        fi
                        if [[ \$(echo "\$x > 1.5"|bc -l) -eq 1 ]];
                        then    echo dope >> above-1.5-.list
                        fi
                        if [[ \$(echo "\$x > 1.6"|bc -l) -eq 1 ]];
                        then    echo dope >> above-1.6-.list
                        fi
                        if [[ \$(echo "\$x > 1.7"|bc -l) -eq 1 ]];
                        then    echo dope >> above-1.7-.list
                        fi
                        if [[ \$(echo "\$x > 1.8"|bc -l) -eq 1 ]];
                        then    echo dope >> above-1.8-.list
                        fi
                        if [[ \$(echo "\$x > 1.9"|bc -l) -eq 1 ]];
                        then    echo dope >> above-1.9-.list
                        fi
                        if [[ \$(echo "\$x > 2.0"|bc -l) -eq 1 ]];
                        then    echo dope >> above-2.0-.list
                        fi
                done
                echo "Entropy,Peaks_above" >> ${params.projtag}_ASV_entropy_breakdown.csv
                for z in above*.list;
                do      entrop=\$(echo \$z | awk -F "-" '{print \$2}')
                        echo ""\$entrop", "\$(wc -l \$z | awk '{print \$1}')"" >> ${params.projtag}_ASV_entropy_breakdown.csv
                done
                rm above*
                mv ${params.projtag}_ASVs_Aligned_informativeonly.fasta-ENTROPY ./tmp2.tsv
                cat tmp2.tsv | tr "\\t" "," > tmp.csv
                rm tmp2.tsv
                echo "Base_position,Shannons_Entropy" >> ${params.projtag}_ASVs_Aligned_informativeonly.fasta-ENTROPY
                cat tmp.csv >> ${params.projtag}_ASVs_Aligned_informativeonly.fasta-ENTROPY
                rm tmp.csv
                """
            }

        } else {
            asv_entropy = Channel.empty()
            asv_entro_csv = Channel.empty()
        }

        if (params.aminoMED) {

                process AminoType_Shannon_Entropy_Analysis_step1 {

                    label 'norm_cpus'

                    //publishDir "${params.workingdir}/${params.outdir}/DataCheck/ClusteringTest/Aminoacid/ShannonEntropy", mode: "copy", overwrite: true

                    conda (params.condaActivate ? "-c bioconda -c conda-forge muscle=5.1" : null)

                    container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/muscle:5.1--h9f5acd7" : "quay.io/biocontainers/muscle:5.1--h9f5acd7")

                    input:
                        file(aminos) from amino_med

                    output:
                        file("*_AminoTypes_muscleAlign.fasta") into amino_shannon_trim1
                        file("*.efa")

                    script:
                        """
                        pre=\$(echo ${aminos} | awk -F ".fasta" '{print \$1}' )

                        if [[ ${params.srep} == "true" && ${params.ensemble} == "false" ]];
                        then  if [[ \$( grep -c ">" ${aminos}) -lt 300 ]]
                              then    comm="align"
                              else    comm="super5"
                              fi
                              muscle -"\$comm" ${aminos} -perm ${params.perm} -perturb ${params.pert} -output ${params.projtag}_AminoTypes_muscleAlign.fasta -threads ${task.cpus} -quiet
                              echo "single replicate alignment chosen; look over muscle5 documentation to learn about ensemble alignment approach" >note.efa
                        elif [[ ${params.srep} == "false" && ${params.ensemble} == "true" ]];
                        then  muscle -align ${aminos} -${params.fied} -output \${pre}_muscle.efa -threads ${task.cpus} -quiet
                              muscle -maxcc \${pre}_muscle.efa -output \${pre}_muscle_raw_ALN.fasta
                        else  if [[ \$( grep -c ">" ${aminos}) -lt 300 ]]
                                then    comm="align"
                                else    comm="super5"
                                fi
                                muscle -"\$comm" ${aminos} -output \${pre}_muscle_raw_ALN.fasta -threads ${task.cpus} -quiet
                                echo "single replicate alignment chosen; look over muscle5 documentation to learn about ensemble alignment approach" >note.efa
                        fi
                        """
                }

                process AminoType_Shannon_Entropy_Analysis_step2 {

                  label 'norm_cpus'
                  conda (params.condaActivate ? "bioconda::trimal=1.4.1=h9f5acd7_6" : null)

                  container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/trimal:1.4.1--h9f5acd7_6" : "quay.io/biocontainers/trimal:1.4.1--h9f5acd7_6")

                  //publishDir "${params.workingdir}/${params.outdir}/DataCheck/ClusteringTest/Aminoacid/ShannonEntropy", mode: "copy", overwrite: true

                  conda (params.condaActivate ? "bioconda::trimal=1.4.1=h9f5acd7_6" : null)

                  container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/trimal:1.4.1--h9f5acd7_6" : "quay.io/biocontainers/trimal:1.4.1--h9f5acd7_6")

                  input:
                      file(aminos) from amino_shannon_trim1

                  output:
                      file("*_AminoTypes_muscleAligned.fasta") into amino_shannon_trim2

                  script:
                    """
                    #trimming
                    trimal -in ${aminos} -out ${params.projtag}_AminoTypes_muscleAligned.fasta  -keepheader -fasta -automated1
                    """
                }

                process AminoType_Shannon_Entropy_Analysis_step3 {

                    label 'norm_cpus'

                    publishDir "${params.workingdir}/${params.outdir}/DataCheck/ClusteringTest/Aminoacid/ShannonEntropy", mode: "copy", overwrite: true, pattern: '*{.csv}'

                    conda (params.condaActivate ? "${params.vampdir}/bin/yamls/oligotyping.yml" : null)

                    container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/oligotyping:2.1--py27_0" : "quay.io/biocontainers/oligotyping:2.1--py27_0")

                    input:
                        file(aminos) from amino_shannon_trim2

                    output:
                        file("*AminoType_entropy_breakdown.csv") into amino_entro_csv
                        file ("*Aligned_informativeonly.fasta-ENTROPY") into amino_entropy
                        file("*AminoTypes*") into aminos

                    script:
                        """
                        # CHECK
                        o-trim-uninformative-columns-from-alignment ${aminos}
                        mv ${params.projtag}_AminoTypes_muscleAligned.fasta-TRIMMED ./${params.projtag}_AminoTypes_Aligned_informativeonly.fasta
                        #entropy analysis
                        entropy-analysis ${params.projtag}_AminoTypes_Aligned_informativeonly.fasta --amino-acid-sequences
                        #summarize entropy peaks
                        awk '{print \$2}' ${params.projtag}_AminoTypes_Aligned_informativeonly.fasta-ENTROPY >> tmp_value.list
                        for x in \$(cat tmp_value.list)
                        do      echo "\$x"
                                if [[ \$(echo "\$x > 0.0"|bc -l) -eq 1 ]];
                                then    echo dope >> above-0.0-.list
                                fi
                                if [[ \$(echo "\$x > 0.1"|bc -l) -eq 1 ]];
                                then    echo dope >> above-0.1-.list
                                fi
                                if [[ \$(echo "\$x > 0.2"|bc -l) -eq 1 ]];
                                then    echo dope >> above-0.2-.list
                                fi
                                if [[ \$(echo "\$x > 0.3"|bc -l) -eq 1 ]];
                                then    echo dope >> above-0.3-.list
                                fi
                                if [[ \$(echo "\$x > 0.4"|bc -l) -eq 1 ]];
                                then    echo dope >> above-0.4-.list
                                fi
                                if [[ \$(echo "\$x > 0.5"|bc -l) -eq 1 ]];
                                then    echo dope >> above-0.5-.list
                                fi
                                if [[ \$(echo "\$x > 0.6"|bc -l) -eq 1 ]];
                                then    echo dope >> above-0.6-.list
                                fi
                                if [[ \$(echo "\$x > 0.7"|bc -l) -eq 1 ]];
                                then    echo dope >> above-0.7-.list
                                fi
                                if [[ \$(echo "\$x > 0.8"|bc -l) -eq 1 ]];
                                then    echo dope >> above-0.8-.list
                                fi
                                if [[ \$(echo "\$x > 0.9"|bc -l) -eq 1 ]];
                                then    echo dope >> above-0.9-.list
                                fi
                                if [[ \$(echo "\$x > 1.0"|bc -l) -eq 1 ]];
                                then    echo dope >> above-1.0-.list
                                fi
                                if [[ \$(echo "\$x > 1.1"|bc -l) -eq 1 ]];
                                then    echo dope >> above-1.1-.list
                                fi
                                if [[ \$(echo "\$x > 1.2"|bc -l) -eq 1 ]];
                                then    echo dope >> above-1.2-.list
                                fi
                                if [[ \$(echo "\$x > 1.3"|bc -l) -eq 1 ]];
                                then    echo dope >> above-1.3-.list
                                fi
                                if [[ \$(echo "\$x > 1.4"|bc -l) -eq 1 ]];
                                then    echo dope >> above-1.4-.list
                                fi
                                if [[ \$(echo "\$x > 1.5"|bc -l) -eq 1 ]];
                                then    echo dope >> above-1.5-.list
                                fi
                                if [[ \$(echo "\$x > 1.6"|bc -l) -eq 1 ]];
                                then    echo dope >> above-1.6-.list
                                fi
                                if [[ \$(echo "\$x > 1.7"|bc -l) -eq 1 ]];
                                then    echo dope >> above-1.7-.list
                                fi
                                if [[ \$(echo "\$x > 1.8"|bc -l) -eq 1 ]];
                                then    echo dope >> above-1.8-.list
                                fi
                                if [[ \$(echo "\$x > 1.9"|bc -l) -eq 1 ]];
                                then    echo dope >> above-1.9-.list
                                fi
                                if [[ \$(echo "\$x > 2.0"|bc -l) -eq 1 ]];
                                then    echo dope >> above-2.0-.list
                                fi
                        done
                        echo "Entropy,Peaks_above" >> ${params.projtag}_AminoType_entropy_breakdown.csv
                        for z in above*.list;
                        do      entrop=\$(echo \$z | awk -F "-" '{print \$2}')
                                echo ""\$entrop", "\$(wc -l \$z | awk '{print \$1}')"" >> ${params.projtag}_AminoType_entropy_breakdown.csv
                        done
                        rm above*
                        mv ${params.projtag}_AminoTypes_Aligned_informativeonly.fasta-ENTROPY ./tmp2.tsv
                        cat tmp2.tsv | tr "\t" "," > tmp.csv
                        rm tmp2.tsv
                        echo "Base_position,Shannons_Entropy" >> ${params.projtag}_AminoTypes_Aligned_informativeonly.fasta-ENTROPY
                        cat tmp.csv >> ${params.projtag}_AminoTypes_Aligned_informativeonly.fasta-ENTROPY
                        rm tmp.csv
                        """
                }

        } else {
            amino_entro_csv = Channel.empty()
            amino_entropy = Channel.empty()
        }

        if (!params.skipReadProcessing || !params.skipMerging ) {

            process combine_csv_DC {

                conda (params.condaActivate ? null : null)

                container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? null : null)

                input:
                    file(csv) from fastp_csv_in1
                        .collect()

                output:
                    file("final_reads_stats.csv") into fastp_csv_dc

                script:
                    """
                    cat ${csv} >all_reads_stats.csv
                    head -n1 all_reads_stats.csv >tmp.names.csv
                    cat all_reads_stats.csv | grep -v ""Sample,Total_"" >tmp.reads.stats.csv
                    cat tmp.names.csv tmp.reads.stats.csv >final_reads_stats.csv
                    rm tmp.names.csv tmp.reads.stats.csv
                    """
            }

        } else {

            process skip_combine_csv_DC {
                output:
                    file("filter_reads.txt") into fastp_csv_dc

                    conda (params.condaActivate ? null : null)

                    container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? null : null)

                script:
                    """
                    echo "Read processing steps skipped." >filter_reads.txt
                    """
            }
        }

        report_dc_in = Channel.create()
        fastp_csv_dc.mix( reads_per_sample_preFilt, reads_per_sample_postFilt, prefilt_basefreq, postFilt_basefreq, prefilt_qualityscore, postFilt_qualityscore, prefilt_gccontent, postFilt_gccontent, prefilt_averagequality, postFilt_averagequality, prefilt_length, postFilt_length, number_per_percent_nucl_plot, number_per_percent_prot_plot, amino_entro_csv, amino_entropy, asv_entro_csv, asv_entropy, asvpdm, aminopdm).into(report_dc_in)

        process Report_DataCheck {

            label 'norm_cpus'

            publishDir "${params.workingdir}/${params.outdir}/DataCheck/Report", mode: "copy", overwrite: true, pattern: '*.{html}'

            conda (params.condaActivate ? "${params.vampdir}/bin/yamls/R.yml" : null)

            container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/mulled-v2-90fb71778d446b480b1050805be9ef346794df6b:6d9281817ba120685d81edae0370730c1ff554cc-0" : "quay.io/biocontainers/mulled-v2-90fb71778d446b480b1050805be9ef346794df6b:6d9281817ba120685d81edae0370730c1ff554cc-0")

            input:
                file(files) from report_dc_in
                    .collect()

            output:
                file("*.html") into datacheckreport

            script:
                """
                cp ${params.vampdir}/bin/vAMPirus_DC_Report.Rmd .
                cp ${params.vampdir}/example_data/conf/vamplogo.png .
                Rscript -e "rmarkdown::render('vAMPirus_DC_Report.Rmd',output_file='${params.projtag}_DataCheck_Report.html')" ${params.projtag} \
                ${params.skipReadProcessing} \
                ${params.skipMerging} \
                ${params.skipAdapterRemoval} \
                ${params.asvMED} \
                ${params.aminoMED}
                """
        }

    } else if (params.Analyze) {

        if (params.ncASV) {

            reads_vsearch5_ch
    	       .into{ asv_file_for_ncasvs; nuclFastas_forDiamond_asv_ch; nuclFastas_forCounts_asv_ch; nuclFastas_forphylogeny_asv; nuclFastas_forMatrix_asv_ch}

            process NucleotideBased_ASV_clustering {

                label 'norm_cpus'

                tag "${mtag}"

                publishDir "${params.workingdir}/${params.outdir}/Analyze/Clustering/ncASV", mode: "copy", overwrite: true, pattern: '*ncASV*.fasta'

                conda (params.condaActivate ? "-c bioconda -c conda-forge vsearch=2.21.1=hf1761c0_1" : null)

                container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/vsearch:2.21.1--hf1761c0_1" : "quay.io/biocontainers/vsearch:2.21.1--hf1761c0_1")

                input:
                    each x from 1..nnuc
                    file(fasta) from asv_file_for_ncasvs

                output:
                    tuple nid, file("*_ncASV*.fasta") into ( nuclFastas_forphylogeny_ncasv, nuclFastas_forphylogeny_ncasv2, nuclFastas_forDiamond_ncasv_ch, nuclFastas_forCounts_ncasv_ch, nuclFastas_forMatrix_ncasv_ch )

                script:
                    nid=slist.get(x-1)
                    mtag="ID=" + slist.get(x-1)
                    """
                    vsearch --cluster_fast ${fasta} --centroids ${params.projtag}_ncASV${nid}.fasta --threads ${task.cpus} --relabel ncASV --id .${nid}
                    """
            }

            if (!params.skipTaxonomy) {

              if (params.dbtype == "NCBI") {

                process ncASV_Taxonomy_Inference_NCBI { // edit !!!!!!!! //CHECK rename_seq.py in container !!!!!!!!!!!

                    label 'high_cpus'

                    tag "${mtag}"

                    publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ncASV/Taxonomy", mode: "copy", overwrite: true, pattern: '*ncASV*.{fasta,csv,tsv}'
                    publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ncASV/Taxonomy/DiamondOutput", mode: "copy", overwrite: true, pattern: '*ncASV*dmd.out'

                    conda (params.condaActivate ? "-c conda-forge bioconda::diamond=2.0.15=hb97b32f_1" : null)

                    container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/diamond:2.0.15--hb97b32f_1" : "quay.io/biocontainers/diamond:2.0.15--hb97b32f_1")

                    input:
                        tuple nid, file(asvs) from nuclFastas_forDiamond_ncasv_ch

                    output:
                        file("*.fasta") into tax_labeled_fasta_ncasv
                        tuple file("*_phyloformat.csv"), file("*summaryTable.tsv"), file("*dmd.out") into summary_diamond_ncasv
                        tuple nid, file("*ncASV*summary_for_plot.csv") into taxplot_ncasv
                        tuple nid, file("*_quick_Taxbreakdown.csv") into tax_table_ncasv
                        tuple nid, file ("*_quicker_taxbreakdown.csv") into tax_nodCol_ncasv

                    script:
                        mtag="ID=" + nid
                        """
                        cp ${params.vampdir}/bin/rename_seq.py .
                        virdb=${params.dbdir}/${params.dbname}
                        if [[ ${params.measurement} == "bitscore" ]]
                        then    measure="--min-score ${params.bitscore}"
                        elif    [[ ${params.measurement} == "evalue" ]]
                        then    measure="-e ${params.evalue}"
                        else    measure="--min-score ${params.bitscore}"
                        fi
                        grep ">" \${virdb} > headers.list
                        headers="headers.list"
                        name=\$( echo ${asvs} | awk -F ".fasta" '{print \$1}')
                        if [[ ${params.ncbitax} == "true" ]]
                        then   diamond blastx -q ${asvs} -d \${virdb} -p ${task.cpus} --id ${params.minID} -l ${params.minaln} \${measure} --${params.sensitivity} -o "\$name"_dmd.out -f 6 qseqid qlen sseqid qstart qend qseq sseq length qframe evalue bitscore pident btop staxids sskingdoms skingdoms sphylums --max-target-seqs 1 --max-hsps 1
                        else   diamond blastx -q ${asvs} -d \${virdb} -p ${task.cpus} --id ${params.minID} -l ${params.minaln} \${measure} --${params.sensitivity} -o "\$name"_dmd.out -f 6 qseqid qlen sseqid qstart qend qseq sseq length qframe evalue bitscore pident btop --max-target-seqs 1 --max-hsps 1
                        fi
                        echo "Preparing lists to generate summary .csv's"
                        echo "[Best hit accession number]" > access.list
                        echo "[e-value]" > evalue.list
                        echo "[Bitscore]" > bit.list
                        echo "[Percent ID (aa)]" > pid.list
                        echo "[Organism ID]" > "\$name"_virus.list
                        echo "[Gene]" > "\$name"_genes.list
                        echo "[ncASV#]" > otu.list
                        echo "[Sequence length]" > length.list
                        grep ">" ${asvs} | awk -F ">" '{print \$2}' > seqids.lst
                        if [[ ${params.lca} == "T" ]]
                        then  grep -w "LCA" ${params.dbanno}/*.txt > lcainfo.list
                              echo "[Taxonomic classification from RVDB annotations]" > lca_classification.list
                        else
                              echo "[Taxonomic classification from RVDB annotations]" > lca_classification.list
                        fi
                        if [[ ${params.ncbitax} == "true" ]]
                        then echo "[NCBI Taxonomy ID],[Taxonomic classification from NCBI]" > ncbi_classification.list
                        fi
                        echo "extracting genes and names"
                        touch new_"\$name"_asvnames.txt
                        for s in \$(cat seqids.lst);do
                            echo "Checking for \$s hit in diamond output"
                            if [[ "\$(grep -wc "\$s" "\$name"_dmd.out)" -eq 1 ]];then
                                        echo "Yep, there was a hit for \$s"
                                        echo "Extracting the information now:"
                                        acc=\$(grep -w "\$s" "\$name"_dmd.out | awk '{print \$3}')
                                        echo "\$s" >> otu.list
                                        echo "\$acc" >> access.list
                                        line="\$(grep -w "\$s" "\$name"_dmd.out)"
                                        echo "\$line" | awk '{print \$10}' >> evalue.list
                                        echo "\$line" | awk '{print \$11}' >> bit.list
                                        echo "\$line" | awk '{print \$12}' >> pid.list
                                        echo "\$line" | awk '{print \$2}' >> length.list
                                        echo "Extracting virus and gene ID for \$s now"
                                        gene=\$(grep -w "\$acc" "\$headers" | awk -F "." '{ print \$2 }' | awk -F "[" '{ print \$1 }' | awk -F " " '{print substr(\$0, index(\$0,\$2))}' | sed 's/ /_/g') &&
                                        echo "\$gene" | sed 's/_/ /g' >> "\$name"_genes.list
                                        virus=\$(grep -w "\$acc" "\$headers" | awk -F "[" '{ print \$2 }' | awk -F "]" '{ print \$1 }'| sed 's/ /_/g')
                                        echo "\$virus" | sed 's/_/ /g' >> "\$name"_virus.list
                                        echo ">"\${s}"_"\$virus"_"\$gene"" >> new_"\$name"_asvnames.txt
                                        if [[ "${params.lca}" == "T" ]]
                                        then    if [[ \$(grep -w "\$acc" ${params.dbanno}/*.txt | wc -l) -eq 1 ]]
                                                then  group=\$(grep -w "\$acc" ${params.dbanno}/*.txt | awk -F ":" '{print \$1}')
                                                      lcla=\$(grep -w "\$group" lcainfo.list | awk -F "\t" '{print \$2}')
                                                      echo "\$lcla" >> lca_classification.list
                                                else  echo "Viruses" >> lca_classification.list
                                                fi
                                        fi
                                        if [[ ${params.ncbitax} == "true" ]]
                                        then  echo "\$line" | awk -F "\t" '{print \$14","\$16"::"\$18"::"\$17}' >> ncbi_classification.list
                                        fi
                                        echo "\$s done."
                            else
                                        echo "Ugh, there was no hit for \$s .."
                                        echo "We still love \$s though and we will add it to the final fasta file"
                                        echo "\$s" >> otu.list
                                        echo "NO_HIT" >> access.list
                                        echo "NO_HIT" >> "\$name"_genes.list
                                        echo "NO_HIT" >> "\$name"_virus.list
                                        echo "NO_HIT" >> evalue.list
                                        echo "NO_HIT" >> bit.list
                                        echo "NO_HIT" >> pid.list
                                        echo "NO_HIT" >> length.list
                                        virus="NO"
                                        gene="HIT"
                                        echo ">\${s}_"\$virus"_"\$gene"" >> new_"\$name"_asvnames.txt
                                        if [[ "${params.lca}" == "T" ]]
                                        then    echo "N/A" >> lca_classification.list
                                        fi
                                        if [[ "${params.ncbitax}" == "true" ]]
                                        then  echo "N/A" >> ncbi_classification.list
                                        fi
                                        echo "\$s done."
                           fi
                        done
                        echo "Now editing "\$name" fasta headers"
                        ###### rename_seq.py
                        ./rename_seq.py ${asvs} new_"\$name"_asvnames.txt "\$name"_TaxonomyLabels.fasta
                        awk 'BEGIN {RS=">";FS="\\n";OFS=""} NR>1 {print ">"\$1; \$1=""; print}' "\$name"_TaxonomyLabels.fasta >"\$name"_tmpssasv.fasta
                        echo "[Sequence header]" > newnames.list
                        cat new_"\$name"_asvnames.txt >> newnames.list
                        touch sequence.list
                        echo "     " > sequence.list
                        grep -v ">" "\$name"_tmpssasv.fasta >> sequence.list
                        rm "\$name"_tmpssasv.fasta
                        if [[ "${params.lca}" == "T" && "${params.ncbitax}" == "true" ]]
                        then
                              paste -d "," sequence.list "\$name"_virus.list "\$name"_genes.list otu.list lca_classification.list ncbi_classification.list newnames.list length.list bit.list evalue.list pid.list access.list >> "\$name"_phyloformat.csv
                              paste -d"\t" otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list ncbi_classification.list sequence.list length.list bit.list evalue.list pid.list >> "\$name"_summaryTable.tsv
                              paste -d"," otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list ncbi_classification.list >> \${name}_quick_Taxbreakdown.csv
                        elif [[ "${params.lca}" == "T" && "${params.ncbitax}" != "true" ]]
                        then
                              paste -d "," sequence.list "\$name"_virus.list "\$name"_genes.list otu.list lca_classification.list newnames.list length.list bit.list evalue.list pid.list access.list >> "\$name"_phyloformat.csv
                              paste -d"\t" otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list sequence.list length.list bit.list evalue.list pid.list >> "\$name"_summaryTable.tsv
                              paste -d"," otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list >> \${name}_quick_Taxbreakdown.csv
                        elif [[ "${params.ncbitax}" == "true" && "${params.lca}" != "T" ]]
                        then
                              paste -d "," sequence.list "\$name"_virus.list "\$name"_genes.list otu.list ncbi_classification.list newnames.list length.list bit.list evalue.list pid.list access.list >> "\$name"_phyloformat.csv
                              paste -d"\t" otu.list access.list "\$name"_virus.list "\$name"_genes.list ncbi_classification.list sequence.list length.list bit.list evalue.list pid.list >> "\$name"_summaryTable.tsv
                              paste -d"," otu.list access.list "\$name"_virus.list "\$name"_genes.list ncbi_classification.list >> \${name}_quick_Taxbreakdown.csv
                        else
                              paste -d "," sequence.list "\$name"_virus.list "\$name"_genes.list otu.list newnames.list length.list bit.list evalue.list pid.list access.list >> "\$name"_phyloformat.csv
                              paste -d"\t" otu.list access.list "\$name"_virus.list "\$name"_genes.list sequence.list length.list bit.list evalue.list pid.list >> "\$name"_summaryTable.tsv
                              echo "skipped" >> \${name}_quick_Taxbreakdown.csv
                        fi
                        for x in *phyloformat.csv;do
                            echo "\$x"
                            lin=\$(( \$(wc -l \$x | awk '{print \$1}')-1))
                            tail -"\$lin" \$x | awk -F "," '{print \$2}' > tmpcol.list;
                            sed 's/ /_/g' tmpcol.list > tmp2col.list;
                            cat tmp2col.list | sort | uniq -c | sort -nr | awk '{print \$2","\$1}' > \${name}_summary_for_plot.csv;
                            rm tmpcol.list tmp2col.list
                        done
                        awk -F "," '{print \$1","\$3"("\$2")"}' \${name}_quick_Taxbreakdown.csv >> \${name}_quicker_taxbreakdown.csv
                        rm evalue.list sequence.list bit.list pid.list length.list seqids.lst otu.list *asvnames.txt "\$name"_virus.list "\$name"_genes.list newnames.list access.list headers.list
                        """
                      }

            } else if (params.dbtype== "RVDB") {

              process ncASV_Taxonomy_Inference_RVDB { // edit !!!!!!! //CHECK rename_seq.py in container !!!!!!!!!!!

                  label 'high_cpus'

                  tag "${mtag}"

                  publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ncASV/Taxonomy", mode: "copy", overwrite: true, pattern: '*ncASV*.{fasta,csv,tsv}'
                  publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ncASV/Taxonomy/DiamondOutput", mode: "copy", overwrite: true, pattern: '*ncASV*dmd.out'

                  conda (params.condaActivate ? "-c conda-forge bioconda::diamond=2.0.15=hb97b32f_1" : null)

                  container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/diamond:2.0.15--hb97b32f_1" : "quay.io/biocontainers/diamond:2.0.15--hb97b32f_1")

                  input:
                      tuple nid, file(asvs) from nuclFastas_forDiamond_ncasv_ch

                  output:
                      file("*.fasta") into tax_labeled_fasta_ncasv
                      tuple file("*_phyloformat.csv"), file("*summaryTable.tsv"), file("*dmd.out") into summary_diamond_ncasv
                      tuple nid, file("*ncASV*summary_for_plot.csv") into taxplot_ncasv
                      tuple nid, file("*_quick_Taxbreakdown.csv") into tax_table_ncasv
                      tuple nid, file ("*_quicker_taxbreakdown.csv") into tax_nodCol_ncasv

                  script:
                      mtag="ID=" + nid
                      """
                      cp ${params.vampdir}/bin/rename_seq.py .
                      virdb=${params.dbdir}/${params.dbname}
                      if [[ ${params.measurement} == "bitscore" ]]
                      then    measure="--min-score ${params.bitscore}"
                      elif    [[ ${params.measurement} == "evalue" ]]
                      then    measure="-e ${params.evalue}"
                      else    measure="--min-score ${params.bitscore}"
                      fi
                      grep ">" \${virdb} > headers.list
                      headers="headers.list"
                      name=\$( echo ${asvs} | awk -F ".fasta" '{print \$1}')
                      diamond blastx -q ${asvs} -d \${virdb} -p ${task.cpus} --id ${params.minID} -l ${params.minaln} \${measure} --${params.sensitivity} -o "\$name"_dmd.out -f 6 qseqid qlen sseqid qstart qend qseq sseq length qframe evalue bitscore pident btop --max-target-seqs 1 --max-hsps 1
                      echo "Preparing lists to generate summary .csv's"
                      echo "[Best hit accession number]" > access.list
                      echo "[e-value]" > evalue.list
                      echo "[Bitscore]" > bit.list
                      echo "[Percent ID (aa)]" > pid.list
                      echo "[Organism ID]" > "\$name"_virus.list
                      echo "[Gene]" > "\$name"_genes.list
                      echo "[ncASV#]" > otu.list
                      echo "[Sequence length]" > length.list
                      grep ">" ${asvs} | awk -F ">" '{print \$2}' > seqids.lst
                      if [[ ${params.lca} == "T" ]]
                      then  grep -w "LCA" ${params.dbanno}/*.txt > lcainfo.list
                            echo "[Taxonomic classification from RVDB annotations]" > lca_classification.list
                      else  echo "skipped" >> \${name}_quick_Taxbreakdown.csv
                            echo "[Taxonomic classification from RVDB annotations]" > lca_classification.list
                      fi
                      echo "extracting genes and names"
                      touch new_"\$name"_asvnames.txt
                      for s in \$(cat seqids.lst);do
                          echo "Using RVDB headers."
                          if [[ "\$(grep -wc "\$s" "\$name"_dmd.out)" -eq 1 ]];then
                              echo "Yep, there was a hit for \$s"
                              echo "Extracting the information now:"
                              acc=\$(grep -w "\$s" "\$name"_dmd.out | awk '{print \$3}' | awk -F "|" '{print \$3}')
                              echo "\$s" >> otu.list
                              echo "\$acc" >> access.list
                              line="\$(grep -w "\$s" "\$name"_dmd.out)"
                              echo "\$line" | awk '{print \$10}' >> evalue.list
                              echo "\$line" | awk '{print \$11}' >> bit.list
                              echo "\$line" | awk '{print \$12}' >> pid.list
                              echo "\$line" | awk '{print \$2}' >> length.list
                              echo "Extracting virus and gene ID for \$s now"
                              gene=\$(grep -w "\$acc" "\$headers" | awk -F "|" '{ print \$6 }' | awk -F "[" '{ print \$1 }' | sed 's/ /_/g') &&
                              echo "\$gene" | sed 's/_/ /g' >> "\$name"_genes.list
                              virus=\$(grep -w "\$acc" "\$headers" | awk -F "|" '{ print \$6 }' | awk -F "[" '{ print \$2 }' | awk -F "]" '{print \$1}' | sed 's/ /_/g') &&
                              echo "\$virus" | sed 's/_/ /g' >> "\$name"_virus.list
                              echo ">\${s}_"\$virus"_"\$gene"" >> new_"\$name"_asvnames.txt
                              if [[ "${params.lca}" == "T" ]]
                              then    if [[ \$(grep -w "\$acc" ${params.dbanno}/*.txt | wc -l) -eq 1 ]]
                                      then  group=\$(grep -w "\$acc" ${params.dbanno}/*.txt | awk -F ":" '{print \$1}')
                                            lcla=\$(grep -w "\$group" lcainfo.list | awk -F "\t" '{print \$2}')
                                            echo "\$lcla" >> lca_classification.list
                                      else  echo "Viruses" >> lca_classification.list
                                      fi
                              fi
                              echo "\$s done."
                          else
                              echo "Ugh, there was no hit for \$s .."
                              echo "We still love \$s though and we will add it to the final fasta file"
                              echo "\$s" >> otu.list
                              echo "NO_HIT" >> access.list
                              echo "NO_HIT" >> "\$name"_genes.list
                              echo "NO_HIT" >> "\$name"_virus.list
                              echo "NO_HIT" >> evalue.list
                              echo "NO_HIT" >> bit.list
                              echo "NO_HIT" >> pid.list
                              echo "NO_HIT" >> length.list
                              virus="NO"
                              gene="HIT"
                              echo ">\${s}_"\$virus"_"\$gene"" >> new_"\$name"_asvnames.txt
                              if [[ "${params.lca}" == "T" ]]
                              then    echo "N/A" >> lca_classification.list
                              fi
                              echo "\$s done."
                          fi
                      echo "Done with \$s"
                      done
                      echo "Now editing "\$name" fasta headers"
                      ###### rename_seq.py
                      ./rename_seq.py ${asvs} new_"\$name"_asvnames.txt "\$name"_TaxonomyLabels.fasta
                      awk 'BEGIN {RS=">";FS="\\n";OFS=""} NR>1 {print ">"\$1; \$1=""; print}' "\$name"_TaxonomyLabels.fasta >"\$name"_tmpssasv.fasta
                      echo "[Sequence header]" > newnames.list
                      cat new_"\$name"_asvnames.txt >> newnames.list
                      touch sequence.list
                      echo "     " > sequence.list
                      grep -v ">" "\$name"_tmpssasv.fasta >> sequence.list
                      rm "\$name"_tmpssasv.fasta
                      if [[ "${params.lca}" == "T" ]]
                      then  paste -d "," sequence.list "\$name"_virus.list "\$name"_genes.list otu.list lca_classification.list newnames.list length.list bit.list evalue.list pid.list access.list >> "\$name"_phyloformat.csv
                            paste -d"\t" otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list sequence.list length.list bit.list evalue.list pid.list >> "\$name"_summaryTable.tsv
                            paste -d"," otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list >> \${name}_quick_Taxbreakdown.csv
                      else  paste -d "," sequence.list "\$name"_virus.list "\$name"_genes.list otu.list newnames.list length.list bit.list evalue.list pid.list access.list >> "\$name"_phyloformat.csv
                            paste -d"\t" otu.list access.list "\$name"_virus.list "\$name"_genes.list sequence.list length.list bit.list evalue.list pid.list >> "\$name"_summaryTable.tsv
                      fi
                      for x in *phyloformat.csv;do
                                echo "\$x"
                                lin=\$(( \$(wc -l \$x | awk '{print \$1}')-1))
                                tail -"\$lin" \$x | awk -F "," '{print \$2}' > tmpcol.list;
                                sed 's/ /_/g' tmpcol.list > tmp2col.list;
                                cat tmp2col.list | sort | uniq -c | sort -nr | awk '{print \$2","\$1}' > \${name}_summary_for_plot.csv;
                                rm tmpcol.list tmp2col.list
                      done
                      awk -F "," '{print \$1","\$3"("\$2")"}' \${name}_quick_Taxbreakdown.csv >> \${name}_quicker_taxbreakdown.csv
                      rm evalue.list sequence.list bit.list pid.list length.list seqids.lst otu.list *asvnames.txt "\$name"_virus.list "\$name"_genes.list newnames.list access.list headers.list
                      """
                }
            }
          } else {

              process skipncASVtaxonomy {

                  input:
                      tuple nid, file(asvs) from nuclFastas_forDiamond_ncasv_ch

                  output:
                      tuple nid, file("skipncASVtaxonomy1.txt") into ( taxplot_ncasv )
                      tuple nid, file("skipncASVtaxonomy2.txt") into ( tax_table_ncasv )
                      tuple nid, file("skipncASVtaxonomy3.txt") into ( tax_nodCol_ncasv )

                  script:
                      """
                      echo "Skipped" >skipncASVtaxonomy1.txt
                      echo "Skipped" >skipncASVtaxonomy2.txt
                      echo "Skipped" >skipncASVtaxonomy3.txt
                      """
              }
        }

            process Generate_ncASV_Counts_Table {

                label 'norm_cpus'

                tag "${mtag}"

                publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ncASV/Counts", mode: "copy", overwrite: true, pattern: '*ncASV*.{biome,csv}'

                conda (params.condaActivate ? "-c bioconda -c conda-forge vsearch=2.21.1=hf1761c0_1" : null)

                container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/vsearch:2.21.1--hf1761c0_1" : "quay.io/biocontainers/vsearch:2.21.1--hf1761c0_1")

                input:
                    tuple nid, file(notus) from nuclFastas_forCounts_ncasv_ch
                    file(merged) from nuclCounts_mergedreads_ncasv_ch

                output:
                    tuple file("*_counts.csv"), file("*_counts.biome") into counts_vsearch_ncasv
                    tuple nid, file("*ncASV*counts.csv") into notu_counts_plots

                script:
                    mtag="ID=" + nid
                    """
                    name=\$( echo ${notus} | awk -F ".fasta" '{print \$1}')
                    vsearch --usearch_global ${merged} --db ${notus} --id .${nid} --threads ${task.cpus} --otutabout \${name}_counts.txt --biomout \${name}_counts.biome
                    cat \${name}_counts.txt | tr "\t" "," >\${name}_count.csv
                    sed 's/#OTU ID/OTU_ID/g' \${name}_count.csv >\${name}_counts.csv
                    rm \${name}_count.csv
                    """
            }

            process Generate_ncASV_Matrices {
                label 'low_cpus'

                tag "${mtag}"

                publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ncASV/Matrices", mode: "copy", overwrite: true, pattern: '*ncASV*PercentID.matrix'

                conda (params.condaActivate ? "-c bioconda -c conda-forge clustalo=1.2.4=h87f3376_5" : null)

                container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/clustalo:1.2.4--h87f3376_5" : "quay.io/biocontainers/clustalo:1.2.4--h87f3376_5")

                input:
                    tuple nid, file(asvs) from nuclFastas_forMatrix_ncasv_ch

                output:
                    file("*.matrix") into clustmatrices_ncasv
                    tuple nid, file("*ncASV*PercentID.matrix") into notu_heatmap

                script:
                    mtag="ID=" + nid
                    """
                    name=\$( echo ${asvs}| awk -F ".fasta" '{print \$1}')
                    clustalo -i ${asvs} --distmat-out=\${name}_PairwiseDistance.matrix --full --force --threads=${task.cpus}
                    clustalo -i ${asvs} --distmat-out=\${name}_PercentIDq.matrix --percent-id --full --force --threads=${task.cpus}
                    cat \${name}_PercentIDq.matrix | tr " " ","| sed 's/,,/,/g' | grep "," >\${name}_PercentID.matrix
                    rm \${name}_PercentIDq.matrix
                    """
                }

                if (!params.skipPhylogeny) {

                    process ncASV_Phylogeny_step1 {

                          label 'norm_cpus'

                          tag "${mtag}"

                          publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ncASV/Phylogeny/Alignment", mode: "copy", overwrite: true

                          conda (params.condaActivate ? "-c bioconda -c conda-forge muscle=5.1" : null)

                          container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/muscle:5.1--h9f5acd7" : "quay.io/biocontainers/muscle:5.1--h9f5acd7")

                          input:
                              tuple nid, file(asvs) from nuclFastas_forphylogeny_ncasv

                          output:
                              tuple nid, file("*_ALN.fasta") into ncalign1
                              tuple nid, file("*.efa")

                          script:
                              mtag="ID=" + nid
                              """
                              pre=\$(echo ${asvs} | awk -F ".fasta" '{print \$1}' )

                              if [[ ${params.srep} == "true" && ${params.ensemble} == "false" ]];
                              then  if [[ \$( grep -c ">" ${asvs}) -lt 300 ]]
                                    then    comm="align"
                                    else    comm="super5"
                                    fi
                                    muscle -"\$comm" ${asvs} -perm ${params.perm} -perturb ${params.pert} -output \${pre}_muscle_raw_ALN.fasta -threads ${task.cpus} -quiet
                                    echo "single replicate alignment chosen; look over muscle5 documentation to learn about ensemble alignment approach" >note.efa
                              elif [[ ${params.srep} == "false" && ${params.ensemble} == "true" ]];
                              then  muscle -align ${asvs} -${params.fied} -output \${pre}_muscle.efa -threads ${task.cpus} -quiet
                                    muscle -maxcc \${pre}_muscle.efa -output \${pre}_muscle_raw_ALN.fasta
                              else  if [[ \$( grep -c ">" ${asvs}) -lt 300 ]]
                                      then    comm="align"
                                      else    comm="super5"
                                      fi
                                      muscle -"\$comm" ${asvs} -output \${pre}_muscle_raw_ALN.fasta -threads ${task.cpus} -quiet
                                      echo "single replicate alignment chosen; look over muscle5 documentation to learn about ensemble alignment approach" >note.efa
                              fi
                              """
                    }

                    process ncASV_Phylogeny_step2 {

                        label 'norm_cpus'

                        tag "${mtag}"

                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ncASV/Phylogeny/Alignment", mode: "copy", overwrite: true,  pattern: '*ncASV*aln.*'

                        conda (params.condaActivate ? "-c conda-forge bioconda::trimal=1.4.1=h9f5acd7_6" : null)

                        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/trimal:1.4.1--h9f5acd7_6" : "quay.io/biocontainers/trimal:1.4.1--h9f5acd7_6")

                        input:
                            tuple nid, file(asvs) from ncalign1

                        output:
                            tuple nid, file("*_aln.html") into align_results_ncasv2
                            tuple nid, file("*_aln.fasta") into ncalign2

                        script:
                            mtag="ID=" + nid
                            """
                            pre=\$(echo ${asvs} | awk -F "_muscle" '{print \$1}' )
                            trimal -in ${asvs} -out \${pre}_trimal_aln.fasta -keepheader -fasta -automated1 -htmlout \${pre}_aln.html
                            """
                    }

                    process ncASV_Phylogeny_step3 {

                          label 'norm_cpus'

                          tag "${mtag}"

                          publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ncASV/Phylogeny/Alignment", mode: "copy", overwrite: true

                          conda (params.condaActivate ? "${params.vampdir}/bin/yamls/oligotyping.yml" : null)

                          container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/oligotyping:2.1--py27_0" : "quay.io/biocontainers/oligotyping:2.1--py27_0")

                          input:
                              tuple nid, file(asvs) from ncalign2

                          output:
                              tuple nid, file("*_Aligned_informativeonly.fasta") into ncalign3

                          script:
                              mtag="ID=" + nid
                              """
                              pre=\$(echo ${asvs} | awk -F "_trimal" '{print \$1}' )
                              o-trim-uninformative-columns-from-alignment ${asvs}
                              mv ${asvs}-TRIMMED ./\${pre}_Aligned_informativeonly.fasta
                              """
                    }

                    process ncASV_Phylogeny_step4 {

                        label 'norm_cpus'

                        tag "${mtag}"

                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ncASV/Phylogeny/ModelTest", mode: "copy", overwrite: true, pattern: '*ncASV*mt*'

                        conda (params.condaActivate ? "-c conda-forge bioconda::modeltest-ng=0.1.7=h5c6ebe3_0" : null)

                        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/modeltest-ng:0.1.7--h5c6ebe3_0" : "quay.io/biocontainers/modeltest-ng:0.1.7--h5c6ebe3_0")

                        input:
                            tuple nid, file(asvs) from ncalign3

                        output:
                            tuple nid, file("*mt*") into align_results_ncasv4
                            tuple nid, file("*mt.out") into ncalign4

                        script:
                            mtag="ID=" + nid
                            """
                            pre=\$(echo ${asvs} | awk -F "_Aligned" '{print \$1}' )
                            # Nucleotide_ModelTest
                            modeltest-ng -i ${asvs} -p ${task.cpus} -o \${pre}_mt -d nt -s 203 --disable-checkpoint
                            """
                    }

                    process ncASV_Phylogeny_step5 {

                          label 'norm_cpus'

                          tag "${mtag}"

                          publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ncASV/Phylogeny/IQ-TREE", mode: "copy", overwrite: true, pattern: '*ncASV*iq*'

                          conda (params.condaActivate ? "-c conda-forge bioconda::iqtree=2.2.0.3=hb97b32f_1" : null)

                          container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/iqtree:2.2.0.3--hb97b32f_1" : "quay.io/biocontainers/iqtree:2.2.0.3--hb97b32f_1")

                          input:
                              tuple nid, file(asvs) from nuclFastas_forphylogeny_ncasv2
                              tuple nid, file(mtout) from ncalign4

                          output:
                              tuple nid, file("*iq*") into align_results_ncasv5
                              tuple nid, file("*iq.treefile") into nucl_phyl_plot_ncasv

                          script:
                              mtag="ID=" + nid
                              """
                              #grabbing best models from modeltestng
                              modbic=\$(grep "iqtree" ${mtout} | head -1 | awk -F "-m " '{print \$2}')
                              modaic=\$(grep "iqtree" ${mtout} | head -2 | tail -1 | awk -F "-m " '{print \$2}')
                              modaicc=\$(grep "iqtree" ${mtout} | tail -1 | awk -F "-m " '{print \$2}')
                              if [[ ${params.crit} == "BIC" ]]
                              then  mod="\$modbic"
                              elif  [[ ${params.crit} == "AIC" ]]
                              then  mod="\$modaic"
                              elif  [[ ${params.crit} == "AICc" ]]
                              then  mod="\$modaicc"
                              fi
                              # grab prefix
                              pre=\$(echo ${asvs} | awk -F "_aligned" '{print \$1}' )
                              # Nucleotide_Phylogeny
                              if [ "${params.iqCustomnt}" != "" ];then
                                  iqtree -s ${asvs} --prefix \${pre}_iq --redo -T auto ${params.iqCustomnt}
                              elif [[ "${params.ModelTnt}" != "false" && "${params.nonparametric}" != "false" ]];then
                                  iqtree -s ${asvs} --prefix \${pre}_iq -m \${mod} --redo -nt auto -b ${params.boots}
                              elif [[ "${params.ModelTnt}" != "false" && "${params.parametric}" != "false" ]];then
                                  iqtree -s ${asvs} --prefix \${pre}_iq -m \${mod} --redo -nt auto -bb ${params.boots} -bnni
                              elif [ "${params.nonparametric}" != "false" ];then
                                  iqtree -s ${asvs} --prefix \${pre}_iq -m MFP -madd --redo -nt auto -b ${params.boots}
                              elif [ "${params.parametric}" != "false" ];then
                                  iqtree -s ${asvs} --prefix \${pre}_iq -m MFP -madd --redo -nt auto -bb ${params.boots} -bnni
                              else
                                  iqtree -s ${asvs} --prefix \${pre}_iq -m MFP -madd --redo -nt auto -bb ${params.boots} -bnni
                              fi
                              """
                    }

                } else {

                    process skipncASVphylogeny {
                        input:
                            tuple nid, file(asvs) from nuclFastas_forphylogeny_ncasv

                        output:
                            tuple nid, file("skipncASVphylogeny.txt") into nucl_phyl_plot_ncasv

                        script:
                            """
                            echo "Skipped" >skipncASVphylogeny.txt
                            """
                    }
            }

        } else {
            reads_vsearch5_ch
    	       .into{ nuclFastas_forDiamond_asv_ch; nuclFastas_forCounts_asv_ch; nuclFastas_forphylogeny_asv; nuclFastas_forMatrix_asv_ch}
        }

        if (!params.skipTaxonomy) {

          if (params.dbtype == "NCBI") {

                process ASV_Taxonomy_Inference_NCBI { // edit //CHECK rename_seq.py in container !!!!!!!!!!!

                    label 'high_cpus'

                    publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ASVs/Taxonomy", mode: "copy", overwrite: true, pattern: '*_ASV*.{fasta,csv,tsv}'
                    publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ASVs/Taxonomy/DiamondOutput", mode: "copy", overwrite: true, pattern: '*_ASV*dmd.out'

                    conda (params.condaActivate ? "-c conda-forge bioconda::diamond=2.0.15=hb97b32f_1" : null)

                    container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/diamond:2.0.15--hb97b32f_1" : "quay.io/biocontainers/diamond:2.0.15--hb97b32f_1")

                    input:
                        file(asvs) from nuclFastas_forDiamond_asv_ch

                    output:
                        file("*.fasta") into tax_labeled_fasta_asv
                        tuple file("*_phyloformat.csv"), file("*summaryTable.tsv"), file("*dmd.out") into summary_diamond_asv
                        file("*_ASV*_summary_for_plot.csv") into taxplot_asv
                        file("*_quick_Taxbreakdown.csv") into tax_table_asv
                        file ("*_quicker_taxbreakdown.csv") into tax_nodCol_asv

                    script:
                        """
                        cp ${params.vampdir}/bin/rename_seq.py .
                        virdb=${params.dbdir}/${params.dbname}
                        if [[ ${params.measurement} == "bitscore" ]]
                        then    measure="--min-score ${params.bitscore}"
                        elif    [[ ${params.measurement} == "evalue" ]]
                        then    measure="-e ${params.evalue}"
                        else    measure="--min-score ${params.bitscore}"
                        fi
                        grep ">" \${virdb} > headers.list
                        headers="headers.list"
                        name=\$( echo ${asvs} | awk -F ".fasta" '{print \$1}')
                        if [[ ${params.ncbitax} == "true" ]]
                        then   diamond blastx -q ${asvs} -d \${virdb} -p ${task.cpus} --id ${params.minID} -l ${params.minaln} \${measure} --${params.sensitivity} -o "\$name"_dmd.out -f 6 qseqid qlen sseqid qstart qend qseq sseq length qframe evalue bitscore pident btop staxids sskingdoms skingdoms sphylums --max-target-seqs 1 --max-hsps 1
                        else   diamond blastx -q ${asvs} -d \${virdb} -p ${task.cpus} --id ${params.minID} -l ${params.minaln} \${measure} --${params.sensitivity} -o "\$name"_dmd.out -f 6 qseqid qlen sseqid qstart qend qseq sseq length qframe evalue bitscore pident btop --max-target-seqs 1 --max-hsps 1
                        fi
                        echo "Preparing lists to generate summary .csv's"
                        echo "[Best hit accession number]" > access.list
                        echo "[e-value]" > evalue.list
                        echo "[Bitscore]" > bit.list
                        echo "[Percent ID (aa)]" > pid.list
                        echo "[Organism ID]" > "\$name"_virus.list
                        echo "[Gene]" > "\$name"_genes.list
                        echo "[ASV#]" > otu.list
                        echo "[Sequence length]" > length.list
                        grep ">" ${asvs} | awk -F ">" '{print \$2}' > seqids.lst
                        if [[ ${params.lca} == "T" ]]
                        then  grep -w "LCA" ${params.dbanno}/*.txt > lcainfo.list
                              echo "[Taxonomic classification from RVDB annotations]" > lca_classification.list
                        else
                              echo "[Taxonomic classification from RVDB annotations]" > lca_classification.list
                        fi
                        if [[ ${params.ncbitax} == "true" ]]
                        then echo "[NCBI Taxonomy ID],[Taxonomic classification from NCBI]" > ncbi_classification.list
                        fi
                        echo "extracting genes and names"
                        touch new_"\$name"_asvnames.txt
                        for s in \$(cat seqids.lst);do
                            echo "Checking for \$s hit in diamond output"
                            if [[ "\$(grep -wc "\$s" "\$name"_dmd.out)" -eq 1 ]];then
                                        echo "Yep, there was a hit for \$s"
                                        echo "Extracting the information now:"
                                        acc=\$(grep -w "\$s" "\$name"_dmd.out | awk '{print \$3}')
                                        echo "\$s" >> otu.list
                                        echo "\$acc" >> access.list
                                        line="\$(grep -w "\$s" "\$name"_dmd.out)"
                                        echo "\$line" | awk '{print \$10}' >> evalue.list
                                        echo "\$line" | awk '{print \$11}' >> bit.list
                                        echo "\$line" | awk '{print \$12}' >> pid.list
                                        echo "\$line" | awk '{print \$2}' >> length.list
                                        echo "Extracting virus and gene ID for \$s now"
                                        gene=\$(grep -w "\$acc" "\$headers" | awk -F "." '{ print \$2 }' | awk -F "[" '{ print \$1 }' | awk -F " " '{print substr(\$0, index(\$0,\$2))}' | sed 's/ /_/g') &&
                                        echo "\$gene" | sed 's/_/ /g' >> "\$name"_genes.list
                                        virus=\$(grep -w "\$acc" "\$headers" | awk -F "[" '{ print \$2 }' | awk -F "]" '{ print \$1 }'| sed 's/ /_/g')
                                        echo "\$virus" | sed 's/_/ /g' >> "\$name"_virus.list
                                        echo ">"\${s}"_"\$virus"_"\$gene"" >> new_"\$name"_asvnames.txt
                                        if [[ "${params.lca}" == "T" ]]
                                        then    if [[ \$(grep -w "\$acc" ${params.dbanno}/*.txt | wc -l) -eq 1 ]]
                                                then  group=\$(grep -w "\$acc" ${params.dbanno}/*.txt | awk -F ":" '{print \$1}')
                                                      lcla=\$(grep -w "\$group" lcainfo.list | awk -F "\t" '{print \$2}')
                                                      echo "\$lcla" >> lca_classification.list
                                                else  echo "Viruses" >> lca_classification.list
                                                fi
                                        fi
                                        if [[ ${params.ncbitax} == "true" ]]
                                        then  echo "\$line" | awk -F "\t" '{print \$14","\$16"::"\$18"::"\$17}' >> ncbi_classification.list
                                        fi
                                        echo "\$s done."
                            else
                                        echo "Ugh, there was no hit for \$s .."
                                        echo "We still love \$s though and we will add it to the final fasta file"
                                        echo "\$s" >> otu.list
                                        echo "NO_HIT" >> access.list
                                        echo "NO_HIT" >> "\$name"_genes.list
                                        echo "NO_HIT" >> "\$name"_virus.list
                                        echo "NO_HIT" >> evalue.list
                                        echo "NO_HIT" >> bit.list
                                        echo "NO_HIT" >> pid.list
                                        echo "NO_HIT" >> length.list
                                        virus="NO"
                                        gene="HIT"
                                        echo ">\${s}_"\$virus"_"\$gene"" >> new_"\$name"_asvnames.txt
                                        if [[ "${params.lca}" == "T" ]]
                                        then    echo "N/A" >> lca_classification.list
                                        fi
                                        if [[ "${params.ncbitax}" == "true" ]]
                                        then  echo "N/A" >> ncbi_classification.list
                                        fi
                                        echo "\$s done."
                           fi
                        done
                        echo "Now editing "\$name" fasta headers"
                        ###### rename_seq.py
                        ./rename_seq.py ${asvs} new_"\$name"_asvnames.txt "\$name"_TaxonomyLabels.fasta
                        awk 'BEGIN {RS=">";FS="\\n";OFS=""} NR>1 {print ">"\$1; \$1=""; print}' "\$name"_TaxonomyLabels.fasta >"\$name"_tmpssasv.fasta
                        echo "[Sequence header]" > newnames.list
                        cat new_"\$name"_asvnames.txt >> newnames.list
                        touch sequence.list
                        echo "     " > sequence.list
                        grep -v ">" "\$name"_tmpssasv.fasta >> sequence.list
                        rm "\$name"_tmpssasv.fasta
                        if [[ "${params.lca}" == "T" && "${params.ncbitax}" == "true" ]]
                        then
                              paste -d "," sequence.list "\$name"_virus.list "\$name"_genes.list otu.list lca_classification.list ncbi_classification.list newnames.list length.list bit.list evalue.list pid.list access.list >> "\$name"_phyloformat.csv
                              paste -d"\t" otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list ncbi_classification.list sequence.list length.list bit.list evalue.list pid.list >> "\$name"_summaryTable.tsv
                              paste -d"," otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list ncbi_classification.list >> \${name}_quick_Taxbreakdown.csv
                        elif [[ "${params.lca}" == "T" && "${params.ncbitax}" != "true" ]]
                        then
                              paste -d "," sequence.list "\$name"_virus.list "\$name"_genes.list otu.list lca_classification.list newnames.list length.list bit.list evalue.list pid.list access.list >> "\$name"_phyloformat.csv
                              paste -d"\t" otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list sequence.list length.list bit.list evalue.list pid.list >> "\$name"_summaryTable.tsv
                              paste -d"," otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list >> \${name}_quick_Taxbreakdown.csv
                        elif [[ "${params.ncbitax}" == "true" && "${params.lca}" != "T" ]]
                        then
                              paste -d "," sequence.list "\$name"_virus.list "\$name"_genes.list otu.list ncbi_classification.list newnames.list length.list bit.list evalue.list pid.list access.list >> "\$name"_phyloformat.csv
                              paste -d"\t" otu.list access.list "\$name"_virus.list "\$name"_genes.list ncbi_classification.list sequence.list length.list bit.list evalue.list pid.list >> "\$name"_summaryTable.tsv
                              paste -d"," otu.list access.list "\$name"_virus.list "\$name"_genes.list ncbi_classification.list >> \${name}_quick_Taxbreakdown.csv
                        else
                              paste -d "," sequence.list "\$name"_virus.list "\$name"_genes.list otu.list newnames.list length.list bit.list evalue.list pid.list access.list >> "\$name"_phyloformat.csv
                              paste -d"\t" otu.list access.list "\$name"_virus.list "\$name"_genes.list sequence.list length.list bit.list evalue.list pid.list >> "\$name"_summaryTable.tsv
                              echo "skipped" >> \${name}_quick_Taxbreakdown.csv
                        fi
                        for x in *phyloformat.csv;do
                            echo "\$x"
                            lin=\$(( \$(wc -l \$x | awk '{print \$1}')-1))
                            tail -"\$lin" \$x | awk -F "," '{print \$2}' > tmpcol.list;
                            sed 's/ /_/g' tmpcol.list > tmp2col.list;
                            cat tmp2col.list | sort | uniq -c | sort -nr | awk '{print \$2","\$1}' > \${name}_summary_for_plot.csv;
                            rm tmpcol.list tmp2col.list
                        done
                        awk -F "," '{print \$1","\$3"("\$2")"}' \${name}_quick_Taxbreakdown.csv >> \${name}_quicker_taxbreakdown.csv
                        rm evalue.list sequence.list bit.list pid.list length.list seqids.lst otu.list *asvnames.txt "\$name"_virus.list "\$name"_genes.list newnames.list access.list headers.list
                        """
                      }
                } else if (params.dbtype== "RVDB") {

                  process ASV_Taxonomy_Inference_RVDB { // edit //CHECK rename_seq.py in container !!!!!!!!!!!

                      label 'high_cpus'

                      publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ASVs/Taxonomy", mode: "copy", overwrite: true, pattern: '*_ASV*.{fasta,csv,tsv}'
                      publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ASVs/Taxonomy/DiamondOutput", mode: "copy", overwrite: true, pattern: '*_ASV*dmd.out'

                      conda (params.condaActivate ? "-c conda-forge bioconda::diamond=2.0.15=hb97b32f_1" : null)

                      container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/diamond:2.0.15--hb97b32f_1" : "quay.io/biocontainers/diamond:2.0.15--hb97b32f_1")

                      input:
                          file(asvs) from nuclFastas_forDiamond_asv_ch

                      output:
                          file("*.fasta") into tax_labeled_fasta_asv
                          tuple file("*_phyloformat.csv"), file("*summaryTable.tsv"), file("*dmd.out") into summary_diamond_asv
                          file("*_ASV*_summary_for_plot.csv") into taxplot_asv
                          file("*_quick_Taxbreakdown.csv") into tax_table_asv
                          file ("*_quicker_taxbreakdown.csv") into tax_nodCol_asv

                      script:
                          """
                          cp ${params.vampdir}/bin/rename_seq.py .
                          virdb=${params.dbdir}/${params.dbname}
                          grep ">" \${virdb} > headers.list
                          if [[ ${params.measurement} == "bitscore" ]]
                          then    measure="--min-score ${params.bitscore}"
                          elif    [[ ${params.measurement} == "evalue" ]]
                          then    measure="-e ${params.evalue}"
                          else    measure="--min-score ${params.bitscore}"
                          fi
                          headers="headers.list"
                          name=\$( echo ${asvs} | awk -F ".fasta" '{print \$1}')
                          diamond blastx -q ${asvs} -d \${virdb} -p ${task.cpus} --id ${params.minID} -l ${params.minaln} \${measure} --${params.sensitivity} -o "\$name"_dmd.out -f 6 qseqid qlen sseqid qstart qend qseq sseq length qframe evalue bitscore pident btop --max-target-seqs 1 --max-hsps 1
                          echo "Preparing lists to generate summary .csv's"
                          echo "[Best hit accession number]" > access.list
                          echo "[e-value]" > evalue.list
                          echo "[Bitscore]" > bit.list
                          echo "[Percent ID (aa)]" > pid.list
                          echo "[Organism ID]" > "\$name"_virus.list
                          echo "[Gene]" > "\$name"_genes.list
                          echo "[ASV#]" > otu.list
                          echo "[Sequence length]" > length.list
                          grep ">" ${asvs} | awk -F ">" '{print \$2}' > seqids.lst
                          if [[ ${params.lca} == "T" ]]
                          then  grep -w "LCA" ${params.dbanno}/*.txt > lcainfo.list
                                echo "[Taxonomic classification from RVDB annotations]" > lca_classification.list
                          else  echo "skipped" >> \${name}_quick_Taxbreakdown.csv
                                echo "[Taxonomic classification from RVDB annotations]" > lca_classification.list
                          fi
                          echo "extracting genes and names"
                          touch new_"\$name"_asvnames.txt
                          for s in \$(cat seqids.lst);do
                              echo "Using RVDB headers."
                              if [[ "\$(grep -wc "\$s" "\$name"_dmd.out)" -eq 1 ]];then
                                  echo "Yep, there was a hit for \$s"
                                  echo "Extracting the information now:"
                                  acc=\$(grep -w "\$s" "\$name"_dmd.out | awk '{print \$3}' | awk -F "|" '{print \$3}')
                                  echo "\$s" >> otu.list
                                  echo "\$acc" >> access.list
                                  line="\$(grep -w "\$s" "\$name"_dmd.out)"
                                  echo "\$line" | awk '{print \$10}' >> evalue.list
                                  echo "\$line" | awk '{print \$11}' >> bit.list
                                  echo "\$line" | awk '{print \$12}' >> pid.list
                                  echo "\$line" | awk '{print \$2}' >> length.list
                                  echo "Extracting virus and gene ID for \$s now"
                                  gene=\$(grep -w "\$acc" "\$headers" | awk -F "|" '{ print \$6 }' | awk -F "[" '{ print \$1 }' | sed 's/ /_/g') &&
                                  echo "\$gene" | sed 's/_/ /g' >> "\$name"_genes.list
                                  virus=\$(grep -w "\$acc" "\$headers" | awk -F "|" '{ print \$6 }' | awk -F "[" '{ print \$2 }' | awk -F "]" '{print \$1}' | sed 's/ /_/g') &&
                                  echo "\$virus" | sed 's/_/ /g' >> "\$name"_virus.list
                                  echo ">\${s}_"\$virus"_"\$gene"" >> new_"\$name"_asvnames.txt
                                  if [[ "${params.lca}" == "T" ]]
                                  then    if [[ \$(grep -w "\$acc" ${params.dbanno}/*.txt | wc -l) -eq 1 ]]
                                          then  group=\$(grep -w "\$acc" ${params.dbanno}/*.txt | awk -F ":" '{print \$1}')
                                                lcla=\$(grep -w "\$group" lcainfo.list | awk -F "\t" '{print \$2}')
                                                echo "\$lcla" >> lca_classification.list
                                          else  echo "Viruses" >> lca_classification.list
                                          fi
                                  fi
                                  echo "\$s done."
                              else
                                  echo "Ugh, there was no hit for \$s .."
                                  echo "We still love \$s though and we will add it to the final fasta file"
                                  echo "\$s" >> otu.list
                                  echo "NO_HIT" >> access.list
                                  echo "NO_HIT" >> "\$name"_genes.list
                                  echo "NO_HIT" >> "\$name"_virus.list
                                  echo "NO_HIT" >> evalue.list
                                  echo "NO_HIT" >> bit.list
                                  echo "NO_HIT" >> pid.list
                                  echo "NO_HIT" >> length.list
                                  virus="NO"
                                  gene="HIT"
                                  echo ">\${s}_"\$virus"_"\$gene"" >> new_"\$name"_asvnames.txt
                                  if [[ "${params.lca}" == "T" ]]
                                  then    echo "N/A" >> lca_classification.list
                                  fi
                                  echo "\$s done."
                              fi
                          echo "Done with \$s"
                          done
                          echo "Now editing "\$name" fasta headers"
                          ###### rename_seq.py
                          ./rename_seq.py ${asvs} new_"\$name"_asvnames.txt "\$name"_TaxonomyLabels.fasta
                          awk 'BEGIN {RS=">";FS="\\n";OFS=""} NR>1 {print ">"\$1; \$1=""; print}' "\$name"_TaxonomyLabels.fasta >"\$name"_tmpssasv.fasta
                          echo "[Sequence header]" > newnames.list
                          cat new_"\$name"_asvnames.txt >> newnames.list
                          touch sequence.list
                          echo "     " > sequence.list
                          grep -v ">" "\$name"_tmpssasv.fasta >> sequence.list
                          rm "\$name"_tmpssasv.fasta
                          if [[ "${params.lca}" == "T" ]]
                          then  paste -d "," sequence.list "\$name"_virus.list "\$name"_genes.list otu.list lca_classification.list newnames.list length.list bit.list evalue.list pid.list access.list >> "\$name"_phyloformat.csv
                                paste -d"\t" otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list sequence.list length.list bit.list evalue.list pid.list >> "\$name"_summaryTable.tsv
                                paste -d"," otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list >> \${name}_quick_Taxbreakdown.csv
                          else  paste -d "," sequence.list "\$name"_virus.list "\$name"_genes.list otu.list newnames.list length.list bit.list evalue.list pid.list access.list >> "\$name"_phyloformat.csv
                                paste -d"\t" otu.list access.list "\$name"_virus.list "\$name"_genes.list sequence.list length.list bit.list evalue.list pid.list >> "\$name"_summaryTable.tsv
                          fi
                          for x in *phyloformat.csv;do
                                    echo "\$x"
                                    lin=\$(( \$(wc -l \$x | awk '{print \$1}')-1))
                                    tail -"\$lin" \$x | awk -F "," '{print \$2}' > tmpcol.list;
                                    sed 's/ /_/g' tmpcol.list > tmp2col.list;
                                    cat tmp2col.list | sort | uniq -c | sort -nr | awk '{print \$2","\$1}' > \${name}_summary_for_plot.csv;
                                    rm tmpcol.list tmp2col.list
                          done
                          awk -F "," '{print \$1","\$3"("\$2")"}' \${name}_quick_Taxbreakdown.csv >> \${name}_quicker_taxbreakdown.csv
                          rm evalue.list sequence.list bit.list pid.list length.list seqids.lst otu.list *asvnames.txt "\$name"_virus.list "\$name"_genes.list newnames.list access.list headers.list
                          """
                }
            }
        } else {
            taxplot_asv = Channel.value('skipping')
            tax_table_asv = Channel.value('skipping')
            tax_nodCol_asv = Channel.value('skipping')
       }

       process Generate_ASV_Counts_Tables {

            label 'norm_cpus'

            publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ASVs/Counts", mode: "copy", overwrite: true, pattern: '*ASV*.{biome,csv}'

            conda (params.condaActivate ? "-c bioconda -c conda-forge vsearch=2.21.1=hf1761c0_1" : null)

            container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/vsearch:2.21.1--hf1761c0_1" : "quay.io/biocontainers/vsearch:2.21.1--hf1761c0_1")

            input:
                file(asvs) from nuclFastas_forCounts_asv_ch
                file(merged) from nuclCounts_mergedreads_asv_ch

            output:
                tuple file("*_counts.csv"), file("*_counts.biome") into counts_vsearch_asv
                file("*_ASV*counts.csv") into (asv_counts_plots, asvcount_med, asvcount_phylogr, asv_deseq)

            script:
                """
                name=\$( echo ${asvs} | awk -F ".fasta" '{print \$1}' | sed 's/ASVs/ASV/g')
                vsearch --usearch_global ${merged} --db ${asvs} --id ${params.asvcountID} --threads ${task.cpus} --otutabout "\$name"_counts.txt --biomout "\$name"_counts.biome
                cat \${name}_counts.txt | tr "\t" "," >\${name}_count.csv
                sed 's/#OTU ID/OTU_ID/g' \${name}_count.csv >\${name}_counts.csv
                rm \${name}_count.csv
                """
        }

        process Generate_ASV_Matrices {

            label 'low_cpus'

            publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ASVs/Matrices", mode: "copy", overwrite: true, pattern: '*ASV*PercentID.matrix'

            conda (params.condaActivate ? "-c bioconda -c conda-forge clustalo=1.2.4=h87f3376_5" : null)

            container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/clustalo:1.2.4--h87f3376_5" : "quay.io/biocontainers/clustalo:1.2.4--h87f3376_5")

            input:
                file(reads) from nuclFastas_forMatrix_asv_ch

            output:
                file("*.matrix") into clustmatrices_asv
                file("*_ASV*PercentID.matrix") into asv_heatmap

            script:
                // remove if statement later (no fin)
                """
                for filename in ${reads};do
                    if [ `echo \${filename} | grep -c "ncASV"` -eq 1 ];then
                        ident=\$( echo \${filename} | awk -F "ncASV" '{print \$2}' | awk -F ".fasta" '{print \$1}')
                        name=\$( echo \${filename}| awk -F ".fasta" '{print \$1}')
                        clustalo -i \${filename} --distmat-out=\${name}_PairwiseDistance.matrix --full --force --threads=${task.cpus}
                        clustalo -i \${filename} --distmat-out=\${name}_PercentIDq.matrix --percent-id --full --force --threads=${task.cpus}
                        for x in *q.matrix;do
                            pre=\$(echo "\$x" | awk -F "q.matrix" '{print \$1}')
                            ya=\$(wc -l \$x | awk '{print \$1}')
                            echo "\$((\$ya-1))"
                            tail -"\$((\$ya-1))" \$x > \${pre}z.matrix
                            rm \$x
                            cat \${pre}z.matrix | sed 's/ /,/g' | sed -E 's/(,*),/,/g' >\${pre}.matrix
                            rm \${pre}z.matrix
                        done
                    else
                        name=\$( echo \${filename} | awk -F ".fasta" '{print \$1}')
                        clustalo -i \${filename} --distmat-out=\${name}_PairwiseDistance.matrix --full --force --threads=${task.cpus}
                        clustalo -i \${filename} --distmat-out=\${name}_PercentIDq.matrix --percent-id --full --force --threads=${task.cpus}
                        for x in *q.matrix;do
                            pre=\$(echo "\$x" | awk -F "q.matrix" '{print \$1}')
                            ya=\$(wc -l \$x | awk '{print \$1}')
                            echo "\$((\$ya-1))"
                            tail -"\$((\$ya-1))" \$x > \${pre}z.matrix
                            rm \$x
                            cat \${pre}z.matrix | sed 's/ /,/g' | sed -E 's/(,*),/,/g' >\${pre}.matrix
                            rm \${pre}z.matrix
                        done
                    fi
                done
                """
            }

            if (!params.skipPhylogeny || params.asvMED || params.asvTClust) {

                process ASV_pre_Phylogeny_step1 {

                      label 'norm_cpus'

                      publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ASVs/Phylogeny/Alignment", mode: "copy", overwrite: true

                      conda (params.condaActivate ? "-c bioconda -c conda-forge muscle=5.1" : null)

                      container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/muscle:5.1--h9f5acd7" : "quay.io/biocontainers/muscle:5.1--h9f5acd7")

                      input:
                          file(asvs) from nuclFastas_forphylogeny_asv

                      output:
                          file("*_ALN.fasta") into asv_align1
                          file("*.efa")

                      script:
                          """
                          pre=\$(echo ${asvs} | awk -F ".fasta" '{print \$1}' )

                          if [[ ${params.srep} == "true" && ${params.ensemble} == "false" ]];
                          then  if [[ \$( grep -c ">" ${asvs}) -lt 300 ]]
                                then    comm="align"
                                else    comm="super5"
                                fi
                                muscle -"\$comm" ${asvs} -perm ${params.perm} -perturb ${params.pert} -output \${pre}_muscle_raw_ALN.fasta -threads ${task.cpus} -quiet
                                echo "single replicate alignment chosen; look over muscle5 documentation to learn about ensemble alignment approach" >note.efa
                          elif [[ ${params.srep} == "false" && ${params.ensemble} == "true" ]];
                          then  muscle -align ${asvs} -${params.fied} -output \${pre}_muscle.efa -threads ${task.cpus} -quiet

                          else  if [[ \$( grep -c ">" ${asvs}) -lt 300 ]]
                                then    comm="align"
                                else    comm="super5"
                                fi
                                muscle -"\$comm" ${asvs} -output \${pre}_muscle_raw_ALN.fasta -threads ${task.cpus} -quiet
                                echo "single replicate alignment chosen; look over muscle5 documentation to learn about ensemble alignment approach" >note.efa
                          fi
                          """
                  }

                  process ASV_pre_Phylogeny_step2 {

                        label 'norm_cpus'

                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ASVs/Phylogeny/Alignment", mode: "copy", overwrite: true

                        conda (params.condaActivate ? "-c bioconda -c conda-forge trimal=1.4.1=h9f5acd7_6" : null)

                        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/trimal:1.4.1--h9f5acd7_6" : "quay.io/biocontainers/trimal:1.4.1--h9f5acd7_6")

                        input:
                            file(align) from asv_align1

                        output:
                            file("*_aln.html") into trimalhtml2
                            file("*_aln.fasta") into asv_align2

                        script:
                            """
                            pre=\$(echo ${align} | awk -F "_muscle" '{print \$1}' )
                            trimal -in ${align} -out \${pre}_trimal_aln.fasta -keepheader -fasta -automated1 -htmlout \${pre}_trimal_aln.html
                            """
                    }

                    process ASV_pre_Phylogeny_step3 {

                          label 'norm_cpus'

                          publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ASVs/Phylogeny/Alignment", mode: "copy", overwrite: true

                          conda (params.condaActivate ? "${params.vampdir}/bin/yamls/oligotyping.yml" : null)

                          container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/oligotyping:2.1--py27_0" : "quay.io/biocontainers/oligotyping:2.1--py27_0")

                          input:
                              file(align) from asv_align2

                          output:
                              file("*.fasta") into asv_align3_mt, asv_align3_iq, asv_align3_med

                          script:
                              """
                              o-trim-uninformative-columns-from-alignment ${align}
                              pre=\$(echo ${align} | awk -F "_trimal" '{print \$1}' )
                              mv ${align}-TRIMMED ./\${pre}_Aligned_informativeonly.fasta
                              """
                      }

                      process ASV_pre_Phylogeny_step4 {

                            label 'norm_cpus'

                            publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ASVs/Phylogeny/ModelTest", mode: "copy", overwrite: true, pattern: '*mt*'

                            conda (params.condaActivate ? "-c conda-forge bioconda::modeltest-ng=0.1.7=h5c6ebe3_0" : null)

                            container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/modeltest-ng:0.1.7--h5c6ebe3_0" : "quay.io/biocontainers/modeltest-ng:0.1.7--h5c6ebe3_0")

                            input:
                                file(align) from asv_align3_mt

                            output:
                                file("*mt.out") into asv_align4
                                file("*mt.*") into mtres

                            script:
                                """
                                # Nucleotide_ModelTest
                                pre=\$(echo ${align} | awk -F "_Aligned_informativeonly.fasta" '{print \$1}' )
                                modeltest-ng -i ${align} -p ${task.cpus} -o \${pre}_mt -d nt -s 203 --disable-checkpoint
                                """
                        }

                    }

                    if (!params.skipPhylogeny || params.asvTClust) {

                        process ASV_Phylogeny_step5 {

                              label 'norm_cpus'

                              publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ASVs/Phylogeny/IQ-TREE", mode: "copy", overwrite: true, pattern: '*iq*'

                              conda (params.condaActivate ? "-c conda-forge bioconda::iqtree=2.2.0.3=hb97b32f_1" : null)

                              container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/iqtree:2.2.0.3--hb97b32f_1" : "quay.io/biocontainers/iqtree:2.2.0.3--hb97b32f_1")

                              input:
                                  file(align) from asv_align3_iq
                                  file(mtout) from asv_align4

                              output:
                                  file("*iq*") into asv_align5
                                  file("*iq.treefile") into (nucl_phyl_plot_asv, asv_treeclust)

                              script:
                                  """
                                  #grabbing best models from modeltestng
                                  modbic=\$(grep "iqtree" ${mtout} | head -1 | awk -F "-m " '{print \$2}')
                                  modaic=\$(grep "iqtree" ${mtout} | head -2 | tail -1 | awk -F "-m " '{print \$2}')
                                  modaicc=\$(grep "iqtree" ${mtout} | tail -1 | awk -F "-m " '{print \$2}')
                                  if [[ ${params.crit} == "BIC" ]]
                                  then  mod="\$modbic"
                                  elif  [[ ${params.crit} == "AIC" ]]
                                  then  mod="\$modaic"
                                  elif  [[ ${params.crit} == "AICc" ]]
                                  then  mod="\$modaicc"
                                  fi
                                  # grab prefix
                                  pre=\$(echo ${align} | awk -F "_trimal_Aligned_informativeonly.fasta" '{print \$1}' )
                                  # Nucleotide_Phylogeny
                                  if [ "${params.iqCustomnt}" != "" ];then
                                      iqtree -s ${align} --prefix \${pre}_iq --redo -T auto ${params.iqCustomnt}
                                  elif [[ "${params.ModelTnt}" != "false" && "${params.nonparametric}" != "false" ]];then
                                      iqtree -s ${align} --prefix \${pre}_iq -m \${mod} --redo -nt auto -b ${params.boots}
                                  elif [[ "${params.ModelTnt}" != "false" && "${params.parametric}" != "false" ]];then
                                      iqtree -s ${align} --prefix \${pre}_iq -m \${mod} --redo -nt auto -bb ${params.boots} -bnni
                                  elif [ "${params.nonparametric}" != "false" ];then
                                      iqtree -s ${align} --prefix \${pre}_iq -m MFP -madd --redo -nt auto -b ${params.boots}
                                  elif [ "${params.parametric}" != "false" ];then
                                      iqtree -s ${align} --prefix \${pre}_iq -m MFP -madd --redo -nt auto -bb ${params.boots} -bnni
                                  else
                                      iqtree -s ${align} --prefix \${pre}_iq -m MFP -madd --redo -nt auto -bb ${params.boots} -bnni
                                  fi
                                  """
                          }
            } else {
                nucl_phyl_plot_asv = Channel.value('skipping')
                asv_treeclust = Channel.value('skipping')
            }

            if (!params.skipPhylogeny || params.asvTClust) {

                process ASV_PhyloClustering {

                      label 'norm_cpus'

                      publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ASVs/TreeClustering", mode: "copy", overwrite: true

                      conda (params.condaActivate ? "${params.vampdir}/bin/yamls/treecluster.yml" : null)

                      container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/treecluster:1.0.3--pyh3252c3a_0" : "quay.io/biocontainers/treecluster:1.0.3--pyh3252c3a_0")

                      input:
                         file(tree) from asv_treeclust
                         file(counts) from asvcount_phylogr
                      output:
                          file("*treeclustering*.out") into asvtreeclustering_res
                          file("${params.projtag}_ASV_phylogroup.csv") into asv_phylogroupcsv
                          file("${params.projtag}_ASV_phyloGroupingcounts.csv") into (asv_phylogroupingcsv, asv_phylogroupingcsv2)

                      script:
                          """
                          TreeCluster.py -i ${tree} ${params.asvTCopp} >${params.projtag}_ASV_treeclustering.out
                          TreeCluster.py -i ${tree} ${params.asvTCopp} >${params.projtag}_ASV_treeclustering_verbose.out
                          #create headless treeclustering.out
                          tail -n +2 ${params.projtag}_ASV_treeclustering.out | sed 's/-1/0X/g' > headless.treeout
                          #extracting singletons
                          grep -w "0X" headless.treeout > single.out
                          awk -F "\\t" '{print \$1}' single.out > single.list
                          #assigning groupID to singletons
                          for x in \$(awk -F "\\t" '{print \$1}' single.out); do echo ""\$x"XX" >> single.groups;done

                          #summarizing clustering results
                          awk -F "\\t" '{print \$2}' headless.treeout | grep -v "0X" | sort | uniq > grp.list
                          cat single.groups >> group.list
                          cat grp.list >> group.list
                          echo "Sequence_ID,phyloGroup_ID" >${params.projtag}_ASV_phylogroup.csv

                          for x in \$(seq "\$(wc -l group.list | awk '{print \$1}')");
                          do      echo "phyloGroup"\$x"" >> grup.list
                          done
                          paste -d "," grup.list group.list > groups.csv
                          rm grup.list group.list
                          awk -F "," '{print \$1}' ${counts} | sed '1d' > asv.list
                          for z in \$(cat asv.list);
                          do      if [[ \$(grep -wc "\$z" single.list) -ge 1 ]]
                                  then
                                          group=\$(grep -w ""\$z"XX" groups.csv | awk -F "," '{print \$1}')
                                  else
                                          grp=\$(grep -w "\$z" headless.treeout | awk -F "\\t" '{print \$2}')
                                          group=\$(grep -w "\$grp" groups.csv | awk -F "," '{print \$1}')
                                  fi
                                  echo ""\$z","\$group"" >>${params.projtag}_ASV_phylogroup.csv
                          done
                          awk -F "," '{print \$2}' ${params.projtag}_ASV_phylogroup.csv > groups.list
                          paste -d"," groups.list ${counts} >${params.projtag}_ASV_phyloGroupingcounts.csv
                          """
                }

            } else {
                asv_phylogroupcsv = Channel.value('skipping')
                asv_phylogroupingcsv = Channel.value('skipping')
                asv_phylogroupingcsv2 = Channel.value('skipping')
            }

            if (params.asvMED) {

                process ASV_Minimum_Entropy_Decomposition_step1 {

                    label 'low_cpus'

                    publishDir "${params.workingdir}/${params.outdir}/Analyze/Clustering/ASVs/MED", mode: "copy", overwrite: true, pattern: '*.{fasta,csv}'
                    publishDir "${params.workingdir}/${params.outdir}/Analyze/Clustering/ASVs/MED/unique", mode: "copy", overwrite: true, pattern: '*_unique'

                    conda (params.condaActivate ? "${params.vampdir}/bin/yamls/oligotyping.yml" : null)

                    container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/oligotyping:2.1--py27_0" : "quay.io/biocontainers/oligotyping:2.1--py27_0")

                    input:
                      file(align) from asv_align3_med

                    output:
                      file("*_unique") into uniqformed
                      file("OLIGO-REPRESENTATIVES.fasta") into oligorep

                    script:
                        """
                        pre=\$(echo ${align} | awk -F "_Aligned" '{print \$1}' )
                        #entopy analysis
                        entropy-analysis ${align}
                        #Decomposition
                        if [[ \$(echo ${params.asvC} | grep -c ",") -ge 1 || "${params.asvSingle}" == "false" ]]
                        then
                              tag=\$(echo ${params.asvC} | sed 's/,/_/g')
                              oligotype ${align} \${pre}_Aligned_informativeonly.fasta-ENTROPY -o ${params.projtag}_asvMED_"\$tag" -M 1 -C ${params.asvC} -N ${task.cpus} --skip-check-input --no-figures --skip-gen-html
                        elif [[ "${params.asvSingle}" == "true" ]]
                        then
                              tag="${params.asvC}"
                              oligotype ${align} \${pre}_Aligned_informativeonly.fasta-ENTROPY -o ${params.projtag}_asvMED_"\$tag" -M 1 -C ${params.asvC} -N ${task.cpus} --skip-check-input --no-figures --skip-gen-html
                        else
                              oligotype ${align} \${pre}_Aligned_informativeonly.fasta-ENTROPY -o ${params.projtag}_asvMED_${params.asvC} -M 1 -c ${params.asvC} -N ${task.cpus} --skip-check-input --no-figures --skip-gen-html
                        fi
                        ###if statement makes sense now to me, need to continue working on the remaining parts -- the fasta files needed for the next process -- looks like this might work - 10/23/22
                        #generatemaps
                        mv ./${params.projtag}_asvMED_*/OLIGO-REPRESENTATIVES.fasta .
                        mv ./${params.projtag}_asvMED_*/OLIGO-REPRESENTATIVES/*_unique .
                        """
                }

                process ASV_Minimum_Entropy_Decomposition_step2 {

                    label 'low_cpus'

                    publishDir "${params.workingdir}/${params.outdir}/Analyze/Clustering/ASVs/MED", mode: "copy", overwrite: true

                    conda (params.condaActivate ? "bioconda::seqtk=1.3=h7132678_4" : null)

                    container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/seqtk:1.3--h7132678_4" : "quay.io/biocontainers/seqtk:1.3--h7132678_4")

                    input:
                      file(asvs) from asv_for_med
                      file(oligo) from oligorep
                    output:
                      file("*_ASV_Grouping.csv") into (asvgroupscsv, asvgroupscsv2)
                      file("${params.projtag}_ASV_group_reps_aligned.fasta") into (groupreps, groupreps2)

                    script:
                        """
                        #copy _unique files here so we can werk with them
                        cp ${params.workingdir}/${params.outdir}/Analyze/Clustering/ASVs/MED/unique .
                        #generatemaps
                        echo "ASV,GroupID,IDPattern"
                        j=1
                        for x in *unique;
                        do      gid=\$(echo \$x | awk -F "_" '{print \$1}')
                                uni=\$(echo \$x | awk -F ""\${gid}"_" '{print \$2}' | awk -F "_uni" '{print \$1}')
                                grep ">"  "\$gid"_"\$uni" | awk -F ">" '{print \$2}' > asv.list
                                seqtk subseq ${asvs} asv.list > Group"\${j}"_sequences.fasta
                                for z in \$( cat asv.list)
                                do      echo ""\$z",Group"\$j","\$uni"" >> ${params.projtag}_ASV_Grouping.csv

                                done
                                rm asv.list
                                echo ">Group\${j}" >> ${params.projtag}_ASV_group_reps_aligned.fasta
                                echo "\$uni" > group.list
                                seqtk subseq ${oligo} group.list > group.fasta
                                tail -1 group.fasta >> ${params.projtag}_ASV_group_reps_aligned.fasta
                                mv "\$gid"_"\$uni" ./Group"\$j"_"\$uni"_aligned.fasta
                                mv "\$gid"_"\$uni"_unique ./Group"\$j"_"\$uni"_unqiues_aligned.fasta
                                #rm "\$gid"*.cPickle
                                j=\$((\$j+1))
                        done
                        #mv ${params.projtag}_ASV_Grouping.csv ../../
                        #mv ${params.projtag}_ASV_group_reps_aligned.fasta ../../
                        #cd ..
                        """
              }

              process ASV_MED_Reps_ModelTesting {

                  label 'low_cpus'

                  publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ASVs/MED/Phylogeny/ModelTest", mode: "copy", overwrite: true, pattern: '*ASV*mt*'

                  conda (params.condaActivate ? "bioconda::modeltest-ng=0.1.7=h5c6ebe3_0" : null)

                  container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/modeltest-ng:0.1.7--h5c6ebe3_0" : "quay.io/biocontainers/modeltest-ng:0.1.7--h5c6ebe3_0")

                  input:
                    file(reps) from groupreps

                  output:
                    file("*mt.out") into asvmedrep_mtout
                    file("*mt.*") into mtres2

                  script:
                      """
                      modeltest-ng -i ${reps} -p ${task.cpus} -o ${params.projtag}_ASV_MEDGroup_Reps_mt -d nt -s 203 --disable-checkpoint
                      """
                  }


                process ASV_MED_Reps_Phylogeny {

                    label 'low_cpus'

                    publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ASVs/MED/Phylogeny/IQ-TREE", mode: "copy", overwrite: true, pattern: '*ASV*iq*'

                    conda (params.condaActivate ? "bioconda::iqtree=2.2.0.3=hb97b32f_1" : null)

                    container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/iqtree:2.2.0.3--hb97b32f_1" : "quay.io/biocontainers/iqtree:2.2.0.3--hb97b32f_1")

                    input:
                      file(reps) from groupreps2
                      file(mtout) from asvmedrep_mtout

                    output:
                      file("*_ASV_Group_Reps*") into align_results_asvmed2
                      file("*iq.treefile") into asv_group_rep_tree

                    script:
                        """
                        #grabbing best models from modeltestng
                        modbic=\$(grep "iqtree" ${mtout} | head -1 | awk -F "-m " '{print \$2}')
                        modaic=\$(grep "iqtree" ${mtout} | head -2 | tail -1 | awk -F "-m " '{print \$2}')
                        modaicc=\$(grep "iqtree" ${mtout} | tail -1 | awk -F "-m " '{print \$2}')
                        if [[ ${params.crit} == "BIC" ]]
                        then  mod="\$modbic"
                        elif  [[ ${params.crit} == "AIC" ]]
                        then  mod="\$modaic"
                        elif  [[ ${params.crit} == "AICc" ]]
                        then  mod="\$modaicc"
                        fi
                        # Protein_Phylogeny
                        if [ "${params.iqCustomaa}" != "" ];then
                            iqtree -s ${reps} --prefix ${params.projtag}_ASV_Group_Reps_iq --redo -T auto ${params.iqCustomaa}

                        elif [[ "${params.ModelTaa}" != "false" && "${params.nonparametric}" != "false" ]];then
                            mod=\$(tail -12 ${reps}.log | head -1 | awk '{print \$6}')
                            iqtree -s ${reps} --prefix ${params.projtag}_ASV_Group_Reps_iq -m \${mod} --redo -nt auto -b ${params.boots}

                        elif [[ "${params.ModelTaa}" != "false" && "${params.parametric}" != "false" ]];then
                            mod=\$(tail -12 ${reps}.log | head -1 | awk '{print \$6}')
                            iqtree -s ${reps} --prefix ${params.projtag}_ASV_Group_Reps_iq -m \${mod} --redo -nt auto -bb ${params.boots} -bnni

                        elif [ "${params.nonparametric}" != "false" ];then
                            iqtree -s ${reps} --prefix ${params.projtag}_ASV_Group_Reps_iq -m MFP -madd --redo -nt auto -b ${params.boots}

                        elif [ "${params.parametric}" != "false" ];then
                            iqtree -s ${reps} --prefix ${params.projtag}_ASV_Group_Reps_iq -m MFP -madd --redo -nt auto -bb ${params.boots} -bnni

                        else
                            iqtree -s ${reps} --prefix ${params.projtag}_ASV_Group_Reps_iq -m MFP -madd --redo -nt auto -bb ${params.boots} -bnni
                        fi
                        """
                    }

              process Adding_ASV_MED_Info {

                  label 'low_cpus'

                  publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ASVs/MED/", mode: "copy", overwrite: true

                  input:
                      file(counts) from asvcount_med
                      file(map) from asvgroupscsv

                  output:
                      file("${params.projtag}_ASV_Groupingcounts.csv") into (asvgroupcounts, asvgroupcounts2)

                  script:
                      """
                      awk -F "," '{print \$1}' ${counts} | sed '1d' > asv.list
                      echo "GroupID" >> group.list
                      for x in \$(cat asv.list);
                      do    group=\$(grep -w \$x ${map} | awk -F "," '{print \$2}')
                            echo "\$group" >> group.list
                      done
                      paste -d',' group.list ${counts} > ${params.projtag}_ASV_Groupingcounts.csv
                      """
                  }
            } else {
                asvgroupscsv = Channel.value('skipping')
                asvgroupscsv2 = Channel.value('skipping')
                asv_group_rep_tree = Channel.value('skipping')
                asvgroupcounts = Channel.value('skipping')
                asvgroupcounts2 = Channel.value('skipping')
            }

            if (!params.skipAminoTyping) {

                process Translate_For_AminoTyping {

                  label 'low_cpus'

                  publishDir "${params.workingdir}/${params.outdir}/Analyze/Clustering/AminoTypes/Translation", mode: "copy", overwrite: true

                  conda (params.condaActivate ? null : null)

                  container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? null : null)

                  input:
                      file(fasta) from asvsforAminotyping

                  output:
                      file("${params.projtag}_all_translations.fasta") into amintypegen
                      file("${params.projtag}_translation_report") into proteinstage_vap_report

                  script:
                      """
                      ${tools}/virtualribosomev2/dna2pep.py ${fasta} -r all -x -o none --fasta ${params.projtag}_all_translaton.fasta --report ${params.projtag}_translation_report
                      awk '/^>/ { print (NR==1 ? "" : RS) \$0; next } { printf "%s", \$0 } END { printf RS }' ${params.projtag}_all_translaton.fasta > ${params.projtag}_all_translations.fasta
                      rm ${params.projtag}_all_translaton.fasta
                      """
                }

                process Generate_AminoTypes { //CHECK rename_seq.py in container !!!!!!!!!!!

                    label 'norm_cpus'

                    publishDir "${params.workingdir}/${params.outdir}/Analyze/Clustering/AminoTypes/SummaryFiles", mode: "copy", overwrite: true, pattern: '*.{clstr,csv,gc}'
                    publishDir "${params.workingdir}/${params.outdir}/Analyze/Clustering/AminoTypes/Problematic", mode: "copy", overwrite: true, pattern: '*problematic*.{fasta}'
                    publishDir "${params.workingdir}/${params.outdir}/Analyze/Clustering/AminoTypes", mode: "copy", overwrite: true, pattern: '*AminoTypes_noTaxonomy.{fasta}'

                    conda (params.condaActivate ? "-c bioconda -c conda-forge vsearch=2.21.1 cd-hit=4.8.1 seqtk=1.3 bbmap=39.01" : null)

                    container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/mulled-v2-089fb9f3537921c3d6dbcc7521fbc33d82301df5:1e1ccff83e5d9864e7f3c008bd4ece458ffbdb8d-0" : "quay.io/biocontainers/mulled-v2-089fb9f3537921c3d6dbcc7521fbc33d82301df5:1e1ccff83e5d9864e7f3c008bd4ece458ffbdb8d-0")

                    input:
                        file(prot) from amintypegen
                        file(asvs) from asvaminocheck

                    output:
                        tuple file("*.fasta"), file("${params.projtag}_AminoTypes.clstr"), file("${params.projtag}_clustered.gc") into ( supplementalfiles )
                        file("${params.projtag}_AminoTypes_noTaxonomy.fasta") into ( aminotypesCounts, aminotypesMafft, aminotypesClustal, aminotypesBlast, aminotypesEmboss, aminos_for_med )
                        file("${params.projtag}_AminoType_summary_map.csv") into aminomapmed

                    script:
                        """
                        set +e
                        cp ${params.vampdir}/bin/rename_seq.py .
                        ####awk 'BEGIN{RS=">";ORS=""}length(\$2)>=${params.minAA}{print ">"\$0}' ${prot} >${params.projtag}_filtered_translations.fasta
                        awk -v RS='>[^\n]+\n' 'length() >= ${params.minAA} {printf "%s", prt \$0} {prt = RT}' ${prot} >${params.projtag}_filtered_translations.fasta
                        ####awk 'BEGIN{RS=">";ORS=""}length(\$2)<${params.minAA}{print ">"\$0}' ${prot} >${params.projtag}_problematic_translations.fasta
                        awk -v RS='>[^\n]+\n' 'length() < ${params.minAA} {printf "%s", prt \$0} {prt = RT}' ${prot} >${params.projtag}_problematic_translations.fasta
                        if [ `wc -l ${params.projtag}_problematic_translations.fasta | awk '{print \$1}'` -gt 1 ];then
                            grep ">" ${params.projtag}_problematic_translations.fasta | awk -F ">" '{print \$2}' > problem_tmp.list
                            seqtk subseq ${asvs} problem_tmp.list > ${params.projtag}_problematic_nucleotides.fasta
                        else
                            rm ${params.projtag}_problematic_translations.fasta
                        fi
                        cd-hit -i ${params.projtag}_filtered_translations.fasta -c 1.0 -o ${params.projtag}_unlabeled_types.fasta
                        sed 's/>Cluster />Cluster_/g' ${params.projtag}_unlabeled_types.fasta.clstr >${params.projtag}_AminoTypes.clstr
                        grep ">Cluster_" ${params.projtag}_AminoTypes.clstr >tmpclusters.list
                        grep -w "*" ${params.projtag}_AminoTypes.clstr | awk '{print \$3}' | awk -F "." '{print \$1}' >tmphead.list
                        grep -w "*" ${params.projtag}_AminoTypes.clstr | awk '{print \$2}' | awk -F "," '{print \$1}' >tmplen.list
                        paste -d"," tmpclusters.list tmphead.list >tmp.info.csv
                        grep ">" ${params.projtag}_unlabeled_types.fasta >lala.list
                        j=1
                        for x in \$(cat lala.list);do
                            echo ">${params.projtag}_AminoType\${j}" >>${params.projtag}_aminoheaders.list
                            echo "\${x},>${params.projtag}_AminoType\${j}" >>tmpaminotype.info.csv
                            j=\$(( \${j}+1 ))
                        done
                        rm lala.list
                        awk -F "," '{print \$2}' tmp.info.csv >>tmporder.list
                        for x in \$(cat tmporder.list);do
                         	grep -w "\$x" tmpaminotype.info.csv | awk -F "," '{print \$2}' >>tmpder.list
                        done
                        paste -d "," tmpclusters.list tmplen.list tmphead.list tmpder.list >${params.projtag}_AminoType_summary_map.csv
                        rm tmp*
                        ./rename_seq.py ${params.projtag}_unlabeled_types.fasta ${params.projtag}_aminoheaders.list ${params.projtag}_AminoTypes_noTaxonomy.fasta
                        stats.sh in=${params.projtag}_AminoTypes_noTaxonomy.fasta gc=${params.projtag}_clustered.gc gcformat=4
                        """
        	    }

                process Generate_AminoType_Matrices {

                    label 'low_cpus'

                    publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/AminoTypes/Matrices", mode: "copy", overwrite: true

                    conda (params.condaActivate ? "-c bioconda -c conda-forge clustalo=1.2.4=h87f3376_5" : null)

                    container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/clustalo:1.2.4--h87f3376_5" : "quay.io/biocontainers/clustalo:1.2.4--h87f3376_5")

                    input:
                        file(prot) from aminotypesClustal

                    output:
                        file("*.matrix") into proclustmatrices
                        file("*PercentID.matrix") into aminotype_heatmap

                    script:
                        """
                        name=\$( echo ${prot} | awk -F "_noTax" '{print \$1}')
                        clustalo -i ${prot} --distmat-out=\${name}_PairwiseDistanceq.matrix --full --force --threads=${task.cpus}
                        clustalo -i ${prot} --distmat-out=\${name}_PercentIDq.matrix --percent-id --full --force --threads=${task.cpus}
                        for x in *q.matrix;do
                            pre=\$(echo "\$x" | awk -F "q.matrix" '{print \$1}')
                            ya=\$(wc -l \$x | awk '{print \$1}')
                            echo "\$((\$ya-1))"
                            tail -"\$(( \$ya-1))" \$x > \${pre}z.matrix
                            rm \$x
                            cat \${pre}z.matrix | sed 's/ /,/g' | sed -E 's/(,*),/,/g' >\${pre}.matrix
                            rm \${pre}z.matrix
                        done
                        """
                }

                if (!params.skipEMBOSS) {

                    process AminoType_EMBOSS_Analyses {

                        label 'low_cpus'

                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/AminoTypes/EMBOSS/2dStructure", mode: "copy", overwrite: true, pattern: '*.{garnier}'
                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/AminoTypes/EMBOSS/HydrophobicMoment", mode: "copy", overwrite: true, pattern: '*HydrophobicMoments*'
                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/AminoTypes/EMBOSS/IsoelectricPoint", mode: "copy", overwrite: true, pattern: '*IsoelectricPoint.{iep,svg}'
                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/AminoTypes/EMBOSS/ProteinProperties", mode: "copy", overwrite: true, pattern: '*.{pepstats,pepinfo}'
                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/AminoTypes/EMBOSS/ProteinProperties/Plots", mode: "copy", overwrite: true, pattern: '*PropertiesPlot.{svg}'
                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/AminoTypes/EMBOSS/2dStructure/Plots", mode: "copy", overwrite: true, pattern: '*Helical*.{svg}'

                        conda (params.condaActivate ? "-c bioconda -c conda-forge emboss=6.6.0=h86d058a_5" : null)

                        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/emboss:6.6.0--h86d058a_5" : "quay.io/biocontainers/emboss:6.6.0--h86d058a_5")

                        input:
                            file(prot) from aminotypesEmboss

                        output:
                            tuple file("*.garnier"), file("*HydrophobicMoments.svg"), file("*IsoelectricPoint*"), file("*.pepstats"), file("*PropertiesPlot*"), file("*Helical*")  into amino_emboss

                        script:
                            """
                            name=\$( echo ${prot} | awk -F ".fasta" '{print \$1}')
                            garnier -sequence ${prot} -outfile \${name}_2dStructures.garnier
                            hmoment -seqall ${prot} -graph svg -plot -double
                            hmoment -seqall ${prot} -double -outfile ${prot}_HydrophobicMoments
                            mv hmoment.svg ./"\${name}"_HydrophobicMoments.svg
                            iep -sequence ${prot} -graph svg -plot -outfile "\${name}"_IsoelectricPoint.iep
                            mv iep.svg ./"\${name}"_IsoelectricPoint.svg
                            pepstats -sequence ${prot} -outfile \${name}_ProteinProperties.pepstats
                            grep ">" ${prot} | awk -F ">" '{print \$2}' > tmpsequence.list
                            for x in \$(cat tmpsequence.list);do
                                grep -A 1 ""\$x"" ${prot} > tmp2.fasta
                                len=\$(tail -1 tmp2.fasta | awk '{print length}')
                                pepinfo -sequence tmp2.fasta -graph svg -outfile "\$x"_PropertiesPlot.pepinfo
                                mv pepinfo.svg ./"\$x"_PropertiesPlot.svg
                                cat "\$x"_PropertiesPlot.pepinfo >> "\${name}"_PropertiesPlot.pepinfo
                                rm "\$x"_PropertiesPlot.pepinfo
                                pepnet -sask -sequence tmp2.fasta -graph svg -sbegin1 1 -send1 \$len
                                mv pepnet.svg ./"\$x"_HelicalNet.svg
                                pepwheel -sequence tmp2.fasta -graph svg -sbegin1 1 -send1 \$len
                                mv pepwheel.svg ./"\$x"_HelicalWheel.svg
                                rm tmp2.fasta
                            done
                            rm tmpsequence.list
                            """
                    }
                }

                if (!params.skipTaxonomy) {

                  if (params.dbtype == "NCBI") {

                    process AminoType_Taxonomy_Inference_NCBI {//CHECK rename_seq.py in container !!!!!!!!!!!

                        label 'high_cpus'

                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/AminoTypes/Taxonomy", mode: "copy", overwrite: true, pattern: '*.{csv,tsv}'
                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/AminoTypes/Taxonomy", mode: "copy", overwrite: true, pattern: '*TaxonomyLabels.fasta'

                        conda (params.condaActivate ? "-c conda-forge bioconda::diamond=2.0.15=hb97b32f_1" : null)

                        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/diamond:2.0.15--hb97b32f_1" : "quay.io/biocontainers/diamond:2.0.15--hb97b32f_1")

                        input:
                            file(asvs) from aminotypesBlast

                        output:
                            tuple file("*_phyloformat.csv"), file("*_summaryTable.tsv"), file("*dmd.out") into summary_AA_diamond
                            file("*_summary_for_plot.csv") into taxplot2
                            file("*TaxonomyLabels.fasta") into tax_labeled_fasta2
                            file("*_quick_Taxbreakdown.csv") into tax_table_amino
                            file ("*_quicker_taxbreakdown.csv") into tax_nodCol_amino

                        script:
                            """
                            cp ${params.vampdir}/bin/rename_seq.py .
                            virdb=${params.dbdir}/${params.dbname}
                            if [[ ${params.measurement} == "bitscore" ]]
                            then    measure="--min-score ${params.bitscore}"
                            elif    [[ ${params.measurement} == "evalue" ]]
                            then    measure="-e ${params.evalue}"
                            else    measure="--min-score ${params.bitscore}"
                            fi
                            grep ">" \${virdb} > headers.list
                            headers="headers.list"
                            name=\$( echo ${asvs} | awk -F ".fasta" '{print \$1}')
                            if [[ ${params.ncbitax} == "true" ]]
                            then   diamond blastp -q ${asvs} -d \${virdb} -p ${task.cpus} --id ${params.minID} -l ${params.minaln} \${measure} --${params.sensitivity} -o "\$name"_dmd.out -f 6 qseqid qlen sseqid qstart qend qseq sseq length qframe evalue bitscore pident btop staxids sskingdoms skingdoms sphylums --max-target-seqs 1 --max-hsps 1
                            else   diamond blastp -q ${asvs} -d \${virdb} -p ${task.cpus} --id ${params.minID} -l ${params.minaln} \${measure} --${params.sensitivity} -o "\$name"_dmd.out -f 6 qseqid qlen sseqid qstart qend qseq sseq length qframe evalue bitscore pident btop --max-target-seqs 1 --max-hsps 1
                            fi
                            echo "Preparing lists to generate summary .csv's"
                            echo "[Best hit accession number]" > access.list
                            echo "[e-value]" > evalue.list
                            echo "[Bitscore]" > bit.list
                            echo "[Percent ID (aa)]" > pid.list
                            echo "[Organism ID]" > "\$name"_virus.list
                            echo "[Gene]" > "\$name"_genes.list
                            echo "[AminoType#]" > otu.list
                            echo "[Sequence length]" > length.list
                            grep ">" ${asvs} | awk -F ">" '{print \$2}' > seqids.lst
                            if [[ ${params.lca} == "T" ]]
                            then  grep -w "LCA" ${params.dbanno}/*.txt > lcainfo.list
                                  echo "[Taxonomic classification from RVDB annotations]" > lca_classification.list
                            else
                                  echo "[Taxonomic classification from RVDB annotations]" > lca_classification.list
                            fi
                            if [[ ${params.ncbitax} == "true" ]]
                            then echo "[NCBI Taxonomy ID],[Taxonomic classification from NCBI]" > ncbi_classification.list
                            fi
                            echo "extracting genes and names"
                            touch new_"\$name"_asvnames.txt
                            for s in \$(cat seqids.lst);do
                                echo "Checking for \$s hit in diamond output"
                                if [[ "\$(grep -wc "\$s" "\$name"_dmd.out)" -eq 1 ]];then
                                            echo "Yep, there was a hit for \$s"
                                            echo "Extracting the information now:"
                                            acc=\$(grep -w "\$s" "\$name"_dmd.out | awk '{print \$3}')
                                            echo "\$s" >> otu.list
                                            echo "\$acc" >> access.list
                                            line="\$(grep -w "\$s" "\$name"_dmd.out)"
                                            echo "\$line" | awk '{print \$10}' >> evalue.list
                                            echo "\$line" | awk '{print \$11}' >> bit.list
                                            echo "\$line" | awk '{print \$12}' >> pid.list
                                            echo "\$line" | awk '{print \$2}' >> length.list
                                            echo "Extracting virus and gene ID for \$s now"
                                            gene=\$(grep -w "\$acc" "\$headers" | awk -F "." '{ print \$2 }' | awk -F "[" '{ print \$1 }' | awk -F " " '{print substr(\$0, index(\$0,\$2))}' | sed 's/ /_/g') &&
                                            echo "\$gene" | sed 's/_/ /g' >> "\$name"_genes.list
                                            virus=\$(grep -w "\$acc" "\$headers" | awk -F "[" '{ print \$2 }' | awk -F "]" '{ print \$1 }'| sed 's/ /_/g')
                                            echo "\$virus" | sed 's/_/ /g' >> "\$name"_virus.list
                                            echo ">"\${s}"_"\$virus"_"\$gene"" >> new_"\$name"_asvnames.txt
                                            if [[ "${params.lca}" == "T" ]]
                                            then    if [[ \$(grep -w "\$acc" ${params.dbanno}/*.txt | wc -l) -eq 1 ]]
                                                    then  group=\$(grep -w "\$acc" ${params.dbanno}/*.txt | awk -F ":" '{print \$1}')
                                                          lcla=\$(grep -w "\$group" lcainfo.list | awk -F "\\t" '{print \$2}')
                                                          echo "\$lcla" >> lca_classification.list
                                                    else  echo "Viruses" >> lca_classification.list
                                                    fi
                                            fi
                                            if [[ ${params.ncbitax} == "true" ]]
                                            then  echo "\$line" | awk -F "\t" '{print \$14","\$16"::"\$18"::"\$17}' >> ncbi_classification.list
                                            fi
                                            echo "\$s done."
                                else
                                            echo "Ugh, there was no hit for \$s .."
                                            echo "We still love \$s though and we will add it to the final fasta file"
                                            echo "\$s" >> otu.list
                                            echo "NO_HIT" >> access.list
                                            echo "NO_HIT" >> "\$name"_genes.list
                                            echo "NO_HIT" >> "\$name"_virus.list
                                            echo "NO_HIT" >> evalue.list
                                            echo "NO_HIT" >> bit.list
                                            echo "NO_HIT" >> pid.list
                                            echo "NO_HIT" >> length.list
                                            virus="NO"
                                            gene="HIT"
                                            echo ">\${s}_"\$virus"_"\$gene"" >> new_"\$name"_asvnames.txt
                                            if [[ "${params.lca}" == "T" ]]
                                            then    echo "N/A" >> lca_classification.list
                                            fi
                                            if [[ "${params.ncbitax}" == "true" ]]
                                            then  echo "N/A" >> ncbi_classification.list
                                            fi
                                            echo "\$s done."
                               fi
                            done
                            echo "Now editing "\$name" fasta headers"
                            ###### rename_seq.py
                            ./rename_seq.py ${asvs} new_"\$name"_asvnames.txt "\$name"_TaxonomyLabels.fasta
                            awk 'BEGIN {RS=">";FS="\\n";OFS=""} NR>1 {print ">"\$1; \$1=""; print}' "\$name"_TaxonomyLabels.fasta >"\$name"_tmpssasv.fasta
                            echo "[Sequence header]" > newnames.list
                            cat new_"\$name"_asvnames.txt >> newnames.list
                            touch sequence.list
                            echo "     " > sequence.list
                            grep -v ">" "\$name"_tmpssasv.fasta >> sequence.list
                            rm "\$name"_tmpssasv.fasta
                            if [[ "${params.lca}" == "T" && "${params.ncbitax}" == "true" ]]
                            then
                                  paste -d "," sequence.list "\$name"_virus.list "\$name"_genes.list otu.list lca_classification.list ncbi_classification.list newnames.list length.list bit.list evalue.list pid.list access.list >> "\$name"_phyloformat.csv
                                  paste -d"\t" otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list ncbi_classification.list sequence.list length.list bit.list evalue.list pid.list >> "\$name"_summaryTable.tsv
                                  paste -d"," otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list ncbi_classification.list >> \${name}_quick_Taxbreakdown.csv
                            elif [[ "${params.lca}" == "T" && "${params.ncbitax}" != "true" ]]
                            then
                                  paste -d "," sequence.list "\$name"_virus.list "\$name"_genes.list otu.list lca_classification.list newnames.list length.list bit.list evalue.list pid.list access.list >> "\$name"_phyloformat.csv
                                  paste -d"\t" otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list sequence.list length.list bit.list evalue.list pid.list >> "\$name"_summaryTable.tsv
                                  paste -d"," otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list >> \${name}_quick_Taxbreakdown.csv
                            elif [[ "${params.ncbitax}" == "true" && "${params.lca}" != "T" ]]
                            then
                                  paste -d "," sequence.list "\$name"_virus.list "\$name"_genes.list otu.list ncbi_classification.list newnames.list length.list bit.list evalue.list pid.list access.list >> "\$name"_phyloformat.csv
                                  paste -d"\t" otu.list access.list "\$name"_virus.list "\$name"_genes.list ncbi_classification.list sequence.list length.list bit.list evalue.list pid.list >> "\$name"_summaryTable.tsv
                                  paste -d"," otu.list access.list "\$name"_virus.list "\$name"_genes.list ncbi_classification.list >> \${name}_quick_Taxbreakdown.csv
                            else
                                  paste -d "," sequence.list "\$name"_virus.list "\$name"_genes.list otu.list newnames.list length.list bit.list evalue.list pid.list access.list >> "\$name"_phyloformat.csv
                                  paste -d"\t" otu.list access.list "\$name"_virus.list "\$name"_genes.list sequence.list length.list bit.list evalue.list pid.list >> "\$name"_summaryTable.tsv
                                  echo "skipped" >> \${name}_quick_Taxbreakdown.csv
                            fi
                            for x in *phyloformat.csv;do
                                echo "\$x"
                                lin=\$(( \$(wc -l \$x | awk '{print \$1}')-1))
                                tail -"\$lin" \$x | awk -F "," '{print \$2}' > tmpcol.list;
                                sed 's/ /_/g' tmpcol.list > tmp2col.list;
                                cat tmp2col.list | sort | uniq -c | sort -nr | awk '{print \$2","\$1}' > \${name}_summary_for_plot.csv;
                                rm tmpcol.list tmp2col.list
                            done
                            awk -F "," '{print \$1","\$3"("\$2")"}' \${name}_quick_Taxbreakdown.csv >> \${name}_quicker_taxbreakdown.csv
                            rm evalue.list sequence.list bit.list pid.list length.list seqids.lst otu.list *asvnames.txt "\$name"_virus.list "\$name"_genes.list newnames.list access.list headers.list
                            """
                        }
                    } else if (params.dbtype== "RVDB") {

                      process AminoType_Taxonomy_Inference_RVDB {//CHECK rename_seq.py in container !!!!!!!!!!!

                          label 'high_cpus'

                          publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/AminoTypes/Taxonomy", mode: "copy", overwrite: true, pattern: '*.{csv,tsv}'
                          publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/AminoTypes/Taxonomy", mode: "copy", overwrite: true, pattern: '*TaxonomyLabels.fasta'

                          conda (params.condaActivate ? "-c conda-forge bioconda::diamond=2.0.15=hb97b32f_1" : null)

                          container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/diamond:2.0.15--hb97b32f_1" : "quay.io/biocontainers/diamond:2.0.15--hb97b32f_1")

                          input:
                              file(asvs) from aminotypesBlast

                          output:
                              tuple file("*_phyloformat.csv"), file("*_summaryTable.tsv"), file("*dmd.out") into summary_AA_diamond
                              file("*_summary_for_plot.csv") into taxplot2
                              file("*TaxonomyLabels.fasta") into tax_labeled_fasta2
                              file("*_quick_Taxbreakdown.csv") into tax_table_amino
                              file ("*_quicker_taxbreakdown.csv") into tax_nodCol_amino

                          script:
                              """
                              cp ${params.vampdir}/bin/rename_seq.py .
                              virdb=${params.dbdir}/${params.dbname}
                              if [[ ${params.measurement} == "bitscore" ]]
                              then    measure="--min-score ${params.bitscore}"
                              elif    [[ ${params.measurement} == "evalue" ]]
                              then    measure="-e ${params.evalue}"
                              else    measure="--min-score ${params.bitscore}"
                              fi
                              grep ">" \${virdb} > headers.list
                              headers="headers.list"
                              name=\$( echo ${asvs} | awk -F ".fasta" '{print \$1}')
                              diamond blastp -q ${asvs} -d \${virdb} -p ${task.cpus} --id ${params.minID} -l ${params.minaln} \${measure} --${params.sensitivity} -o "\$name"_dmd.out -f 6 qseqid qlen sseqid qstart qend qseq sseq length qframe evalue bitscore pident btop --max-target-seqs 1 --max-hsps 1
                              echo "Preparing lists to generate summary .csv's"
                              echo "[Best hit accession number]" > access.list
                              echo "[e-value]" > evalue.list
                              echo "[Bitscore]" > bit.list
                              echo "[Percent ID (aa)]" > pid.list
                              echo "[Organism ID]" > "\$name"_virus.list
                              echo "[Gene]" > "\$name"_genes.list
                              echo "[AminoType#]" > otu.list
                              echo "[Sequence length]" > length.list
                              grep ">" ${asvs} | awk -F ">" '{print \$2}' > seqids.lst
                              if [[ ${params.lca} == "T" ]]
                              then  grep -w "LCA" ${params.dbanno}/*.txt > lcainfo.list
                                    echo "[Taxonomic classification from RVDB annotations]" > lca_classification.list
                              else  echo "skipped" >> \${name}_quick_Taxbreakdown.csv
                                    echo "[Taxonomic classification from RVDB annotations]" > lca_classification.list
                              fi
                              echo "extracting genes and names"
                              touch new_"\$name"_asvnames.txt
                              for s in \$(cat seqids.lst);do
                                  echo "Using RVDB headers."
                                  if [[ "\$(grep -wc "\$s" "\$name"_dmd.out)" -eq 1 ]];then
                                      echo "Yep, there was a hit for \$s"
                                      echo "Extracting the information now:"
                                      acc=\$(grep -w "\$s" "\$name"_dmd.out | awk '{print \$3}' | awk -F "|" '{print \$3}')
                                      echo "\$s" >> otu.list
                                      echo "\$acc" >> access.list
                                      line="\$(grep -w "\$s" "\$name"_dmd.out)"
                                      echo "\$line" | awk '{print \$10}' >> evalue.list
                                      echo "\$line" | awk '{print \$11}' >> bit.list
                                      echo "\$line" | awk '{print \$12}' >> pid.list
                                      echo "\$line" | awk '{print \$2}' >> length.list
                                      echo "Extracting virus and gene ID for \$s now"
                                      gene=\$(grep -w "\$acc" "\$headers" | awk -F "|" '{ print \$6 }' | awk -F "[" '{ print \$1 }' | sed 's/ /_/g') &&
                                      echo "\$gene" | sed 's/_/ /g' >> "\$name"_genes.list
                                      virus=\$(grep -w "\$acc" "\$headers" | awk -F "|" '{ print \$6 }' | awk -F "[" '{ print \$2 }' | awk -F "]" '{print \$1}' | sed 's/ /_/g') &&
                                      echo "\$virus" | sed 's/_/ /g' >> "\$name"_virus.list
                                      echo ">\${s}_"\$virus"_"\$gene"" >> new_"\$name"_asvnames.txt
                                      if [[ "${params.lca}" == "T" ]]
                                      then    if [[ \$(grep -w "\$acc" ${params.dbanno}/*.txt | wc -l) -eq 1 ]]
                                              then  group=\$(grep -w "\$acc" ${params.dbanno}/*.txt | awk -F ":" '{print \$1}')
                                                    lcla=\$(grep -w "\$group" lcainfo.list | awk -F "\t" '{print \$2}')
                                                    echo "\$lcla" >> lca_classification.list
                                              else  echo "Viruses" >> lca_classification.list
                                              fi
                                      fi
                                      echo "\$s done."
                                  else
                                      echo "Ugh, there was no hit for \$s .."
                                      echo "We still love \$s though and we will add it to the final fasta file"
                                      echo "\$s" >> otu.list
                                      echo "NO_HIT" >> access.list
                                      echo "NO_HIT" >> "\$name"_genes.list
                                      echo "NO_HIT" >> "\$name"_virus.list
                                      echo "NO_HIT" >> evalue.list
                                      echo "NO_HIT" >> bit.list
                                      echo "NO_HIT" >> pid.list
                                      echo "NO_HIT" >> length.list
                                      virus="NO"
                                      gene="HIT"
                                      echo ">\${s}_"\$virus"_"\$gene"" >> new_"\$name"_asvnames.txt
                                      if [[ "${params.lca}" == "T" ]]
                                      then    echo "N/A" >> lca_classification.list
                                      fi
                                      echo "\$s done."
                                  fi
                              echo "Done with \$s"
                              done
                              echo "Now editing "\$name" fasta headers"
                              ###### rename_seq.py
                              ./rename_seq.py ${asvs} new_"\$name"_asvnames.txt "\$name"_TaxonomyLabels.fasta
                              awk 'BEGIN {RS=">";FS="\\n";OFS=""} NR>1 {print ">"\$1; \$1=""; print}' "\$name"_TaxonomyLabels.fasta >"\$name"_tmpssasv.fasta
                              echo "[Sequence header]" > newnames.list
                              cat new_"\$name"_asvnames.txt >> newnames.list
                              touch sequence.list
                              echo "     " > sequence.list
                              grep -v ">" "\$name"_tmpssasv.fasta >> sequence.list
                              rm "\$name"_tmpssasv.fasta
                              if [[ "${params.lca}" == "T" ]]
                              then  paste -d "," sequence.list "\$name"_virus.list "\$name"_genes.list otu.list lca_classification.list newnames.list length.list bit.list evalue.list pid.list access.list >> "\$name"_phyloformat.csv
                                    paste -d"\t" otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list sequence.list length.list bit.list evalue.list pid.list >> "\$name"_summaryTable.tsv
                                    paste -d"," otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list >> \${name}_quick_Taxbreakdown.csv
                              else  paste -d "," sequence.list "\$name"_virus.list "\$name"_genes.list otu.list newnames.list length.list bit.list evalue.list pid.list access.list >> "\$name"_phyloformat.csv
                                    paste -d"\t" otu.list access.list "\$name"_virus.list "\$name"_genes.list sequence.list length.list bit.list evalue.list pid.list >> "\$name"_summaryTable.tsv
                              fi
                              for x in *phyloformat.csv;do
                                        echo "\$x"
                                        lin=\$(( \$(wc -l \$x | awk '{print \$1}')-1))
                                        tail -"\$lin" \$x | awk -F "," '{print \$2}' > tmpcol.list;
                                        sed 's/ /_/g' tmpcol.list > tmp2col.list;
                                        cat tmp2col.list | sort | uniq -c | sort -nr | awk '{print \$2","\$1}' > \${name}_summary_for_plot.csv;
                                        rm tmpcol.list tmp2col.list
                              done
                              awk -F "," '{print \$1","\$3"("\$2")"}' \${name}_quick_Taxbreakdown.csv >> \${name}_quicker_taxbreakdown.csv
                              rm evalue.list sequence.list bit.list pid.list length.list seqids.lst otu.list *asvnames.txt "\$name"_virus.list "\$name"_genes.list newnames.list access.list headers.list
                              """
                    }
                }
             } else {
                 taxplot2 = Channel.value('skipping')
                 tax_table_amino = Channel.value('skipping')
                 tax_nodCol_amino = Channel.value('skipping')
             }


                if (!params.skipPhylogeny || params.aminoMED || params.aminoTClust) {

                    process AminoType_pre_Phylogeny_step1 {

                        label 'norm_cpus'

                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/AminoTypes/Phylogeny/Alignment", mode: "copy", overwrite: true

                        conda (params.condaActivate ? "-c bioconda -c conda-forge muscle=5.1" : null)

                        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/muscle:5.1--h9f5acd7" : "quay.io/biocontainers/muscle:5.1--h9f5acd7")
                        input:
                            file(amino) from aminotypesMafft

                        output:
                            file("*_ALN.fasta") into amino_align1
                            file("*.efa")

                        script:
                            """
                            pre=\$(echo ${amino} | awk -F ".fasta" '{print \$1}' )

                            if [[ ${params.srep} == "true" && ${params.ensemble} == "false" ]];
                            then  if [[ \$( grep -c ">" ${amino}) -lt 300 ]]
                                  then    comm="align"
                                  else    comm="super5"
                                  fi
                                  muscle -"\$comm" ${amino} -perm ${params.perm} -perturb ${params.pert} -output \${pre}_muscle_raw_ALN.fasta -threads ${task.cpus} -quiet
                                  echo "single replicate alignment chosen; look over muscle5 documentation to learn about ensemble alignment approach" >note.efa
                            elif [[ ${params.srep} == "false" && ${params.ensemble} == "true" ]];
                            then  muscle -align ${amino} -${params.fied} -output \${pre}_muscle.efa -threads ${task.cpus} -quiet
                                  muscle -maxcc \${pre}_muscle.efa -output \${pre}_muscle_raw_ALN.fasta
                            else  if [[ \$( grep -c ">" ${amino}) -lt 300 ]]
                                  then    comm="align"
                                  else    comm="super5"
                                  fi
                                  muscle -"\$comm" ${amino} -output \${pre}_muscle_raw_ALN.fasta -threads ${task.cpus} -quiet
                                  echo "single replicate alignment chosen; look over muscle5 documentation to learn about ensemble alignment approach" >note.efa
                            fi
                            """
                    }

                    process AminoType_pre_Phylogeny_step2 {

                        label 'norm_cpus'

                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/AminoTypes/Phylogeny/Alignment", mode: "copy", overwrite: true

                        conda (params.condaActivate ? "-c bioconda -c conda-forge trimal=1.4.1" : null)

                        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/trimal:1.4.1--h9f5acd7_6" : "quay.io/biocontainers/trimal:1.4.1--h9f5acd7_6")

                        input:
                            file(align) from amino_align1

                        output:
                            file("*_aln.html") into trimalhtml
                            file("*_aln.fasta") into amino_align2

                        script:
                            """
                            pre=\$(echo ${align} | awk -F "_muscle" '{print \$1}' )
                            trimal -in ${align} -out \${pre}_trimal_aln.fasta -keepheader -fasta -automated1 -htmlout \${pre}_trimal_aln.html
                            """
                    }

                    process AminoType_pre_Phylogeny_step3 {

                        label 'norm_cpus'

                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/AminoTypes/Phylogeny/Alignment", mode: "copy", overwrite: true

                        conda (params.condaActivate ? "${params.vampdir}/bin/yamls/oligotyping.yml" : null)

                        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/oligotyping:2.1--py27_0" : "quay.io/biocontainers/oligotyping:2.1--py27_0")

                        input:
                            file(align) from amino_align2

                        output:
                            file("*.fasta") into amino_align3_mt, amino_align3_iq, amino_align3_med

                        script:
                            """
                            o-trim-uninformative-columns-from-alignment ${align}
                            pre=\$(echo ${align} | awk -F "_trimal" '{print \$1}' )
                            mv ${align}-TRIMMED ./\${pre}_Aligned_informativeonly.fasta
                            """
                    }

                    process AminoType_pre_Phylogeny_step4 {

                          label 'norm_cpus'

                          publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/AminoTypes/Phylogeny/ModelTest", mode: "copy", overwrite: true, pattern: '*mt*'

                          conda (params.condaActivate ? "-c conda-forge bioconda::modeltest-ng=0.1.7=h5c6ebe3_0" : null)

                          container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/modeltest-ng:0.1.7--h5c6ebe3_0" : "quay.io/biocontainers/modeltest-ng:0.1.7--h5c6ebe3_0")

                          input:
                              file(align) from amino_align3_mt

                          output:
                              file("*mt.out") into amino_align4
                              file("*mt.*") into aminomtres

                          script:
                              """
                              # Nucleotide_ModelTest
                              pre=\$(echo ${align} | awk -F "_Aligned_informativeonly.fasta" '{print \$1}' )
                              modeltest-ng -i ${align} -p ${task.cpus} -o \${pre}_mt -d aa -s 203 --disable-checkpoint
                              """
                      }

                  }

                  if (!params.skipPhylogeny || params.asvTClust) {

                      process AminoType_Phylogeny_step5 {

                            label 'norm_cpus'

                            publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ASVs/Phylogeny/IQ-TREE", mode: "copy", overwrite: true, pattern: '*iq*'

                            conda (params.condaActivate ? "-c conda-forge bioconda::iqtree=2.2.0.3=hb97b32f_1" : null)

                            container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/iqtree:2.2.0.3--hb97b32f_1" : "quay.io/biocontainers/iqtree:2.2.0.3--hb97b32f_1")

                            input:
                                file(align) from amino_align3_iq
                                file(mtout) from amino_align4

                            output:
                                file("*iq*") into amino_align5
                                file("*iq.treefile") into (amino_rax_plot, amino_repphy, amino_treeclust)

                            script:
                                """
                                #grabbing best models from modeltestng
                                modbic=\$(grep "iqtree" ${mtout} | head -1 | awk -F "-m " '{print \$2}')
                                modaic=\$(grep "iqtree" ${mtout} | head -2 | tail -1 | awk -F "-m " '{print \$2}')
                                modaicc=\$(grep "iqtree" ${mtout} | tail -1 | awk -F "-m " '{print \$2}')
                                if [[ ${params.crit} == "BIC" ]]
                                then  mod="\$modbic"
                                elif  [[ ${params.crit} == "AIC" ]]
                                then  mod="\$modaic"
                                elif  [[ ${params.crit} == "AICc" ]]
                                then  mod="\$modaicc"
                                fi
                                # grab prefix
                                pre=\$(echo ${align} | awk -F "_Aligned_informativeonly.fasta" '{print \$1}' )
                                # Nucleotide_Phylogeny
                                if [ "${params.iqCustomnt}" != "" ];then
                                    iqtree -s ${align} --prefix \${pre}_iq --redo -T auto ${params.iqCustomnt}
                                elif [[ "${params.ModelTnt}" != "false" && "${params.nonparametric}" != "false" ]];then
                                    iqtree -s ${align} --prefix \${pre}_iq -m \${mod} --redo -nt auto -b ${params.boots}
                                elif [[ "${params.ModelTnt}" != "false" && "${params.parametric}" != "false" ]];then
                                    iqtree -s ${align} --prefix \${pre}_iq -m \${mod} --redo -nt auto -bb ${params.boots} -bnni
                                elif [ "${params.nonparametric}" != "false" ];then
                                    iqtree -s ${align} --prefix \${pre}_iq -m MFP -madd --redo -nt auto -b ${params.boots}
                                elif [ "${params.parametric}" != "false" ];then
                                    iqtree -s ${align} --prefix \${pre}_iq -m MFP -madd --redo -nt auto -bb ${params.boots} -bnni
                                else
                                    iqtree -s ${align} --prefix \${pre}_iq -m MFP -madd --redo -nt auto -bb ${params.boots} -bnni
                                fi
                                """
                        }
              } else {
                  amino_rax_plot = Channel.value('skipping')
                  amino_rephy = Channel.value('skipping')
                  amino_treeclust = Channel.value('skipping')
              }

                process Generate_AminoTypes_Counts_Table {

                    label 'high_cpus'

                    publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/AminoTypes/Counts", mode: "copy", overwrite: true

                    conda (params.condaActivate ? "-c conda-forge bioconda::diamond=2.0.15=hb97b32f_1" : null)

                    container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/diamond:2.0.15--hb97b32f_1" : "quay.io/biocontainers/diamond:2.0.15--hb97b32f_1")

                    input:
                        file(fasta) from aminotypesCounts
                        file(merged) from mergeforprotcounts
                        file(samplist) from samplelist

                    output:
                        tuple file("*_AminoType_counts.csv"), file("*dmd.out") into counts_summary
                        file("*_AminoType_counts.csv") into (aminocounts_plot, aminocountmed, amino_countphylo, aminocounts_deseq)

                    script:
                        """
                        set +e
                        diamond makedb --in ${fasta} --db ${fasta}
                        diamond blastx -q ${merged} -d ${fasta} -p ${task.cpus} --min-score ${params.ProtCountsBit} --id ${params.ProtCountID} -l ${params.ProtCountsLength} --${params.sensitivity} -o ${params.projtag}_protCounts_dmd.out -f 6 qseqid qlen sseqid qstart qend qseq sseq length qframe evalue bitscore pident btop --max-target-seqs 1 --max-hsps 1 --max-hsps 1
                        echo "OTU_ID" >tmp.col1.txt
                        echo "Generating sample id list"
                        grep ">" ${fasta} | awk -F ">" '{print \$2}' | sort | uniq > otuid.list
                        cat otuid.list >> tmp.col1.txt
                        echo "Beginning them counts tho my g"
                        for y in \$( cat ${samplist} );do
                            echo "Starting with \$y now ..."
                            grep -w "\$y" ${params.projtag}_protCounts_dmd.out > tmp."\$y".out
                            echo "Isolated hits"
                            echo "Created uniq subject id list"
                            echo "\$y" > "\$y"_col.txt
                            echo "Starting my counts"
                            for z in \$(cat otuid.list);do
                                echo "Counting \$z hits"
                	            echo "grep -wc "\$z" >> "\$y"_col.txt"
                	            grep -wc "\$z" tmp."\$y".out >> "\$y"_col.txt
                	            echo "\$z counted"
                            done
                       done
                       paste -d "," tmp.col1.txt *col.txt > ${params.projtag}_AminoType_counts.csv
                       rm tmp*
                       rm *col.txt
                       """
                }
            }

            if (params.aminoTClust) {

                process AminoType_PhyloClustering {

                      label 'norm_cpus'

                      publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/AminoTypes/TreeCluster", mode: "copy", overwrite: true

                      conda (params.condaActivate ? "${params.vampdir}/bin/yamls/treecluster.yml" : null)

                      container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/treecluster:1.0.3--pyh3252c3a_0" : "quay.io/biocontainers/treecluster:1.0.3--pyh3252c3a_0")

                      input:
                         file(tree) from amino_treeclust
                         file(counts) from amino_countphylo
                      output:
                         file("*treeclustering*.out") into aminotreeclustering_res
                         file("${params.projtag}_amino_phylogroup.csv") into (amino_phylogroupcsv, amino_phylogroupcsv2)
                         file("${params.projtag}_amino_phyloGroupingcounts.csv") into (amino_phylogroupingcsv, amino_phylogroupingcsv2)

                      script:
                          """
                          TreeCluster.py -i ${tree} ${params.asvTCopp} > ${params.projtag}_AminoType_treeclustering.out
                          TreeCluster.py -i ${tree} ${params.asvTCopp} > ${params.projtag}_AminoType_treeclustering_verbose.out
                          #create headless treeclustering.out
                          tail -n +2 ${params.projtag}_AminoType_treeclustering.out | sed 's/-1/0X/g' > headless.treeout
                          #extracting singletons
                          grep -w "0X" headless.treeout > single.out
                          awk -F "\\t" '{print \$1}' single.out > single.list
                          #assigning groupID to singletons
                          for x in \$(awk -F "\\t" '{print \$1}' single.out); do echo ""\$x"XX" >> single.groups;done
                          #summarizing clustering results
                          awk -F "\\t" '{print \$2}' headless.treeout | grep -v "0X" | sort | uniq > grp.list
                          cat single.groups >> group.list
                          cat grp.list >> group.list
                          echo "Sequence_ID,phyloGroup_ID" > ${params.projtag}_amino_phylogroup.csv

                          for x in \$(seq "\$(wc -l group.list | awk '{print \$1}')");
                          do      echo "phyloGroup"\$x"" >> grup.list
                          done
                          paste -d "," grup.list group.list > groups.csv
                          rm grup.list group.list
                          awk -F "," '{print \$1}' ${counts} | sed '1d' > asv.list
                          for z in \$(cat asv.list);
                          do      if [[ \$(grep -wc "\$z" single.list) -ge 1 ]]
                                  then
                                          group=\$(grep -w ""\$z"XX" groups.csv | awk -F "," '{print \$1}')
                                  else
                                          grp=\$(grep -w "\$z" headless.treeout | awk -F "\\t" '{print \$2}')
                                          group=\$(grep -w "\$grp" groups.csv | awk -F "," '{print \$1}')
                                  fi
                                  echo ""\$z","\$group"" >> ${params.projtag}_amino_phylogroup.csv
                          done
                          awk -F "," '{print \$2}' ${params.projtag}_amino_phylogroup.csv > groups.list
                          paste -d"," groups.list ${counts} > ${params.projtag}_amino_phyloGroupingcounts.csv
                         """
                     }
            } else {
                amino_phylogroupcsv = Channel.value('skipping')
                amino_phylogroupcsv2 = Channel.value('skipping')
                amino_phylogroupingcsv = Channel.value('skipping')
            }

              if (params.aminoMED) {

                    process AminoType_Minimum_Entropy_Decomposition_step1 {

                        label 'low_cpus'

                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Clustering/AminoTypes/MED", mode: "copy", overwrite: true, pattern: '*.{fasta,csv}'
                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Clustering/AminoTypes/MED", mode: "copy", overwrite: true, pattern: '*_unique'

                        conda (params.condaActivate ? "${params.vampdir}/bin/yamls/oligotyping.yml" : null)

                        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/oligotyping:2.1--py27_0" : "quay.io/biocontainers/oligotyping:2.1--py27_0")

                        input:
                            file(align) from amino_align3_med

                        output:

                            file("*_unique") into uniqformedamino
                            file("OLIGO-REPRESENTATIVES.fasta") into oligorepamino

                        script:
                            """
                            pre=\$(echo ${align} | awk -F "_Aligned_informativeonly.fasta" '{print \$1}' )
                            #entopy analysis
                            entropy-analysis ${align}
                            #Decomposition
                            if [[ \$(echo ${params.aminoC} | grep -c ",") -gt 0 ]]
                            then
                                  tag=\$(echo ${params.aminoC} | sed 's/,/_/g')
                                  oligotype ${align} \${pre}_Aligned_informativeonly.fasta-ENTROPY -o ${params.projtag}_AminoTypeMED_"\$tag" -M 1 -C ${params.aminoC} -N ${task.cpus} --skip-check-input --no-figures --skip-gen-html
                            elif [[ "${params.aminoSingle}" == "true" ]]
                            then
                                  tag="${params.aminoC}"
                                  oligotype ${align} \${pre}_Aligned_informativeonly.fasta-ENTROPY -o ${params.projtag}_AminoTypeMED_"\$tag" -M 1 -C ${params.aminoC} -N ${task.cpus} --skip-check-input --no-figures --skip-gen-html
                            else
                                  oligotype ${align} \${pre}_Aligned_informativeonly.fasta-ENTROPY -o ${params.projtag}_AminoTypeMED_${params.aminoC} -M 1 -c ${params.aminoC} -N ${task.cpus} --skip-check-input --no-figures --skip-gen-html
                            fi
                            #generatemaps
                            mv ./${params.projtag}_AminoTypeMED_*/OLIGO-REPRESENTATIVES.fasta .
                            mv ./${params.projtag}_AminoTypeMED_*/OLIGO-REPRESENTATIVES/*_unique .
                            """
                    }

                    process AminoType_Minimum_Entropy_Decomposition_step2 {

                        label 'low_cpus'

                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Clustering/AminoTypes/MED", mode: "copy", overwrite: true

                        conda (params.condaActivate ? "bioconda::seqtk=1.3=h7132678_4" : null)

                        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/seqtk:1.3--h7132678_4" : "quay.io/biocontainers/seqtk:1.3--h7132678_4")

                        input:
                            file(aminos) from aminos_for_med
                            file(oligo) from oligorepamino

                        output:
                            file("*_AminoType_Grouping.csv") into (atygroupscsv, atygroupscsv2)
                            file("${params.projtag}_AminoType_group_reps_aligned.fasta") into (grouprepsatalign, grouprepsatalign2)

                        script:
                            """
                            #copy _unique files here so we can werk with them
                            cp ${params.workingdir}/${params.outdir}/Analyze/Clustering/AminoTypes/MED/unique .
                            #generatemaps
                            echo "AminoType,Group,IDPattern"
                            j=1
                            for x in *_unique;
                            do      gid=\$(echo \$x | awk -F "_" '{print \$1}')
                                    uni=\$(echo \$x | awk -F ""\${gid}"_" '{print \$2}' | awk -F "_uni" '{print \$1}')
                                    grep ">"  "\$gid"_"\$uni" | awk -F ">" '{print \$2}' > asv.list
                                    seqtk subseq ../../${aminos} asv.list > Group"\${j}"_sequences.fasta
                                    for z in \$( cat asv.list)
                                    do      echo ""\$z",Group"\$j","\$uni"" >> ${params.projtag}_AminoType_Grouping.csv

                                    done
                                    rm asv.list
                                    echo ">Group\${j}" >> ${params.projtag}_AminoType_group_reps_aligned.fasta
                                    echo "\$uni" > group.list
                                    seqtk subseq ${oligo} group.list > group.fasta
                                    tail -1 group.fasta >> ${params.projtag}_AminoType_group_reps_aligned.fasta
                                    mv "\$gid"_"\$uni" ./Group"\$j"_"\$uni"_aligned.fasta
                                    mv "\$gid"_"\$uni"_unique ./Group"\$j"_"\$uni"_unqiues_aligned.fasta
                                    rm "\$gid"*.cPickle
                                    j=\$((\$j+1))
                            done
                            #mv ${params.projtag}_AminoType_Grouping.csv ../../
                            #mv ${params.projtag}_AminoType_group_reps_aligned.fasta ../../
                            #cd ..
                            """
                    }

                    if (!params.skipPhylogeny) {

                        process AminoType_MED_Reps_ModelTesting {

                            label 'low_cpus'

                            publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/AminoTypes/MED/Phylogeny/ModelTest", mode: "copy", overwrite: true, pattern: '*ASV*mt*'

                            conda (params.condaActivate ? "bioconda::modeltest-ng=0.1.7=h5c6ebe3_0" : null)

                            container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/modeltest-ng:0.1.7--h5c6ebe3_0" : "quay.io/biocontainers/modeltest-ng:0.1.7--h5c6ebe3_0")

                            input:
                                file(reps) from grouprepsatalign

                            output:
                                file("*mt.out") into aminomedrep_mtout
                                file("*mt.*") into mtresamino

                            script:
                                """
                                # Protein_ModelTest
                                modeltest-ng -i ${reps} -p ${task.cpus} -o ${params.projtag}_AminoType_MEDGroup_Reps_mt -d aa -s 203 --disable-checkpoint
                                """
                        }

                        process AminoType_MED_Reps_phylogeny {

                            label 'low_cpus'

                            publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/AminoTypes/MED/Phylogeny/IQ-TREE", mode: "copy", overwrite: true, pattern: '*iq*'

                            conda (params.condaActivate ? "bioconda::iqtree=2.2.0.3=hb97b32f_1" : null)

                            container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/iqtree:2.2.0.3--hb97b32f_1" : "quay.io/biocontainers/iqtree:2.2.0.3-hb97b32f_1")

                            input:
                                file(reps) from grouprepsatalign2
                                file(mtout) from aminomedrep_mtout

                            output:
                                file("*iq.treefile") into amino_group_rep_tree
                                file("*_AminoType_Group_Reps*") into align_results_asvmed

                            script:
                                """
                                #grabbing best models from modeltestng
                                modbic=\$(grep "iqtree" ${mtout} | head -1 | awk -F "-m " '{print \$2}')
                                modaic=\$(grep "iqtree" ${mtout} | head -2 | tail -1 | awk -F "-m " '{print \$2}')
                                modaicc=\$(grep "iqtree" ${mtout} | tail -1 | awk -F "-m " '{print \$2}')
                                if [[ ${params.crit} == "BIC" ]]
                                then  mod="\$modbic"
                                elif  [[ ${params.crit} == "AIC" ]]
                                then  mod="\$modaic"
                                elif  [[ ${params.crit} == "AICc" ]]
                                then  mod="\$modaicc"
                                fi
                                # Protein_Phylogeny
                                if [ "${params.iqCustomaa}" != "" ];then
                                    iqtree -s ${reps} --prefix ${params.projtag}_AminoType_Group_Reps_iq --redo -T auto ${params.iqCustomaa}

                                elif [[ "${params.ModelTaa}" != "false" && "${params.nonparametric}" != "false" ]];then
                                    mod=\$(tail -12 ${reps}.log | head -1 | awk '{print \$6}')
                                    iqtree -s ${reps} --prefix ${params.projtag}_AminoType_Group_Reps_iq -m \${mod} --redo -nt auto -b ${params.boots}

                                elif [[ "${params.ModelTaa}" != "false" && "${params.parametric}" != "false" ]];then
                                    mod=\$(tail -12 ${reps}.log | head -1 | awk '{print \$6}')
                                    iqtree -s ${reps} --prefix ${params.projtag}_AminoType_Group_Reps_iq -m \${mod} --redo -nt auto -bb ${params.boots} -bnni

                                elif [ "${params.nonparametric}" != "false" ];then
                                    iqtree -s ${reps} --prefix ${params.projtag}_AminoType_Group_Reps_iq -m MFP -madd --redo -nt auto -b ${params.boots}

                                elif [ "${params.parametric}" != "false" ];then
                                    iqtree -s ${reps} --prefix ${params.projtag}_AminoType_Group_Reps_iq -m MFP -madd --redo -nt auto -bb ${params.boots} -bnni

                                else
                                    iqtree -s ${reps} --prefix ${params.projtag}_AminoType_Group_Reps_iq -m MFP -madd --redo -nt auto -bb ${params.boots} -bnni
                                fi
                                """
                            }
                        }

                      process Adding_AminoType_MED_Info {

                          label 'low_cpus'

                          publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/AminoTypes/MED/", mode: "copy", overwrite: true

                          input:
                              file(counts) from aminocountmed
                              file(tree) from amino_repphy
                              file(map) from atygroupscsv

                          output:
                              file("${params.projtag}_AminoType_Groupingcounts.csv") into (amino_groupcounts, amino_groupcounts2)

                          script:
                              """
                              awk -F "," '{print \$1}' ${counts} | sed '1d' > amino.list
                              echo "GroupID" >> group.list
                              for x in \$(cat amino.list);
                              do    group=\$(grep -w \$x ${map} | awk -F "," '{print \$2}')
                                    echo "\$group" >> group.list
                              done
                              paste -d',' group.list ${counts} > ${params.projtag}_AminoType_Groupingcounts.csv
                              """
                    }
        } else {
            atygroupscsv = Channel.value('skipping')
            atygroupscsv2 = Channel.value('skipping')
            amino_group_rep_tree = Channel.value('skipping')
            amino_groupcounts = Channel.value('skipping')
            amino_groupcounts2 = Channel.value('skipping')
        }

            if (params.pcASV) {        // ASV_nucl -> ASV_aa -> clusteraa by %id with ch-hit -> extract representative nucl sequences to generate new OTU file

                process Translation_For_pcASV_Generation {

                      label 'low_cpus'

                      publishDir "${params.workingdir}/${params.outdir}/Analyze/Clustering/pcASV/Translation", mode: "copy", overwrite: true, pattern: '*_ASV_translations*'

                      conda (params.condaActivate ? null : null)

                      container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? null : null)

                      input:
                          file(fasta) from nucl2aa

                      output:
                          file("*ASV*translations.fasta") into clustering_aa
                          file("*_ASV_translations_report") into reportaa_VR
                          file("*_ASV_nucleotide.fasta") into asvfastaforaaclust

                      script:
                          """
                          ${tools}/virtualribosomev2/dna2pep.py ${fasta} -r all -x -o none --fasta ${params.projtag}_ASV_translation.fasta --report ${params.projtag}_ASV_translations_report
                          awk '/^>/ { print (NR==1 ? "" : RS) \$0; next } { printf "%s", \$0 } END { printf RS }' ${params.projtag}_ASV_translation.fasta > ${params.projtag}_ASV_translations.fasta
                          cp ${fasta} ${params.projtag}_ASV_nucleotide.fasta
                          """
                }

                process Generate_pcASVs { //CHECK rename_seq.py in container !!!!!!!!!!!

                    label 'norm_cpus'

                    tag "${mtag}"

                    publishDir "${params.workingdir}/${params.outdir}/Analyze/Clustering/pcASV", mode: "copy", overwrite: true, pattern: '*pcASV*.{fasta}'
                    publishDir "${params.workingdir}/${params.outdir}/Analyze/Clustering/pcASV/SummaryFiles", mode: "copy", overwrite: true, pattern: '*.{clstr,csv,gc}'
                    publishDir "${params.workingdir}/${params.outdir}/Analyze/Clustering/pcASV/Problematic", mode: "copy", overwrite: true, pattern: '*problem*.{fasta}'

                    conda (params.condaActivate ? "-c bioconda -c conda-forge vsearch=2.21.1 cd-hit=4.8.1 seqtk=1.3 bbmap=39.01" : null)

                    container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/mulled-v2-089fb9f3537921c3d6dbcc7521fbc33d82301df5:1e1ccff83e5d9864e7f3c008bd4ece458ffbdb8d-0" : "quay.io/biocontainers/mulled-v2-089fb9f3537921c3d6dbcc7521fbc33d82301df5:1e1ccff83e5d9864e7f3c008bd4ece458ffbdb8d-0")

                    input:
                        each x from 1..naa
                        file(fasta) from clustering_aa
                        file(asvs) from asvfastaforaaclust

                    output:
                        tuple nid, file("${params.projtag}_nucleotide_pcASV*.fasta") into ( pcASV_ntDiamond_ch, pcASV_nt_counts_ch, pcASV_ntmatrix_ch, pcASV_ntmuscle_ch, pcASV_ntmuscle_ch2 )
                        tuple nid, file("*_aminoacid_pcASV*_noTaxonomy.fasta") into ( pcASV_aaMatrix_ch, pcASV_aaDiamond_ch, pcASV_aaMafft_ch, pcASV_aaMafft_ch2, pcASV_aaCounts_ch, pcASVEMBOSS )
                        tuple nid, file("*.fasta"), file("*.clstr"), file("*.csv"), file("*.gc") into ( pcASVsupplementalfiles )

                    script:
                        // add awk script to count seqs
                        nid=slist2.get(x-1)
                        mtag="ID=" + slist2.get(x-1)
                        """
                        set +e
                        cp ${params.vampdir}/bin/rename_seq.py .
                        ####awk 'BEGIN{RS=">";ORS=""}length(\$2)>=${params.minAA}{print ">"\$0}' ${fasta} > ${params.projtag}_filtered_proteins.fasta
                        awk -v RS='>[^\n]+\n' 'length() >= ${params.minAA} {printf "%s", prt \$0} {prt = RT}' ${fasta} > ${params.projtag}_filtered_proteins.fasta
                        cd-hit -i ${params.projtag}_filtered_proteins.fasta -c .${nid} -o ${params.projtag}_pcASV${nid}.fasta
                        sed 's/>Cluster />Cluster_/g' ${params.projtag}_pcASV${nid}.fasta.clstr >${params.projtag}_pcASV${nid}.clstr
                        grep ">Cluster_" ${params.projtag}_pcASV${nid}.clstr >temporaryclusters.list
                        y=\$(grep -c ">Cluster_" ${params.projtag}_pcASV${nid}.clstr)
                        echo ">Cluster_"\${y}"" >> ${params.projtag}_pcASV${nid}.clstr
                        t=1
                        b=1
                        for x in \$(cat temporaryclusters.list);do
                            echo "Extracting \$x"
                            name="\$( echo \$x | awk -F ">" '{print \$2}')"
                            clust="pcASV"\${t}""
                            echo "\${name}"
                            awk '/^>'\${name}'\$/,/^>Cluster_'\${b}'\$/' ${params.projtag}_pcASV${nid}.clstr > "\${name}"_"\${clust}"_tmp.list
                            t=\$(( \${t}+1 ))
                            b=\$(( \${b}+1 ))
                        done

        		        ls *_tmp.list
                        u=1
                        for x in *_tmp.list;do
                            name="\$(echo \$x | awk -F "_p" '{print \$1}')"
                            echo "\${name}"
                            cluster="\$(echo \$x | awk -F "_" '{print \$3}')"
                            echo "\${cluster}"
                            grep "ASV" \$x | awk -F ", " '{print \$2}' | awk -F "_" '{print \$1}' | awk -F ">" '{print \$2}' > \${name}_\${cluster}_seqs_tmps.list
                            seqtk subseq ${asvs} \${name}_\${cluster}_seqs_tmps.list > \${name}_\${cluster}_nucleotide_sequences.fasta
                            vsearch --cluster_fast \${name}_\${cluster}_nucleotide_sequences.fasta --id 0.2 --centroids \${name}_\${cluster}_centroids.fasta
                            grep ">" \${name}_\${cluster}_centroids.fasta >> \${name}_\${cluster}_tmp_centroids.list
                            for y in \$( cat \${name}_\${cluster}_tmp_centroids.list );do
                                echo ">\${cluster}_type"\$u"" >> \${name}_\${cluster}_tmp_centroid.newheaders
                                u=\$(( \${u}+1 ))
                            done
                            u=1
                            ./rename_seq.py \${name}_\${cluster}_centroids.fasta \${name}_\${cluster}_tmp_centroid.newheaders \${cluster}_types_labeled.fasta
                        done
                        cat *_types_labeled.fasta >> ${params.projtag}_nucleotide_pcASV${nid}_noTaxonomy.fasta
                        grep -w "*" ${params.projtag}_pcASV${nid}.clstr | awk '{print \$3}' | awk -F "." '{print \$1}' >tmphead.list
                        grep -w "*" ${params.projtag}_pcASV${nid}.clstr | awk '{print \$2}' | awk -F "," '{print \$1}' >tmplen.list
                        paste -d"," temporaryclusters.list tmphead.list >tmp.info.csv
                        grep ">" ${params.projtag}_pcASV${nid}.fasta >lala.list
                        j=1
                        for x in \$(cat lala.list);do
                            echo ">${params.projtag}_pcASV\${j}" >>${params.projtag}_aminoheaders.list
                            echo "\${x},>${params.projtag}_pcASV\${j}" >>tmpaminotype.info.csv
                            j=\$(( \${j}+1 ))
                        done
                        rm lala.list
                        awk -F "," '{print \$2}' tmp.info.csv >>tmporder.list
                        for x in \$(cat tmporder.list);do
                            grep -w "\$x" tmpaminotype.info.csv | awk -F "," '{print \$2}' >>tmpder.list
                        done
                        paste -d "," temporaryclusters.list tmplen.list tmphead.list tmpder.list >${params.projtag}_pcASVCluster${nid}_summary.csv
                        ./rename_seq.py ${params.projtag}_pcASV${nid}.fasta ${params.projtag}_aminoheaders.list ${params.projtag}_aminoacid_pcASV${nid}_noTaxonomy.fasta
                        stats.sh in=${params.projtag}_aminoacid_pcASV${nid}_noTaxonomy.fasta gc=${params.projtag}_pcASV${nid}_aminoacid_clustered.gc gcformat=4
                        stats.sh in=${params.projtag}_nucleotide_pcASV${nid}_noTaxonomy.fasta gc=${params.projtag}_pcASV${nid}_nucleotide_clustered.gc gcformat=4
                        ####awk 'BEGIN{RS=">";ORS=""}length(\$2)<${params.minAA}{print ">"\$0}' ${fasta} >${params.projtag}_pcASV${nid}_problematic_translations.fasta
                        awk -v RS='>[^\n]+\n' 'length() < ${params.minAA} {printf "%s", prt \$0} {prt = RT}' ${fasta} >${params.projtag}_pcASV${nid}_problematic_translations.fasta
                        if [ `wc -l ${params.projtag}_pcASV${nid}_problematic_translations.fasta | awk '{print \$1}'` -gt 1 ];then
                            grep ">" ${params.projtag}_pcASV${nid}_problematic_translations.fasta | awk -F ">" '{print \$2}' > problem_tmp.list
                            seqtk subseq ${asvs} problem_tmp.list > ${params.projtag}_pcASV${nid}_problematic_nucleotides.fasta
                        else
                           rm ${params.projtag}_pcASV${nid}_problematic_translations.fasta
                        fi
                        rm *.list
                        rm Cluster*
                        rm *types*
                        rm *tmp*
                        rm ${params.projtag}_pcASV${nid}.fast*
                        """
                    }

                if (!params.skipTaxonomy) {

                  if (params.dbtype == "NCBI") { //CHECK rename_seq.py in container !!!!!!!!!!!

                    process pcASV_Nucleotide_Taxonomy_Inference_NCBI {

                        label 'high_cpus'

                        tag "${mtag}"

                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Nucleotide/Taxonomy/SummaryFiles", mode: "copy", overwrite: true, pattern: '*.{csv,tsv}'
                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Nucleotide/Taxonomy/DiamondOutput", mode: "copy", overwrite: true, pattern: '*dmd.{out}'
                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Nucleotide/Taxonomy", mode: "copy", overwrite: true, pattern: '*.{fasta}'

                        conda (params.condaActivate ? "-c conda-forge bioconda::diamond=2.0.15=hb97b32f_1" : null)

                        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/diamond:2.0.15--hb97b32f_1" : "quay.io/biocontainers/diamond:2.0.15--hb97b32f_1")

                        input:
                            tuple nid, file(asvs) from pcASV_ntDiamond_ch

                        output:
                            file("*.fasta") into ( pcASV_labeled )
                            tuple file("*_phyloformat.csv"), file("*_summaryTable.tsv"), file("*dmd.out") into summary_AAdiamond
                            tuple nid, file("*_summary_for_plot.csv") into taxplot3
                            tuple nid, file("*_quick_Taxbreakdown.csv") into tax_table_pcasvnt
                            tuple nid, file ("*_quicker_taxbreakdown.csv") into tax_nodCol_pcasvnt

                        script:
                            mtag="ID=" + nid
                            """
                            set +e
                            cp ${params.vampdir}/bin/rename_seq.py .
                            virdb=${params.dbdir}/${params.dbname}
                            if [[ ${params.measurement} == "bitscore" ]]
                            then    measure="--min-score ${params.bitscore}"
                            elif    [[ ${params.measurement} == "evalue" ]]
                            then    measure="-e ${params.evalue}"
                            else    measure="--min-score ${params.bitscore}"
                            fi
                            grep ">" \${virdb} > headers.list
                            headers="headers.list"
                            name=\$( echo ${asvs} | awk -F ".fasta" '{print \$1}')
                            if [[ ${params.ncbitax} == "true" ]]
                            then   diamond blastx -q ${asvs} -d \${virdb} -p ${task.cpus} --id ${params.minID} -l ${params.minaln} \${measure} --${params.sensitivity} -o "\$name"_dmd.out -f 6 qseqid qlen sseqid qstart qend qseq sseq length qframe evalue bitscore pident btop staxids sskingdoms skingdoms sphylums --max-target-seqs 1 --max-hsps 1
                            else   diamond blastx -q ${asvs} -d \${virdb} -p ${task.cpus} --id ${params.minID} -l ${params.minaln} \${measure} --${params.sensitivity} -o "\$name"_dmd.out -f 6 qseqid qlen sseqid qstart qend qseq sseq length qframe evalue bitscore pident btop --max-target-seqs 1 --max-hsps 1
                            fi
                            echo "Preparing lists to generate summary .csv's"
                            echo "[Best hit accession number]" > access.list
                            echo "[e-value]" > evalue.list
                            echo "[Bitscore]" > bit.list
                            echo "[Percent ID (aa)]" > pid.list
                            echo "[Organism ID]" > "\$name"_virus.list
                            echo "[Gene]" > "\$name"_genes.list
                            echo "[pcASV#]" > otu.list
                            echo "[Sequence length]" > length.list
                            grep ">" ${asvs} | awk -F ">" '{print \$2}' > seqids.lst
                            if [[ ${params.lca} == "T" ]]
                            then  grep -w "LCA" ${params.dbanno}/*.txt > lcainfo.list
                                  echo "[Taxonomic classification from RVDB annotations]" > lca_classification.list
                            else
                                  echo "[Taxonomic classification from RVDB annotations]" > lca_classification.list
                            fi
                            if [[ ${params.ncbitax} == "true" ]]
                            then echo "[NCBI Taxonomy ID],[Taxonomic classification from NCBI]" > ncbi_classification.list
                            fi
                            echo "extracting genes and names"
                            touch new_"\$name"_asvnames.txt
                            for s in \$(cat seqids.lst);do
                                echo "Checking for \$s hit in diamond output"
                                if [[ "\$(grep -wc "\$s" "\$name"_dmd.out)" -eq 1 ]];then
                                            echo "Yep, there was a hit for \$s"
                                            echo "Extracting the information now:"
                                            acc=\$(grep -w "\$s" "\$name"_dmd.out | awk '{print \$3}')
                                            echo "\$s" >> otu.list
                                            echo "\$acc" >> access.list
                                            line="\$(grep -w "\$s" "\$name"_dmd.out)"
                                            echo "\$line" | awk '{print \$10}' >> evalue.list
                                            echo "\$line" | awk '{print \$11}' >> bit.list
                                            echo "\$line" | awk '{print \$12}' >> pid.list
                                            echo "\$line" | awk '{print \$2}' >> length.list
                                            echo "Extracting virus and gene ID for \$s now"
                                            gene=\$(grep -w "\$acc" "\$headers" | awk -F "." '{ print \$2 }' | awk -F "[" '{ print \$1 }' | awk -F " " '{print substr(\$0, index(\$0,\$2))}' | sed 's/ /_/g') &&
                                            echo "\$gene" | sed 's/_/ /g' >> "\$name"_genes.list
                                            virus=\$(grep -w "\$acc" "\$headers" | awk -F "[" '{ print \$2 }' | awk -F "]" '{ print \$1 }'| sed 's/ /_/g')
                                            echo "\$virus" | sed 's/_/ /g' >> "\$name"_virus.list
                                            echo ">"\${s}"_"\$virus"_"\$gene"" >> new_"\$name"_asvnames.txt
                                            if [[ "${params.lca}" == "T" ]]
                                            then    if [[ \$(grep -w "\$acc" ${params.dbanno}/*.txt | wc -l) -eq 1 ]]
                                                    then  group=\$(grep -w "\$acc" ${params.dbanno}/*.txt | awk -F ":" '{print \$1}')
                                                          lcla=\$(grep -w "\$group" lcainfo.list | awk -F "\t" '{print \$2}')
                                                          echo "\$lcla" >> lca_classification.list
                                                    else  echo "Viruses" >> lca_classification.list
                                                    fi
                                            fi
                                            if [[ ${params.ncbitax} == "true" ]]
                                            then  echo "\$line" | awk -F "\t" '{print \$14","\$16"::"\$18"::"\$17}' >> ncbi_classification.list
                                            fi
                                            echo "\$s done."
                                else
                                            echo "Ugh, there was no hit for \$s .."
                                            echo "We still love \$s though and we will add it to the final fasta file"
                                            echo "\$s" >> otu.list
                                            echo "NO_HIT" >> access.list
                                            echo "NO_HIT" >> "\$name"_genes.list
                                            echo "NO_HIT" >> "\$name"_virus.list
                                            echo "NO_HIT" >> evalue.list
                                            echo "NO_HIT" >> bit.list
                                            echo "NO_HIT" >> pid.list
                                            echo "NO_HIT" >> length.list
                                            virus="NO"
                                            gene="HIT"
                                            echo ">\${s}_"\$virus"_"\$gene"" >> new_"\$name"_asvnames.txt
                                            if [[ "${params.lca}" == "T" ]]
                                            then    echo "N/A" >> lca_classification.list
                                            fi
                                            if [[ "${params.ncbitax}" == "true" ]]
                                            then  echo "N/A" >> ncbi_classification.list
                                            fi
                                            echo "\$s done."
                               fi
                            done
                            echo "Now editing "\$name" fasta headers"
                            ###### rename_seq.py
                            ./rename_seq.py ${asvs} new_"\$name"_asvnames.txt "\$name"_TaxonomyLabels.fasta
                            awk 'BEGIN {RS=">";FS="\\n";OFS=""} NR>1 {print ">"\$1; \$1=""; print}' "\$name"_TaxonomyLabels.fasta >"\$name"_tmpssasv.fasta
                            echo "[Sequence header]" > newnames.list
                            cat new_"\$name"_asvnames.txt >> newnames.list
                            touch sequence.list
                            echo "     " > sequence.list
                            grep -v ">" "\$name"_tmpssasv.fasta >> sequence.list
                            rm "\$name"_tmpssasv.fasta
                            if [[ "${params.lca}" == "T" && "${params.ncbitax}" == "true" ]]
                            then
                                  paste -d "," sequence.list "\$name"_virus.list "\$name"_genes.list otu.list lca_classification.list ncbi_classification.list newnames.list length.list bit.list evalue.list pid.list access.list >> "\$name"_phyloformat.csv
                                  paste -d"\t" otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list ncbi_classification.list sequence.list length.list bit.list evalue.list pid.list >> "\$name"_summaryTable.tsv
                                  paste -d"," otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list ncbi_classification.list >> \${name}_quick_Taxbreakdown.csv
                            elif [[ "${params.lca}" == "T" && "${params.ncbitax}" != "true" ]]
                            then
                                  paste -d "," sequence.list "\$name"_virus.list "\$name"_genes.list otu.list lca_classification.list newnames.list length.list bit.list evalue.list pid.list access.list >> "\$name"_phyloformat.csv
                                  paste -d"\t" otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list sequence.list length.list bit.list evalue.list pid.list >> "\$name"_summaryTable.tsv
                                  paste -d"," otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list >> \${name}_quick_Taxbreakdown.csv
                            elif [[ "${params.ncbitax}" == "true" && "${params.lca}" != "T" ]]
                            then
                                  paste -d "," sequence.list "\$name"_virus.list "\$name"_genes.list otu.list ncbi_classification.list newnames.list length.list bit.list evalue.list pid.list access.list >> "\$name"_phyloformat.csv
                                  paste -d"\t" otu.list access.list "\$name"_virus.list "\$name"_genes.list ncbi_classification.list sequence.list length.list bit.list evalue.list pid.list >> "\$name"_summaryTable.tsv
                                  paste -d"," otu.list access.list "\$name"_virus.list "\$name"_genes.list ncbi_classification.list >> \${name}_quick_Taxbreakdown.csv
                            else
                                  paste -d "," sequence.list "\$name"_virus.list "\$name"_genes.list otu.list newnames.list length.list bit.list evalue.list pid.list access.list >> "\$name"_phyloformat.csv
                                  paste -d"\t" otu.list access.list "\$name"_virus.list "\$name"_genes.list sequence.list length.list bit.list evalue.list pid.list >> "\$name"_summaryTable.tsv
                                  echo "skipped" >> \${name}_quick_Taxbreakdown.csv
                            fi
                            for x in *phyloformat.csv;do
                                echo "\$x"
                                lin=\$(( \$(wc -l \$x | awk '{print \$1}')-1))
                                tail -"\$lin" \$x | awk -F "," '{print \$2}' > tmpcol.list;
                                sed 's/ /_/g' tmpcol.list > tmp2col.list;
                                cat tmp2col.list | sort | uniq -c | sort -nr | awk '{print \$2","\$1}' > \${name}_summary_for_plot.csv;
                                rm tmpcol.list tmp2col.list
                            done
                            awk -F "," '{print \$1","\$3"("\$2")"}' \${name}_quick_Taxbreakdown.csv >> \${name}_quicker_taxbreakdown.csv
                            rm evalue.list sequence.list bit.list pid.list length.list seqids.lst otu.list *asvnames.txt "\$name"_virus.list "\$name"_genes.list newnames.list access.list headers.list
                            """
                        }
                    } else if (params.dbtype== "RVDB") {

                      process pcASV_Nucleotide_Taxonomy_Inference_RVDB { //CHECK rename_seq.py in container !!!!!!!!!!!

                          label 'high_cpus'

                          tag "${mtag}"

                          publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Nucleotide/Taxonomy/SummaryFiles", mode: "copy", overwrite: true, pattern: '*.{csv,tsv}'
                          publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Nucleotide/Taxonomy/DiamondOutput", mode: "copy", overwrite: true, pattern: '*dmd.{out}'
                          publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Nucleotide/Taxonomy", mode: "copy", overwrite: true, pattern: '*.{fasta}'

                          conda (params.condaActivate ? "-c conda-forge bioconda::diamond=2.0.15=hb97b32f_1" : null)

                          container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/diamond:2.0.15--hb97b32f_1" : "quay.io/biocontainers/diamond:2.0.15--hb97b32f_1")

                          input:
                              tuple nid, file(asvs) from pcASV_ntDiamond_ch

                          output:
                              file("*.fasta") into ( pcASV_labeled )
                              tuple file("*_phyloformat.csv"), file("*_summaryTable.tsv"), file("*dmd.out") into summary_AAdiamond
                              tuple nid, file("*_summary_for_plot.csv") into taxplot3
                              tuple nid, file("*_quick_Taxbreakdown.csv") into tax_table_pcasvnt
                              tuple nid, file ("*_quicker_taxbreakdown.csv") into tax_nodCol_pcasvnt

                          script:
                              mtag="ID=" + nid
                              """
                              set +e
                              cp ${params.vampdir}/bin/rename_seq.py .
                              virdb=${params.dbdir}/${params.dbname}
                              if [[ ${params.measurement} == "bitscore" ]]
                              then    measure="--min-score ${params.bitscore}"
                              elif    [[ ${params.measurement} == "evalue" ]]
                              then    measure="-e ${params.evalue}"
                              else    measure="--min-score ${params.bitscore}"
                              fi
                              grep ">" \${virdb} > headers.list
                              headers="headers.list"
                              name=\$( echo ${asvs} | awk -F ".fasta" '{print \$1}')
                              diamond blastx -q ${asvs} -d \${virdb} -p ${task.cpus} --id ${params.minID} -l ${params.minaln} \${measure} --${params.sensitivity} -o "\$name"_dmd.out -f 6 qseqid qlen sseqid qstart qend qseq sseq length qframe evalue bitscore pident btop --max-target-seqs 1 --max-hsps 1
                              echo "Preparing lists to generate summary .csv's"
                              echo "[Best hit accession number]" > access.list
                              echo "[e-value]" > evalue.list
                              echo "[Bitscore]" > bit.list
                              echo "[Percent ID (aa)]" > pid.list
                              echo "[Organism ID]" > "\$name"_virus.list
                              echo "[Gene]" > "\$name"_genes.list
                              echo "[pcASV#]" > otu.list
                              echo "[Sequence length]" > length.list
                              grep ">" ${asvs} | awk -F ">" '{print \$2}' > seqids.lst
                              if [[ ${params.lca} == "T" ]]
                              then  grep -w "LCA" ${params.dbanno}/*.txt > lcainfo.list
                                    echo "[Taxonomic classification from RVDB annotations]" > lca_classification.list
                              else  echo "skipped" >> \${name}_quick_Taxbreakdown.csv
                                    echo "[Taxonomic classification from RVDB annotations]" > lca_classification.list
                              fi
                              echo "extracting genes and names"
                              touch new_"\$name"_asvnames.txt
                              for s in \$(cat seqids.lst);do
                                  echo "Using RVDB headers."
                                  if [[ "\$(grep -wc "\$s" "\$name"_dmd.out)" -eq 1 ]];then
                                      echo "Yep, there was a hit for \$s"
                                      echo "Extracting the information now:"
                                      acc=\$(grep -w "\$s" "\$name"_dmd.out | awk '{print \$3}' | awk -F "|" '{print \$3}')
                                      echo "\$s" >> otu.list
                                      echo "\$acc" >> access.list
                                      line="\$(grep -w "\$s" "\$name"_dmd.out)"
                                      echo "\$line" | awk '{print \$10}' >> evalue.list
                                      echo "\$line" | awk '{print \$11}' >> bit.list
                                      echo "\$line" | awk '{print \$12}' >> pid.list
                                      echo "\$line" | awk '{print \$2}' >> length.list
                                      echo "Extracting virus and gene ID for \$s now"
                                      gene=\$(grep -w "\$acc" "\$headers" | awk -F "|" '{ print \$6 }' | awk -F "[" '{ print \$1 }' | sed 's/ /_/g') &&
                                      echo "\$gene" | sed 's/_/ /g' >> "\$name"_genes.list
                                      virus=\$(grep -w "\$acc" "\$headers" | awk -F "|" '{ print \$6 }' | awk -F "[" '{ print \$2 }' | awk -F "]" '{print \$1}' | sed 's/ /_/g') &&
                                      echo "\$virus" | sed 's/_/ /g' >> "\$name"_virus.list
                                      echo ">\${s}_"\$virus"_"\$gene"" >> new_"\$name"_asvnames.txt
                                      if [[ "${params.lca}" == "T" ]]
                                      then    if [[ \$(grep -w "\$acc" ${params.dbanno}/*.txt | wc -l) -eq 1 ]]
                                              then  group=\$(grep -w "\$acc" ${params.dbanno}/*.txt | awk -F ":" '{print \$1}')
                                                    lcla=\$(grep -w "\$group" lcainfo.list | awk -F "\t" '{print \$2}')
                                                    echo "\$lcla" >> lca_classification.list
                                              else  echo "Viruses" >> lca_classification.list
                                              fi
                                      fi
                                      echo "\$s done."
                                  else
                                      echo "Ugh, there was no hit for \$s .."
                                      echo "We still love \$s though and we will add it to the final fasta file"
                                      echo "\$s" >> otu.list
                                      echo "NO_HIT" >> access.list
                                      echo "NO_HIT" >> "\$name"_genes.list
                                      echo "NO_HIT" >> "\$name"_virus.list
                                      echo "NO_HIT" >> evalue.list
                                      echo "NO_HIT" >> bit.list
                                      echo "NO_HIT" >> pid.list
                                      echo "NO_HIT" >> length.list
                                      virus="NO"
                                      gene="HIT"
                                      echo ">\${s}_"\$virus"_"\$gene"" >> new_"\$name"_asvnames.txt
                                      if [[ "${params.lca}" == "T" ]]
                                      then    echo "N/A" >> lca_classification.list
                                      fi
                                      echo "\$s done."
                                  fi
                              echo "Done with \$s"
                              done
                              echo "Now editing "\$name" fasta headers"
                              ###### rename_seq.py
                              ./rename_seq.py ${asvs} new_"\$name"_asvnames.txt "\$name"_TaxonomyLabels.fasta
                              awk 'BEGIN {RS=">";FS="\\n";OFS=""} NR>1 {print ">"\$1; \$1=""; print}' "\$name"_TaxonomyLabels.fasta >"\$name"_tmpssasv.fasta
                              echo "[Sequence header]" > newnames.list
                              cat new_"\$name"_asvnames.txt >> newnames.list
                              touch sequence.list
                              echo "     " > sequence.list
                              grep -v ">" "\$name"_tmpssasv.fasta >> sequence.list
                              rm "\$name"_tmpssasv.fasta
                              if [[ "${params.lca}" == "T" ]]
                              then  paste -d "," sequence.list "\$name"_virus.list "\$name"_genes.list otu.list lca_classification.list newnames.list length.list bit.list evalue.list pid.list access.list >> "\$name"_phyloformat.csv
                                    paste -d"\t" otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list sequence.list length.list bit.list evalue.list pid.list >> "\$name"_summaryTable.tsv
                                    paste -d"," otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list >> \${name}_quick_Taxbreakdown.csv
                              else  paste -d "," sequence.list "\$name"_virus.list "\$name"_genes.list otu.list newnames.list length.list bit.list evalue.list pid.list access.list >> "\$name"_phyloformat.csv
                                    paste -d"\t" otu.list access.list "\$name"_virus.list "\$name"_genes.list sequence.list length.list bit.list evalue.list pid.list >> "\$name"_summaryTable.tsv
                              fi
                              for x in *phyloformat.csv;do
                                        echo "\$x"
                                        lin=\$(( \$(wc -l \$x | awk '{print \$1}')-1))
                                        tail -"\$lin" \$x | awk -F "," '{print \$2}' > tmpcol.list;
                                        sed 's/ /_/g' tmpcol.list > tmp2col.list;
                                        cat tmp2col.list | sort | uniq -c | sort -nr | awk '{print \$2","\$1}' > \${name}_summary_for_plot.csv;
                                        rm tmpcol.list tmp2col.list
                              done
                              awk -F "," '{print \$1","\$3"("\$2")"}' \${name}_quick_Taxbreakdown.csv >> \${name}_quicker_taxbreakdown.csv
                              rm evalue.list sequence.list bit.list pid.list length.list seqids.lst otu.list *asvnames.txt "\$name"_virus.list "\$name"_genes.list newnames.list access.list headers.list
                              """
                    }
                }
             } else {

                 process skippcASVnuctaxonomy {

                     input:
                         tuple nid, file(asvs) from pcASV_ntDiamond_ch

                     output:
                         tuple nid, file("skipncASVnubtaxonomy1.txt") into ( taxplot3 )
                         tuple nid, file("skipncASVnubtaxonomy2.txt") into ( tax_table_pcasvnt )
                         tuple nid, file("skipncASVnubtaxonomy3.txt") into ( tax_nodCol_pcasvnt )

                     script:
                         """
                         echo "Skipped" >skipncASVnubtaxonomy1.txt
                         echo "Skipped" >skipncASVnubtaxonomy2.txt
                         echo "Skipped" >skipncASVnubtaxonomy3.txt
                         """
                 }
             }

                process Generate_Nucleotide_pcASV_Counts {

                    label 'norm_cpus'

                    tag "${mtag}"

                    publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Nucleotide/Counts", mode: "copy", overwrite: true, pattern: '*.{biome,csv,txt}'

                    conda (params.condaActivate ? "-c bioconda -c conda-forge vsearch=2.21.1=hf1761c0_1" : null)

                    container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/vsearch:2.21.1--hf1761c0_1" : "quay.io/biocontainers/vsearch:2.21.1--hf1761c0_1")

                    input:
                        tuple nid, file(potus) from pcASV_nt_counts_ch
                        file(merged) from pcASV_mergedreads_ch

                    output:
                        tuple file("*_counts.txt"), file("*_counts.biome") into pcASVcounts_vsearch
                        tuple nid, file("*.csv") into potu_Ncounts_for_report

                    script:
                        mtag="ID=" + nid
                        """
                    	name=\$( echo ${potus} | awk -F ".fasta" '{print \$1}')
                    	vsearch --usearch_global ${merged} --db ${potus} --id .${params.npcasvcountID} --threads ${task.cpus} --otutabout \${name}_counts.txt --biomout \${name}_counts.biome
                    	cat \${name}_counts.txt | tr "\t" "," >\${name}_count.csv
                    	sed 's/#OTU ID/OTU_ID/g' \${name}_count.csv >\${name}_counts.csv
                    	rm \${name}_count.csv
                        """
                }

                process Generate_pcASV_Nucleotide_Matrix {

                    label 'low_cpus'

                    tag "${mtag}"

                    publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Nucleotide/Matrix", mode: "copy", overwrite: true

                    conda (params.condaActivate ? "-c bioconda -c conda-forge clustalo=1.2.4=h87f3376_5" : null)

                    container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/clustalo:1.2.4--h87f3376_5" : "quay.io/biocontainers/clustalo:1.2.4--h87f3376_5")

                    input:
                        tuple nid, file(potus) from pcASV_ntmatrix_ch

                    output:
                        file("*.matrix") into pcASVclustmatrices
                        tuple nid, file("*PercentID.matrix") into potu_nucl_heatmap

                    script:
                        //check --percent-id second clustalo
                        mtag="ID=" + nid
                        """
                        name=\$( echo ${potus} | awk -F ".fasta" '{print \$1}')
                        clustalo -i ${potus} --distmat-out=\${name}_PairwiseDistanceq.matrix --full --force --threads=${task.cpus}
                        clustalo -i ${potus} --distmat-out=\${name}_PercentIDq.matrix --percent-id --full --force --threads=${task.cpus}
                        cat \${name}_PercentIDq.matrix | tr " " "," | sed 's/,,/,/g' | grep "," >\${name}_PercentID.matrix
                        rm \${name}_PercentIDq.matrix
                        """
                }

                if (!params.skipPhylogeny) {

                    process pcASV_Nucleotide_Phylogeny_step1 {

                        label 'norm_cpus'

                        tag "${mtag}"

                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Nucleotide/Phylogeny/Alignment", mode: "copy", overwrite: true, pattern: '*aln.*'

                        conda (params.condaActivate ? "-c bioconda -c conda-forge muscle=5.1" : null)

                        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/muscle:5.1--h9f5acd7" : "quay.io/biocontainers/muscle:5.1--h9f5acd7")

                        input:
                            tuple nid, file(pcASVn) from pcASV_ntmuscle_ch

                        output:
                            tuple nid, file("*_ALN.fasta") into potu_align1
                            file("*.efa")

                        script:
                            mtag="ID=" + nid
                            """
                            pre=\$(echo ${pcASVn} | awk -F ".fasta" '{print \$1}' )

                            if [[ ${params.srep} == "true" && ${params.ensemble} == "false" ]];
                            then  if [[ \$( grep -c ">" ${pcASVn}) -lt 300 ]]
                                  then    comm="align"
                                  else    comm="super5"
                                  fi
                                  muscle -"\$comm" ${pcASVn} -perm ${params.perm} -perturb ${params.pert} -output \${pre}_muscle_raw_ALN.fasta -threads ${task.cpus} -quiet
                                  echo "single replicate alignment chosen; look over muscle5 documentation to learn about ensemble alignment approach" >note.efa
                            elif [[ ${params.srep} == "false" && ${params.ensemble} == "true" ]];
                            then  muscle -align ${pcASVn} -${params.fied} -output \${pre}_muscle.efa -threads ${task.cpus} -quiet
                                  muscle -maxcc \${pre}_muscle.efa -output \${pre}_muscle_raw_ALN.fasta
                            else  if [[ \$( grep -c ">" ${pcASVn}) -lt 300 ]]
                                  then    comm="align"
                                  else    comm="super5"
                                  fi
                                  muscle -"\$comm" ${pcASVn} -output \${pre}_muscle_raw_ALN.fasta -threads ${task.cpus} -quiet
                                  echo "single replicate alignment chosen; look over muscle5 documentation to learn about ensemble alignment approach" >note.efa
                            fi
                            """
                    }


                    process pcASV_Nucleotide_Phylogeny_step2 {

                        label 'norm_cpus'

                        tag "${mtag}"

                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Nucleotide/Phylogeny/Alignment", mode: "copy", overwrite: true

                        conda (params.condaActivate ? "-c conda-forge bioconda::trimal=1.4.1=h9f5acd7_6" : null)

                        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/trimal:1.4.1--h9f5acd7_6" : "quay.io/biocontainers/trimal:1.4.1--h9f5acd7_6")

                        input:
                            tuple nid, file(pcASVn) from potu_align1

                        output:
                            tuple nid, file("*_aln.html") into pcASV_nucleotide_phylogeny_results1
                            tuple nid, file("*_aln.fasta") into potu_align2

                        script:
                            mtag="ID=" + nid
                            """
                            pre=\$( echo ${pcASVn} | awk -F "_muscle" '{print \$1}' )
                            trimal -in ${pcASVn} -out \${pre}_trimal_aln.fasta -keepheader -fasta -automated1 -htmlout \${pre}_aln.html
                            """
                    }

                    process pcASV_Nucleotide_Phylogeny_step3 {

                        label 'norm_cpus'

                        tag "${mtag}"

                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Nucleotide/Phylogeny/Alignment", mode: "copy", overwrite: true

                        conda (params.condaActivate ? "${params.vampdir}/bin/yamls/oligotyping.yml" : null)

                        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/oligotyping:2.1--py27_0" : "quay.io/biocontainers/oligotyping:2.1--py27_0")

                        input:
                            tuple nid, file(pcASVn) from potu_align2

                        output:
                            tuple nid, file("*_Aligned_informativeonly.fasta") into potu_align3

                        script:
                            mtag="ID=" + nid
                            """
                            pre=\$( echo ${pcASVn} | awk -F "_trimal" '{print \$1}' )
                            o-trim-uninformative-columns-from-alignment ${pcASVn}
                            mv ${pcASVn}-TRIMMED ./\${pre}_Aligned_informativeonly.fasta
                            """
                    }

                    process pcASV_Nucleotide_Phylogeny_step4 {

                        label 'norm_cpus'

                        tag "${mtag}"

                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Nucleotide/Phylogeny/ModelTest", mode: "copy", overwrite: true, pattern: '*iq*'

                        conda (params.condaActivate ? "-c conda-forge bioconda::modeltest-ng=0.1.7=h5c6ebe3_0" : null)

                        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/modeltest-ng:0.1.7--h5c6ebe3_0" : "quay.io/biocontainers/modeltest-ng:0.1.7--h5c6ebe3_0")

                        input:
                            tuple nid, file(pcASVn) from potu_align3

                        output:
                            tuple nid, file("*mt*") into pcASV_nucleotide_phylogeny_results2
                            tuple nid, file("*mt.out") into pcASVmtout
                        script:
                            mtag="ID=" + nid
                            """
                            pre=\$( echo ${pcASVn} | awk -F "_Aligned" '{print \$1}' )
                            # pcASV_Nucleotide_ModelTest
                            modeltest-ng -i ${pcASVn} -p ${task.cpus} -o \${pre}_mt -d nt -s 203 --disable-checkpoint
                            """
                    }

                    process pcASV_Nucleotide_Phylogeny_step5 {

                        label 'norm_cpus'

                        tag "${mtag}"

                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Nucleotide/Phylogeny/IQ-TREE", mode: "copy", overwrite: true, pattern: '*iq*'

                        conda (params.condaActivate ? "-c conda-forge bioconda::iqtree=2.2.0.3=hb97b32f_1" : null)

                        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/iqtree:2.2.0.3--hb97b32f_1" : "quay.io/biocontainers/iqtree:2.2.0.3--hb97b32f_1")

                        input:
                            tuple nid, file(pcASVn) from pcASV_ntmuscle_ch2
                            tuple nid, file(mtout) from pcASVmtout

                        output:
                            tuple nid, file("*.tree"), file("*.log"), file("*iq*") into pcASV_nucleotide_phylogeny_results3
                            tuple nid, file("*iq.treefile") into potu_Ntree_plot

                        script:
                            mtag="ID=" + nid
                            """
                            #grabbing best models from modeltestng
                            modbic=\$(grep "iqtree" ${mtout} | head -1 | awk -F "-m " '{print \$2}')
                            modaic=\$(grep "iqtree" ${mtout} | head -2 | tail -1 | awk -F "-m " '{print \$2}')
                            modaicc=\$(grep "iqtree" ${mtout} | tail -1 | awk -F "-m " '{print \$2}')
                            if [[ ${params.crit} == "BIC" ]]
                            then  mod="\$modbic"
                            elif  [[ ${params.crit} == "AIC" ]]
                            then  mod="\$modaic"
                            elif  [[ ${params.crit} == "AICc" ]]
                            then  mod="\$modaicc"
                            fi
                            # grab prefix
                            pre=\$( echo ${pcASVn} | awk -F "_Aligned" '{print \$1}' )
                            # pcASV_Nucleotide_Phylogeny
                            if [ "${params.iqCustomnt}" != "" ];then
                                iqtree -s ${pcASVn} --prefix \${pre}_noTaxonomy_iq --redo -T auto ${params.iqCustomnt}
                            elif [[ "${params.ModelTnt}" != "false" && "${params.nonparametric}" != "false" ]];then
                                iqtree -s ${pcASVn} --prefix \${pre}_noTaxonomy_iq -m \${mod} --redo-nt auto -b ${params.boots}
                            elif [[ "${params.ModelTnt}" != "false" && "${params.parametric}" != "false" ]];then
                                iqtree -s ${pcASVn} --prefix \${pre}_noTaxonomy_iq -m \${mod} --redo -nt auto -bb ${params.boots} -bnni
                            elif [ "${params.nonparametric}" != "false" ];then
                                iqtree -s ${pcASVn} --prefix \${pre}_noTaxonomy_iq -m MFP -madd --redo -nt auto -b ${params.boots}
                            elif [ "${params.parametric}" != "false" ];then
                                iqtree -s ${pcASVn} --prefix \${pre}_noTaxonomy_iq -m MFP -madd --redo -nt auto -bb ${params.boots} -bnni
                            else
                                iqtree -s ${pcASVn} --prefix \${pre}_noTaxonomy_iq -m MFP -madd --redo -nt auto -bb ${params.boots} -bnni
                            fi
                            """
                    }

                } else {

                    process skippcASVnucphylogeny {

                        input:
                            tuple nid, file(prot) from pcASV_ntmuscle_ch

                        output:
                            tuple nid, file("skipncASVnucphy.txt") into ( potu_Ntree_plot )

                        script:
                            """
                            echo "Skipped" >skipncASVnucphy.txt
                            """
                    }
                }

                process pcASV_AminoAcid_Matrix {

                    label 'low_cpus'

                    tag "${mtag}"

                    publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Aminoacid/Matrix", mode: "copy", overwrite: true

                    conda (params.condaActivate ? "-c bioconda -c conda-forge clustalo=1.2.4=h87f3376_5" : null)

                    container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/clustalo:1.2.4--h87f3376_5" : "quay.io/biocontainers/clustalo:1.2.4--h87f3376_5")

                    input:
                        tuple nid, file(prot) from pcASV_aaMatrix_ch

                    output:
                        file("*.matrix") into pcASVaaMatrix
                        tuple nid, file("*PercentID.matrix") into potu_aa_heatmap

                    script:
                        mtag="ID=" + nid
                        """
                        name=\$( echo ${prot} | awk -F ".fasta" '{print \$1}')
                        clustalo -i ${prot} --distmat-out=\${name}_PairwiseDistanceq.matrix --full --force --threads=${task.cpus}
                        clustalo -i ${prot} --distmat-out=\${name}_PercentIDq.matrix --percent-id --full --force --threads=${task.cpus}
                        cat \${name}_PercentIDq.matrix | tr " " "," | sed 's/,,/,/g' | grep "," >\${name}_PercentID.matrix
                        rm \${name}_PercentIDq.matrix
                        """
                }

                if (!params.skipEMBOSS) {

                    process pcASV_EMBOSS_Analyses {

                        label 'low_cpus'

                        tag "${mtag}"

                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Aminoacid/EMBOSS/2dStructure", mode: "copy", overwrite: true, pattern: '*.{garnier}'
                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Aminoacid/EMBOSS/HydrophobicMoment", mode: "copy", overwrite: true, pattern: '*HydrophobicMoments*'
                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Aminoacid/EMBOSS/IsoelectricPoint", mode: "copy", overwrite: true, pattern: '*IsoelectricPoint.{iep,svg}'
                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Aminoacid/EMBOSS/ProteinProperties", mode: "copy", overwrite: true, pattern: '*.{pepstats,pepinfo}'
                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Aminoacid/EMBOSS/ProteinProperties/Plots", mode: "copy", overwrite: true, pattern: '*PropertiesPlot.{svg}'
                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Aminoacid/EMBOSS/2dStructure/Plots", mode: "copy", overwrite: true, pattern: '*Helical*.{svg}'

                        conda (params.condaActivate ? "-c bioconda -c conda-forge emboss=6.6.0=h86d058a_5" : null)

                        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/emboss:6.6.0--h86d058a_5" : "quay.io/biocontainers/emboss:6.6.0--h86d058a_5")

                        input:
                            tuple nid, file(prot) from pcASVEMBOSS

                        output:
                            tuple file("*.garnier"), file("*HydrophobicMoments.svg"), file("*IsoelectricPoint*"), file("*.pepstats"), file("*PropertiesPlot*"), file("*Helical*")  into pcASV_emboss

                        script:
                            // check do I need for loop
                            mtag="ID=" + nid
                            """
                            name=\$( echo ${prot} | awk -F ".fasta" '{print \$1}')
                            garnier -sequence ${prot} -outfile \${name}_2dStructures.garnier
                            hmoment -seqall ${prot} -graph svg -plot -double
                            hmoment -seqall ${prot} -double -outfile ${prot}_HydrophobicMoments
                            mv hmoment.svg ./"\${name}"_HydrophobicMoments.svg
                            iep -sequence ${prot} -graph svg -plot -outfile "\${name}"_IsoelectricPoint.iep
                            mv iep.svg ./"\${name}"_IsoelectricPoint.svg
                            pepstats -sequence ${prot} -outfile \${name}_ProteinProperties.pepstats
                            grep ">" ${prot} | awk -F ">" '{print \$2}' > tmpsequence.list
                            for x in \$(cat tmpsequence.list);do
                                grep -A 1 ""\$x"" ${prot} > tmp2.fasta
                                len=\$(tail -1 tmp2.fasta | awk '{print length}')
                                pepinfo -sequence tmp2.fasta -graph svg -outfile "\$x"_PropertiesPlot.pepinfo
                                mv pepinfo.svg ./"\$x"_PropertiesPlot.svg
                                cat "\$x"_PropertiesPlot.pepinfo >> "\${name}"_PropertiesPlot.pepinfo
                                rm "\$x"_PropertiesPlot.pepinfo
                                pepnet -sask -sequence tmp2.fasta -graph svg -sbegin1 1 -send1 \$len
                                mv pepnet.svg ./"\$x"_HelicalNet.svg
                                pepwheel -sequence tmp2.fasta -graph svg -sbegin1 1 -send1 \$len
                                mv pepwheel.svg ./"\$x"_HelicalWheel.svg
                                rm tmp2.fasta
                            done
                            rm tmpsequence.list
                            """
                        }
                }

                if (!params.skipTaxonomy) {

                  if (params.dbtype == "NCBI") {

                    process pcASV_AminoAcid_Taxonomy_Inference_NCBI { //CHECK rename_seq.py in container !!!!!!!!!!!

                        label 'high_cpus'

                        tag "${mtag}"

                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Aminoacid/Taxonomy/SummaryFiles", mode: "copy", overwrite: true, pattern: '*.{csv,tsv}'
                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Aminoacid/Taxonomy/DiamondOutput", mode: "copy", overwrite: true, pattern: '*dmd.{out}'
                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Aminoacid/Taxonomy", mode: "copy", overwrite: true, pattern: '*.{fasta}'

                        conda (params.condaActivate ? "-c conda-forge bioconda::diamond=2.0.15=hb97b32f_1" : null)

                        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/diamond:2.0.15--hb97b32f_1" : "quay.io/biocontainers/diamond:2.0.15--hb97b32f_1")

                        input:
                            tuple nid, file(asvs) from pcASV_aaDiamond_ch

                        output:
                            file("*.fasta") into ( pcASV_labeledAA )
                            tuple file("*phyloformat.csv"), file("*summaryTable.tsv"), file("*dmd.out") into summary_potuaadiamond
                            tuple nid, file("*_summary_for_plot.csv") into taxplot4
                            tuple nid, file("*_quick_Taxbreakdown.csv") into tax_table_pcasvaa
                            tuple nid, file ("*_quicker_taxbreakdown.csv") into tax_nodCol_pcasvaa

                        script:
                            mtag="ID=" + nid
                            """
                            cp ${params.vampdir}/bin/rename_seq.py .
                            virdb=${params.dbdir}/${params.dbname}
                            if [[ ${params.measurement} == "bitscore" ]]
                            then    measure="--min-score ${params.bitscore}"
                            elif    [[ ${params.measurement} == "evalue" ]]
                            then    measure="-e ${params.evalue}"
                            else    measure="--min-score ${params.bitscore}"
                            fi
                            grep ">" \${virdb} > headers.list
                            headers="headers.list"
                            name=\$( echo ${asvs} | awk -F ".fasta" '{print \$1}')
                            if [[ ${params.ncbitax} == "true" ]]
                            then   diamond blastp -q ${asvs} -d \${virdb} -p ${task.cpus} --id ${params.minID} -l ${params.minaln} \${measure} --${params.sensitivity} -o "\$name"_dmd.out -f 6 qseqid qlen sseqid qstart qend qseq sseq length qframe evalue bitscore pident btop staxids sskingdoms skingdoms sphylums --max-target-seqs 1 --max-hsps 1
                            else   diamond blastp -q ${asvs} -d \${virdb} -p ${task.cpus} --id ${params.minID} -l ${params.minaln} \${measure} --${params.sensitivity} -o "\$name"_dmd.out -f 6 qseqid qlen sseqid qstart qend qseq sseq length qframe evalue bitscore pident btop --max-target-seqs 1 --max-hsps 1
                            fi
                            echo "Preparing lists to generate summary .csv's"
                            echo "[Best hit accession number]" > access.list
                            echo "[e-value]" > evalue.list
                            echo "[Bitscore]" > bit.list
                            echo "[Percent ID (aa)]" > pid.list
                            echo "[Organism ID]" > "\$name"_virus.list
                            echo "[Gene]" > "\$name"_genes.list
                            echo "[pcASV#]" > otu.list
                            echo "[Sequence length]" > length.list
                            grep ">" ${asvs} | awk -F ">" '{print \$2}' > seqids.lst
                            if [[ ${params.lca} == "T" ]]
                            then  grep -w "LCA" ${params.dbanno}/*.txt > lcainfo.list
                                  echo "[Taxonomic classification from RVDB annotations]" > lca_classification.list
                            else  echo "skipped" >> \${name}_quick_Taxbreakdown.csv
                                  echo "[Taxonomic classification from RVDB annotations]" > lca_classification.list
                            fi
                            if [[ ${params.ncbitax} == "true" ]]
                            then echo "[NCBI Taxonomy ID],[Taxonomic classification from NCBI]" > ncbi_classification.list
                            fi
                            echo "extracting genes and names"
                            touch new_"\$name"_asvnames.txt
                            for s in \$(cat seqids.lst);do
                                echo "Checking for \$s hit in diamond output"
                                if [[ "\$(grep -wc "\$s" "\$name"_dmd.out)" -eq 1 ]];then
                                            echo "Yep, there was a hit for \$s"
                                            echo "Extracting the information now:"
                                            acc=\$(grep -w "\$s" "\$name"_dmd.out | awk '{print \$3}')
                                            echo "\$s" >> otu.list
                                            echo "\$acc" >> access.list
                                            line="\$(grep -w "\$s" "\$name"_dmd.out)"
                                            echo "\$line" | awk '{print \$10}' >> evalue.list
                                            echo "\$line" | awk '{print \$11}' >> bit.list
                                            echo "\$line" | awk '{print \$12}' >> pid.list
                                            echo "\$line" | awk '{print \$2}' >> length.list
                                            echo "Extracting virus and gene ID for \$s now"
                                            gene=\$(grep -w "\$acc" "\$headers" | awk -F "." '{ print \$2 }' | awk -F "[" '{ print \$1 }' | awk -F " " '{print substr(\$0, index(\$0,\$2))}' | sed 's/ /_/g') &&
                                            echo "\$gene" | sed 's/_/ /g' >> "\$name"_genes.list
                                            virus=\$(grep -w "\$acc" "\$headers" | awk -F "[" '{ print \$2 }' | awk -F "]" '{ print \$1 }'| sed 's/ /_/g')
                                            echo "\$virus" | sed 's/_/ /g' >> "\$name"_virus.list
                                            echo ">"\${s}"_"\$virus"_"\$gene"" >> new_"\$name"_asvnames.txt
                                            if [[ "${params.lca}" == "T" ]]
                                            then    if [[ \$(grep -w "\$acc" ${params.dbanno}/*.txt | wc -l) -eq 1 ]]
                                                    then  group=\$(grep -w "\$acc" ${params.dbanno}/*.txt | awk -F ":" '{print \$1}')
                                                          lcla=\$(grep -w "\$group" lcainfo.list | awk -F "\t" '{print \$2}')
                                                          echo "\$lcla" >> lca_classification.list
                                                    else  echo "Viruses" >> lca_classification.list
                                                    fi
                                            fi
                                            if [[ ${params.ncbitax} == "true" ]]
                                            then  echo "\$line" | awk -F "\t" '{print \$14","\$16"::"\$18"::"\$17}' >> ncbi_classification.list
                                            fi
                                            echo "\$s done."
                                else
                                            echo "Ugh, there was no hit for \$s .."
                                            echo "We still love \$s though and we will add it to the final fasta file"
                                            echo "\$s" >> otu.list
                                            echo "NO_HIT" >> access.list
                                            echo "NO_HIT" >> "\$name"_genes.list
                                            echo "NO_HIT" >> "\$name"_virus.list
                                            echo "NO_HIT" >> evalue.list
                                            echo "NO_HIT" >> bit.list
                                            echo "NO_HIT" >> pid.list
                                            echo "NO_HIT" >> length.list
                                            virus="NO"
                                            gene="HIT"
                                            echo ">\${s}_"\$virus"_"\$gene"" >> new_"\$name"_asvnames.txt
                                            if [[ "${params.lca}" == "T" ]]
                                            then    echo "N/A" >> lca_classification.list
                                            fi
                                            if [[ "${params.ncbitax}" == "true" ]]
                                            then  echo "N/A" >> ncbi_classification.list
                                            fi
                                            echo "\$s done."
                               fi
                            done
                            echo "Now editing "\$name" fasta headers"
                            ###### rename_seq.py
                            ./rename_seq.py ${asvs} new_"\$name"_asvnames.txt "\$name"_TaxonomyLabels.fasta
                            awk 'BEGIN {RS=">";FS="\\n";OFS=""} NR>1 {print ">"\$1; \$1=""; print}' "\$name"_TaxonomyLabels.fasta >"\$name"_tmpssasv.fasta
                            echo "[Sequence header]" > newnames.list
                            cat new_"\$name"_asvnames.txt >> newnames.list
                            touch sequence.list
                            echo "     " > sequence.list
                            grep -v ">" "\$name"_tmpssasv.fasta >> sequence.list
                            rm "\$name"_tmpssasv.fasta
                            if [[ "${params.lca}" == "T" && "${params.ncbitax}" == "true" ]]
                            then
                                  paste -d "," sequence.list "\$name"_virus.list "\$name"_genes.list otu.list lca_classification.list ncbi_classification.list newnames.list length.list bit.list evalue.list pid.list access.list >> "\$name"_phyloformat.csv
                                  paste -d"\t" otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list ncbi_classification.list sequence.list length.list bit.list evalue.list pid.list >> "\$name"_summaryTable.tsv
                                  paste -d"," otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list ncbi_classification.list >> \${name}_quick_Taxbreakdown.csv
                            elif [[ "${params.lca}" == "T" && "${params.ncbitax}" != "true" ]]
                            then
                                  paste -d "," sequence.list "\$name"_virus.list "\$name"_genes.list otu.list lca_classification.list newnames.list length.list bit.list evalue.list pid.list access.list >> "\$name"_phyloformat.csv
                                  paste -d"\t" otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list sequence.list length.list bit.list evalue.list pid.list >> "\$name"_summaryTable.tsv
                                  paste -d"," otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list >> \${name}_quick_Taxbreakdown.csv
                            elif [[ "${params.ncbitax}" == "true" && "${params.lca}" != "T" ]]
                            then
                                  paste -d "," sequence.list "\$name"_virus.list "\$name"_genes.list otu.list ncbi_classification.list newnames.list length.list bit.list evalue.list pid.list access.list >> "\$name"_phyloformat.csv
                                  paste -d"\t" otu.list access.list "\$name"_virus.list "\$name"_genes.list ncbi_classification.list sequence.list length.list bit.list evalue.list pid.list >> "\$name"_summaryTable.tsv
                                  paste -d"," otu.list access.list "\$name"_virus.list "\$name"_genes.list ncbi_classification.list >> \${name}_quick_Taxbreakdown.csv
                            else
                                  paste -d "," sequence.list "\$name"_virus.list "\$name"_genes.list otu.list newnames.list length.list bit.list evalue.list pid.list access.list >> "\$name"_phyloformat.csv
                                  paste -d"\t" otu.list access.list "\$name"_virus.list "\$name"_genes.list sequence.list length.list bit.list evalue.list pid.list >> "\$name"_summaryTable.tsv
                                  echo "skipped" >> \${name}_quick_Taxbreakdown.csv
                            fi
                            for x in *phyloformat.csv;do
                                echo "\$x"
                                lin=\$(( \$(wc -l \$x | awk '{print \$1}')-1))
                                tail -"\$lin" \$x | awk -F "," '{print \$2}' > tmpcol.list;
                                sed 's/ /_/g' tmpcol.list > tmp2col.list;
                                cat tmp2col.list | sort | uniq -c | sort -nr | awk '{print \$2","\$1}' > \${name}_summary_for_plot.csv;
                                rm tmpcol.list tmp2col.list
                            done
                            awk -F "," '{print \$1","\$3"("\$2")"}' \${name}_quick_Taxbreakdown.csv >> \${name}_quicker_taxbreakdown.csv
                            rm evalue.list sequence.list bit.list pid.list length.list seqids.lst otu.list *asvnames.txt "\$name"_virus.list "\$name"_genes.list newnames.list access.list headers.list
                            """
                        }
                    } else if (params.dbtype== "RVDB") {

                      process pcASV_AminoAcid_Taxonomy_Inference_RVDB { //CHECK rename_seq.py in container !!!!!!!!!!!

                          label 'high_cpus'

                          tag "${mtag}"

                          publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Aminoacid/Taxonomy/SummaryFiles", mode: "copy", overwrite: true, pattern: '*.{csv,tsv}'
                          publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Aminoacid/Taxonomy/DiamondOutput", mode: "copy", overwrite: true, pattern: '*dmd.{out}'
                          publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Aminoacid/Taxonomy", mode: "copy", overwrite: true, pattern: '*.{fasta}'

                          conda (params.condaActivate ? "-c conda-forge bioconda::diamond=2.0.15=hb97b32f_1" : null)

                          container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/diamond:2.0.15--hb97b32f_1" : "quay.io/biocontainers/diamond:2.0.15--hb97b32f_1")

                          input:
                              tuple nid, file(asvs) from pcASV_aaDiamond_ch

                          output:
                              file("*.fasta") into ( pcASV_labeledAA )
                              tuple file("*phyloformat.csv"), file("*summaryTable.tsv"), file("*dmd.out") into summary_potuaadiamond
                              tuple nid, file("*_summary_for_plot.csv") into taxplot4
                              tuple nid, file("*_quick_Taxbreakdown.csv") into tax_table_pcasvaa
                              tuple nid, file ("*_quicker_taxbreakdown.csv") into tax_nodCol_pcasvaa

                          script:
                              mtag="ID=" + nid
                              """
                              cp ${params.vampdir}/bin/rename_seq.py .
                              virdb=${params.dbdir}/${params.dbname}
                              if [[ ${params.measurement} == "bitscore" ]]
                              then    measure="--min-score ${params.bitscore}"
                              elif    [[ ${params.measurement} == "evalue" ]]
                              then    measure="-e ${params.evalue}"
                              else    measure="--min-score ${params.bitscore}"
                              fi
                              grep ">" \${virdb} > headers.list
                              headers="headers.list"
                              name=\$( echo ${asvs} | awk -F ".fasta" '{print \$1}')
                              diamond blastp -q ${asvs} -d \${virdb} -p ${task.cpus} --id ${params.minID} -l ${params.minaln} \${measure} --${params.sensitivity} -o "\$name"_dmd.out -f 6 qseqid qlen sseqid qstart qend qseq sseq length qframe evalue bitscore pident btop --max-target-seqs 1 --max-hsps 1
                              echo "Preparing lists to generate summary .csv's"
                              echo "[Best hit accession number]" > access.list
                              echo "[e-value]" > evalue.list
                              echo "[Bitscore]" > bit.list
                              echo "[Percent ID (aa)]" > pid.list
                              echo "[Organism ID]" > "\$name"_virus.list
                              echo "[Gene]" > "\$name"_genes.list
                              echo "[pcASV#]" > otu.list
                              echo "[Sequence length]" > length.list
                              grep ">" ${asvs} | awk -F ">" '{print \$2}' > seqids.lst
                              if [[ ${params.lca} == "T" ]]
                              then  grep -w "LCA" ${params.dbanno}/*.txt > lcainfo.list
                                    echo "[Taxonomic classification from RVDB annotations]" > lca_classification.list
                              else  echo "skipped" >> \${name}_quick_Taxbreakdown.csv
                                    echo "[Taxonomic classification from RVDB annotations]" > lca_classification.list
                              fi
                              echo "extracting genes and names"
                              touch new_"\$name"_asvnames.txt
                              for s in \$(cat seqids.lst);do
                                  echo "Using RVDB headers."
                                  if [[ "\$(grep -wc "\$s" "\$name"_dmd.out)" -eq 1 ]];then
                                      echo "Yep, there was a hit for \$s"
                                      echo "Extracting the information now:"
                                      acc=\$(grep -w "\$s" "\$name"_dmd.out | awk '{print \$3}' | awk -F "|" '{print \$3}')
                                      echo "\$s" >> otu.list
                                      echo "\$acc" >> access.list
                                      line="\$(grep -w "\$s" "\$name"_dmd.out)"
                                      echo "\$line" | awk '{print \$10}' >> evalue.list
                                      echo "\$line" | awk '{print \$11}' >> bit.list
                                      echo "\$line" | awk '{print \$12}' >> pid.list
                                      echo "\$line" | awk '{print \$2}' >> length.list
                                      echo "Extracting virus and gene ID for \$s now"
                                      gene=\$(grep -w "\$acc" "\$headers" | awk -F "|" '{ print \$6 }' | awk -F "[" '{ print \$1 }' | sed 's/ /_/g') &&
                                      echo "\$gene" | sed 's/_/ /g' >> "\$name"_genes.list
                                      virus=\$(grep -w "\$acc" "\$headers" | awk -F "|" '{ print \$6 }' | awk -F "[" '{ print \$2 }' | awk -F "]" '{print \$1}' | sed 's/ /_/g') &&
                                      echo "\$virus" | sed 's/_/ /g' >> "\$name"_virus.list
                                      echo ">\${s}_"\$virus"_"\$gene"" >> new_"\$name"_asvnames.txt
                                      if [[ "${params.lca}" == "T" ]]
                                      then    if [[ \$(grep -w "\$acc" ${params.dbanno}/*.txt | wc -l) -eq 1 ]]
                                              then  group=\$(grep -w "\$acc" ${params.dbanno}/*.txt | awk -F ":" '{print \$1}')
                                                    lcla=\$(grep -w "\$group" lcainfo.list | awk -F "\t" '{print \$2}')
                                                    echo "\$lcla" >> lca_classification.list
                                              else  echo "Viruses" >> lca_classification.list
                                              fi
                                      fi
                                      echo "\$s done."
                                  else
                                      echo "Ugh, there was no hit for \$s .."
                                      echo "We still love \$s though and we will add it to the final fasta file"
                                      echo "\$s" >> otu.list
                                      echo "NO_HIT" >> access.list
                                      echo "NO_HIT" >> "\$name"_genes.list
                                      echo "NO_HIT" >> "\$name"_virus.list
                                      echo "NO_HIT" >> evalue.list
                                      echo "NO_HIT" >> bit.list
                                      echo "NO_HIT" >> pid.list
                                      echo "NO_HIT" >> length.list
                                      virus="NO"
                                      gene="HIT"
                                      echo ">\${s}_"\$virus"_"\$gene"" >> new_"\$name"_asvnames.txt
                                      if [[ "${params.lca}" == "T" ]]
                                      then    echo "N/A" >> lca_classification.list
                                      fi
                                      echo "\$s done."
                                  fi
                              echo "Done with \$s"
                              done
                              echo "Now editing "\$name" fasta headers"
                              ###### rename_seq.py
                              ./rename_seq.py ${asvs} new_"\$name"_asvnames.txt "\$name"_TaxonomyLabels.fasta
                              awk 'BEGIN {RS=">";FS="\\n";OFS=""} NR>1 {print ">"\$1; \$1=""; print}' "\$name"_TaxonomyLabels.fasta >"\$name"_tmpssasv.fasta
                              echo "[Sequence header]" > newnames.list
                              cat new_"\$name"_asvnames.txt >> newnames.list
                              touch sequence.list
                              echo "     " > sequence.list
                              grep -v ">" "\$name"_tmpssasv.fasta >> sequence.list
                              rm "\$name"_tmpssasv.fasta
                              if [[ "${params.lca}" == "T" ]]
                              then  paste -d "," sequence.list "\$name"_virus.list "\$name"_genes.list otu.list lca_classification.list newnames.list length.list bit.list evalue.list pid.list access.list >> "\$name"_phyloformat.csv
                                    paste -d"\t" otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list sequence.list length.list bit.list evalue.list pid.list >> "\$name"_summaryTable.tsv
                                    paste -d"," otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list >> \${name}_quick_Taxbreakdown.csv
                              else  paste -d "," sequence.list "\$name"_virus.list "\$name"_genes.list otu.list newnames.list length.list bit.list evalue.list pid.list access.list >> "\$name"_phyloformat.csv
                                    paste -d"\t" otu.list access.list "\$name"_virus.list "\$name"_genes.list sequence.list length.list bit.list evalue.list pid.list >> "\$name"_summaryTable.tsv
                              fi
                              for x in *phyloformat.csv;do
                                        echo "\$x"
                                        lin=\$(( \$(wc -l \$x | awk '{print \$1}')-1))
                                        tail -"\$lin" \$x | awk -F "," '{print \$2}' > tmpcol.list;
                                        sed 's/ /_/g' tmpcol.list > tmp2col.list;
                                        cat tmp2col.list | sort | uniq -c | sort -nr | awk '{print \$2","\$1}' > \${name}_summary_for_plot.csv;
                                        rm tmpcol.list tmp2col.list
                              done
                              awk -F "," '{print \$1","\$3"("\$2")"}' \${name}_quick_Taxbreakdown.csv >> \${name}_quicker_taxbreakdown.csv
                              rm evalue.list sequence.list bit.list pid.list length.list seqids.lst otu.list *asvnames.txt "\$name"_virus.list "\$name"_genes.list newnames.list access.list headers.list
                              """
                    }
                }
            } else {

                process skippcASVprotTaxonomy {

                    input:
                        tuple nid, file(asvs) from pcASV_aaDiamond_ch

                    output:
                        tuple nid, file("skipncASVprottax1.txt") into taxplot4
                        tuple nid, file("skipncASVprottax2.txt") into tax_table_pcasvaa
                        tuple nid, file("skipncASVprottax3.txt") into tax_nodCol_pcasvaa

                    script:
                        """
                        echo "Skipped" >skipncASVprottax1.txt
                        echo "Skipped" >skipncASVprottax2.txt
                        echo "Skipped" >skipncASVprottax3.txt
                        """
                }
            }

                if (!params.skipPhylogeny) {

                    process pcASV_Protein_Phylogeny_step1 {

                        label 'norm_cpus'

                        tag "${mtag}"

                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Aminoacid/Phylogeny/Alignment", mode: "copy", overwrite: true

                        conda (params.condaActivate ? "-c bioconda -c conda-forge muscle=5.1" : null)

                        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/muscle:5.1--h9f5acd7" : "quay.io/biocontainers/muscle:5.1--h9f5acd7")

            	        input:
                            tuple nid, file(pcASV) from pcASV_aaMafft_ch

                        output:
                            tuple nid, file("*_ALN.fasta") into pcASV_align1
                            tuple nid, file("*.efa")

                        script:
                        mtag="ID=" + nid
                            """
                            pre=\$(echo ${pcASV} | awk -F ".fasta" '{print \$1}' )

                            if [[ ${params.srep} == "true" && ${params.ensemble} == "false" ]];
                            then  if [[ \$( grep -c ">" ${pcASV}) -lt 300 ]]
                                  then    comm="align"
                                  else    comm="super5"
                                  fi
                                  muscle -"\$comm" ${pcASV} -perm ${params.perm} -perturb ${params.pert} -output \${pre}_muscle_raw_ALN.fasta -threads ${task.cpus} -quiet
                                  echo "single replicate alignment chosen; look over muscle5 documentation to learn about ensemble alignment approach" >note.efa
                            elif [[ ${params.srep} == "false" && ${params.ensemble} == "true" ]];
                            then  muscle -align ${pcASV} -${params.fied} -output \${pre}_muscle.efa -threads ${task.cpus} -quiet
                                  muscle -maxcc \${pre}_muscle.efa -output \${pre}_muscle_raw_ALN.fasta
                            else  if [[ \$( grep -c ">" ${pcASV}) -lt 300 ]]
                                  then    comm="align"
                                  else    comm="super5"
                                  fi
                                  muscle -"\$comm" ${pcASV} -output \${pre}_muscle_raw_ALN.fasta -threads ${task.cpus} -quiet
                                  echo "single replicate alignment chosen; look over muscle5 documentation to learn about ensemble alignment approach" >note.efa
                            fi
                            """
                    }

                    process pcASV_Protein_Phylogeny_step2 {

                        label 'norm_cpus'

                        tag "${mtag}"

                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Aminoacid/Phylogeny/Alignment", mode: "copy", overwrite: true, pattern: '*aln.*'
                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Aminoacid/Phylogeny/Modeltest", mode: "copy", overwrite: true, pattern: '*mt*'
                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Aminoacid/Phylogeny/IQ-TREE", mode: "copy", overwrite: true, pattern: '*iq*'

                        conda (params.condaActivate ? "-c conda-forge bioconda::trimal=1.4.1=h9f5acd7_6" : null)

                        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/trimal:1.4.1--h9f5acd7_6" : "quay.io/biocontainers/trimal:1.4.1--h9f5acd7_6")

                        input:
                            tuple nid, file(pcASV) from pcASV_align1

                        output:
                            tuple nid, file("*_aln.html") into pcASV_protein_phylogeny_results2
                            tuple nid, file("*_aln.fasta") into pcASV_align2

                        script:
                            mtag="ID=" + nid
                            """
                            pre=\$( echo ${pcASV} | awk -F "_muscle" '{print \$1}' )
                            trimal -in ${pcASV} -out \${pre}_trimal_aln.fasta -keepheader -fasta -automated1 -htmlout \${pre}_aln.html
                            """
                    }

                    process pcASV_Protein_Phylogeny_step3 {

                        label 'norm_cpus'

                        tag "${mtag}"

                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Aminoacid/Phylogeny/Alignment", mode: "copy", overwrite: true

                        conda (params.condaActivate ? "${params.vampdir}/bin/yamls/oligotyping.yml" : null)

                        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/oligotyping:2.1--py27_0" : "quay.io/biocontainers/oligotyping:2.1--py27_0")

                        input:
                            tuple nid, file(pcASV) from pcASV_align2

                        output:
                            tuple nid, file("*_Aligned_informativeonly.fasta") into pcASV_align3

                        script:
                            mtag="ID=" + nid
                            """
                            pre=\$( echo ${pcASV} | awk -F "_trimal" '{print \$1}' )
                            o-trim-uninformative-columns-from-alignment ${pcASV}
                            mv ${pcASV}-TRIMMED ./\${pre}_Aligned_informativeonly.fasta
                            """
                    }

                    process pcASV_Protein_Phylogeny_step4 {

                        label 'norm_cpus'

                        tag "${mtag}"

                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Aminoacid/Phylogeny/Modeltest", mode: "copy", overwrite: true, pattern: '*mt*'

                        conda (params.condaActivate ? "-c conda-forge bioconda::modeltest-ng=0.1.7=h5c6ebe3_0" : null)

                        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/modeltest-ng:0.1.7--h5c6ebe3_0" : "quay.io/biocontainers/modeltest-ng:0.1.7--h5c6ebe3_0")

                        input:
                            tuple nid, file(pcASV) from pcASV_align3

                        output:
                            tuple nid, file("*mt*") into pcASV_protein_phylogeny_results3
                            tuple nid, file("*mt.out") into potu_mtout

                        script:
                            mtag="ID=" + nid
                            """
                            pre=\$( echo ${pcASV} | awk -F "_Aligned" '{print \$1}' )
                            # pcASV_Protein_ModelTest
                            modeltest-ng -i ${pcASV} -p ${task.cpus} -o \${pre}_mt -d aa -s 203 --disable-checkpoint
                            """
                    }

                    process pcASV_Protein_Phylogeny_step5 {

                        label 'norm_cpus'

                        tag "${mtag}"

                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Aminoacid/Phylogeny/Alignment", mode: "copy", overwrite: true, pattern: '*aln.*'
                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Aminoacid/Phylogeny/Modeltest", mode: "copy", overwrite: true, pattern: '*mt*'
                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Aminoacid/Phylogeny/IQ-TREE", mode: "copy", overwrite: true, pattern: '*iq*'

                        conda (params.condaActivate ? "-c conda-forge bioconda::iqtree=2.2.0.3=hb97b32f_1" : null)

                        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/iqtree:2.2.0.3--hb97b32f_1" : "quay.io/biocontainers/iqtree:2.2.0.3--hb97b32f_1")

                        input:
                            tuple nid, file(prot) from pcASV_aaMafft_ch2
                            tuple nid, file(mtout) from potu_mtout

                        output:
                            tuple nid, file("*iq*") into pcASV_protein_phylogeny_results
                            tuple nid, file("*iq.treefile") into potu_Atree_plot

                        script:
                            mtag="ID=" + nid
                            """
                            #grabbing best models from modeltestng
                            modbic=\$(grep "iqtree" ${mtout} | head -1 | awk -F "-m " '{print \$2}')
                            modaic=\$(grep "iqtree" ${mtout} | head -2 | tail -1 | awk -F "-m " '{print \$2}')
                            modaicc=\$(grep "iqtree" ${mtout} | tail -1 | awk -F "-m " '{print \$2}')
                            if [[ ${params.crit} == "BIC" ]]
                            then  mod="\$modbic"
                            elif  [[ ${params.crit} == "AIC" ]]
                            then  mod="\$modaic"
                            elif  [[ ${params.crit} == "AICc" ]]
                            then  mod="\$modaicc"
                            fi
                            # grab prefix
                            pre=\$( echo ${prot} | awk -F "_Aligned" '{print \$1}' )
                            # pcASV_Protein_Phylogeny
                            if [ "${params.iqCustomaa}" != "" ];then
                                iqtree -s ${prot} --prefix \${pre}_iq --redo -T auto ${params.iqCustomaa}
                            elif [[ "${params.ModelTaa}" != "false" && "${params.nonparametric}" != "false" ]];then
                                iqtree -s ${prot} --prefix \${pre}_iq -m \${mod} --redo  -nt auto -b ${params.boots}
                            elif [[ "${params.ModelTaa}" != "false" && "${params.parametric}" != "false" ]];then
                                iqtree -s ${prot} --prefix \${pre}_iq -m \${mod} --redo -nt auto -bb ${params.boots} -bnni
                            elif [ "${params.nonparametric}" != "false" ];then
                                iqtree -s ${prot} --prefix \${pre}_iq -m MFP -madd --redo -nt auto -b ${params.boots}
                            elif [ "${params.parametric}" != "false" ];then
                                iqtree -s ${prot} --prefix \${pre}_iq -m MFP -madd --redo -nt auto -bb ${params.boots} -bnni
                            else
                                iqtree -s ${prot} --prefix \${pre}_iq -m MFP -madd --redo -nt auto -bb ${params.boots} -bnni
                            fi
                            """
                    }

                } else {

                    process skippcASVprotPhylogeny {

                        input:
                            tuple nid, file(prot) from pcASV_aaMafft_ch

                        output:
                            tuple nid, file("skippcASVprotPhylogeny.txt") into ( potu_Atree_plot )

                        script:
                            """
                            echo "Skipped" >skippcASVprotPhylogeny.txt
                            """
                    }

                }

                process Generate_pcASV_Protein_Counts {

                    label 'high_cpus'

                    tag "${mtag}"

                    publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Aminoacid/Counts", mode: "copy", overwrite: true

                    conda (params.condaActivate ? "-c conda-forge bioconda::diamond=2.0.15=hb97b32f_1" : null)

                    container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/diamond:2.0.15--hb97b32f_1" : "quay.io/biocontainers/diamond:2.0.15--hb97b32f_1")

                    input:
                        tuple nid, file(fasta) from pcASV_aaCounts_ch
                        file(merged) from mergeforpcASVaacounts
                        file(samplist) from samplistpotu

                    output:
                        tuple file("*_counts.csv"), file("*dmd.out") into potuaacounts_summary
                        tuple nid, file("*counts.csv") into potu_Acounts

                    script:
                        // check do I need for loop
                        mtag="ID=" + nid
                        """
                        set +e
                        potu="\$( echo ${fasta} | awk -F "_" '{print \$3}')"
                        diamond makedb --in ${fasta} --db ${fasta}
                        diamond blastx -q ${merged} -d ${fasta} -p ${task.cpus} --min-score ${params.ProtCountsBit} --id ${nid} -l ${params.ProtCountsLength} --${params.sensitivity} -o ${params.projtag}_\${potu}_Counts_dmd.out -f 6 qseqid qlen sseqid qstart qend qseq sseq length qframe evalue bitscore pident btop --max-target-seqs 1 --max-hsps 1 --max-hsps 1
                        echo "OTU_ID" >tmp.col1.txt
                        echo "Generating sample id list"
                        grep ">" ${fasta} | awk -F ">" '{print \$2}' | sort | uniq > otuid.list
                        cat otuid.list >> tmp.col1.txt
                        echo "Beginning them counts tho my g"
                        for y in \$( cat ${samplist} );do
                            echo "Starting with \$y now ..."
                            grep -w "\$y" ${params.projtag}_\${potu}_Counts_dmd.out > tmp."\$y".out
                            echo "Isolated hits"
                            echo "Created uniq subject id list"
                            echo "\$y" > "\$y"_col.txt
                            echo "Starting my counts"
                            for z in \$(cat otuid.list);do
                                echo "Counting \$z hits"
                	            echo "grep -wc "\$z" >> "\$y"_col.txt"
                	            grep -wc "\$z" tmp."\$y".out >> "\$y"_col.txt
            		            echo "\$z counted"
            	            done
                       done
                       paste -d "," tmp.col1.txt *col.txt > ${params.projtag}_aminoacid_\${potu}_noTaxonomy_counts.csv
                       rm tmp*
                       rm *col.txt
                       """
                   }
            }

            if (!params.skipReport) {

                if (!params.skipAdapterRemoval || !params.skipReadProcessing || !params.skipMerging) {

                    process combine_csv {

                        conda (params.condaActivate ? null : null)

                        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? null : null)

                        input:
                            file(csv) from fastp_csv_in2
                                .collect()

                        output:
                            file("final_reads_stats.csv") into fastp_csv_in

                        script:
                            """
                            cat ${csv} >all_reads_stats.csv
                            head -n1 all_reads_stats.csv >tmp.names.csv
                            cat all_reads_stats.csv | grep -v ""Sample,Total_"" >tmp.reads.stats.csv
                            cat tmp.names.csv tmp.reads.stats.csv >final_reads_stats.csv
                            rm tmp.names.csv tmp.reads.stats.csv
                            """
                    }
                } else {

                    process skip_combine_csv {

                        output:
                            file("filter_reads.txt") into fastp_csv_in

                        script:
                            """
                            echo "Read processing steps skipped." >filter_reads.txt
                            """
                    }
                }

                report_asv = Channel.create()
                asv_counts_plots.mix(taxplot_asv, asv_heatmap, nucl_phyl_plot_asv, asvgroupscsv, asvgroupcounts, asv_group_rep_tree, tax_table_asv, tax_nodCol_asv, asv_phylogroupcsv, asv_phylogroupingcsv).flatten().buffer(size:11).dump(tag:'asv').into(report_asv)

                if (params.ncASV) {
                    report_ncasv = Channel.create()
                    notu_counts_plots.mix(taxplot_ncasv, notu_heatmap, nucl_phyl_plot_ncasv, tax_table_ncasv, tax_nodCol_ncasv).groupTuple(by:0, size:6).dump(tag:'ncasv').into(report_ncasv)
                } else {
                    report_ncasv = Channel.empty()
                }

                if (params.pcASV) {
                    report_pcasv_aa = Channel.create()
                    potu_Acounts.mix(taxplot4, potu_aa_heatmap, potu_Atree_plot, tax_table_pcasvaa, tax_nodCol_pcasvaa).groupTuple(by:0, size:6).dump(tag:'pcasv1').into(report_pcasv_aa)

                    report_pcasv_nucl = Channel.create()
                    potu_Ncounts_for_report.mix(taxplot3, potu_nucl_heatmap, potu_Ntree_plot, tax_table_pcasvnt, tax_nodCol_pcasvnt).groupTuple(by:0, size:6).dump(tag:'pcasv2').into(report_pcasv_nucl)

                } else {
                     report_pcasv_aa = Channel.empty()
                     report_pcasv_nucl = Channel.empty()
                }

                if (!params.skipAminoTyping) {
                    report_aminotypes = Channel.create()
                    aminocounts_plot.mix(taxplot2, aminotype_heatmap, amino_rax_plot, atygroupscsv, amino_group_rep_tree, amino_groupcounts, tax_table_amino, tax_nodCol_amino, amino_phylogroupcsv, amino_phylogroupingcsv).flatten().buffer(size:11).dump(tag:'amino').into(report_aminotypes)

                } else {
                    report_aminotypes = Channel.empty()
                }

                report_all_ch = Channel.create()
                report_asv.mix(report_ncasv, report_pcasv_aa, report_pcasv_nucl, report_aminotypes).map{it.flatten()}.dump(tag:'report').into(report_all_ch)

                process Report {

                    label 'norm_cpus'

                    publishDir "${params.workingdir}/${params.outdir}/Analyze/FinalReports", mode: "copy", overwrite: true

                    conda (params.condaActivate ? "${params.vampdir}/bin/yamls/R.yml" : null)

                    container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/mulled-v2-90fb71778d446b480b1050805be9ef346794df6b:6d9281817ba120685d81edae0370730c1ff554cc-0" : "quay.io/biocontainers/mulled-v2-90fb71778d446b480b1050805be9ef346794df6b:6d9281817ba120685d81edae0370730c1ff554cc-0")

                    input:
                        file(csv) from fastp_csv_in
                        file(files) from report_all_ch

                    output:
                        file("*.html") into report_all_out

                    script:
                        """
                        name=\$( ls *_counts.csv | awk -F "_counts.csv" '{print \$1}')
                        type=\$( ls *_counts.csv | awk -F "${params.projtag}" '{print \$2}' | awk -F "_" '{print \$2}'  )
                        cp ${params.vampdir}/bin/vAMPirus_Report.Rmd .
                        cp ${params.vampdir}/example_data/conf/vamplogo.png .
                        Rscript -e "rmarkdown::render('vAMPirus_Report.Rmd',output_file='\${name}_Report.html')" \${name} \
                        ${params.skipReadProcessing} \
                        ${params.skipMerging} \
                        ${params.skipAdapterRemoval} \
                        ${params.skipTaxonomy} \
                        ${params.skipPhylogeny} \
                        ${params.trymax} \
                        ${params.stats} \
                        ${params.metadata} \
                        ${params.minimumCounts} \
                        ${params.asvMED} \
                        ${params.aminoMED} \
                        \${type} \
                        ${params.nodeCol} \
                        ${params.asvTClust} \
                        ${params.aminoTClust} \
                        """
                }
            }
        }

} else {
    println("\n\t\033[0;31mMandatory argument not specified. For more info use `nextflow run vAMPirus.nf --help`\n\033[0m")
    exit 0
}
if (params.DataCheck) {
    workflow.onComplete {
        log.info ( workflow.success ? \
            "---------------------------------------------------------------------------------" \
            + "\n\033[0;32mDone! Open the following reports in your browser\033[0m" \
            + "\n\033[0;32mPipeline performance report: ${params.workingdir}/${params.outdir}/${params.tracedir}/vampirus_report.html\033[0m" \
            + "\n\033[0;32mvAMPirus --DataCheck interactive report: ${params.workingdir}/${params.outdir}/DataCheck/Report/*.hmtl\033[0m" \
            : \
            "---------------------------------------------------------------------------------" \
            + "\n\033[0;31mSomething went wrong. Check error message below and/or log files.\033[0m" )
    }
} else if (params.Analyze) {
    workflow.onComplete {
        log.info ( workflow.success ? \
            "---------------------------------------------------------------------------------" \
            + "\n\033[0;32mDone! Open the following reports in your browser\033[0m" \
            + "\n\033[0;32mPipeline performance report: ${params.workingdir}/${params.outdir}/${params.tracedir}/vampirus_report.html\033[0m" \
            + "\n\033[0;32mvAMPirus --Analyze interactive report: ${params.workingdir}/${params.outdir}/Analyze/*.hmtl\033[0m" \
            : \
            "---------------------------------------------------------------------------------" \
            + "\n\033[0;31mSomething went wrong. Check error message below and/or log files.\033[0m" )
    }
}
