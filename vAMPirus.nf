muscle#!/usr/bin/env nextflow

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


            --Primer removal--

                General primer removal parameters

                    --primerLength                  Use this option to set the max primer length to restrict bbduk.sh primer trimming to the first x number of bases

                    --maxkmer                       Maximum kmer length for bbduk.sh to use for primer detection and removal (must be shorter than your primer length; default = 13)

                    --minkmer                       Minimum kmer length for primer removal (default = 3)

                    --minilen                       Minimum read length after adapter and primer removal (default = 200)

                Single primer set removal-

                    --GlobTrim                      Set this option to perform global trimming to reads to remove primer sequences. Example usage "--GlobTrim #basesfromforward,#basesfromreverse"

                    --fwd                           Forward primer sequence for reads to be detected and removed from reads (must specify reverse sequence if providing forward)

                    --rev                           Reverse primer sequence for reads to be detected and removed from reads (must specify forward sequence if providing reverse)

                Multiple primer set removal-

                    --multi                         Use this option to signal multiple primer sequence removal within the specified pipeline

                    --primers                       Use this option to set the path to a fasta file with all of the primer sequences to be detected and removed from reads


            --Amplicon Sequence Variant (ASV) genration and clustering--

                    --alpha                         Alpha value for denoising - the higher the alpha the higher the chance of false positives in ASV generation (1 or 2)

                    --minSize                       Minimum size or representation for sequence to be considered in ASV generation

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


              --Minimum Entropy Decomposition arguments--

                      --asvMED                        Set this option to perform Minimum Entropy Decomposition on ASV sequences, see manual for more information. You will need to set a value for --asvC to perform this analysis

                      --amino_med                     Set this option to perform Minimum Entropy Decomposition on AminoType sequences, see manual for more information. You will need to set a value for --aminoC to perform this analysis

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


              --Primer removal--

                  General primer removal parameters

                      --primerLength                  Use this option to set the max primer length to restrict bbduk.sh primer trimming to the first x number of bases

                      --maxkmer                       Maximum kmer length for bbduk.sh to use for primer detection and removal (must be shorter than your primer length; default = 13)

                      --minkmer                       Minimum kmer length for primer removal (default = 3)

                      --minilen                       Minimum read length after adapter and primer removal (default = 200)

                  Single primer set removal-

                      --GlobTrim                      Set this option to perform global trimming to reads to remove primer sequences. Example usage "--GlobTrim #basesfromforward,#basesfromreverse"

                      --fwd                           Forward primer sequence for reads to be detected and removed from reads (must specify reverse sequence if providing forward)

                      --rev                           Reverse primer sequence for reads to be detected and removed from reads (must specify forward sequence if providing reverse)

                  Multiple primer set removal-

                      --multi                         Use this option to signal multiple primer sequence removal within the specified pipeline

                      --primers                       Use this option to set the path to a fasta file with all of the primer sequences to be detected and removed from reads


              --Amplicon Sequence Variant (ASV) genration and clustering--

                      --alpha                         Alpha value for denoising - the higher the alpha the higher the chance of false positives in ASV generation (1 or 2)

                      --minSize                       Minimum size or representation for sequence to be considered in ASV generation

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

        Example 1. Launching the vAMPirus DataCheck pipeline using conda

            nextflow run vAMPirus.nf -c vampirus.config -profile conda --DataCheck

        Example 2. Launching the vAMPirus DataCheck pipeline using Singularity and multiple primer removal with the path to the fasta file with the primer sequences set in the launch command

            nextflow run vAMPirus.nf -c vampirus.config -profile singularity --DataCheck --multi --primers /PATH/TO/PRIMERs.fa

        Example 3. Launching the vAMPirus DataCheck pipeline with primer removal by global trimming of 20 bp from forward reads and 26 bp from reverse reads

            nextflow run vAMPirus.nf -c vampirus.config -profile conda --DataCheck --GlobTrim 20,26


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
                                        Project name:                ${params.projtag}
                                        Working directory:           ${params.workingdir}
                                        Results directory:           ${params.outdir}
                                        Database directory:          ${params.dbdir}
                                        Database name:               ${params.dbname}
                                        Metadata file:               ${params.metadata}
        """.stripIndent()

if (params.readsTest) {
    println("\n\tRunning vAMPirus with TEST dataset\n")
    Channel
        .fromFilePairs(params.readsTest)
        .ifEmpty{ exit 1, "params.readTest was empty - no input files supplied" }
        .into{ reads_ch; reads_qc_ch; reads_processing }
} else {
    println("\n\tEverything ready for launch.\n")
    Channel
        .fromFilePairs("${params.reads}", checkIfExists: true)
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

    println("\n\tRunning vAMPirus Analyze pipeline - This might take a while, check out Nextflow tower (tower.nf) to remotely monitor the run.\n")

    if (!params.skipTaxonomy) {

        process Database_Check {
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
                                if [[ ${ncbitax} == "true" && ${dbtype} == "NCBI" ]]
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
                            diamond makedb --in ${params.dbdir}/${params.dbname} -d ${params.dbdir}/${params.dbname}
                            export virdb=${params.dbdir}/${params.dbname}
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

            process QualityCheck_1 {

                label 'low_cpus'

                tag "${sample_id}"

                publishDir "${params.workingdir}/${params.outdir}/ReadProcessing/FastQC/PreClean", mode: "copy", overwrite: true

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

        if (!params.skipAdapterRemoval) {

            process Adapter_Removal {

                label 'norm_cpus'

                tag "${sample_id}"

                publishDir "${params.workingdir}/${params.outdir}/ReadProcessing/AdapterRemoval", mode: "copy", overwrite: true, pattern: "*.filter.fq"
                publishDir "${params.workingdir}/${params.outdir}/ReadProcessing/AdapterRemoval/fastpOut", mode: "copy", overwrite: true, pattern: "*.fastp.{json,html}"

                input:
                    tuple sample_id, file(reads) from reads_ch

                output:
                    tuple sample_id, file("*.fastp.{json,html}") into fastp_results
                    tuple sample_id, file("*.filter.fq") into reads_fastp_ch
                    file("*.csv") into ( fastp_csv_in1, fastp_csv_in2 )

                script:
                    """
                    echo ${sample_id}

                    fastp -i ${reads[0]} -I ${reads[1]} -o left-${sample_id}.filter.fq -O right-${sample_id}.filter.fq --detect_adapter_for_pe \
                    --average_qual 25 -c --overrepresentation_analysis --html ${sample_id}.fastp.html --json ${sample_id}.fastp.json --thread ${task.cpus} \
                    --report_title ${sample_id}

                    bash get_readstats.sh ${sample_id}.fastp.json
                    """
                }
        } else {
            reads_ch
                .set{ reads_fastp_ch }
            fastp_results = Channel.empty()

        }

        if (!params.skipPrimerRemoval) {

            process Primer_Removal {

                label 'norm_cpus'

                tag "${sample_id}"

                publishDir "${params.workingdir}/${params.outdir}/ReadProcessing/PrimerRemoval", mode: "copy", overwrite: true

                input:
                    tuple sample_id, file(reads) from reads_fastp_ch

                output:
                    tuple sample_id, file("*bbduk*.fastq.gz") into ( reads_bbduk_ch, readsforqc2 )

                script:
                    // check if we need to check this outside processes
                    if ( params.fwd == "" && params.rev == "" && !params.multi) {
                        """
                        bbduk.sh in1=${reads[0]} out=${sample_id}_bb_R1.fastq.gz ftl=${params.defaultFwdTrim} t=${task.cpus}
                        bbduk.sh in=${reads[1]} out=${sample_id}_bb_R2.fastq.gz ftl=${params.defaultRevTrim} t=${task.cpus}
        		            repair.sh in1=${sample_id}_bb_R1.fastq.gz in2=${sample_id}_bb_R2.fastq.gz out1=${sample_id}_bbduk_R1.fastq.gz out2=${sample_id}_bbduk_R2.fastq.gz outs=sing.fq repair
                        """
                    } else if ( params.GlobTrim && !params.GlobTrim == "" ) {
                        """
                        FTRIM=\$( echo ${GlobTrim} | cut -f 1 -d "," )
                        RTRIM=\$( echo ${GlobTrim} | cut -f 2 -d "," )
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
        } else {
            reads_fastp_ch
                .set{ reads_bbduk_ch }

        }

        if (!params.skipFastQC && !params.skipPrimerRemoval) {

            process QualityCheck_2 {

                label 'low_cpus'

                tag "${sample_id}"

                publishDir "${params.workingdir}/${params.outdir}/ReadProcessing/FastQC/PostClean", mode: "copy", overwrite: true

                input:
                    tuple sample_id, file(reads) from readsforqc2

                output:
                    tuple sample_id, file("*_fastqc.{zip,html}") into fastqc2_results

                script:
                    """
                    fastqc --quiet --threads ${task.cpus} ${reads}
                    """
            }
        }
    } else {
        reads_ch
            .set{ reads_bbduk_ch }

    }

    if (!params.skipMerging) {

        process Read_Merging {

            label 'norm_cpus'

            tag "${sample_id}"

            publishDir "${params.workingdir}/${params.outdir}/ReadProcessing/ReadMerging/Individual", mode: "copy", overwrite: true, pattern: "*mergedclean.fastq"
            publishDir "${params.workingdir}/${params.outdir}/ReadProcessing/ReadMerging/Individual/notmerged", mode: "copy", overwrite: true, pattern: "*notmerged*.fastq"

            input:
                tuple sample_id, file(reads) from reads_bbduk_ch

            output:
                file("*_mergedclean.fastq") into reads_vsearch1_ch
                file("*.name") into names
                file("*notmerged*.fastq") into notmerged

            script:
                """
                vsearch --fastq_mergepairs ${reads[0]} --reverse ${reads[1]} --threads ${task.cpus} --fastqout ${sample_id}_mergedclean.fastq --fastqout_notmerged_fwd ${sample_id}_notmerged_fwd.fastq --fastqout_notmerged_rev ${sample_id}_notmerged_rev.fastq  --fastq_maxdiffs ${params.diffs} --fastq_maxns ${params.maxn} --fastq_allowmergestagger --fastq_maxee ${params.maxEE} --relabel ${sample_id}.
                echo ${sample_id} > ${sample_id}.name
                """

        }

    } else {
        reads_bbduk_ch
          .set{ reads_vsearch1_ch }
    }


    process Filtering_Prep1 {

        label 'low_cpus'

        publishDir "${params.workingdir}/${params.outdir}/ReadProcessing/ReadMerging/LengthFiltering", mode: "copy", overwrite: true

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
            #bbduk.sh in=${reads} bhist=${params.projtag}_all_merged_preFilt_preClean_baseFrequency_hist.txt qhist=${params.projtag}_all_merged_preFilt_preClean_qualityScore_hist.txt gchist=${params.projtag}_all_merged_preFilt_preClean_gcContent_hist.txt aqhist=${params.projtag}_all_merged_preFilt_preClean_averageQuality_hist.txt lhist=${params.projtag}_all_merged_preFilt__preClean_length_hist.txt gcbins=auto
            #fastp -i ${reads} -o ${params.projtag}_merged_preFilt_clean.fastq -b ${params.maxLen} -l ${params.minLen} --thread ${task.cpus} -n 1
            #reformat.sh in=${params.projtag}_merged_preFilt_clean.fastq out=${params.projtag}_merged_preFilt_clean.fasta t=${task.cpus}
            #bbduk.sh in=${params.projtag}_merged_preFilt_clean.fastq out=${params.projtag}_merged_clean_Lengthfiltered.fastq minlength=${params.maxLen} maxlength=${params.maxLen} t=${task.cpus}
            #bbduk.sh in=${params.projtag}_merged_clean_Lengthfiltered.fastq bhist=${params.projtag}_all_merged_postFilt_baseFrequency_hist.txt qhist=${params.projtag}_all_merged_postFilt_qualityScore_hist.txt gchist=${params.projtag}_all_merged_postFilt_gcContent_hist.txt aqhist=${params.projtag}_all_merged_postFilt_averageQuaulity_hist.txt lhist=${params.projtag}_all_merged_postFilt_length_hist.txt gcbins=auto

            # from DC
            bbduk.sh in=${reads} bhist=${params.projtag}_all_merged_preFilt_preClean_baseFrequency_hist.txt qhist=${params.projtag}_all_merged_preFilt_preClean_qualityScore_hist.txt gchist=${params.projtag}_all_merged_preFilt_preClean_gcContent_hist.txt aqhist=${params.projtag}_all_merged_preFilt_preClean_averageQuality_hist.txt lhist=${params.projtag}_all_merged_preFilt_preClean_length_hist.txt gcbins=auto
            for x in *preFilt*hist.txt;do
                pre=\$(echo \$x | awk -F ".txt" '{print \$1}')
                cat \$x | tr "\t" "," > \${pre}.csv
                rm \$x
            done
            reformat.sh in=${reads} out=${params.projtag}_preFilt_preclean.fasta t=${task.cpus}
            echo "sample,reads" >> reads_per_sample_preFilt_preClean.csv
            grep ">" ${params.projtag}_preFilt_preclean.fasta | awk -F ">" '{print \$2}' | awk -F "." '{print \$1}' | sort --parallel=${task.cpus} | uniq -c | sort -brg --parallel=${task.cpus} | awk '{print \$2","\$1}' >> reads_per_sample_preFilt_preClean.csv
            rm ${params.projtag}_preFilt_preclean.fasta
            fastp -i ${reads} -o ${params.projtag}_merged_preFilt_clean.fastq -b ${params.maxLen} -l ${params.minLen} --thread ${task.cpus} -n 1
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

        input:
            file(reads) from reads_vsearch2_ch

        output:
            file("*unique_sequences.fasta") into reads_vsearch3_ch

        script:
            """
            vsearch --derep_fulllength ${reads} --sizeout --relabel_keep --output ${params.projtag}_unique_sequences.fasta
            """
    }

    process Identify_ASVs {

        label 'norm_cpus'

        publishDir "${params.workingdir}/${params.outdir}/Analyze/Clustering/ASVs/ChimeraCheck", mode: "copy", overwrite: true

        input:
            file(reads) from reads_vsearch3_ch

        output:
            file("*notChecked.fasta") into reads_vsearch4_ch

        script:
            """
            vsearch --cluster_unoise ${reads} --unoise_alpha ${params.alpha} --relabel ASV --centroids ${params.projtag}_notChecked.fasta --minsize ${params.minSize}
            """
    }

    process Chimera_Check {

        label 'low_cpus'

        publishDir "${params.workingdir}/${params.outdir}/Analyze/Clustering/ASVs", mode: "copy", overwrite: true

        input:
            file(fasta) from reads_vsearch4_ch

        output:
            file("*ASVs.fasta") into ( reads_vsearch5_ch, asv_med, nucl2aa, asvsforAminotyping, asvfastaforcounts, asvaminocheck )

        script:
            """
	        vsearch --uchime3_denovo ${fasta} --relabel ASV --nonchimeras ${params.projtag}_ASVs.fasta
            """
    }

    // UNTIL HERE DEFAULT

    if (params.DataCheck) {

        process NucleotideBased_ASV_clustering_DC {

        label 'norm_cpus'

        publishDir "${params.workingdir}/${params.outdir}/DataCheck/Clustering/Nucleotide", mode: "copy", overwrite: true, pattern: '*{.csv}'

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

       process Translation_For_ProteinBased_Clustering_DC {

           label 'norm_cpus'

           publishDir "${params.workingdir}/${params.outdir}/DataCheck/Clustering/Aminoacid/translation", mode: "copy", overwrite: true

           input:
                file(fasta) from nucl2aa

            output:
                file("*ASVprotforclust.fasta") into clustering_aa
                file("*_translation_report") into reportaa_VR
                file("*_ASV_all.fasta") into asvfastaforaaclust

            script:
                """
                ${tools}/virtualribosomev2/dna2pep.py ${fasta} -r all -x -o none --fasta ${params.projtag}_ASVprotforclust.fasta --report ${params.projtag}_translation_report
                cp ${fasta} ${params.projtag}_ASV_all.fasta
                """
       }

        process Protein_clustering_DC {

            label 'norm_cpus'

            publishDir "${params.workingdir}/${params.outdir}/DataCheck/Clustering/Aminoacid", mode: "copy", overwrite: true, pattern: '*{.csv}'

            input:
                file(fasta) from clustering_aa
                file(asvs) from asvfastaforaaclust

            output:
                file("number_per_percentage_prot.csv") into number_per_percent_prot_plot
                file("*aminoacid_pcASV1.0_noTaxonomy.fasta") into amino_med
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
                    awk 'BEGIN{RS=">";ORS=""}length(\$2)>="${params.minAA}"{print ">"\$0}' ${fasta} > ${params.projtag}_filtered_proteins.fasta
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
                    awk 'BEGIN{RS=">";ORS=""}length(\$2)<"${params.minAA}"{print ">"\$0}' ${fasta} >${params.projtag}_pcASV\${id}_problematic_translations.fasta
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

        process ASV_Shannons_Entropy_Analysis {

          label 'norm_cpus'

          publishDir "${params.workingdir}/${params.outdir}/DataCheck/Clustering/Nucleotide/ShannonEntropy", mode: "copy", overwrite: true

          input:
              file(asvs) from asv_med

          output:

              file("*_ASV_entropy_breakdown.csv") into asv_entro_csv
              file("*Aligned_informativeonly.fasta-ENTROPY") into asv_entropy
              file("*ASV*") into entrop

          script:
          """
            set +e
            #alignment
            ${tools}/muscle5.0.1278_linux64 -in ${asvs} -out ${params.projtag}_ASVs_muscleAlign.fasta -thread ${task.cpus} -quiet
            #mafft --thread ${task.cpus} --maxiterate 15000 --auto ${asvs} > ${params.projtag}_ASVs_muscleAlign.fasta
            #trimming
            trimal -in ${params.projtag}_ASVs_mafftAlign.fasta -out ${params.projtag}_ASVs_muscleAligned.fasta  -keepheader -fasta -automated1
            rm ${params.projtag}_ASVs_muscleAlign.fasta
            o-trim-uninformative-columns-from-alignment ${params.projtag}_ASVs_muscleAligned.fasta
            mv ${params.projtag}_ASVs_mafftAligned.fasta-TRIMMED ./${params.projtag}_ASVs_Aligned_informativeonly.fasta
            #entopy analysis
            entropy-analysis ${params.projtag}_ASVs_Aligned_informativeonly.fasta
            #summarize entropy peaks
            awk '{print \$2}' ${params.projtag}_ASVs_Aligned_informativeonly.fasta-ENTROPY >> tmp_value.list
            for x in \$(cat tmp_value.list)
            do      echo "\$x"
                    if [[ \$(echo ""\$x" > 0.0"|bc -l) -eq 1 ]];
                    then    echo dope >> above-0.0-.list
                    fi
                    if [[ \$(echo ""\$x" > 0.1"|bc -l) -eq 1 ]];
                    then    echo dope >> above-0.1-.list
                    fi
                    if [[ \$(echo ""\$x" > 0.2"|bc -l) -eq 1 ]];
                    then    echo dope >> above-0.2-.list
                    fi
                    if [[ \$(echo ""\$x" > 0.3"|bc -l) -eq 1 ]];
                    then    echo dope >> above-0.3-.list
                    fi
                    if [[ \$(echo ""\$x" > 0.4"|bc -l) -eq 1 ]];
                    then    echo dope >> above-0.4-.list
                    fi
                    if [[ \$(echo ""\$x" > 0.5"|bc -l) -eq 1 ]];
                    then    echo dope >> above-0.5-.list
                    fi
                    if [[ \$(echo ""\$x" > 0.6"|bc -l) -eq 1 ]];
                    then    echo dope >> above-0.6-.list
                    fi
                    if [[ \$(echo ""\$x" > 0.7"|bc -l) -eq 1 ]];
                    then    echo dope >> above-0.7-.list
                    fi
                    if [[ \$(echo ""\$x" > 0.8"|bc -l) -eq 1 ]];
                    then    echo dope >> above-0.8-.list
                    fi
                    if [[ \$(echo ""\$x" > 0.9"|bc -l) -eq 1 ]];
                    then    echo dope >> above-0.9-.list
                    fi
                    if [[ \$(echo ""\$x" > 1.0"|bc -l) -eq 1 ]];
                    then    echo dope >> above-1.0-.list
                    fi
                    if [[ \$(echo ""\$x" > 1.5"|bc -l) -eq 1 ]];
                    then    echo dope >> above-1.5-.list
                    fi
            done
            echo "Entropy,Peaks_above" >> ${params.projtag}_ASV_entropy_breakdown.csv
            for z in above*.list;
            do      entrop=\$(echo \$z | awk -F "-" '{print \$2}')
                    echo ""\$entrop", "\$(wc -l \$z | awk '{print \$1}')"" >> ${params.projtag}_ASV_entropy_breakdown.csv
            done
            rm above*
            mv ${params.projtag}_ASVs_Aligned_informativeonly.fasta-ENTROPY ./tmp.fasta
            echo "Base_position  Shannons_Entropy" >> ${params.projtag}_ASVs_Aligned_informativeonly.fasta-ENTROPY
            cat tmp.fasta >> ${params.projtag}_ASVs_Aligned_informativeonly.fasta-ENTROPY
            rm tmp.fasta

          """

        }

        process AminoType_Shannons_Entropy_Analysis {

          label 'norm_cpus'

          publishDir "${params.workingdir}/${params.outdir}/DataCheck/Clustering/Aminoacid/ShannonEntropy", mode: "copy", overwrite: true, pattern: '*{.csv}'

          input:
              file(aminos) from amino_med

          output:
              file("*AminoType_entropy_breakdown.csv") into amino_entro_csv
              file ("*Aligned_informativeonly.fasta-ENTROPY") into amino_entropy
              file("*AminoTypes*") into aminos

          script:
            """
            #alignment
            if [[ $(grep -c ">" ${aminos}) -gt 499 ]]; then algo="super5"; else algo="mpc"; fi
            ${tools}/muscle5.0.1278_linux64 -"\${algo}" ${aminos} -out ${params.projtag}_AminoTypes_muscleAlign.fasta -thread ${task.cpus} -quiet
            #mafft --thread ${task.cpus} --maxiterate 15000 --auto ${aminos} > ${params.projtag}_AminoTypes_muscleAlign.fasta
            #trimming
            trimal -in ${params.projtag}_AminoTypes_muscleAlign.fasta -out ${params.projtag}_AminoTypes_muscleAligned.fasta  -keepheader -fasta -automated1
            rm ${params.projtag}_AminoTypes_muscleAlign.fasta
            o-trim-uninformative-columns-from-alignment ${params.projtag}_AminoTypes_muscleAligned.fasta
            mv ${params.projtag}_AminoTypes_muscleAligned.fasta-TRIMMED ./${params.projtag}_AminoTypes_Aligned_informativeonly.fasta
            #entropy analysis
            entropy-analysis ${params.projtag}_AminoTypes_Aligned_informativeonly.fasta --amino-acid-sequences
            #summarize entropy peaks
            awk '{print \$2}' ${params.projtag}_AminoTypes_Aligned_informativeonly.fasta-ENTROPY >> tmp_value.list
            for x in \$(cat tmp_value.list)
            do      echo "\$x"
                    if [[ \$(echo ""\$x" > 0.0"|bc -l) -eq 1 ]];
                    then    echo dope >> above-0.0-.list
                    fi
                    if [[ \$(echo ""\$x" > 0.1"|bc -l) -eq 1 ]];
                    then    echo dope >> above-0.1-.list
                    fi
                    if [[ \$(echo ""\$x" > 0.2"|bc -l) -eq 1 ]];
                    then    echo dope >> above-0.2-.list
                    fi
                    if [[ \$(echo ""\$x" > 0.3"|bc -l) -eq 1 ]];
                    then    echo dope >> above-0.3-.list
                    fi
                    if [[ \$(echo ""\$x" > 0.4"|bc -l) -eq 1 ]];
                    then    echo dope >> above-0.4-.list
                    fi
                    if [[ \$(echo ""\$x" > 0.5"|bc -l) -eq 1 ]];
                    then    echo dope >> above-0.5-.list
                    fi
                    if [[ \$(echo ""\$x" > 0.6"|bc -l) -eq 1 ]];
                    then    echo dope >> above-0.6-.list
                    fi
                    if [[ \$(echo ""\$x" > 0.7"|bc -l) -eq 1 ]];
                    then    echo dope >> above-0.7-.list
                    fi
                    if [[ \$(echo ""\$x" > 0.8"|bc -l) -eq 1 ]];
                    then    echo dope >> above-0.8-.list
                    fi
                    if [[ \$(echo ""\$x" > 0.9"|bc -l) -eq 1 ]];
                    then    echo dope >> above-0.9-.list
                    fi
                    if [[ \$(echo ""\$x" > 1.0"|bc -l) -eq 1 ]];
                    then    echo dope >> above-1.0-.list
                    fi
                    if [[ \$(echo ""\$x" > 1.5"|bc -l) -eq 1 ]];
                    then    echo dope >> above-1.5-.list
                    fi
            done
            echo "Entropy,Peaks_above" >> ${params.projtag}_AminoType_entropy_breakdown.csv
            for z in above*.list;
            do      entrop=\$(echo \$z | awk -F "-" '{print \$2}')
                    echo ""\$entrop", "\$(wc -l \$z | awk '{print \$1}')"" >> ${params.projtag}_AminoType_entropy_breakdown.csv
            done
            rm above*
            mv ${params.projtag}_AminoTypes_Aligned_informativeonly.fasta-ENTROPY ./tmp.fasta
            echo "Base_position  Shannons_Entropy" >> ${params.projtag}_AminoTypes_Aligned_informativeonly.fasta-ENTROPY
            cat tmp.fasta >> ${params.projtag}_AminoTypes_Aligned_informativeonly.fasta-ENTROPY
            rm tmp.fasta
            """
        }

        if (!params.skipReadProcessing || !params.skipMerging ) {

            process combine_csv_DC {

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

                script:
                    """
                    echo "Read processing steps skipped." >filter_reads.txt
                    """
            }
        }

        report_dc_in = Channel.create()
        fastp_csv_dc.mix( reads_per_sample_preFilt, reads_per_sample_postFilt, prefilt_basefreq, postFilt_basefreq, prefilt_qualityscore, postFilt_qualityscore, prefilt_gccontent, postFilt_gccontent, prefilt_averagequality, postFilt_averagequality, prefilt_length, postFilt_length, number_per_percent_nucl_plot, number_per_percent_prot_plot, amino_entro_csv, amino_entropy, asv_entro_csv, asv_entropy).into(report_dc_in)

        process Report_DataCheck {

            label 'norm_cpus'

            publishDir "${params.workingdir}/${params.outdir}/DataCheck/Report", mode: "copy", overwrite: true, pattern: '*.{html}'

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
                ${params.skipAdapterRemoval}
                """
        }

    } else if (params.Analyze) {

        if (params.ncASV) {

            reads_vsearch5_ch
    	       .into{ asv_file_for_ncasvs; nuclFastas_forDiamond_asv_ch; nuclFastas_forCounts_asv_ch; nuclFastas_forphylogeny_asv; nuclFastas_forMatrix_asv_ch; asv_for_med }

            process NucleotideBased_ASV_clustering {

                label 'norm_cpus'

                tag "${mtag}"

                publishDir "${params.workingdir}/${params.outdir}/Analyze/Clustering/ncASV", mode: "copy", overwrite: true, pattern: '*ncASV*.fasta'

                input:
                    each x from 1..nnuc
                    file(fasta) from asv_file_for_ncasvs

                output:
                    tuple nid, file("*_ncASV*.fasta") into ( nuclFastas_forphylogeny_ncasv, nuclFastas_forDiamond_ncasv_ch, nuclFastas_forCounts_ncasv_ch, nuclFastas_forMatrix_ncasv_ch )

                script:
                    nid=slist.get(x-1)
                    mtag="ID=" + slist.get(x-1)
                    """
                    vsearch --cluster_fast ${fasta} --centroids ${params.projtag}_ncASV${nid}.fasta --threads ${task.cpus} --relabel ncASV --id .${nid}
                    """
            }

            if (!params.skipTaxonomy) {

              if (params.dbtype == "NCBI") {

                process ncASV_Taxonomy_Inference_NCBI { /////// editttt

                    label 'high_cpus'

                    tag "${mtag}"

                    publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ncASV/Taxonomy", mode: "copy", overwrite: true, pattern: '*ncASV*.{fasta,csv,tsv}'
                    publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ncASV/Taxonomy/DiamondOutput", mode: "copy", overwrite: true, pattern: '*ncASV*dmd.out'

                    input:
                        tuple nid, file(asvs) from nuclFastas_forDiamond_ncasv_ch

                    output:
                        file("*.fasta") into tax_labeled_fasta_ncasv
                        tuple file("*_phyloformat.csv"), file("*summaryTable.tsv"), file("*dmd.out") into summary_diamond_ncasv
                        tuple nid, file("*ncASV*summary_for_plot.csv") into taxplot_ncasv
                        tuple nid, file("*_quick_Taxbreakdown.csv") into tax_table_ncasv

                    script:
                        mtag="ID=" + nid
                        """
                        cp ${params.vampdir}/bin/rename_seq.py .
                        virdb=${params.dbdir}/${params.dbname}
                        grep ">" \${virdb} > headers.list
                        headers="headers.list"
                        name=\$( echo ${asvs} | awk -F ".fasta" '{print \$1}')
                        if [[ ${ncbitax} == "true" ]]
                        then   diamond blastx -q ${asvs} -d \${virdb} -p ${task.cpus} --id ${params.minID} -l ${params.minaln} --min-score ${params.bitscore} --${sensitivity} -o "\$name"_dmd.out -f 6 qseqid qlen sseqid qstart qend qseq sseq length qframe evalue bitscore pident btop staxids sskingdoms skingdoms sphylums --max-target-seqs 1 --max-hsps 1
                        else   diamond blastx -q ${asvs} -d \${virdb} -p ${task.cpus} --id ${params.minID} -l ${params.minaln} --min-score ${params.bitscore} --${sensitivity} -o "\$name"_dmd.out -f 6 qseqid qlen sseqid qstart qend qseq sseq length qframe evalue bitscore pident btop --max-target-seqs 1 --max-hsps 1
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
                        if [[ ${ncbitax} == "true" ]]
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
                                        gene=\$(grep -w "\$acc" "\$headers" | awk -F "." '{ print \$2 }' | awk -F "[" '{ print \$1 }' | awk -F " " print substr(\$0, index(\$0,\$2)) | sed 's/ /_/g') &&
                                        echo "\$gene" | sed 's/_/ /g' >> "\$name"_genes.list
                                        virus=\$(grep -w "\$acc" "\$headers" | awk -F "[" '{ print \$2 }' | awk -F "]" '{ print \$1 }'| sed 's/ /_/g')
                                        echo "\$virus" | sed 's/_/ /g' >> "\$name"_virus.list
                                        echo ">"\${s}"_"\$virus"_"\$gene"" >> new_"\$name"_asvnames.txt
                                        if [[ "${params.lca}" == "T" ]]
                                        then    if [[ \$(grep -w "\$acc" ${params.dbanno}/*.txt | wc -l) -eq 1 ]]
                                                then  group=\$(grep -w "\$acc" ${params.dbanno}/*.txt | awk -F ":" '{print \$1}')
                                                      lcla=\$(grep -w "\$group" lcainfo.list | awk -F "\t" '{print \$2}')
                                                      echo "\$lcla" >> lca_classification.list
                                                else  echo "Viruses::unclassified" >> lca_classification.list
                                                fi
                                        fi
                                        if [[ ${ncbitax} == "true" ]]
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
                                        if [[ "${ncbitax}" == "true" ]]
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
                        if [[ "${params.lca}" == "T" && "${ncbitax}" == "true" ]]
                        then
                              paste -d "," sequence.list "\$name"_virus.list "\$name"_genes.list otu.list lca_classification.list ncbi_classification.list newnames.list length.list bit.list evalue.list pid.list access.list >> "\$name"_phyloformat.csv
                              paste -d"\t" otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list ncbi_classification.list sequence.list length.list bit.list evalue.list pid.list >> "\$name"_summaryTable.tsv
                              paste -d"," otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list ncbi_classification.list >> \${name}_quick_Taxbreakdown.csv
                        elif [[ "${params.lca}" == "T" && "${ncbitax}" != "true" ]]
                        then
                              paste -d "," sequence.list "\$name"_virus.list "\$name"_genes.list otu.list lca_classification.list newnames.list length.list bit.list evalue.list pid.list access.list >> "\$name"_phyloformat.csv
                              paste -d"\t" otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list sequence.list length.list bit.list evalue.list pid.list >> "\$name"_summaryTable.tsv
                              paste -d"," otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list >> \${name}_quick_Taxbreakdown.csv
                        elif [[ "${ncbitax}" == "true" && "${params.lca}" != "T"]]
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
                        rm evalue.list sequence.list bit.list pid.list length.list seqids.lst otu.list *asvnames.txt "\$name"_virus.list "\$name"_genes.list newnames.list access.list headers.list
                        """
                      }
            } else if (params.dbtype== "RVDB") {

              process ncASV_Taxonomy_Inference_RVDB { /////// editttt

                  label 'high_cpus'

                  tag "${mtag}"

                  publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ncASV/Taxonomy", mode: "copy", overwrite: true, pattern: '*ncASV*.{fasta,csv,tsv}'
                  publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ncASV/Taxonomy/DiamondOutput", mode: "copy", overwrite: true, pattern: '*ncASV*dmd.out'

                  input:
                      tuple nid, file(asvs) from nuclFastas_forDiamond_ncasv_ch

                  output:
                      file("*.fasta") into tax_labeled_fasta_ncasv
                      tuple file("*_phyloformat.csv"), file("*summaryTable.tsv"), file("*dmd.out") into summary_diamond_ncasv
                      tuple nid, file("*ncASV*summary_for_plot.csv") into taxplot_ncasv
                      tuple nid, file("*_quick_Taxbreakdown.csv") into tax_table_ncasv

                  script:
                      mtag="ID=" + nid
                      """
                      cp ${params.vampdir}/bin/rename_seq.py .
                      virdb=${params.dbdir}/${params.dbname}
                      grep ">" \${virdb} > headers.list
                      headers="headers.list"
                      name=\$( echo ${asvs} | awk -F ".fasta" '{print \$1}')
                      diamond blastx -q ${asvs} -d \${virdb} -p ${task.cpus} --id ${params.minID} -l ${params.minaln} --min-score ${params.bitscore} --${sensitivity} -o "\$name"_dmd.out -f 6 qseqid qlen sseqid qstart qend qseq sseq length qframe evalue bitscore pident btop --max-target-seqs 1 --max-hsps 1
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
                                      else  echo "Viruses::unclassified" >> lca_classification.list
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
                      rm evalue.list sequence.list bit.list pid.list length.list seqids.lst otu.list *asvnames.txt "\$name"_virus.list "\$name"_genes.list newnames.list access.list headers.list
                      """
                }
            }
          }


            process Generate_ncASV_Counts_Table {

                label 'norm_cpus'

                tag "${mtag}"

                publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ASVs/Counts", mode: "copy", overwrite: true, pattern: '*_ASV*.{biome,csv}'
                publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ncASV/Counts", mode: "copy", overwrite: true, pattern: '*ncASV*.{biome,csv}'

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
                    cat \${name}_PercentIDq.matrix | tr " " ","  | grep "," >\${name}_PercentID.matrix
                    rm \${name}_PercentIDq.matrix
                    """
                }

                if (!params.skipPhylogeny) {

                    process ncASV_Phylogeny {

                          label 'norm_cpus'

                          tag "${mtag}"

                          publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ncASV/Phylogeny/Alignment", mode: "copy", overwrite: true,  pattern: '*ncASV*aln.*'
                          publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ncASV/Phylogeny/ModelTest", mode: "copy", overwrite: true, pattern: '*ncASV*mt*'
                          publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ncASV/Phylogeny/IQ-TREE", mode: "copy", overwrite: true, pattern: '*ncASV*iq*'

                          input:
                              tuple nid, file(asvs) from nuclFastas_forphylogeny_ncasv

                          output:
                              tuple nid, file("*_aln.fasta"), file("*_aln.html"), file("*.tree"), file("*.log"), file("*iq*"), file("*mt*") into align_results_ncasv
                              tuple nid, file("*iq.treefile") into nucl_phyl_plot_ncasv

                          script:
                              mtag="ID=" + nid
                              """
                              pre=\$(echo ${asvs} | awk -F ".fasta" '{print \$1}' )
                              ${tools}/muscle5.0.1278_linux64 -in ${asvs} -out \${pre}_ALN.fasta -thread ${task.cpus} -quiet
                              #mafft --thread ${task.cpus} --maxiterate 15000 --auto ${asvs} >\${pre}_ALN.fasta
                              trimal -in \${pre}_ALN.fasta -out \${pre}_aln.fasta -keepheader -fasta -automated1 -htmlout \${pre}_aln.html
                              o-trim-uninformative-columns-from-alignment \${pre}_aln.fasta
                              mv \${pre}_aln.fasta-TRIMMED ./\${pre}_Aligned_informativeonly.fasta
                              # Nucleotide_ModelTest
                              modeltest-ng -i \${pre}_Aligned_informativeonly.fasta -p ${task.cpus} -o \${pre}_mt -d nt -s 203 --disable-checkpoint
                              # Nucleotide_Phylogeny
                              if [ "${params.iqCustomnt}" != "" ];then
                                  iqtree -s \${pre}_Aligned_informativeonly.fasta --prefix \${pre}_iq --redo -t \${pre}_mt.tree -T auto ${params.iqCustomnt}
                              elif [[ "${params.ModelTnt}" != "false" && "${params.nonparametric}" != "false" ]];then
                                  mod=\$(tail -12 \${pre}_Aligned_informativeonly.fasta.log | head -1 | awk '{print \$6}')
                                  iqtree -s \${pre}_Aligned_informativeonly.fasta --prefix \${pre}_iq -m \${mod} --redo -t \${pre}_mt.tree -nt auto -b ${params.boots}
                              elif [[ "${params.ModelTnt}" != "false" && "${params.parametric}" != "false" ]];then
                                  mod=\$(tail -12 \${pre}_Aligned_informativeonly.fasta.log | head -1 | awk '{print \$6}')
                                  iqtree -s \${pre}_Aligned_informativeonly.fasta --prefix \${pre}_iq -m \${mod} --redo -t \${pre}_mt.tree -nt auto -bb ${params.boots} -bnni
                              elif [ "${params.nonparametric}" != "false" ];then
                                  iqtree -s \${pre}_Aligned_informativeonly.fasta --prefix \${pre}_iq -m MFP --redo -t \${pre}_mt.tree -nt auto -b ${params.boots}
                              elif [ "${params.parametric}" != "false" ];then
                                  iqtree -s \${pre}_Aligned_informativeonly.fasta --prefix \${pre}_iq -m MFP --redo -t \${pre}_mt.tree -nt auto -bb ${params.boots} -bnni
                              else
                                  iqtree -s \${pre}_Aligned_informativeonly.fasta --prefix \${pre}_iq -m MFP --redo -t \${pre}_mt.tree -nt auto -bb ${params.boots} -bnni
                              fi
                              """
                      }

                }

        } else {
            reads_vsearch5_ch
    	       .into{ nuclFastas_forDiamond_asv_ch; nuclFastas_forCounts_asv_ch; nuclFastas_forphylogeny_asv; nuclFastas_forMatrix_asv_ch; asv_for_med }
        }

        if (!params.skipTaxonomy) {

          if (params.dbtype == "NCBI") {

                process ASV_Taxonomy_Inference_NCBI { /////// editttt

                    label 'high_cpus'

                    publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ASVs/Taxonomy", mode: "copy", overwrite: true, pattern: '*_ASV*.{fasta,csv,tsv}'
                    publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ASVs/Taxonomy/DiamondOutput", mode: "copy", overwrite: true, pattern: '*_ASV*dmd.out'

                    input:
                        file(asvs) from nuclFastas_forDiamond_asv_ch

                    output:
                        file("*.fasta") into tax_labeled_fasta_asv
                        tuple file("*_phyloformat.csv"), file("*summaryTable.tsv"), file("*dmd.out") into summary_diamond_asv
                        file("*_ASV*_summary_for_plot.csv") into taxplot_asv
                        file("*_quick_Taxbreakdown.csv") into tax_table_asv

                    script:
                        """
                        cp ${params.vampdir}/bin/rename_seq.py .
                        virdb=${params.dbdir}/${params.dbname}
                        grep ">" \${virdb} > headers.list
                        headers="headers.list"
                        name=\$( echo ${asvs} | awk -F ".fasta" '{print \$1}')
                        if [[ ${ncbitax} == "true" ]]
                        then   diamond blastx -q ${asvs} -d \${virdb} -p ${task.cpus} --id ${params.minID} -l ${params.minaln} --min-score ${params.bitscore} --${sensitivity} -o "\$name"_dmd.out -f 6 qseqid qlen sseqid qstart qend qseq sseq length qframe evalue bitscore pident btop staxids sskingdoms skingdoms sphylums --max-target-seqs 1 --max-hsps 1
                        else   diamond blastx -q ${asvs} -d \${virdb} -p ${task.cpus} --id ${params.minID} -l ${params.minaln} --min-score ${params.bitscore} --${sensitivity} -o "\$name"_dmd.out -f 6 qseqid qlen sseqid qstart qend qseq sseq length qframe evalue bitscore pident btop --max-target-seqs 1 --max-hsps 1
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
                        if [[ ${ncbitax} == "true" ]]
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
                                        gene=\$(grep -w "\$acc" "\$headers" | awk -F "." '{ print \$2 }' | awk -F "[" '{ print \$1 }' | awk -F " " print substr(\$0, index(\$0,\$2)) | sed 's/ /_/g') &&
                                        echo "\$gene" | sed 's/_/ /g' >> "\$name"_genes.list
                                        virus=\$(grep -w "\$acc" "\$headers" | awk -F "[" '{ print \$2 }' | awk -F "]" '{ print \$1 }'| sed 's/ /_/g')
                                        echo "\$virus" | sed 's/_/ /g' >> "\$name"_virus.list
                                        echo ">"\${s}"_"\$virus"_"\$gene"" >> new_"\$name"_asvnames.txt
                                        if [[ "${params.lca}" == "T" ]]
                                        then    if [[ \$(grep -w "\$acc" ${params.dbanno}/*.txt | wc -l) -eq 1 ]]
                                                then  group=\$(grep -w "\$acc" ${params.dbanno}/*.txt | awk -F ":" '{print \$1}')
                                                      lcla=\$(grep -w "\$group" lcainfo.list | awk -F "\t" '{print \$2}')
                                                      echo "\$lcla" >> lca_classification.list
                                                else  echo "Viruses::unclassified" >> lca_classification.list
                                                fi
                                        fi
                                        if [[ ${ncbitax} == "true" ]]
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
                                        if [[ "${ncbitax}" == "true" ]]
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
                        if [[ "${params.lca}" == "T" && "${ncbitax}" == "true" ]]
                        then
                              paste -d "," sequence.list "\$name"_virus.list "\$name"_genes.list otu.list lca_classification.list ncbi_classification.list newnames.list length.list bit.list evalue.list pid.list access.list >> "\$name"_phyloformat.csv
                              paste -d"\t" otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list ncbi_classification.list sequence.list length.list bit.list evalue.list pid.list >> "\$name"_summaryTable.tsv
                              paste -d"," otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list ncbi_classification.list >> \${name}_quick_Taxbreakdown.csv
                        elif [[ "${params.lca}" == "T" && "${ncbitax}" != "true" ]]
                        then
                              paste -d "," sequence.list "\$name"_virus.list "\$name"_genes.list otu.list lca_classification.list newnames.list length.list bit.list evalue.list pid.list access.list >> "\$name"_phyloformat.csv
                              paste -d"\t" otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list sequence.list length.list bit.list evalue.list pid.list >> "\$name"_summaryTable.tsv
                              paste -d"," otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list >> \${name}_quick_Taxbreakdown.csv
                        elif [[ "${ncbitax}" == "true" && "${params.lca}" != "T"]]
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
                        rm evalue.list sequence.list bit.list pid.list length.list seqids.lst otu.list *asvnames.txt "\$name"_virus.list "\$name"_genes.list newnames.list access.list headers.list
                        """
                      }
                } else if (params.dbtype== "RVDB") {

                  process ASV_Taxonomy_Inference_RVDB { /////// editttt

                      label 'high_cpus'

                      publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ASVs/Taxonomy", mode: "copy", overwrite: true, pattern: '*_ASV*.{fasta,csv,tsv}'
                      publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ASVs/Taxonomy/DiamondOutput", mode: "copy", overwrite: true, pattern: '*_ASV*dmd.out'

                      input:
                          file(asvs) from nuclFastas_forDiamond_asv_ch

                      output:
                          file("*.fasta") into tax_labeled_fasta_asv
                          tuple file("*_phyloformat.csv"), file("*summaryTable.tsv"), file("*dmd.out") into summary_diamond_asv
                          file("*_ASV*_summary_for_plot.csv") into taxplot_asv
                          file("*_quick_Taxbreakdown.csv") into tax_table_asv

                      script:
                          """
                          cp ${params.vampdir}/bin/rename_seq.py .
                          virdb=${params.dbdir}/${params.dbname}
                          grep ">" \${virdb} > headers.list
                          headers="headers.list"
                          name=\$( echo ${asvs} | awk -F ".fasta" '{print \$1}')
                          diamond blastx -q ${asvs} -d \${virdb} -p ${task.cpus} --id ${params.minID} -l ${params.minaln} --min-score ${params.bitscore} --${sensitivity} -o "\$name"_dmd.out -f 6 qseqid qlen sseqid qstart qend qseq sseq length qframe evalue bitscore pident btop --max-target-seqs 1 --max-hsps 1
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
                                          else  echo "Viruses::unclassified" >> lca_classification.list
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
                          rm evalue.list sequence.list bit.list pid.list length.list seqids.lst otu.list *asvnames.txt "\$name"_virus.list "\$name"_genes.list newnames.list access.list headers.list
                          """
                }
            }
        }

        process Generate_ASV_Counts_Tables {

            label 'norm_cpus'

            publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ASVs/Counts", mode: "copy", overwrite: true, pattern: '*ASV*.{biome,csv}'

            input:
                file(asvs) from nuclFastas_forCounts_asv_ch
                file(merged) from nuclCounts_mergedreads_asv_ch

            output:
                tuple file("*_counts.csv"), file("*_counts.biome") into counts_vsearch_asv
                file("*_ASV*counts.csv") into (asv_counts_plots, asvcount_med)

            script:
                """
                name=\$( echo ${asvs} | awk -F ".fasta" '{print \$1}')
                vsearch --usearch_global ${merged} --db ${asvs} --id .${params.asvcountID} --threads ${task.cpus} --otutabout "\$name"_counts.txt --biomout "\$name"_counts.biome
                cat \${name}_counts.txt | tr "\t" "," >\${name}_count.csv
                sed 's/#OTU ID/OTU_ID/g' \${name}_count.csv >\${name}_counts.csv
                rm \${name}_count.csv
                """
        }

        process Generate_ASV_Matrices {

            label 'low_cpus'

            publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ASVs/Matrices", mode: "copy", overwrite: true, pattern: '*ASV*PercentID.matrix'

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

            if (!params.skipPhylogeny) { // need to edit paths

                process ASV_Phylogeny {

                      label 'norm_cpus'

                      publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ASVs/Phylogeny/Alignment", mode: "copy", overwrite: true,  pattern: '*ASV*aln.*'
                      publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ASVs/Phylogeny/ModelTest", mode: "copy", overwrite: true, pattern: '*ASV*mt*'
                      publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ASVs/Phylogeny/IQ-TREE", mode: "copy", overwrite: true, pattern: '*ASV*iq*'

                      input:
                          file(asvs) from nuclFastas_forphylogeny_asv

                      output:
                          tuple file("*_aln.fasta"), file("*_aln.html"), file("*.tree"), file("*.log"), file("*iq*"), file("*mt*") into align_results_asv
                          file("*iq.treefile") into (nucl_phyl_plot_asv, asvphy_med)

                      script:
                          """
                          pre=\$(echo ${asvs} | awk -F ".fasta" '{print \$1}' )
                          ${tools}/muscle5.0.1278_linux64 -in ${asvs} -out \${pre}_ALN.fasta -thread ${task.cpus} -quiet
                          #mafft --thread ${task.cpus} --maxiterate 15000 --auto ${asvs} >\${pre}_ALN.fasta
                          trimal -in \${pre}_ALN.fasta -out \${pre}_aln.fasta -keepheader -fasta -automated1 -htmlout \${pre}_aln.html
                          o-trim-uninformative-columns-from-alignment \${pre}_aln.fasta
                          mv \${pre}_aln.fasta-TRIMMED ./\${pre}_Aligned_informativeonly.fasta
                          # Nucleotide_ModelTest
                          modeltest-ng -i \${pre}_Aligned_informativeonly.fasta -p ${task.cpus} -o \${pre}_mt -d nt -s 203 --disable-checkpoint
                          # Nucleotide_Phylogeny
                          if [ "${params.iqCustomnt}" != "" ];then
                              iqtree -s \${pre}_Aligned_informativeonly.fasta --prefix \${pre}_iq --redo -T auto ${params.iqCustomnt}
                          elif [[ "${params.ModelTnt}" != "false" && "${params.nonparametric}" != "false" ]];then
                              mod=\$(tail -12 \${pre}_Aligned_informativeonly.fasta.log | head -1 | awk '{print \$6}')
                              iqtree -s \${pre}_Aligned_informativeonly.fasta --prefix \${pre}_iq -m \${mod} --redo -nt auto -b ${params.boots}
                          elif [[ "${params.ModelTnt}" != "false" && "${params.parametric}" != "false" ]];then
                              mod=\$(tail -12 \${pre}_Aligned_informativeonly.fasta.log | head -1 | awk '{print \$6}')
                              iqtree -s \${pre}_Aligned_informativeonly.fasta --prefix \${pre}_iq -m \${mod} --redo -nt auto -bb ${params.boots} -bnni
                          elif [ "${params.nonparametric}" != "false" ];then
                              iqtree -s \${pre}_Aligned_informativeonly.fasta --prefix \${pre}_iq -m MFP --redo -nt auto -b ${params.boots}
                          elif [ "${params.parametric}" != "false" ];then
                              iqtree -s \${pre}_Aligned_informativeonly.fasta --prefix \${pre}_iq -m MFP --redo -nt auto -bb ${params.boots} -bnni
                          else
                              iqtree -s \${pre}_Aligned_informativeonly.fasta --prefix \${pre}_iq -m MFP --redo -nt auto -bb ${params.boots} -bnni
                          fi
                          """
                  }
            }

            if (params.asvMED) {

                process ASV_Minimum_Entropy_Decomposition {

                label 'low_cpus'

                publishDir "${params.workingdir}/${params.outdir}/Analyze/Clustering/ASVs/MED", mode: "copy", overwrite: true

                input:
                  file(asvs) from asv_for_med

                output:
                  file("*_ASV_Grouping.csv") into asvgroupscsv
                  file("${params.projtag}_ASV_group_reps_aligned.fasta") into groupreps
                  file("${params.projtag}_asvMED_${params.asvC}")

                script:
                    """
                    #alignment
                    ${tools}/muscle5.0.1278_linux64 -in ${asvs} -out  -thread ${task.cpus} -quiet
                    #mafft --thread ${task.cpus} --maxiterate 15000 --auto ${asvs} > ${params.projtag}_ASVs_muscleAlign.fasta
                    #trimming
                    trimal -in ${params.projtag}_ASVs_muscleAlign.fasta -out ${params.projtag}_ASVs_muscleAligned.fasta  -keepheader -fasta -automated1
                    rm ${params.projtag}_ASVs_muscleAlign.fasta
                    o-trim-uninformative-columns-from-alignment ${params.projtag}_ASVs_muscleAligned.fasta
                    mv ${params.projtag}_ASVs_muscleAligned.fasta-TRIMMED ./${params.projtag}_ASVs_Aligned_informativeonly.fasta
                    #entopy analysis
                    entropy-analysis ${params.projtag}_ASVs_Aligned_informativeonly.fasta
                    #Decomposition
                    if [[ \$(echo ${params.asvC} | grep -c ",") -eq 1 ]]
                    then
                          tag=$(echo ${params.asvC} | sed 's/,/_/g')
                          oligotype ${params.projtag}_ASVs_Aligned_informativeonly.fasta ${params.projtag}_ASVs_Aligned_informativeonly.fasta-ENTROPY -o ${params.projtag}_asvMED_"\$tag -M 1 -C ${params.asvC} -N ${task.cpus} --skip-check-input --no-figures --skip-gen-html
                    else
                          oligotype ${params.projtag}_ASVs_Aligned_informativeonly.fasta ${params.projtag}_ASVs_Aligned_informativeonly.fasta-ENTROPY -o ${params.projtag}_asvMED_${params.asvC} -M 1 -c ${params.asvC} -N ${task.cpus} --skip-check-input --no-figures --skip-gen-html
                    fi
                    #generatemaps
                    cd ./${params.projtag}_asvMED_${params.asvC}/OLIGO-REPRESENTATIVES/
                    echo "ASV,GroupID,IDPattern"
                    j=1
                    for x in *_unique;
                    do      gid=\$(echo \$x | awk -F "_" '{print \$1}')
                            uni=\$(echo \$x | awk -F ""\${gid}"_" '{print \$2}' | awk -F "_uni" '{print \$1}')
                            grep ">"  "\$gid"_"\$uni" | awk -F ">" '{print \$2}' > asv.list
                            seqtk subseq ../../${asvs} asv.list > Group"\${j}"_sequences.fasta
                            for z in \$( cat asv.list)
                            do      echo ""\$z",Group"\$j","\$uni"" >> ${params.projtag}_ASV_Grouping.csv

                            done
                            rm asv.list
                            echo ">Group\${j}" >> ${params.projtag}_ASV_group_reps_aligned.fasta
                            echo "\$uni" > group.list
                            seqtk subseq ../OLIGO-REPRESENTATIVES.fasta group.list > group.fasta
                            tail -1 group.fasta >> ${params.projtag}_ASV_group_reps_aligned.fasta
                            mv "\$gid"_"\$uni" ./Group"\$j"_"\$uni"_aligned.fasta
                            mv "\$gid"_"\$uni"_unique ./Group"\$j"_"\$uni"_unqiues_aligned.fasta
                            rm "\$gid"*.cPickle
                            j=\$((\$j+1))
                    done
                    mv ${params.projtag}_ASV_Grouping.csv ../../
                    mv ${params.projtag}_ASV_group_reps_aligned.fasta ../../
                    cd ..
                    """
              }

              if (!params.skipPhylogeny) {

                process ASV_MED_Reps_phylogeny {

                label 'low_cpus'

                publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ASVs/MED/Phylogeny/ModelTest", mode: "copy", overwrite: true, pattern: '*ASV*mt*'
                publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ASVs/MED/Phylogeny/IQ-TREE", mode: "copy", overwrite: true, pattern: '*ASV*iq*'

                input:
                  file(reps) from groupreps

                output:
                  file("*_ASV_Group_Reps*") into align_results_asvmed
                  file("*iq.treefile") into asv_group_rep_tree

                script:
                    """
                    # Protein_ModelTest
                    modeltest-ng -i ${reps} -p ${task.cpus} -o ${params.projtag}_ASV_Group_Reps_mt -d aa -s 203 --disable-checkpoint

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
                        iqtree -s ${reps} --prefix ${params.projtag}_ASV_Group_Reps_iq -m MFP --redo -nt auto -b ${params.boots}

                    elif [ "${params.parametric}" != "false" ];then
                        iqtree -s ${reps} --prefix ${params.projtag}_ASV_Group_Reps_iq -m MFP --redo -nt auto -bb ${params.boots} -bnni

                    else
                        iqtree -s ${reps} --prefix ${params.projtag}_ASV_Group_Reps_iq -m MFP --redo -nt auto -bb ${params.boots} -bnni
                    fi
                    """
                }
              }
              process Adding_ASV_MED_Info {

              label 'low_cpus'

              publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/ASVs/MED/", mode: "copy", overwrite: true

              input:
                  file(counts) from asvcount_med
                  file(map) from asvgroupscsv

              output:
                  file("${params.projtag}_ASV_Groupingcounts.csv") into asvgroupcounts

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
            }

            if (!params.skipAminoTyping) {

                process Translate_For_AminoTyping {

                  label 'low_cpus'

                  publishDir "${params.workingdir}/${params.outdir}/Analyze/Clustering/AminoTypes/Translation", mode: "copy", overwrite: true

                  input:
                      file(fasta) from asvsforAminotyping

                  output:
                      file("${params.projtag}_all_translations.fasta") into amintypegen
                      file("${params.projtag}_translation_report") into proteinstage_vap_report

                  script:
                      """
                      ${tools}/virtualribosomev2/dna2pep.py ${fasta} -r all -x -o none --fasta ${params.projtag}_all_translations.fasta --report ${params.projtag}_translation_report
                      """
                }

                process Generate_AminoTypes {

                    label 'norm_cpus'

                    publishDir "${params.workingdir}/${params.outdir}/Analyze/Clustering/AminoTypes/SummaryFiles", mode: "copy", overwrite: true, pattern: '*.{clstr,csv,gc}'
                    publishDir "${params.workingdir}/${params.outdir}/Analyze/Clustering/AminoTypes/Problematic", mode: "copy", overwrite: true, pattern: '*problematic*.{fasta}'
                    publishDir "${params.workingdir}/${params.outdir}/Analyze/Clustering/AminoTypes", mode: "copy", overwrite: true, pattern: '*AminoTypes_noTaxonomy.{fasta}'

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
                        awk 'BEGIN{RS=">";ORS=""}length(\$2)>="${params.minAA}"{print ">"\$0}' ${prot} >${params.projtag}_filtered_translations.fasta
                        awk 'BEGIN{RS=">";ORS=""}length(\$2)<"${params.minAA}"{print ">"\$0}' ${prot} >${params.projtag}_problematic_translations.fasta
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
                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/AminoTypes/EMBOSS/HydrophobicMoment", mode: "copy", overwrite: true, pattern: '*HydrophobicMoments.{svg}'
                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/AminoTypes/EMBOSS/IsoelectricPoint", mode: "copy", overwrite: true, pattern: '*IsoelectricPoint.{iep,svg}'
                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/AminoTypes/EMBOSS/ProteinProperties", mode: "copy", overwrite: true, pattern: '*.{pepstats,pepinfo}'
                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/AminoTypes/EMBOSS/ProteinProperties/Plots", mode: "copy", overwrite: true, pattern: '*PropertiesPlot.{svg}'
                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/AminoTypes/EMBOSS/2dStructure/Plots", mode: "copy", overwrite: true, pattern: '*Helical*.{svg}'

                        input:
                            file(prot) from aminotypesEmboss

                        output:
                            tuple file("*.garnier"), file("*HydrophobicMoments.svg"), file("*IsoelectricPoint*"), file("*.pepstats"), file("*PropertiesPlot*"), file("*Helical*")  into amino_emboss

                        script:
                            """
                            name=\$( echo ${prot} | awk -F ".fasta" '{print \$1}')
                            garnier -sequence ${prot} -outfile \${name}_2dStructures.garnier
                            hmoment -seqall ${prot} -graph svg -plot
                            mv hmoment.svg ./"\${name}"_HydrophobicMoments.svg
                            iep -sequence ${prot} -graph svg -plot -outfile "\${name}"_IsoelectricPoint.iep
                            mv iep.svg ./"\${name}"_IsoelectricPoint.svg
                            pepstats -sequence ${prot} -outfile \${name}_ProteinProperties.pepstats
                            grep ">" ${prot} | awk -F ">" '{print \$2}' > tmpsequence.list
                            for x in \$(cat tmpsequence.list);do
                                echo \$x > tmp1.list
                                seqtk subseq ${prot} tmp1.list > tmp2.fasta
                                len=\$(tail -1 tmp2.fasta | awk '{print length}')
                                pepinfo -sequence tmp2.fasta -graph svg -outfile "\$x"_PropertiesPlot.pepinfo
                                mv pepinfo.svg ./"\$x"_PropertiesPlot.svg
                                cat "\$x"_PropertiesPlot.pepinfo >> "\${name}"_PropertiesPlot.pepinfo
                                rm "\$x"_PropertiesPlot.pepinfo
                                pepnet -sask -sequence tmp2.fasta -graph svg -sbegin1 1 -send1 \$len
                                mv pepnet.svg ./"\$x"_HelicalNet.svg
                                pepwheel -sequence tmp2.fasta -graph svg -sbegin1 1 -send1 \$len
                                mv pepwheel.svg ./"\$x"_HelicalWheel.svg
                                rm tmp1.list tmp2.fasta
                            done
                            rm tmpsequence.list
                            """
                    }
                }

                if (!params.skipTaxonomy) {

                  if (params.dbtype == "NCBI") {

                    process AminoType_Taxonomy_Inference_NCBI {

                        label 'high_cpus'

                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/AminoTypes/Taxonomy", mode: "copy", overwrite: true, pattern: '*.{csv,tsv}'
                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/AminoTypes/Taxonomy", mode: "copy", overwrite: true, pattern: '*TaxonomyLabels.fasta'

                        input:
                            file(asvs) from aminotypesBlast

                        output:
                            tuple file("*_phyloformat.csv"), file("*_summaryTable.tsv"), file("*dmd.out") into summary_AA_diamond
                            file("*_summary_for_plot.csv") into taxplot2
                            file("*TaxonomyLabels.fasta") into tax_labeled_fasta2
                            file("*_quick_Taxbreakdown.csv") into tax_table_amino

                        script:
                            """
                            cp ${params.vampdir}/bin/rename_seq.py .
                            virdb=${params.dbdir}/${params.dbname}
                            grep ">" \${virdb} > headers.list
                            headers="headers.list"
                            name=\$( echo ${asvs} | awk -F ".fasta" '{print \$1}')
                            if [[ ${ncbitax} == "true" ]]
                            then   diamond blastp -q ${asvs} -d \${virdb} -p ${task.cpus} --id ${params.minID} -l ${params.minaln} --min-score ${params.bitscore} --${sensitivity} -o "\$name"_dmd.out -f 6 qseqid qlen sseqid qstart qend qseq sseq length qframe evalue bitscore pident btop staxids sskingdoms skingdoms sphylums --max-target-seqs 1 --max-hsps 1
                            else   diamond blastp -q ${asvs} -d \${virdb} -p ${task.cpus} --id ${params.minID} -l ${params.minaln} --min-score ${params.bitscore} --${sensitivity} -o "\$name"_dmd.out -f 6 qseqid qlen sseqid qstart qend qseq sseq length qframe evalue bitscore pident btop --max-target-seqs 1 --max-hsps 1
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
                            if [[ ${ncbitax} == "true" ]]
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
                                            gene=\$(grep -w "\$acc" "\$headers" | awk -F "." '{ print \$2 }' | awk -F "[" '{ print \$1 }' | awk -F " " print substr(\$0, index(\$0,\$2)) | sed 's/ /_/g') &&
                                            echo "\$gene" | sed 's/_/ /g' >> "\$name"_genes.list
                                            virus=\$(grep -w "\$acc" "\$headers" | awk -F "[" '{ print \$2 }' | awk -F "]" '{ print \$1 }'| sed 's/ /_/g')
                                            echo "\$virus" | sed 's/_/ /g' >> "\$name"_virus.list
                                            echo ">"\${s}"_"\$virus"_"\$gene"" >> new_"\$name"_asvnames.txt
                                            if [[ "${params.lca}" == "T" ]]
                                            then    if [[ \$(grep -w "\$acc" ${params.dbanno}/*.txt | wc -l) -eq 1 ]]
                                                    then  group=\$(grep -w "\$acc" ${params.dbanno}/*.txt | awk -F ":" '{print \$1}')
                                                          lcla=\$(grep -w "\$group" lcainfo.list | awk -F "\t" '{print \$2}')
                                                          echo "\$lcla" >> lca_classification.list
                                                    else  echo "Viruses::unclassified" >> lca_classification.list
                                                    fi
                                            fi
                                            if [[ ${ncbitax} == "true" ]]
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
                                            if [[ "${ncbitax}" == "true" ]]
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
                            if [[ "${params.lca}" == "T" && "${ncbitax}" == "true" ]]
                            then
                                  paste -d "," sequence.list "\$name"_virus.list "\$name"_genes.list otu.list lca_classification.list ncbi_classification.list newnames.list length.list bit.list evalue.list pid.list access.list >> "\$name"_phyloformat.csv
                                  paste -d"\t" otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list ncbi_classification.list sequence.list length.list bit.list evalue.list pid.list >> "\$name"_summaryTable.tsv
                                  paste -d"," otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list ncbi_classification.list >> \${name}_quick_Taxbreakdown.csv
                            elif [[ "${params.lca}" == "T" && "${ncbitax}" != "true" ]]
                            then
                                  paste -d "," sequence.list "\$name"_virus.list "\$name"_genes.list otu.list lca_classification.list newnames.list length.list bit.list evalue.list pid.list access.list >> "\$name"_phyloformat.csv
                                  paste -d"\t" otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list sequence.list length.list bit.list evalue.list pid.list >> "\$name"_summaryTable.tsv
                                  paste -d"," otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list >> \${name}_quick_Taxbreakdown.csv
                            elif [[ "${ncbitax}" == "true" && "${params.lca}" != "T"]]
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
                            rm evalue.list sequence.list bit.list pid.list length.list seqids.lst otu.list *asvnames.txt "\$name"_virus.list "\$name"_genes.list newnames.list access.list headers.list
                            """
                        }
                    } else if (params.dbtype== "RVDB") {

                      process AminoType_Taxonomy_Inference_RVDB {

                          label 'high_cpus'

                          publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/AminoTypes/Taxonomy", mode: "copy", overwrite: true, pattern: '*.{csv,tsv}'
                          publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/AminoTypes/Taxonomy", mode: "copy", overwrite: true, pattern: '*TaxonomyLabels.fasta'

                          input:
                              file(asvs) from aminotypesBlast

                          output:
                              tuple file("*_phyloformat.csv"), file("*_summaryTable.tsv"), file("*dmd.out") into summary_AA_diamond
                              file("*_summary_for_plot.csv") into taxplot2
                              file("*TaxonomyLabels.fasta") into tax_labeled_fasta2
                              file("*_quick_Taxbreakdown.csv") into tax_table_amino

                          script:
                              """
                              cp ${params.vampdir}/bin/rename_seq.py .
                              virdb=${params.dbdir}/${params.dbname}
                              grep ">" \${virdb} > headers.list
                              headers="headers.list"
                              name=\$( echo ${asvs} | awk -F ".fasta" '{print \$1}')
                              diamond blastp -q ${asvs} -d \${virdb} -p ${task.cpus} --id ${params.minID} -l ${params.minaln} --min-score ${params.bitscore} --${sensitivity} -o "\$name"_dmd.out -f 6 qseqid qlen sseqid qstart qend qseq sseq length qframe evalue bitscore pident btop --max-target-seqs 1 --max-hsps 1
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
                                              else  echo "Viruses::unclassified" >> lca_classification.list
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
                              rm evalue.list sequence.list bit.list pid.list length.list seqids.lst otu.list *asvnames.txt "\$name"_virus.list "\$name"_genes.list newnames.list access.list headers.list
                              """
                    }
                }
             }

                if (!params.skipPhylogeny) {

                    process AminoType_Phylogeny {

                        label 'norm_cpus'

                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/AminoTypes/Phylogeny/Alignment", mode: "copy", overwrite: true, pattern: '*aln.*'
                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/AminoTypes/Phylogeny/Modeltest", mode: "copy", overwrite: true, pattern: '*mt*'
                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/AminoTypes/Phylogeny/IQ-TREE", mode: "copy", overwrite: true, pattern: '*iq*'

                        input:
                            file(prot) from aminotypesMafft

                        output:
                            tuple file("*_aln.fasta"), file("*_aln.html"), file("*.log"), file("*iq*"), file("*mt*") into alignprot_results
                            file("*iq.treefile") into (amino_rax_plot, amino_repphy)

                        script:
                            """
                            # Protein_Alignment
                            pre=\$(echo ${prot} | awk -F "_noTax" '{print \$1}' )
                            if [[ $(grep -c ">" ${prot}) -gt 499 ]]; then algo="super5"; else algo="mpc"; fi
                            ${tools}/muscle5.0.1278_linux64 -"\${algo}" ${prot} -out \${pre}_ALN.fasta -thread ${task.cpus} -quiet
                            #mafft --thread ${task.cpus} --maxiterate 15000 --auto ${prot} >\${pre}_ALN.fasta
                            trimal -in \${pre}_ALN.fasta -out \${pre}_aln.fasta -keepheader -fasta -automated1 -htmlout \${pre}_aln.html
                            o-trim-uninformative-columns-from-alignment \${pre}_aln.fasta
                            mv \${pre}_aln.fasta-TRIMMED ./\${pre}_Aligned_informativeonly.fasta
                            # Protein_ModelTest
                            modeltest-ng -i \${pre}_Aligned_informativeonly.fasta -p ${task.cpus} -o \${pre}_mt -d aa -s 203 --disable-checkpoint

                            # Protein_Phylogeny
                            if [ "${params.iqCustomaa}" != "" ];then
                                iqtree -s \${pre}_Aligned_informativeonly.fasta --prefix \${pre}_iq --redo -T auto ${params.iqCustomaa}

                            elif [[ "${params.ModelTaa}" != "false" && "${params.nonparametric}" != "false" ]];then
                                mod=\$(tail -12 \${pre}_Aligned_informativeonly.fasta.log | head -1 | awk '{print \$6}')
                                iqtree -s \${pre}_Aligned_informativeonly.fasta --prefix \${pre}_iq -m \${mod} --redo -nt auto -b ${params.boots}

                            elif [[ "${params.ModelTaa}" != "false" && "${params.parametric}" != "false" ]];then
                                mod=\$(tail -12 \${pre}_Aligned_informativeonly.fasta.log | head -1 | awk '{print \$6}')
                                iqtree -s \${pre}_Aligned_informativeonly.fasta --prefix \${pre}_iq -m \${mod} --redo -nt auto -bb ${params.boots} -bnni

                            elif [ "${params.nonparametric}" != "false" ];then
                                iqtree -s \${pre}_Aligned_informativeonly.fasta --prefix \${pre}_iq -m MFP --redo -nt auto -b ${params.boots}

                            elif [ "${params.parametric}" != "false" ];then
                                iqtree -s\${pre}_Aligned_informativeonly.fasta --prefix \${pre}_iq -m MFP --redo -nt auto -bb ${params.boots} -bnni

                            else
                                iqtree -s \${pre}_Aligned_informativeonly.fasta --prefix \${pre}_iq -m MFP --redo -nt auto -bb ${params.boots} -bnni
                            fi
                            """
                    }
                }

                process Generate_AminoTypes_Counts_Table {

                    label 'high_cpus'

                    publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/AminoTypes/Counts", mode: "copy", overwrite: true

                    input:
                        file(fasta) from aminotypesCounts
                        file(merged) from mergeforprotcounts
                        file(samplist) from samplelist

                    output:
                        tuple file("*_AminoType_counts.csv"), file("*dmd.out") into counts_summary
                        file("*_AminoType_counts.csv") into (aminocounts_plot, aminocountmed)

                    script:
                        """
                        set +e
                        diamond makedb --in ${fasta} --db ${fasta}
                        diamond blastx -q ${merged} -d ${fasta} -p ${task.cpus} --min-score ${params.ProtCountsBit} --id ${params.ProtCountID} -l ${params.ProtCountsLength} --${sensitivity} -o ${params.projtag}_protCounts_dmd.out -f 6 qseqid qlen sseqid qstart qend qseq sseq length qframe evalue bitscore pident btop --max-target-seqs 1 --max-hsps 1 --max-hsps 1
                        echo "OTU_ID" >tmp.col1.txt
                        echo "Generating sample id list"
                        grep ">" ${fasta} | awk -F ">" '{print \$2}' | sort | uniq > otuid.list
                        cat otuid.list >> tmp.col1.txt
                        echo "Beginning them counts tho my g"
                        for y in \$( cat ${samplist} );do
                            echo "Starting with \$y now ..."
                            grep "\$y" ${params.projtag}_protCounts_dmd.out > tmp."\$y".out
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

              if (params.aminoMED) {

                    process AminoType_Minimum_Entropy_Decomposition {

                    label 'low_cpus'

                    publishDir "${params.workingdir}/${params.outdir}/Analyze/Clustering/AminoTypes/MED", mode: "copy", overwrite: true

                    input:
                      file(aminos) from aminos_for_med

                    output:
                      file("*_AminoType_Grouping.csv") into atygroupscsv
                      file("${params.projtag}_AminoType_group_reps_aligned.fasta") into atygroupreps

                    script:
                    """
                    #alignment
                    if [[ $(grep -c ">" ${aminos}) -gt 499 ]]; then algo="super5"; else algo="mpc"; fi
                    ${tools}/muscle5.0.1278_linux64 -"\${algo}" ${aminos} -out ${params.projtag}_AminoTypes_muscleAlign.fasta -thread ${task.cpus} -quiet
                    #mafft --thread ${task.cpus} --maxiterate 15000 --auto ${aminos} > ${params.projtag}_AminoTypes_mafftAlign.fasta
                    #trimming
                    trimal -in ${params.projtag}_AminoTypes_muscleAlign.fasta -out ${params.projtag}_AminoTypes_muscleAligned.fasta  -keepheader -fasta -automated1
                    rm ${params.projtag}_AminoTypes_muscleAlign.fasta
                    o-trim-uninformative-columns-from-alignment ${params.projtag}_AminoTypes_muscleAligned.fasta
                    mv ${params.projtag}_AminoTypes_muscleAligned.fasta-TRIMMED ./${params.projtag}_AminoTypes_Aligned_informativeonly.fasta
                    #entopy analysis
                    entropy-analysis ${params.projtag}_AminoTypes_Aligned_informativeonly.fasta
                    #Decomposition
                    if [[ \$(echo ${params.aminoC} | grep -c ",") -eq 1 ]]
                    then
                          tag=$(echo ${params.aminoC} | sed 's/,/_/g')
                          oligotype ${params.projtag}_ASVs_Aligned_informativeonly.fasta ${params.projtag}_ASVs_Aligned_informativeonly.fasta-ENTROPY -o ${params.projtag}_asvMED_"\$tag -M 1 -C ${params.aminoC} -N ${task.cpus} --skip-check-input --no-figures --skip-gen-html
                    else
                          oligotype ${params.projtag}_ASVs_Aligned_informativeonly.fasta ${params.projtag}_ASVs_Aligned_informativeonly.fasta-ENTROPY -o ${params.projtag}_asvMED_${params.aminoC} -M 1 -c ${params.aminoC} -N ${task.cpus} --skip-check-input --no-figures --skip-gen-html
                    fi
                    #generatemaps
                    cd ./${params.projtag}_AminoTypeMED_${params.aminoC}/OLIGO-REPRESENTATIVES/
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
                            seqtk subseq ../OLIGO-REPRESENTATIVES.fasta group.list > group.fasta
                            tail -1 group.fasta >> ${params.projtag}_AminoType_group_reps_aligned.fasta
                            mv "\$gid"_"\$uni" ./Group"\$j"_"\$uni"_aligned.fasta
                            mv "\$gid"_"\$uni"_unique ./Group"\$j"_"\$uni"_unqiues_aligned.fasta
                            rm "\$gid"*.cPickle
                            j=\$((\$j+1))
                    done
                    mv ${params.projtag}_AminoType_Grouping.csv ../../
                    mv ${params.projtag}_AminoType_group_reps_aligned.fasta ../../
                    cd ..

                    """
                    }
                  }
                  if (!params.skipPhylogeny) {

                    process AminoType_MED_Reps_phylogeny {

                    label 'low_cpus'

                    publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/AminoTypes/MED/Phylogeny/Modeltest", mode: "copy", overwrite: true, pattern: '*mt*'
                    publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/AminoTypes/MED/Phylogeny/IQ-TREE", mode: "copy", overwrite: true, pattern: '*iq*'

                    input:
                      file(reps) from atygroupreps

                    output:
                      file("*_AminoType_Group_Reps*") into align_results_aminmed
                      file("*iq.treefile") into amino_group_rep_tree

                    script:
                        """
                        # Protein_ModelTest
                        modeltest-ng -i ${reps} -p ${task.cpus} -o ${params.projtag}_AminoType_Group_Reps_mt -d aa -s 203 --disable-checkpoint

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
                            iqtree -s ${reps} --prefix ${params.projtag}_AminoType_Group_Reps_iq -m MFP --redo -nt auto -b ${params.boots}

                        elif [ "${params.parametric}" != "false" ];then
                            iqtree -s ${reps} --prefix ${params.projtag}_AminoType_Group_Reps_iq -m MFP --redo -nt auto -bb ${params.boots} -bnni

                        else
                            iqtree -s ${reps} --prefix ${params.projtag}_AminoType_Group_Reps_iq -m MFP --redo -nt auto -bb ${params.boots} -bnni
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
                      file("${params.projtag}_AminoType_Groupingcounts.csv") into amino_groupcounts

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
            }

            if (params.pcASV) {        // ASV_nucl -> ASV_aa -> clusteraa by %id with ch-hit -> extract representative nucl sequences to generate new OTU file

                process Translation_For_pcASV_Generation {

                      label 'low_cpus'

                      publishDir "${params.workingdir}/${params.outdir}/Analyze/Clustering/pcASV/Translation", mode: "copy", overwrite: true, pattern: '*_ASV_translations*'

                      input:
                          file(fasta) from nucl2aa

                      output:
                          file("*ASV*translations.fasta") into clustering_aa
                          file("*_ASV_translations_report") into reportaa_VR
                          file("*_ASV_nucleotide.fasta") into asvfastaforaaclust

                      script:
                          """
                          ${tools}/virtualribosomev2/dna2pep.py ${fasta} -r all -x -o none --fasta ${params.projtag}_ASV_translations.fasta --report ${params.projtag}_ASV_translations_report
                          cp ${fasta} ${params.projtag}_ASV_nucleotide.fasta
                          """
                }

                process Generate_pcASVs {

                    label 'norm_cpus'

                    tag "${mtag}"

                    publishDir "${params.workingdir}/${params.outdir}/Analyze/Clustering/pcASV", mode: "copy", overwrite: true, pattern: '*pcASV*.{fasta}'
                    publishDir "${params.workingdir}/${params.outdir}/Analyze/Clustering/pcASV/SummaryFiles", mode: "copy", overwrite: true, pattern: '*.{clstr,csv,gc}'
                    publishDir "${params.workingdir}/${params.outdir}/Analyze/Clustering/pcASV/Problematic", mode: "copy", overwrite: true, pattern: '*problem*.{fasta}'

                    input:
                        each x from 1..naa
                        file(fasta) from clustering_aa
                        file(asvs) from asvfastaforaaclust

                    output:
                        tuple nid, file("${params.projtag}_nucleotide_pcASV*.fasta") into ( pcASV_ntDiamond_ch, pcASV_nt_counts_ch, pcASV_ntmatrix_ch, pcASV_ntmafft_ch )
                        tuple nid, file("*_aminoacid_pcASV*_noTaxonomy.fasta") into ( pcASV_aaMatrix_ch, pcASV_aaDiamond_ch, pcASV_aaMafft_ch, pcASV_aaCounts_ch, pcASVEMBOSS )
                        tuple nid, file("*.fasta"), file("*.clstr"), file("*.csv"), file("*.gc") into ( pcASVsupplementalfiles )

                    script:
                        // add awk script to count seqs
                        nid=slist2.get(x-1)
                        mtag="ID=" + slist2.get(x-1)
                        """
                        set +e
                        cp ${params.vampdir}/bin/rename_seq.py .
                        awk 'BEGIN{RS=">";ORS=""}length(\$2)>="${params.minAA}"{print ">"\$0}' ${fasta} > ${params.projtag}_filtered_proteins.fasta
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
                        awk 'BEGIN{RS=">";ORS=""}length(\$2)<"${params.minAA}"{print ">"\$0}' ${fasta} >${params.projtag}_pcASV${nid}_problematic_translations.fasta
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

                  if (params.dbtype == "NCBI") {

                    process pcASV_Nucleotide_Taxonomy_Inference_NCBI {

                        label 'high_cpus'

                        tag "${mtag}"

                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Nucleotide/Taxonomy/SummaryFiles", mode: "copy", overwrite: true, pattern: '*.{csv,tsv}'
                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Nucleotide/Taxonomy/DiamondOutput", mode: "copy", overwrite: true, pattern: '*dmd.{out}'
                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Nucleotide/Taxonomy", mode: "copy", overwrite: true, pattern: '*.{fasta}'

                        input:
                            tuple nid, file(asvs) from pcASV_ntDiamond_ch

                        output:
                            file("*.fasta") into ( pcASV_labeled )
                            tuple file("*_phyloformat.csv"), file("*_summaryTable.tsv"), file("*dmd.out") into summary_AAdiamond
                            tuple nid, file("*_summary_for_plot.csv") into taxplot3
                            tuple nid, file("*_quick_Taxbreakdown.csv") into tax_table_pcasvnt

                        script:
                            mtag="ID=" + nid
                            """
                            set +e
                            cp ${params.vampdir}/bin/rename_seq.py .
                            virdb=${params.dbdir}/${params.dbname}
                            grep ">" \${virdb} > headers.list
                            headers="headers.list"
                            name=\$( echo ${asvs} | awk -F ".fasta" '{print \$1}')
                            if [[ ${ncbitax} == "true" ]]
                            then   diamond blastx -q ${asvs} -d \${virdb} -p ${task.cpus} --id ${params.minID} -l ${params.minaln} --min-score ${params.bitscore} --${sensitivity} -o "\$name"_dmd.out -f 6 qseqid qlen sseqid qstart qend qseq sseq length qframe evalue bitscore pident btop staxids sskingdoms skingdoms sphylums --max-target-seqs 1 --max-hsps 1
                            else   diamond blastx -q ${asvs} -d \${virdb} -p ${task.cpus} --id ${params.minID} -l ${params.minaln} --min-score ${params.bitscore} --${sensitivity} -o "\$name"_dmd.out -f 6 qseqid qlen sseqid qstart qend qseq sseq length qframe evalue bitscore pident btop --max-target-seqs 1 --max-hsps 1
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
                            if [[ ${ncbitax} == "true" ]]
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
                                            gene=\$(grep -w "\$acc" "\$headers" | awk -F "." '{ print \$2 }' | awk -F "[" '{ print \$1 }' | awk -F " " print substr(\$0, index(\$0,\$2)) | sed 's/ /_/g') &&
                                            echo "\$gene" | sed 's/_/ /g' >> "\$name"_genes.list
                                            virus=\$(grep -w "\$acc" "\$headers" | awk -F "[" '{ print \$2 }' | awk -F "]" '{ print \$1 }'| sed 's/ /_/g')
                                            echo "\$virus" | sed 's/_/ /g' >> "\$name"_virus.list
                                            echo ">"\${s}"_"\$virus"_"\$gene"" >> new_"\$name"_asvnames.txt
                                            if [[ "${params.lca}" == "T" ]]
                                            then    if [[ \$(grep -w "\$acc" ${params.dbanno}/*.txt | wc -l) -eq 1 ]]
                                                    then  group=\$(grep -w "\$acc" ${params.dbanno}/*.txt | awk -F ":" '{print \$1}')
                                                          lcla=\$(grep -w "\$group" lcainfo.list | awk -F "\t" '{print \$2}')
                                                          echo "\$lcla" >> lca_classification.list
                                                    else  echo "Viruses::unclassified" >> lca_classification.list
                                                    fi
                                            fi
                                            if [[ ${ncbitax} == "true" ]]
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
                                            if [[ "${ncbitax}" == "true" ]]
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
                            if [[ "${params.lca}" == "T" && "${ncbitax}" == "true" ]]
                            then
                                  paste -d "," sequence.list "\$name"_virus.list "\$name"_genes.list otu.list lca_classification.list ncbi_classification.list newnames.list length.list bit.list evalue.list pid.list access.list >> "\$name"_phyloformat.csv
                                  paste -d"\t" otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list ncbi_classification.list sequence.list length.list bit.list evalue.list pid.list >> "\$name"_summaryTable.tsv
                                  paste -d"," otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list ncbi_classification.list >> \${name}_quick_Taxbreakdown.csv
                            elif [[ "${params.lca}" == "T" && "${ncbitax}" != "true" ]]
                            then
                                  paste -d "," sequence.list "\$name"_virus.list "\$name"_genes.list otu.list lca_classification.list newnames.list length.list bit.list evalue.list pid.list access.list >> "\$name"_phyloformat.csv
                                  paste -d"\t" otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list sequence.list length.list bit.list evalue.list pid.list >> "\$name"_summaryTable.tsv
                                  paste -d"," otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list >> \${name}_quick_Taxbreakdown.csv
                            elif [[ "${ncbitax}" == "true" && "${params.lca}" != "T"]]
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
                            rm evalue.list sequence.list bit.list pid.list length.list seqids.lst otu.list *asvnames.txt "\$name"_virus.list "\$name"_genes.list newnames.list access.list headers.list
                            """
                        }
                    } else if (params.dbtype== "RVDB") {

                      process pcASV_Nucleotide_Taxonomy_Inference_RVDB {

                          label 'high_cpus'

                          tag "${mtag}"

                          publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Nucleotide/Taxonomy/SummaryFiles", mode: "copy", overwrite: true, pattern: '*.{csv,tsv}'
                          publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Nucleotide/Taxonomy/DiamondOutput", mode: "copy", overwrite: true, pattern: '*dmd.{out}'
                          publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Nucleotide/Taxonomy", mode: "copy", overwrite: true, pattern: '*.{fasta}'

                          input:
                              tuple nid, file(asvs) from pcASV_ntDiamond_ch

                          output:
                              file("*.fasta") into ( pcASV_labeled )
                              tuple file("*_phyloformat.csv"), file("*_summaryTable.tsv"), file("*dmd.out") into summary_AAdiamond
                              tuple nid, file("*_summary_for_plot.csv") into taxplot3
                              tuple nid, file("*_quick_Taxbreakdown.csv") into tax_table_pcasvnt

                          script:
                              mtag="ID=" + nid
                              """
                              set +e
                              cp ${params.vampdir}/bin/rename_seq.py .
                              virdb=${params.dbdir}/${params.dbname}
                              grep ">" \${virdb} > headers.list
                              headers="headers.list"
                              name=\$( echo ${asvs} | awk -F ".fasta" '{print \$1}')
                              diamond blastx -q ${asvs} -d \${virdb} -p ${task.cpus} --id ${params.minID} -l ${params.minaln} --min-score ${params.bitscore} --${sensitivity} -o "\$name"_dmd.out -f 6 qseqid qlen sseqid qstart qend qseq sseq length qframe evalue bitscore pident btop --max-target-seqs 1 --max-hsps 1
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
                                              else  echo "Viruses::unclassified" >> lca_classification.list
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
                              rm evalue.list sequence.list bit.list pid.list length.list seqids.lst otu.list *asvnames.txt "\$name"_virus.list "\$name"_genes.list newnames.list access.list headers.list
                              """
                    }
                }
             }


                process Generate_Nucleotide_pcASV_Counts {

                    label 'norm_cpus'

                    tag "${mtag}"

                    publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Nucleotide/Counts", mode: "copy", overwrite: true, pattern: '*.{biome,csv,txt}'

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
                    	vsearch --usearch_global ${merged} --db ${potus} --id .${nid} --threads ${task.cpus} --otutabout \${name}_counts.txt --biomout \${name}_counts.biome
                    	cat \${name}_counts.txt | tr "\t" "," >\${name}_count.csv
                    	sed 's/#OTU ID/OTU_ID/g' \${name}_count.csv >\${name}_counts.csv
                    	rm \${name}_count.csv
                        """
                }

                process Generate_pcASV_Nucleotide_Matrix {

                    label 'low_cpus'

                    tag "${mtag}"

                    publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Nucleotide/Matrix", mode: "copy", overwrite: true

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
                        cat \${name}_PercentIDq.matrix | tr " " "," | grep "," >\${name}_PercentID.matrix
                        rm \${name}_PercentIDq.matrix
                        """
                }

                if (!params.skipPhylogeny) {

                    process pcASV_Nucleotide_Phylogeny {

                        label 'norm_cpus'

                        tag "${mtag}"

                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Nucleotide/Phylogeny/Alignment", mode: "copy", overwrite: true, pattern: '*aln.*'
                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Nucleotide/Phylogeny/ModelTest", mode: "copy", overwrite: true, pattern: '*mt*'
                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Nucleotide/Phylogeny/IQ-TREE", mode: "copy", overwrite: true, pattern: '*iq*'

                        input:
                            tuple nid, file(prots) from pcASV_ntmafft_ch

                        output:
                            tuple file("*_aln.fasta"), file("*_aln.html"), file("*.tree"), file("*.log"), file("*iq*"), file("*mt*") into pcASV_nucleotide_phylogeny_results
                            tuple nid, file("*iq.treefile") into potu_Ntree_plot

                        script:
                            mtag="ID=" + nid
                            """
                            pre=\$( echo ${prots} | awk -F "_noTax" '{print \$1}' )
                            if [[ $(grep -c ">" ${prots}) -gt 499 ]]; then algo="super5"; else algo="mpc"; fi
                            ${tools}/muscle5.0.1278_linux64 -"\${algo}" ${prots} -out \${pre}_ALN.fasta -thread ${task.cpus} -quiet
                            #mafft --maxiterate 5000 --auto ${reads} >\${pre}_ALN.fasta
                            trimal -in \${pre}_ALN.fasta -out \${pre}_aln.fasta -keepheader -fasta -automated1 -htmlout \${pre}_aln.html
                            o-trim-uninformative-columns-from-alignment \${pre}_aln.fasta
                            mv \${pre}_aln.fasta-TRIMMED ./\${pre}_Aligned_informativeonly.fasta
                            # pcASV_Nucleotide_ModelTest
                            modeltest-ng -i \${pre}_Aligned_informativeonly.fasta -p ${task.cpus} -o \${pre}_noTaxonomy_mt -d nt -s 203 --disable-checkpoint

                            # pcASV_Nucleotide_Phylogeny
                            if [ "${params.iqCustomnt}" != "" ];then
                                iqtree -s \${pre}_Aligned_informativeonly.fasta --prefix \${pre}_noTaxonomy_iq --redo -T auto ${params.iqCustomnt}

                            elif [[ "${params.ModelTnt}" != "false" && "${params.nonparametric}" != "false" ]];then
                                mod=\$(tail -12 \${pre}_Aligned_informativeonly.fasta.log | head -1 | awk '{print \$6}')
                                iqtree -s \${pre}_Aligned_informativeonly.fasta --prefix \${pre}_noTaxonomy_iq -m \${mod} --redo-nt auto -b ${params.boots}

                            elif [[ "${params.ModelTnt}" != "false" && "${params.parametric}" != "false" ]];then
                                mod=\$(tail -12 \${pre}_Aligned_informativeonly.fasta.log | head -1 | awk '{print \$6}')
                                iqtree -s \${pre}_Aligned_informativeonly.fasta --prefix \${pre}_noTaxonomy_iq -m \${mod} --redo -nt auto -bb ${params.boots} -bnni

                            elif [ "${params.nonparametric}" != "false" ];then
                                iqtree -s \${pre}_Aligned_informativeonly.fasta --prefix \${pre}_noTaxonomy_iq -m MFP --redo -nt auto -b ${params.boots}

                            elif [ "${params.parametric}" != "false" ];then
                                iqtree -s \${pre}_Aligned_informativeonly.fasta --prefix \${pre}_noTaxonomy_iq -m MFP --redo -nt auto -bb ${params.boots} -bnni

                            else
                                iqtree -s \${pre}_Aligned_informativeonly.fasta --prefix \${pre}_noTaxonomy_iq -m MFP --redo -nt auto -bb ${params.boots} -bnni
                            fi
                            """
                    }
                }

                process pcASV_AminoAcid_Matrix {

                    label 'low_cpus'

                    tag "${mtag}"

                    publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Aminoacid/Matrix", mode: "copy", overwrite: true

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
                        cat \${name}_PercentIDq.matrix | tr " " "," | grep "," >\${name}_PercentID.matrix
                        rm \${name}_PercentIDq.matrix
                        """
                }

                if (!params.skipEMBOSS) {

                    process pcASV_EMBOSS_Analyses {

                        label 'low_cpus'

                        tag "${mtag}"

                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Aminoacid/EMBOSS/2dStructure", mode: "copy", overwrite: true, pattern: '*.{garnier}'
                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Aminoacid/EMBOSS/HydrophobicMoment", mode: "copy", overwrite: true, pattern: '*HydrophobicMoments.{svg}'
                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Aminoacid/EMBOSS/IsoelectricPoint", mode: "copy", overwrite: true, pattern: '*IsoelectricPoint.{iep,svg}'
                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Aminoacid/EMBOSS/ProteinProperties", mode: "copy", overwrite: true, pattern: '*.{pepstats,pepinfo}'
                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Aminoacid/EMBOSS/ProteinProperties/Plots", mode: "copy", overwrite: true, pattern: '*PropertiesPlot.{svg}'
                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Aminoacid/EMBOSS/2dStructure/Plots", mode: "copy", overwrite: true, pattern: '*Helical*.{svg}'

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
                            hmoment -seqall ${prot} -graph svg -plot
                            mv hmoment.svg ./"\${name}"_HydrophobicMoments.svg
                            iep -sequence ${prot} -graph svg -plot -outfile "\${name}"_IsoelectricPoint.iep
                            mv iep.svg ./"\${name}"_IsoelectricPoint.svg
                            pepstats -sequence ${prot} -outfile \${name}_ProteinProperties.pepstats
                            grep ">" ${prot} | awk -F ">" '{print \$2}' > tmpsequence.list
                            for x in \$(cat tmpsequence.list);do
                                echo \$x > tmp1.list
                                seqtk subseq ${prot} tmp1.list > tmp2.fasta
                                len=\$(tail -1 tmp2.fasta | awk '{print length}')
                                pepinfo -sequence tmp2.fasta -graph svg -outfile "\$x"_PropertiesPlot.pepinfo
                                mv pepinfo.svg ./"\$x"_PropertiesPlot.svg
                                cat "\$x"_PropertiesPlot.pepinfo >> "\${name}"_PropertiesPlot.pepinfo
                                rm "\$x"_PropertiesPlot.pepinfo
                                pepnet -sask -sequence tmp2.fasta -graph svg -sbegin1 1 -send1 \$len
                                mv pepnet.svg ./"\$x"_HelicalNet.svg
                                pepwheel -sequence tmp2.fasta -graph svg -sbegin1 1 -send1 \$len
                                mv pepwheel.svg ./"\$x"_HelicalWheel.svg
                                rm tmp1.list tmp2.fasta
                            done
                            rm tmpsequence.list
                            """
                        }
                }

                if (!params.skipTaxonomy) {

                  if (params.dbtype == "NCBI") {

                    process pcASV_AminoAcid_Taxonomy_Inference_NCBI {

                        label 'high_cpus'

                        tag "${mtag}"

                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Aminoacid/Taxonomy/SummaryFiles", mode: "copy", overwrite: true, pattern: '*.{csv,tsv}'
                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Aminoacid/Taxonomy/DiamondOutput", mode: "copy", overwrite: true, pattern: '*dmd.{out}'
                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Aminoacid/Taxonomy", mode: "copy", overwrite: true, pattern: '*.{fasta}'

                        input:
                            tuple nid, file(asvs) from pcASV_aaDiamond_ch

                        output:
                            file("*.fasta") into ( pcASV_labeledAA )
                            tuple file("*phyloformat.csv"), file("*summaryTable.tsv"), file("*dmd.out") into summary_potuaadiamond
                            tuple nid, file("*_summary_for_plot.csv") into taxplot4
                            tuple nid, file("*_quick_Taxbreakdown.csv") into tax_table_pcasvaa

                        script:
                            mtag="ID=" + nid
                            """
                            cp ${params.vampdir}/bin/rename_seq.py .
                            virdb=${params.dbdir}/${params.dbname}
                            grep ">" \${virdb} > headers.list
                            headers="headers.list"
                            name=\$( echo ${asvs} | awk -F ".fasta" '{print \$1}')
                            if [[ ${ncbitax} == "true" ]]
                            then   diamond blastp -q ${asvs} -d \${virdb} -p ${task.cpus} --id ${params.minID} -l ${params.minaln} --min-score ${params.bitscore} --${sensitivity} -o "\$name"_dmd.out -f 6 qseqid qlen sseqid qstart qend qseq sseq length qframe evalue bitscore pident btop staxids sskingdoms skingdoms sphylums --max-target-seqs 1 --max-hsps 1
                            else   diamond blastp -q ${asvs} -d \${virdb} -p ${task.cpus} --id ${params.minID} -l ${params.minaln} --min-score ${params.bitscore} --${sensitivity} -o "\$name"_dmd.out -f 6 qseqid qlen sseqid qstart qend qseq sseq length qframe evalue bitscore pident btop --max-target-seqs 1 --max-hsps 1
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
                            if [[ ${ncbitax} == "true" ]]
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
                                            gene=\$(grep -w "\$acc" "\$headers" | awk -F "." '{ print \$2 }' | awk -F "[" '{ print \$1 }' | awk -F " " print substr(\$0, index(\$0,\$2)) | sed 's/ /_/g') &&
                                            echo "\$gene" | sed 's/_/ /g' >> "\$name"_genes.list
                                            virus=\$(grep -w "\$acc" "\$headers" | awk -F "[" '{ print \$2 }' | awk -F "]" '{ print \$1 }'| sed 's/ /_/g')
                                            echo "\$virus" | sed 's/_/ /g' >> "\$name"_virus.list
                                            echo ">"\${s}"_"\$virus"_"\$gene"" >> new_"\$name"_asvnames.txt
                                            if [[ "${params.lca}" == "T" ]]
                                            then    if [[ \$(grep -w "\$acc" ${params.dbanno}/*.txt | wc -l) -eq 1 ]]
                                                    then  group=\$(grep -w "\$acc" ${params.dbanno}/*.txt | awk -F ":" '{print \$1}')
                                                          lcla=\$(grep -w "\$group" lcainfo.list | awk -F "\t" '{print \$2}')
                                                          echo "\$lcla" >> lca_classification.list
                                                    else  echo "Viruses::unclassified" >> lca_classification.list
                                                    fi
                                            fi
                                            if [[ ${ncbitax} == "true" ]]
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
                                            if [[ "${ncbitax}" == "true" ]]
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
                            if [[ "${params.lca}" == "T" && "${ncbitax}" == "true" ]]
                            then
                                  paste -d "," sequence.list "\$name"_virus.list "\$name"_genes.list otu.list lca_classification.list ncbi_classification.list newnames.list length.list bit.list evalue.list pid.list access.list >> "\$name"_phyloformat.csv
                                  paste -d"\t" otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list ncbi_classification.list sequence.list length.list bit.list evalue.list pid.list >> "\$name"_summaryTable.tsv
                                  paste -d"," otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list ncbi_classification.list >> \${name}_quick_Taxbreakdown.csv
                            elif [[ "${params.lca}" == "T" && "${ncbitax}" != "true" ]]
                            then
                                  paste -d "," sequence.list "\$name"_virus.list "\$name"_genes.list otu.list lca_classification.list newnames.list length.list bit.list evalue.list pid.list access.list >> "\$name"_phyloformat.csv
                                  paste -d"\t" otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list sequence.list length.list bit.list evalue.list pid.list >> "\$name"_summaryTable.tsv
                                  paste -d"," otu.list access.list "\$name"_virus.list "\$name"_genes.list lca_classification.list >> \${name}_quick_Taxbreakdown.csv
                            elif [[ "${ncbitax}" == "true" && "${params.lca}" != "T"]]
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
                            rm evalue.list sequence.list bit.list pid.list length.list seqids.lst otu.list *asvnames.txt "\$name"_virus.list "\$name"_genes.list newnames.list access.list headers.list
                            """
                        }
                    } else if (params.dbtype== "RVDB") {

                      process pcASV_AminoAcid_Taxonomy_Inference_RVDB {

                          label 'high_cpus'

                          tag "${mtag}"

                          publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Aminoacid/Taxonomy/SummaryFiles", mode: "copy", overwrite: true, pattern: '*.{csv,tsv}'
                          publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Aminoacid/Taxonomy/DiamondOutput", mode: "copy", overwrite: true, pattern: '*dmd.{out}'
                          publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Aminoacid/Taxonomy", mode: "copy", overwrite: true, pattern: '*.{fasta}'

                          input:
                              tuple nid, file(asvs) from pcASV_aaDiamond_ch

                          output:
                              file("*.fasta") into ( pcASV_labeledAA )
                              tuple file("*phyloformat.csv"), file("*summaryTable.tsv"), file("*dmd.out") into summary_potuaadiamond
                              tuple nid, file("*_summary_for_plot.csv") into taxplot4
                              tuple nid, file("*_quick_Taxbreakdown.csv") into tax_table_pcasvaa

                          script:
                              mtag="ID=" + nid
                              """
                              cp ${params.vampdir}/bin/rename_seq.py .
                              virdb=${params.dbdir}/${params.dbname}
                              grep ">" \${virdb} > headers.list
                              headers="headers.list"
                              name=\$( echo ${asvs} | awk -F ".fasta" '{print \$1}')
                              diamond blastp -q ${asvs} -d \${virdb} -p ${task.cpus} --id ${params.minID} -l ${params.minaln} --min-score ${params.bitscore} --${sensitivity} -o "\$name"_dmd.out -f 6 qseqid qlen sseqid qstart qend qseq sseq length qframe evalue bitscore pident btop --max-target-seqs 1 --max-hsps 1
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
                                              else  echo "Viruses::unclassified" >> lca_classification.list
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
                              rm evalue.list sequence.list bit.list pid.list length.list seqids.lst otu.list *asvnames.txt "\$name"_virus.list "\$name"_genes.list newnames.list access.list headers.list
                              """
                    }
                }
            }

                if (!params.skipPhylogeny) {

                    process pcASV_Protein_Phylogeny {

                        label 'norm_cpus'

                        tag "${mtag}"

                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Aminoacid/Phylogeny/Alignment", mode: "copy", overwrite: true, pattern: '*aln.*'
                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Aminoacid/Phylogeny/Modeltest", mode: "copy", overwrite: true, pattern: '*mt*'
                        publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Aminoacid/Phylogeny/IQ-TREE", mode: "copy", overwrite: true, pattern: '*iq*'

            	        input:
                            tuple nid, file(prot) from pcASV_aaMafft_ch

                        output:
                            tuple file("*_aln.fasta"), file("*_aln.html"), file("*.tree"), file("*.log"), file("*iq*"), file("*mt*") into pcASV_protein_phylogeny_results
                            tuple nid, file("*iq.treefile") into potu_Atree_plot

                        script:
                            mtag="ID=" + nid
                            """
                            pre=\$( echo ${prot} | awk -F ".fasta" '{print \$1}' )
                            if [[ $(grep -c ">" ${prots}) -gt 499 ]]; then algo="super5"; else algo="mpc"; fi
                            ${tools}/muscle5.0.1278_linux64 -"\${algo}" ${prots} -out \${pre}_ALN.fasta -thread ${task.cpus} -quiet
                            #mafft --maxiterate 5000 --auto ${prot} >\${pre}_ALN.fasta
                            trimal -in \${pre}_ALN.fasta -out \${pre}_aln.fasta -keepheader -fasta -automated1 -htmlout \${pre}_aln.html
                            o-trim-uninformative-columns-from-alignment \${pre}_aln.fasta
                            mv \${pre}_aln.fasta-TRIMMED ./\${pre}_Aligned_informativeonly.fasta
                            # pcASV_Protein_ModelTest
                            modeltest-ng -i \${pre}_Aligned_informativeonly.fasta -p ${task.cpus} -o \${pre}_mt -d aa -s 203 --disable-checkpoint

                            # pcASV_Protein_Phylogeny
                            if [ "${params.iqCustomaa}" != "" ];then
                                iqtree -s \${pre}_Aligned_informativeonly.fasta --prefix \${pre}_iq --redo -T auto ${params.iqCustomaa}

                            elif [[ "${params.ModelTaa}" != "false" && "${params.nonparametric}" != "false" ]];then
                                mod=\$(tail -12 \${pre}_Aligned_informativeonly.fasta.log | head -1 | awk '{print \$6}')
                                iqtree -s \${pre}_Aligned_informativeonly.fasta --prefix \${pre}_iq -m \${mod} --redo  -nt auto -b ${params.boots}

                            elif [[ "${params.ModelTaa}" != "false" && "${params.parametric}" != "false" ]];then
                                mod=\$(tail -12 \${pre}_Aligned_informativeonly.fasta.log | head -1 | awk '{print \$6}')
                                iqtree -s \${pre}_Aligned_informativeonly.fasta --prefix \${pre}_iq -m \${mod} --redo -nt auto -bb ${params.boots} -bnni

                            elif [ "${params.nonparametric}" != "false" ];then
                                iqtree -s \${pre}_Aligned_informativeonly.fasta --prefix \${pre}_iq -m MFP --redo -nt auto -b ${params.boots}

                            elif [ "${params.parametric}" != "false" ];then
                                iqtree -s \${pre}_Aligned_informativeonly.fasta --prefix \${pre}_iq -m MFP --redo -nt auto -bb ${params.boots} -bnni

                            else
                                iqtree -s \${pre}_Aligned_informativeonly.fasta --prefix \${pre}_iq -m MFP --redo -nt auto -bb ${params.boots} -bnni
                            fi
                            """
                    }
                }

                process Generate_pcASV_Protein_Counts {

                    label 'high_cpus'

                    tag "${mtag}"

                    publishDir "${params.workingdir}/${params.outdir}/Analyze/Analyses/pcASV/Aminoacid/Counts", mode: "copy", overwrite: true

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
                        diamond blastx -q ${merged} -d ${fasta} -p ${task.cpus} --min-score ${params.ProtCountsBit} --id ${params.ProtCountID} -l ${params.ProtCountsLength} --${sensitivity} -o ${params.projtag}_\${potu}_Counts_dmd.out -f 6 qseqid qlen sseqid qstart qend qseq sseq length qframe evalue bitscore pident btop --max-target-seqs 1 --max-hsps 1 --max-hsps 1
                        echo "OTU_ID" >tmp.col1.txt
                        echo "Generating sample id list"
                        grep ">" ${fasta} | awk -F ">" '{print \$2}' | sort | uniq > otuid.list
                        cat otuid.list >> tmp.col1.txt
                        echo "Beginning them counts tho my g"
                        for y in \$( cat ${samplist} );do
                            echo "Starting with \$y now ..."
                            grep "\$y" ${params.projtag}_\${potu}_Counts_dmd.out > tmp."\$y".out
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

                //NEW REPORT !!!!!!!!!!!!!!!!!
                /*Report_ASV
                asv_counts_plots -> ${params.projtag}_ASV_counts.csv
                taxplot1 -> ${params.projtag}_ASV_summary_for_plot.csv
                asv_heatmap -> ${params.projtag}_ASV_PercentID.matrix
                nucl_phyl_plot -> ${params.projtag}_ASV_iq.treefile
                file("*_ASV_Grouping.csv") into asvgroupscsv
                "${params.projtag}_ASV_Groupingcounts.csv") into asvgroupcounts
                *_quick_Taxbreakdown.csv") into tax_table_asv
                \\${params.projtag}_ASV_Group_Reps_iq.treefile
                */

                report_asv = Channel.create()
                asv_counts_plots.mix(taxplot_asv, asv_heatmap, nucl_phyl_plot_asv, asvgroupscsv, asvgroupcounts, asv_group_rep_tree, tax_table_asv).flatten().buffer(size:8).dump(tag:'asv').into(report_asv)

                if (params.ncASV) {
                    report_ncasv = Channel.create()
                    notu_counts_plots.mix(taxplot_ncasv, notu_heatmap, nucl_phyl_plot_ncasv, tax_table_ncasv).groupTuple(by:0, size:5).dump(tag:'ncasv').into(report_ncasv)
                    /*
                    notu_counts_plots -> ${params.projtag}_ncASV${id}_counts.csv
                    taxplot1a -> ${params.projtag}_ncASV${id}_summary_for_plot.csv
                    notu_heatmap -> ${params.projtag}_ncASV${id}_PercentID.matrix
                    nucl_phyl_plot -> ${params.projtag}_ncASV${id}_iq.treefile
                                      ${params.projtag}_ncASV${id}_quick_Taxbreakdown.csv
                    */
                } else {
                    report_ncasv = Channel.empty()
                }

                if (params.pcASV) {
                    report_pcasv_aa = Channel.create()
                    potu_Acounts.mix(taxplot4, potu_aa_heatmap, potu_Atree_plot, tax_table_pcasvaa).groupTuple(by:0, size:5).dump(tag:'pcasv1').into(report_pcasv_aa)
                    /*Report_pcASV_AminoAcid
                    potu_Acounts -> ${params.projtag}_pcASV${id}_noTaxonomy_counts.csv
                    taxplot4 -> ${params.projtag}_aminoacid_pcASV${id}_noTaxonomy_summary_for_plot.csv
                    potu_aa_heatmap -> ${params.projtag}_aminoacid_pcASV${id}_noTaxonomy_PercentID.matrix
                    potu_Atree_plot -> ${params.projtag}_aminoacid_pcASV${id}_noTaxonomy_iq.treefile
                    tax_table_pcasvaa -> ${params.projtag}_aminoacid_pcASV${id}_quick_Taxbreakdown.csv
                    */
                    report_pcasv_nucl = Channel.create()
                    potu_Ncounts_for_report.mix(taxplot3, potu_nucl_heatmap, potu_Ntree_plot, tax_table_pcasvnt).groupTuple(by:0, size:5).dump(tag:'pcasv2').into(report_pcasv_nucl)
                    /*Report_pcASV_Nucleotide
                    potu_Ncounts_for_report -> ${params.projtag}_nucleotide_pcASV${id}_noTaxonomy_counts.csv
                    taxplot3 -> ${params.projtag}_nucleotide_pcASV${id}_noTaxonomy_summary_for_plot.csv
                    potu_nucl_heatmap -> ${params.projtag}_nucleotide_pcASV${id}_noTaxonomy_PercentID.matrix
                    potu_Ntree_plot -> ${params.projtag}_nucleotide_pcASV${id}_noTaxonomy_iq.treefile
                    tax_table_pcasvnt -> ${params.projtag}_nucleotide_pcASV${id}_quick_Taxbreakdown.csv
                    */
                } else {
                     report_pcasv_aa = Channel.empty()
                     report_pcasv_nucl = Channel.empty()
                }

                if (!params.skipAminoTyping) {
                    report_aminotypes = Channel.create()
                    aminocounts_plot.mix(taxplot2, aminotype_heatmap, amino_rax_plot, atygroupscsv, amino_group_rep_tree, amino_groupcounts, tax_table_amino).flatten().buffer(size:8).dump(tag:'amino').into(report_aminotypes)
                    /*
                    Report_AminoTypes
                    aminocounts_plot -> ${params.projtag}_AminoType_counts.csv
                    taxplot2 -> ${params.projtag}_AminoTypes_summary_for_plot.csv
                    aminotype_heatmap -> ${params.projtag}_AminoTypes_PercentID.matrix
                    amino_rax_plot -> ${params.projtag}_AminoTypes_iq.treefile
                    atygroupscsv  -> *_AminoType_Grouping.csv
                    amino_group_rep_tree -> ${params.projtag}_AminoType_Group_Reps_iq.treefile
                                         params.projtag}_AminoType_Groupingcounts.csv") into amino_groupcounts
                    *_quick_Taxbreakdown.csv") into tax_table_amino
                    */
                } else {
                    report_aminotypes = Channel.empty()
                }

                report_all_ch = Channel.create()
                report_asv.mix(report_ncasv, report_pcasv_aa, report_pcasv_nucl, report_aminotypes).map{it.flatten()}.dump(tag:'report').into(report_all_ch)

                process Report {

                    label 'norm_cpus'

                    publishDir "${params.workingdir}/${params.outdir}/Analyze/FinalReports", mode: "copy", overwrite: true

                    input:
                        file(csv) from fastp_csv_in
                        file(files) from report_all_ch

                    output:
                        file("*.html") into report_all_out

                    script:
                        """
                        name=\$( ls *summary_for_plot.csv | awk -F "_summary_for_plot.csv" '{print \$1}')
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
                        ${params.minimumCounts}
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
