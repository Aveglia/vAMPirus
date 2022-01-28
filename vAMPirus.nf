/*
=============================================================================================================================================================
                                                                    Configuration File vAMPirus
=============================================================================================================================================================
                                                                            vAMPirus
                                                           Author: Alex J. Veglia and Ramón Rivera-Vicéns
-------------------------------------------------------------------------------------------------------------------------------------------------------------
*/

params {
// --------------------------  EDIT variables where needed  -------------------------- //

    // Project specific information

        // Project name - Name that will be used as a prefix for naming files by vAMPirus
             projtag="vAMPirusAnalysis"
        // Path to metadata spreadsheet file to be used for plot
             metadata="/PATH/TO/vampirus_meta.csv"
        // reads directory, must specify the path with "*R{1,2}*" or for reads to be properly read by Nextflow
             reads="/PATH/TO/reads/*_R{1,2}*"
        // Single-end data? Make single = true if this is the case.
            single = false
        // PATH to working directory of your choosing, will automatically be set to vAMPirus installation
             workingdir="VAMPDIR"
        // Name of directory created to store output of vAMPirus analyses (Nextflow will create this directory in the working directory)
             outdir="results"

    // Quality filter/trimming options

        // Average read quality - forward or reverse reads will be discarded if average base quality across the read is below the number set below (25 is a good start)
            avQ="25"
        // Maximum number of "N"s acceptable in the forward or reverse reads (default for fastp is 5)
            mN="5"
        // Minmum base quality to be trimmed
            trimq="15"

    // Primer Removal parameters

        // If not specifying primer sequences, forward and reverse reads will be trimmed by number of bases specified using "--GlobTrim #basesfromforward,#basesfromreverse". Use "trim" below for trimming single end reads.
            GlobTrim=""
        // If not specifying primer sequences, single reads will be trimmed by number of bases specified below.
            trim=""
        // Specific primer sequence on forward reads to be removed. Also, if using single-end mode add the primer sequence here and leave the rev="". NOTE - bbduk.sh which is used to trim the primers does not recognize Inosine (I) in the primer sequence, replace "I" with "N" in the sequence. It recognizes all other  IUPAC degenerate base codes.
            fwd=""
        // Reverse primer sequence. Leave blank if using single-end mode. -- NOTE - bbduk.sh which is used to trim the primers does not recognize Inosine (I) in the primer sequence, replace "I" with "N" in the sequence. It recognizes all other  IUPAC degenerate base codes.
            rev=""
        // Path to fasta file with primer sequences to remove (need to specify if using --multi option ). You can use this with single-end as well. -- NOTE - bbduk.sh which is used to trim the primers does not recognize Inosine (I) in the primer sequence, replace "I" with "N" in the sequence. It recognizes all other  IUPAC degenerate base codes.
            primers="/PATH/TO/PRIMERS.fasta"
        // Primer length (default 26)- If trimming primers with the --multi option or by specifying primer sequences above, change this to the length of the longer of the two primer sequences
            primerLength="26"
        // Maximum kmer length for primer removal (must be shorter than your primer length; default = 13)
            maxkmer="13"
        // Minimum kmer length for primer removal (default = 3)
            minkmer="3"
        // Minimum non-merged read length after adapter and primer removal (default = 200)
            minilen="100"

    // Merged read length filtering parameters

        // Minimum merged read length - reads with lengths greater than minLen and below the specified maximum read length will be used for counts only
            minLen="400"
        // Maximum merged read length - reads with length equal to the specified max read length will be used to generate uniques and ASVs (safe to set at expected amplicon size to start)
            maxLen="420"
        // Maximum expected error for vsearch merge command - vsearch discard sequences with more than the speciﬁed number of expected errors
            maxEE="3"
        // Maximum number of non-matching nucleotides allowed in overlap region
            diffs="10"
        // Maximum number of "N"'s in a sequence - if above the specified value, sequence will be discarded (should be similar to what is set for "mN" above in fastp parameters)
            maxn="10"
        // Minimum length of overlap for sequence merging to occur for a pair
            minoverlap="10"

    // ASV generation parameters

        // Alpha value for denoising - the higher the alpha the higher the chance of false positives in ASV generation (1 or 2)
            alpha="1"
        // Minimum size or representation in dataset for sequence to be considered in ASV generation (ex. If set to 4, any unique sequence that is not seen in the data more 3 times is removed)
            minSize="8"

    // ASV filtering parameters - You can set the filtering to run with the command --filter

        // Path to database containing sequences that if ASVs match, are then removed prior to any analyses
            filtDB=""
        // Path to database containing sequences that if ASVs match to, are kept for final ASV file to be used in susequent analyses
            keepDB=""
        // Keep any sequences without hits - for yes, set keepnohit to ="true"
            keepnohit="true"

        //Parameters for diamond command for ASV filtering

            // Set minimum percent amino acid similarity for best hit to be counted in taxonomy assignment
                filtminID="80"
            // Set minimum amino acid alignment length for best hit to be counted in taxonomy assignment
                filtminaln="30"
            // Set sensitivity parameters for DIAMOND aligner (read more here: https://github.com/bbuchfink/diamond/wiki; default = ultra-sensitive)
                filtsensitivity="ultra-sensitive"
            // Set the max e-value for best hit to be recorded
                filtevalue="0.001"

    // ASV clustering parameters

        // Percent similarity to cluster nucleotide ASV sequences (used when --ncASV is set)
            clusterNuclID="85"
        // List of percent similarities to cluster nucleotide ASV sequences - must be separated by  a comma (ex. "95,96")
            clusterNuclIDlist=""
        // Default percent similarity to cluster aminoacid sequences
            clusterAAID="97"
        // List of percent similarities to cluster aminoacid sequences - must be separated by "95,96"
            clusterAAIDlist=""
        // Minimum length of amino acid translation to be considered during protein clustered ASV (pcASV) generation. Recommended to put this at the expected amino acid sequence length based on your maximum read length (e.g. if maxLen="420", then minAA should be 420/3 so 140)
            minAA="140"

    // Minimum Entropy Decomposition (MED) parameters for clustering (https://merenlab.org/2012/05/11/oligotyping-pipeline-explained/)

        // If you plan to do MED on ASVs using the option "--asvMED" you can set here the number of entopy peak positions or a comma seperated list of biologically meaningful positons (e.g. 35,122,21) for oligotyping to take into consideration. If you want to use a single specific position, make "asvSingle="true"".
            asvC=""
            asvSingle="false"
        // If you plan to do MED on ASVs using the option "--aminoMED" you can set here the number of positions for oligotyping to take into consideration. If you want to use a single specific position, make "aminoSingle="true"".
            aminoC=""
            aminoSingle="false"

    // Phylogeny-based ASV/AminoType clustering parameters using the program TreeCluster (https://github.com/niemasd/TreeCluster)

      // Add the "--asvTClust" option to the launch command to run phylogeny-based clustering of ASVs ; Add "--aminoTClust" to launch command for phylogeny-based clustering on AminoTypes
      // NOTE: you can't use "--skipPhylogeny" when doing this form of sequence clustering

        // TreeCluster command options for ASV clustering (--asvTClust) -- (Example: "-option1 A -option2 B -option3 C -option4 D") - See TreeCluster paper and github page to determine the best options (a good start is what is below)
            asvTCopp="-t 0.045 -m max_clade"

        // TreeCluster command options for AminoType clustering (--aminoTClust) -- (Example: "-option1 A -option2 B -option3 C -option4 D") - See TreeCluster paper and github page to determine the best options
            aminoTCopp="-t 0.045 -m max_clade"

    // Counts table generation parameters

        // Percent similarity to use for ASV/cASV counts table generation with vsearch
            asvcountID="97"
        // Parameters for protein counts table generation
            // Minimum Bitscore for counts
                ProtCountsBit="50"
            // Minimum aminoacid sequence similarity for hit to count
                ProtCountID="85"
            // Minimum alignment length for hit to count
                ProtCountsLength="50"

    // Taxonomy inference parameters

        //Parameters for diamond command
            // Set which measurement to use for a minimum threshold in taxonomy inference - must be either "evalue" or "bitscore"
                measurement="bitscore"
            // Set maximum e-value for hits to be counted
                evalue="0.001"
            // Set minimum bitscore for best hit in taxonomy assignment (default = 30)
                bitscore="30"
            // Set minimum percent amino acid similarity for best hit to be counted in taxonomy assignment
                minID="40"
            // Set minimum amino acid alignment length for best hit to be counted in taxonomy assignment
                minaln="30"
            // Set sensitivity parameters for DIAMOND aligner (read more here: https://github.com/bbuchfink/diamond/wiki; default = ultra-sensitive)
                sensitivity="ultra-sensitive"

        // Database information
            // Specify name of database to use for analysis
                dbname="DATABASENAME"
            // Path to Directory where database is being stored - vAMPirus will look here to make sure the database with the name provided above is present and built
                dbdir="DATABASEDIR"
            // Set database type (NCBI or RVDB). Lets vAMPirus know which sequence header format is being used and must be set to NCBI when using RefSeq or Non-Redundant databases. -> dbtype="NCBI" to toggle use of RefSeq header format; set to "RVDB" to signal the use of Reverence Viral DataBase (RVDB) headers (see manual)
                dbtype="TYPE"

        // Classification settings - if planning on inferring LCA from RVDB annotation files OR using NCBI taxonomy files, confirm options below are accurate.
            // Path to directory RVDB hmm annotation .txt file - see manual for information on this. Leave as is if not planning on using RVDB LCA.
                dbanno="DATABASEANNOT"
            // Set lca="T" if you would like to add "Least Common Ancestor" classifications to taxonomy results using information provided by RVDB annotation files (works when using NCBI or RVDB databases) - example: "ASV1, Viruses::Duplodnaviria::Heunggongvirae::Peploviricota::Herviviricetes::Herpesvirales::Herpesviridae::Gammaherpesvirinae::Macavirus"
                lca="LCA"
            // DIAMOND taxonomy inference using NCBI taxmap files (can be downloaded using the startup script using the option -t); set to "true" for this to run (ONLY WORKS WITH dbtype="NCBI")
                ncbitax="false"

    // Phylogeny analysis parameters

        // Color nodes on phylogenetic tree in Analyze report with MED Group information (nodeCol="MED"), taxonomy (nodeCol=TAX) hit, or TreeCluster Group Information (nodeCol=TC). If you would like nodes colored by sequence ID, leave nodeCol="" below.
            nodeCol="empty" //Do not leave this parameter without a value, will cause error.

        // Customs options for IQ-TREE (Example: "-option1 A -option2 B -option3 C -option4 D")
            iqCustomnt=""
            iqCustomaa=""

        // These options below you can set at the command like, for example, to set to use model from ModelTest-NG with parametric bootstrapping --ModelTnt --ModelTaa --parametric

            // Signal for IQ-TREE to use model determined by ModelTest-NG (Default is IQ-TREE will do automatic model testing with ModelFinder Plus)
                ModelTnt=false
                ModelTaa=false
            // Set to have non-parametric bootstrapping IQ-TREE to perform
                parametric=false
                nonparametric=false
            // Number of bootstraps (recommended 1000 for parametric and 100 for non-parametric)
                boots="1000"

    // Stats options
        // Tell vAMPirus to perform statistical analyses by setting "stats = true" below or in the launch command by adding "--stats" to it
             stats = false
        // Minimum number of hit counts for a sample to have to be included in the downstream statistical analyses and report generation
             minimumCounts="1000"
        // Maximum number of iteration performed by metaMDS
             trymax="900"

    // Conda env PATH (added automatically by startup script)
        condaDir="CONDADIR"

/*
// ----------------------------------------------------------------------     STOP     ---------------------------------------------------------------------- //
// -------------------------------------------------------- Do not modify variables below this line. -------------------------------------------------------- //
// ----------------------------------------------------------- Proceed to modify processes at end ----------------------------------------------------------- //
// -----------------------------------------------------------------        If needed       ----------------------------------------------------------------- //
*/

// Path to vAMPirus installation directory, will be filled automatically when startup script is run, otherwise, edit below
    vampdir="VAMPDIR"
// Pipeline options
    help = false
    fullHelp = false
// Manadotory arguments
    Analyze=false
    DataCheck=false
// Non-Mandatory options
    // Cluster nucleotide sequences (ncASVs)
        ncASV = false
    // Cluster by aminoacid translations and generate protein-based OTUs (pcASVs)
        pcASV = false
    // Generate virus types with MED of ASV sequences
        asvMED = false
    // Generate virus types with MED of ASV sequences
        aminoMED = false
    // Filter ASVs
        filter = false
    // Phylogeny-based ASV clustering
        asvTClust = false
    // Phylogeny-based AminoType clustering
        aminoTClust = false
// Skip options
    // Skip all Read Processing steps
        skipReadProcessing = false
    // Skip quality control processes only
        skipFastQC = false
    // Skip adapter removal process only
        skipAdapterRemoval = false
    // Skip primer removal process only
        skipPrimerRemoval = false
    // Skip AminoTyping
        skipAminoTyping = false
    // Skip Taxonomy
        skipTaxonomy = false
    // Skip phylogeny
        skipPhylogeny = false
    // Skip EMBOSS analyses
        skipEMBOSS = false
    // Skip Reports
        skipReport = false
    // Skip Merging steps -> will also skip all read Processing
        skipMerging = false

// Data check parameters
    datacheckntIDlist=".55,.65,.75,.80,.81,.82,.83,.84,.85,.86,.87,.88,.89,.90,.91,.92,.93,.94,.95,.96,.97,.98,.99"
    datacheckaaIDlist=".55,.65,.75,.80,.81,.82,.83,.84,.85,.86,.87,.88,.89,.90,.91,.92,.93,.94,.95,.96,.97,.98,.99,1.0"
// If not specifying primer sequences OR --GlobTrim, forward reads will be trimmed by number of bases specified below
        defaultFwdTrim="20"
// If not specifying primer sequences OR --GlobTrim, reverse reads will be trimmed by number of bases specified below
        defaultRevTrim="26"
// Option for multi-barcoding approach
        multi = false
// Directory for pipeline info
        tracedir="PipelinePerformance"
// Path to logo
        logo="${params.vampdir}/exampledata/conf"
// readTest
        readsTest = false
// These options will chnage how the profiles work.
    // Run with conda installed by the precheck
        myConda = false
        condaActivate = false

    // vAMPirus container with all programs
        oneContainer = false

    // Cache directory for conda and singularity files. Leave in blank if not sure
        envCacheDir = ""

    // Singularity
    // Use singularity image created after pulling from docker and not from Galaxy depot (singularity image ready to use).
        singularity_pull_docker_container = false
        sing = false
}
/*
// ------------------------- Process variables below  ------------------------- //
    Proceed to modify processes if needed. Choose the scheduler and options:
            Executor = SLURM, PBS, etc.
            Cluster Options = Partition, Nodes, Priority, Email, etc.
    If running locally leave the comments (the ""\\") on "executor" and "clusterOptions".
    For more info see the README and/or Nextflow documentation.
*/

process {
    withLabel: low_cpus {
        cpus='1'
        memory='2 GB'
        //executor='slurm'
        //clusterOptions='--cluster=cm2 --partition=cm2_tiny --qos=cm2_tiny --nodes=1'
    }
    withLabel: norm_cpus {
        cpus='2'
        memory='2 GB'
        //executor='slurm'
        //clusterOptions='--cluster=cm2 --partition=cm2_tiny --qos=cm2_tiny --nodes=1'
    }
    withLabel: high_cpus {
        cpus='2'
        memory='2 GB'
        //executor='slurm'
        //clusterOptions='--cluster=cm2 --partition=cm2_tiny --qos=cm2_tiny --nodes=1'
    }
    errorStrategy='finish'
}

// env variables (only for nextflow)
env.tools="${params.vampdir}/bin/"

process.shell = ['/bin/bash', '-euo', 'pipefail']

timeline {
  enabled = true
  file = "${params.workingdir}/${params.outdir}/${params.tracedir}/vampirus_timeline.html"
}
report {
  enabled = true
  file = "${params.workingdir}/${params.outdir}/${params.tracedir}/vampirus_report.html"
}
trace {
  enabled = true
  file = "${params.workingdir}/${params.outdir}/${params.tracedir}/vampirus_trace.txt"
}
dag {
  enabled = true
  file = "${params.workingdir}/${params.outdir}/${params.tracedir}/vampirus_dag.html"
}

// Get PATH for cache environments
params.localCacheDir = (params.envCacheDir ? "${params.envCacheDir}" : "${launchDir}")

profiles {
    conda {
        params.condaActivate = true
        // cache for condaEnv created individually
        conda.cacheDir = "${params.localCacheDir}/condaEnv/"
        process.conda = "${params.condaDir}"
    }
    docker {
        docker.enabled = true
        docker.runOptions = '-u \$(id -u):\$(id -g)'
        process.container="aveglia/vampirus:v2.0.0"
    }
    singularity {
        singularity.enabled = true
        singularity.autoMounts = true
        // cache for images from docker pull
        singularity.cacheDir="${params.localCacheDir}/singularityCache/"
        process.container="aveglia/vampirus:v2.0.0"
        params.sing = true
    }
    podman {
        podman.enabled = true
        process.container="aveglia/vampirus:v2.0.0"
    }
    vampContainer {
        process {
            params.oneContainer = true
            params.VPcontainer="aveglia/vampirus:v2.0.0"
        }
    }
    test {
        includeConfig 'example_data/conf/test.config'
    }
}
manifest {
    name = 'vAMPirus'
    author = 'Alex J. Veglia,Ramón Rivera-Vicéns'
    description = 'Automated virus amplicon sequencing analysis program'
    mainScript = 'vAMPirus.nf'
    nextflowVersion = '>=21.04.1'
    version = '2.0.2'
}
