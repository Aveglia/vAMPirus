#!/usr/bin/env nextflow

/*
========================================================================================
                                vAMPirusv0.1.0
========================================================================================
                       Virus Amplicon Sequencing Analysis Pipeline
                       Author: Alex J. Veglia
                       Version: 0.1.0 (dev) LAST EDIT: 5/30/20
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info """
    vAMPirus: A program for non-rRNA amplicon seq data analysis currenlty tailored towards viruses

    If you have any comments or questions, feel free to contact Alex Veglia at ajv5@rice.edu

    General exicution:

    Move to directory containing raw or cleaned read files and exicute script

    General exicution:

    nextflow run vAMPirusv0.1.0.sh --all (other_options_here) -[hrabBefiIlLnmMpztx] 2>&1 | tee output_vAMP.out

    ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

    Command line options:

        [ -h ]                       	Print help information

        [ -r ]                       	Runs a loop to cut read library file names, does not run analysis when set; BEFORE RUNNING: NEED TO SET COLUMN NUMBER TO CUT TO
                            					FROM THE LEFT - see lines 179.

        [ -a ]	                    	Detect and remove adapter sequences with fastp ( should be coupled with -b or -B options for primer removal).

        [ -b ]                       	Primer detection and removal from reads using BBduk.sh post-adapter contamination removal.

        [ -B <f,r|n>]                	Remove primers by global trimming with BBduk.sh: chopping off specified amount of basepairs: f (forward) and r (reverse), or n (both).

        [ -c ]	                 	    With this option set, ASVs will be culled by most represented sequence length.

        [ -d <1|2> ]          		    Set the db header format: 1 (default) - Reference Viral Database ; 2 - RefSeq format

        [ -e <exmple@email.com> ]		Set email for notification to be sent to after analysis is completed.

        [ -f <1|2> ]          		    Sets alpha parameter for denoising, default is a more strict value of 1. Higher the alpha the higher the chance of false positives.

        [ -g ]                			Set this option to have nucleotide and protein RAxML commands run with the GTR substitution model instead of best model spit out by modeltest

        [ -i <percent id> ]			    Sets OTU clustering percentage used by vsearch, default is $defid; can't be set with -I.

        [ -I ]            				Set this option to run the analysis with several different cluster percentages. If this is set, a single copy of the list file must
        					            be in the wokring directory and end with *id.list.

        [ -j <1|2|3> ]                  This option signals that the user would like to apply a custom RAxML-ng command for:
          								                                1 - Nucleotide RAxML run only; default or for protein tree if set
        								                                2 - Protein RAxML run only; default for nucleotide tree
        								                                3 - Custom RAxML command supplied for both nucleotide and protein tree
        								                                    **NOTE: with this option set, you must have *.raxcom file set up and present in your working directory

        [ -l <value> ]             		Use this option to set the minimum read length to be included in the analysis. Default is 400 bp.

        [ -L <value> ]			        Use this option to set the maximum read length to be included in the analysis. Default is 440 bp.

        [ -m ]                			Determine nucleotide model of substitution by running ModelTest-NG with non-labeled OTU alignments.

        [ -M ]                          Runs nucleotide substitution model analysis and creates a phylogenetic tree with RAxML-ng.

        [ -n <Exmp_proj> ]			    Use this option to set the prefix for all output files produced.

        [ -p <1|2|3|4|5> ]			    Set this option with [1-5] for different combination of protein analysis:
           								                                1 - Generate protein translation of all OTU fastas generated in analysis
        								                                2 - Generate protein translation of all OTU fastas and run model test
        								                                3 - Generate protein translation of all OTU fastas, run model test, and run RAxML
        								                                4 - Generate protein translation of all OTU fastas, run model test, run RAxML and diamond blastp
        								                                5 - Generate protein translation of all OTU fastas and run diamond blastp

        [ -t <value> ]			         Use this option to set the number of threads available for use in the analysis; default is 1.

        [ -x ]				             Set this option to generate percent ID and pairwise distance matrices with all OTU files generated

        [ -y ]				             This option combined with adapter/primer removal options -[ab|B] will limit the script to only read cleaning.

        [ -z ] 				             Set this option to include unclustered ASV file in counts table generation and all extra specified anlyses (e.g. model testing et al.)

        |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    """.stripIndent()
}
def fullHelpMessage() {
    log.info """
    ==============================================================================================================================================================================================
                                                        THIS IS A LONGER HELP WITH USAGE EXAMPLES vAMPirus
    ==============================================================================================================================================================================================

        Steps:
            1- Run the `precheck_TransPi.sh` to install tools, set up the databases and directories for
            2- Run vAMPirusv0.1.0.sh

            Usage:

                nextflow run vAMPirusv0.1.0.sh --all (other_options_here)

        Mandatory arguments (choose one):

                --only-def                     Run read processing steps and generate ASVs, OTU's, OTU counts tables with vsearch, and pair-wise distance and %ID matrices

                --def-ntmt                     Run default analysis plus nucleotide alignment with MAFFT and model testing with ModelTest-NG

                --def-ntmtrx                   Run default analysis plus nucleotide alignment with MAFFT, model testing with ModelTest-NG, and RAxML phylogeny

                --def-ptmt                     Run default analysis plus protein translation and protein model testing

                --def-ptmtrx                   Run default analysis plus protein translation, protein model testing, and RAxML phylogeny

                --def-ptmtrxpc                 Run default analysis plus protein translation, protein model testing, RAxML phylogeny, and generate protein counts tables

                --def-ntptmt                   Run default analysis plus nucleotide and protein model testing (creates translated sequences as well)

                --def-mtrx                     Runs all model testing for nucleotide and protein translation

                --noclust-ntptmt

                --noclust-todos

                --prot-cnts-only               Runs only the protein counts script, need to define directory where protein fastas are

                --con-todos                    Run absolutely everything

                --read-processing              Run only read processing steps

        Other options:

                [ -h ]                         Print help information

                --alpha                        Sets alpha parameter for denoising, default is a more strict value of 1. Higher the alpha the higher the chance of false positives.

                --GlobTrim                     Remove primers by global trimming with BBduk.sh: chopping off specified amount of basepairs: f (forward) and r (reverse), or n (both).

                --refseq <F|T>          	   Set the database header format: F (default) - Reference Viral Database ; T - RefSeq format

                --idlist <.id1,.id2,.id3,..>

                --clustid <.[0-99]>		       Sets OTU clustering percentage used by vsearch, default is $defid; can't be set with -I.


                [ -j <1|2|3> ]                 This option signals that the user would like to apply a custom RAxML-ng command for:
                                                              1 - Nucleotide RAxML run only; default or for protein tree if set
                                                              2 - Protein RAxML run only; default for nucleotide tree
                                                              3 - Custom RAxML command supplied for both nucleotide and protein tree
                                                               **NOTE: with this option set, you must have *.raxcom file set up and present in your working directory

                [ -p ]                         Signals to do prot counts with files from defined directory in config ${params.protdir}

                --min_rd_len <value>           Use this option to set the minimum read length to be included in the analysis. Default is 400 bp.

                --max_rd_len <value>           Use this option to set the maximum read length to be included in the analysis. Default is 440 bp.

                [ -t <value> ]			       Use this option to set the number of threads available for use in the analysis; default is 1.

                [ -z ] 				           Set this option to include unclustered ASV file in counts table generation and all extra specified anlyses (e.g. model testing et al.)

        |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

        #################################################################################################

                                  Various examples on how to deploy vAMPirus

        #################################################################################################

        I. Steps for running on a local cluster and local installation of TransPi

            1- Run the `precheck_TransPi.sh` to install tools, set up the databases and directories for TransPi
            2- Run TransPi

            Usage:

                nextflow run TransPi.nf --all -with-conda ~/anaconda3/envs/TransPi

        #################################################################################################

        II. Steps for running on a local cluster and conda installation by nextflow

            1- Run the `precheck_TransPi.sh` to install tools, set up the databases and directories for TransPi
            2- Run TransPi

            Usage:

                nextflow run TransPi.nf --all -profile conda

            NOTE:
                Not recommended, not all programs are installed by conda. Use if other dependencies are manually installed

        #################################################################################################

        III. Steps for running on docker (in development)

            1- Run the `container_precheck_TransPi.sh` to install tools, set up the databases and directories for TransPi
            2- Run TransPi

            Usage:

                nextflow run TransPi.nf --all -profile docker

            NOTE:
                All necesary tools for running TransPi are pre-installed in the container
                `container_precheck_TransPi.sh` will install only the database used by TransPi

        #################################################################################################

        IV. Steps for running on singualarity (in development)

            1- Run the `container_precheck_TransPi.sh` to install tools, set up the databases and directories for TransPi
            2- Run TransPi

            Usage:

                nextflow run TransPi.nf --all -profile singularity

            NOTE:
                All necesary tools for running TransPi are pre-installed in the container
                `container_precheck_TransPi.sh` will install only the database used by TransPi

        #################################################################################################

        V. Steps for running with multiple profiles (in development - depending on profile selected)

            1- Run the `precheck_TransPi.sh` to install tools, set up the databases and directories for TransPi
            2- Run TransPi

            Usage:

                nextflow run TransPi.nf --all -profile conda,test

            NOTE:
                This will run TransPi using a test dataset and a conda environment created by Nextflow
                To run with containers first run the `container_precheck_TransPi.sh` and use -profile docker,test

        #################################################################################################

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

// This will be printed to the user in each run. here thy can check if the values the selected are fine
log.info """\

        ======================================================
        vAMPirus - Virus Amplicon Sequencing Analysis Pipeline
        ======================================================
        Project name:                ${params.projtag}
        Email:                       ${params.email}
        Working directory:           ${params.mypwd}
        Minimum read length:	     ${params.minLen}
        Maximum read length:         ${params.maxLen}
        OTU cluster percentages:     ${params.clustid}
        Alpha value:                 ${params.alpha}
        Database directory:          ${params.dbdir}
        Database name:               ${params.dbname}
        Database header format:      ${params.dbhead}

        """.stripIndent()

// read files from here. This will be changed!
if (params.readsTest) {
    println("\n\tRunning vAMPirus with TEST dataset\n")
    Channel
        .from(params.readsTest)
        .map{ row -> [ row[0], [ file(row[1][0],checkIfExists: true),file(row[2][0],checkIfExists: true) ] ] }
        .ifEmpty{ exit 1, "params.readsTest was empty - no input files supplied" }
        .into{ reads_ch; reads_qc_ch }
} else {
    println("\n\tRunning vAMPirus with your dataset\n")
    Channel
        .fromFilePairs("${params.reads}", checkIfExists: true)
        .into{ reads_ch; reads_qc_ch }
}

process build_diamond_db {
    script:
        """
        cd ${params.mypwd}
        echo -e "-- Checking if Diamond database folder is present --\\n"
        if [ ! -d DBs/diamonddb_custom/ ];then
            echo -e "-- Folder is not present, creating one and the Diamond database --\\n"
            mkdir -p DBs/diamonddb_custom/
            cd DBs/diamonddb_custom
            cp ${params.dbdir} .
            diamond makedb --in ${params.dbname} -d ${params.dbname}
            export virdb=`pwd`/${params.dbname}
            cd ../
        elif [ -d DBs/diamonddb_custom/ ];then
            echo -e "-- Folder is present. Checking if Diamond database is built --\\n"
            cd DBs/diamonddb_custom
            if [ ! -e ${params.dbname}.dmnd ];then
                echo -e "-- Diamond database not present, creating one --\\n"
                cp ${params.dbdir} .
                diamond makedb --in ${params.dbname} -d ${params.dbname}
                export virdb=`pwd`/${params.dbname}
            elif [ -e ${params.dbname}.dmnd  ];then
                echo -e "-- Diamond database already created --\\n"
                export virdb=`pwd`/${params.dbname}
            fi
            cd ../
        fi
        """
}

if (params.conTodos) {

    println("\n\tRunning vAMPirus \n")

    if (!params.skipFastQC) {

        process fastqc_VAP {

            label 'low_cpus'

            tag "${sample_id}"

            publishDir "${params.mypwd}/${params.outdir}/Read_Processing/fastqc", mode: "copy", overwrite: true

            input:
                tuple sample_id, file(reads) from reads_qc_ch

            output:
                tuple sample_id, file("*_fastqc.{zip,html}") into fastqc_results_OAS

            script:
                """
                fastqc --quiet --threads $task.cpus $reads
                """
        }
    }

    if (!params.skipFastP) {

        process adapterremoval_fastp_VAP {

            label 'low_cpus'

            tag "${sample_id}"

            publishDir "${params.mypwd}/${params.outdir}/Read_Processing/adapter_removal", mode: "copy", overwrite: true, pattern: "*.fastp.{json,html}"

            input:
                tuple sample_id, file(reads) from reads_ch

            output:
                tuple sample_id, file("*.fastp.{json,html}") into fastp_results
                tuple sample_id, file("*.filter.fq") into reads_fastp_ch

            script:
                """
                echo ${sample_id}

                fastp -i ${reads[0]} -I ${reads[1]} -o left-${sample_id}.filter.fq -O right-${sample_id}.filter.fq --detect_adapter_for_pe \
                --average_qual 25 --overrepresentation_analysis --html ${sample_id}.fastp.html --json ${sample_id}.fastp.json --thread ${task.cpus} \
                --report_title ${sample_id}
                """
            }
    } else {
        reads_ch
            .set{ reads_fastp_ch }
        fastp_results = Channel.empty()
    }

    if (!params.skipBBduk) {

        process primerremoval_bbduk_VAP {

            label 'low_cpus'

            tag "${sample_id}"

            publishDir "${params.mypwd}/${params.outdir}/Read_Processing/primer_removal", mode: "copy", overwrite: true

            input:
                tuple sample_id, file(reads) from reads_fastp_ch

            output:
                tuple sample_id, file("*bbduk*.fastq.gz") into ( reads_bbduk_ch, readsforqc2 )

            script:
                // check if we need to check this outside processes
                if ( params.fwd == "" && params.rev == "" ) {
                    """
                    bbduk.sh in1=${reads[0]} out=${sample_id}_bb_R1.fastq.gz ftl=${params.FTRIM} t=${task.cpus}
                    bbduk.sh in=${reads[1]} out=${sample_id}_bb_R2.fastq.gz ftl=${params.RTRIM} t=${task.cpus}
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
                } else {
                    """
                    bbduk.sh in=${reads[0]} in2=${reads[1]} out=${sample_id}_bbduck_R1.fastq.gz out2=${sample_id}_bbduck_R2.fastq.gz literal=${params.fwd},${params.rev} copyundefined=t t=${task.cpus}
                    """
                }
    	  }
    } else {
        reads_fastp_ch
            .set{ reads_bbduk_ch }
    }

    if (!params.skipFastQC) {

        process fastqc2_VAP {

            label 'low_cpus'

            tag "${sample_id}"

            publishDir "${params.mypwd}/${params.outdir}/Read_Processing/fastqc/fastqc2", mode: "copy", overwrite: true

            input:
                tuple sample_id, file(reads) from readsforqc2

            output:
                tuple sample_id, file("*_fastqc.{zip,html}") into fastqc2_results_OAS

            script:
                """
                fastqc --quiet --threads $task.cpus $reads
                """
        }
    }

    process mergereads_vsearch_VAP {

        label 'norm_cpus'

        tag "${sample_id}"

        publishDir "${params.mypwd}/${params.outdir}/Read_Processing/read_merging", mode: "copy", overwrite: true

        input:
            tuple sample_id, file(reads) from reads_bbduk_ch

        output:
            tuple sample_id, file("*_mergered.fq") into reads_vsearch1_ch

        script:
            """
            vsearch --fastq_mergepairs ${reads[0]} --reverse ${reads[1]} --threads ${task.cpus} --fastqout ${sample_id}_mergered.fq --fastq_maxee 1 --relabel ${sample_id}.
            """

    }

    process collect_samples {

        label 'low_cpus'

        publishDir "${params.mypwd}/${params.outdir}/Read_Processing/readmerging", mode: "copy", overwrite: true

        input:
            file(reads) from reads_vsearch1_ch
                .collect()

        output:
            file("all_merged.fq") into collect_samples_ch
            file("*sample_ids.list") into samplelist

        script:
            """
            cat ${reads} >>all_merged.fq
            echo ${reads} | awk -F "_mergered." '{print \$1}' >>${params.projtag}_sample_ids.list
	        """
    }

    process readlength_trimming_VAP {

        label 'norm_cpus'

        publishDir "${params.mypwd}/${params.outdir}/Read_Processing/read_merging/clean", mode: "copy", overwrite: true

        input:
            file(reads) from collect_samples_ch

        output:
            file("all_merged_clean.fasta") into ( mergedforcounts, mergeforprotcounts, mergeforpOTUcounts, mergeforpOTUaacounts )
            file("all_merged_clean_filter.fasta") into ( reads_vsearch2_ch )

        script:
            """
            fastp -i ${reads} -o all_merged_cln.fq -b ${params.maxLen} -l ${params.minLen} --thread ${task.cpus} -n 1

            # reformat.sh from BBtools
            reformat.sh in=all_merged_cln.fq out=all_merged_clean.fasta t=${task.cpus}

            bbduk.sh in=all_merged_clean.fasta out=all_merged_clean_filter.fasta minlength=${params.maxLen} maxlength=${params.maxLen} t=${task.cpus}
	        """
    }

    process uniques_vsearch_VAP {

        label 'norm_cpus'

        publishDir "${params.mypwd}/${params.outdir}/Read_Processing/uniques", mode: "copy", overwrite: true

        input:
            file(reads) from reads_vsearch2_ch

        output:
            file("all_unique_seq.fasta") into reads_vsearch3_ch

        script:
            """
            vsearch --derep_fulllength ${reads} --sizeout --relabel_keep --output all_unique_seq.fasta
            """
    }

    process asvs_vsearch_VAP {

        label 'norm_cpus'

        publishDir "${params.mypwd}/${params.outdir}/Read_Processing/uniques/asvs/wchimeras", mode: "copy", overwrite: true

        input:
            file(reads) from reads_vsearch3_ch

        output:
            file("asv_chim.fasta") into reads_vsearch4_ch

        script:
            """
            vsearch --cluster_unoise ${reads} --unoise_alpha ${params.alpha} --relabel ASV --centroids asv_chim.fasta
            """
    }

    process chimeradetection_vsearch_VAP {

        label 'norm_cpus'

        publishDir "${params.mypwd}/${params.outdir}/Clustering/asvs", mode: "copy", overwrite: true

        input:
            file(reads) from reads_vsearch4_ch

        output:
            file("*ASVs_all.fasta") into ( reads_vsearch5_ch, nucl2aa, asvsforAminotyping, asvfastaforcounts )

        script:
            """
	        vsearch --uchime3_denovo asv_chim.fasta --relabel ASV --nonchimeras ${params.projtag}_ASVs_all.fasta
            """
    }

    if (params.clusterntASV) {

        process nucleotide_cluster_VAP {

            label 'norm_cpus'

            publishDir "${params.mypwd}/${params.outdir}/Clustering/nucleotide_OTU", mode: "copy", overwrite: true, pattern: '*ASVs_all.fasta'

            input:
                file(fasta) from reads_vsearch5_ch

            output:
                tuple file("*_otu*.fasta"), file("*ASV_all.fasta") into ( reads_diamond_ch, fileforcounts )
                file("*_otu*.fasta") into ( fileforMafft, fileforClustal )

            script:
            if (params.idlist) {
                """
                for id in `echo ${params.idlist} | tr "," "\n"`;done
                    vsearch --cluster_fast ${fasta} --centroids ${params.projtag}_otu\${id}.fasta --threads ${task.cpus} --relabel OTU --id \${id}
                done
                """
            } else if (params.sinid) {
                """
                id=${params.sinid}
                vsearch --cluster_fast ${fasta} --centroids ${params.projtag}_otu\${id}.fasta --threads ${task.cpus} --relabel OTU --id \${id}
                """
            }
        }
    } else {
        reads_vsearch5_ch
	       .into{ reads_diamond_ch; fileforcounts; fileforMafft; fileforClustal }
    }

    process diamond_VAP {

        label 'norm_cpus'

        publishDir "${params.mypwd}/${params.outdir}/Taxonomy/", mode: "copy", overwrite: true, pattern: '*.{fasta,csv,tsv}'
        publishDir "${params.mypwd}/${params.outdir}/Taxonomy/diamond_outputs", mode: "copy", overwrite: true, pattern: '*dmd.out'

        input:
            file(reads) from reads_diamond_ch

        output:
            file("*.fasta") into ( tax_labeled_fasta )
            tuple file("*_sumPHYtab.csv"), file("*_summarytable.tsv"), file("*dmd.out") into summary_diamond

        script:
            """
            virdb=${params.mypwd}/DBs/diamonddb_custom/${params.dbname}
            name=\$(ls ${reads} | awk -F ".fasta" '{print \$1}')
            diamond blastx -q ${reads} -d \${virdb} -p ${task.cpus} --min-score 50 --more-sensitive -o "\$name"_dmd.out -f 6 qseqid qlen sseqid qstart qend qseq sseq length qframe evalue bitscore pident btop --max-target-seqs 1
            echo "Preparing lists to generate summary .csv's"
            echo "[Best hit accession number]" > access.list
            echo "[OTU sequence length]" > length.list
            echo "[e-value]" > evalue.list
            echo "[Bitscore]" > bit.list
            echo "[Percent ID (aa)]" > pid.list
            echo "[OTU#]" > otu.list
            echo "[Virus ID]" > "\$name"_virus.list
            echo "[Gene]" > "\$name"_genes.list
            grep ">" ${reads} | awk -F ">" '{print \$2}' > seqids.lst
            echo "extracting genes and names"
            touch new_"\$name"_asvnames.txt
            j=1
            for s in \$(cat seqids.lst);do
                echo "Checking for \$s hit in diamond output"
                    if [[ ${params.refseq} == "T" ]];then
                        echo "RefSeq headers specified"
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
                            echo ">OTU\${j}_"\$virus"_"\$gene"" >> new_"\$name"_asvnames.txt
                            j=\$((\$j+1))
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
                            echo ">OTU\${j}_"\$virus"_"\$gene"" >> new_"\$name"_asvnames.txt
                            j=\$((\$j+1))
                            echo "\$s done."
                    fi
                    else
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
                            echo ">OTU\${j}_"\$virus"_"\$gene"" >> new_"\$name"_asvnames.txt
                            j=\$((\$j+1))
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
                            echo ">OTU\${j}_"\$virus"_"\$gene"" >> new_"\$name"_asvnames.txt
                            j=\$((\$j+1))
                            echo "\$s done."
                        fi
                    fi
            echo "Done with \$s"
            done
            echo "Now editing "\$name" fasta headers"
            ###### rename_seq.py
            rename_seq.py ${reads} new_"\$name"_asvnames.txt "\$name"_fin.fasta
            awk 'BEGIN {RS=">";FS="\\n";OFS=""} NR>1 {print ">"\$1; \$1=""; print}' "\$name"_fin.fasta > "\$name"_tmpssasv.fasta
            echo "[Sequence header]" > newnames.list
            cat new_"\$name"_asvnames.txt >> newnames.list
            touch sequence.list
            echo "     " > sequence.list
            grep -E '^[ACGT]+\$' "\$name"_tmpssasv.fasta >> sequence.list
            rm "\$name"_tmpssasv.fasta
            paste -d "," sequence.list "\$name"_virus.list "\$name"_genes.list otu.list newnames.list length.list bit.list evalue.list pid.list access.list >> "\$name"_sumPHYtab.csv
            paste -d"\t" otu.list access.list "\$name"_virus.list "\$name"_genes.list sequence.list length.list bit.list evalue.list pid.list >> "\$name"_summarytable.tsv
            rm evalue.list ; rm sequence.list ; rm bit.list ; rm pid.list ; rm length.list ;
            rm new_"\$name"_asvnames.txt
            rm "\$name"_virus.list
            rm "\$name"_genes.list
            rm newnames.list
            rm access.list
            echo "Done Assigning Taxonomy To : ${reads} "
            """
    }

    process vsearch_global_VAP {

        label 'norm_cpus'

        publishDir "${params.mypwd}/${params.outdir}/Counts", mode: "copy", overwrite: true, pattern: '*.{biome,csv}'

        input:
            file(reads) from fileforcounts
            file(merged) from mergedforcounts

        output:
            tuple file("*_counts.csv"), file("*_counts.biome") into counts_vsearch

        script:
            """
             if [ `echo ${reads} | grep -c "otu"` -eq 1 ];then
                ident=\$( echo ${reads} | awk -F "otu" '{print \$2}' | awk -F ".fasta" '{print \$1}')
                name=\$( echo ${reads} | awk -F ".fasta" '{print \$1}')
                vsearch --usearch_global all_merged_clean.fasta --db ${reads} --id \$ident --threads ${task.cpus} --otutabout \${name}_counts.txt --biomout \${name}_counts.biome
                cat \${name}_counts.txt | tr "\t" "," >\${name}_counts.csv
            else
                name=\$( echo ${reads} | awk -F ".fasta" '{print \$1}')
                vsearch --usearch_global all_merged_clean.fasta --db ${reads} --id ${params.idForasvCounts} --threads ${task.cpus} --otutabout "\$name"_counts.txt --biomout "\$name"_counts.biome
                cat \${name}_counts.txt | tr "\t" "," >\${name}_counts.csv
            fi
            """
    }

    process clustalo_VAP {

        label 'norm_cpus'

        publishDir "${params.mypwd}/${params.outdir}/Clustering/matrix", mode: "copy", overwrite: true

        input:
            file(reads) from fileforClustal

        output:
            file("*.matrix") into clustmatrices

        script:
            """
            if [ `echo ${reads} | grep -c "fin"` -eq 1 ];then
                dent=\$( echo ${reads} | awk -F "otu" '{print \$2}' | awk -F "_fin" '{print \$1}')
                name=\$( echo ${reads} | awk -F "_fin" '{print \$1}')
                clustalo -i ${reads} --distmat-out=\${name}_dist.matrix --full --force --threads=${task.cpus}
                clustalo -i ${reads} --distmat-out=\${name}_id.matrix --percent-id --full --force --threads=${task.cpus}
            else
                ident=\$( echo ${reads} | awk -F "otu" '{print \$2}' | awk -F ".fasta" '{print \$1}')
                name=\$( echo ${reads} | awk -F ".fasta" '{print \$1}')
                clustalo -i ${reads} --distmat-out=\${name}_dist.matrix --full --force --threads=${task.cpus}
                clustalo -i ${reads} --distmat-out=\${name}_id.matrix --percent-id --full --force --threads=${task.cpus}
            fi
            """
    }

    process mafft_VAP {

        label 'norm_cpus'

        publishDir "${params.mypwd}/${params.outdir}/Analyses/alignment", mode: "copy", overwrite: true

        input:
            file(reads) from fileforMafft

        output:
            file("*_aln.fasta") into ( ntmodeltest, ntrax_align )
            file("\${name}_aln.html") into align_html

        script:
            """
            name=\$( echo ${reads} | awk -F ".fasta" '{print \$1}' )
            mafft --maxiterate 5000 --auto ${reads} >\${name}_ALN.fasta
            trimal -in \${name}_ALN.fasta -out \${name}_aln.fasta -keepheader -fasta -automated1 -htmlout \${name}_aln.html
            """
    }

    process model_test_VAP {

        label 'norm_cpus'

        publishDir "${params.mypwd}/${params.outdir}/Analyses/model_test", mode: "copy", overwrite: true, pattern: '*.{tree,log}'
        publishDir "${params.mypwd}/${params.outdir}/Analyses/model_test/summaryfiles", mode: "copy", overwrite: true, pattern: '*model.summary'

        input:
            file(reads) from ntmodeltest

        output:
            file("*.tree") into tree_rax
            file("*.log") into log_rax
            file("*model.summary") into ntmodelsum

        script:
            """
            name=\$( echo ${reads} | awk -F "_aln" '{print \$1}')
            modeltest-ng -i ${reads} -p ${task.cpus} -d nt -s 203 --disable-checkpoint
            echo "Modeltest analysis complete :)"
            echo "Modeltest-ng results summary:" > \${name}_model.summary
            tail -7 \${name}_aln.fasta.log >> \${name}_model.summary
            """
    }

    process raxml_VAP {

        label 'rax_cpus'

        publishDir "${params.mypwd}/${params.outdir}/Analyses/trees", mode: "copy", overwrite: true

        input:
            file(align) from ntrax_align
            file(tree) from tree_rax
            file(log) from log_rax

        output:
            file("*raxml*") into raxml_summary

        script:
        if (params.ntraxcust) {
            """
            pre=\$(echo ${align} | awk -F "_aln.fasta" '{print \$1}' )
            raxml-ng --all --msa ${align} --prefix \${pre} --threads ${task.cpus} ${params.ntraxcust}
            """
        } else if (params.ntmodelt) {
            """
            pre=\$(echo ${align} | awk -F "_aln.fasta" '{print \$1}' )
            mod=\$(tail -14 ${log} | head -1 | awk -F "--model " '{print \$2}')
            raxml-ng --all --msa ${align} --model \${mod} --seed 2 --blopt nr_safe --threads ${task.cpus} --tree ${tree} --bs-trees autoMRE --prefix \${pre} --redo
            """
        } else {
            """
            pre=\$(echo ${align} | awk -F "_aln.fasta" '{print \$1}' )
            raxml-ng --all --msa ${align} --model GTR --seed 2 --blopt nr_safe --threads ${task.cpus} --tree ${tree} --bs-trees autoMRE --prefix \${pre} --redo
            """
        }
    }

    process proteinstage_VAP {

        label 'norm_cpus'

        conda 'python=2.7'

        publishDir "${params.mypwd}/${params.outdir}/AminoTypes/translation", mode: "copy", overwrite: true, pattern: '*'

        input:
            file(fasta) from asvsforAminotyping

        output:
            file("${params.projtag}_filtered_prots.fasta") into amintypegen
            file("${params.projtag}_allprot.fasta") into allprots
            file("${params.projtag}_translation_report") into proteinstage_vap_report

        script:
            """
            ${tools}/virtualribosomev2/dna2pep.py ${fasta} -r all -x -o none --fasta ${params.projtag}_allprot.fasta --report ${params.projtag}_translation_report
            awk 'BEGIN{RS=">";ORS=""}length(\$2)>${params.minAA}{print ">"\$0}' ${params.projtag}_allprot.fasta >${params.projtag}_filtered_prots.fasta
            """
    }

    process aminotyping {

        label 'norm_cpus'

        publishDir "${params.mypwd}/${params.outdir}/AminoTypes/", mode: "copy", overwrite: true

        input:
            file(prot) from amintypegen
            file(fasta) from allprots

        output:
            tuple file("${params.projtag}_aTypes.fasta"), file("${params.projtag}_aTypes.clstr"), file("${params.projtag}_aminoClus_sum.csv"), file ("${params.projtag}_problematic_translations.fasta"), file("${params.projtag}_problematic_nucleotides.fasta") into ( supplementalfiles )
            file("${params.projtag}_AminoTypes_noTax.fasta") into ( aminotypesCounts, aminotypesMafft, aminotypesClustal, aminotypesBlast )

        script:
            """
            cd-hit -i ${prot} -c 1.0 -o ${params.projtag}_aTypes.fasta
            sed 's/>Cluster />Cluster_/g' ${params.projtag}_aTypes.fasta.clstr >${params.projtag}_aTypes.clstr
            grep ">Cluster_" ${params.projtag}_aTypes.clstr >tmpclusters.list
            grep -w "*" ${params.projtag}_aTypes.clstr | awk '{print \$3}' | awk -F "." '{print \$1}' >tmphead.list
            grep -w "*" ${params.projtag}_aTypes.clstr | awk '{print \$2}' | awk -F "," '{print \$1}' >tmplen.list
            paste -d"," tmpclusters.list tmphead.list >tmp.info.csv
            grep ">" ${params.projtag}_aTypes.fasta >lala.list
            j=1
            for x in \$(cat lala.list);do
                echo ">${params.projtag}_AminoType\${j}" >>${params.projtag}_aminoheaders.list
                echo "\${x},>${params.projtag}_AminoType\${j}" >>tmpaminotype.info.csv
                j=\$(( \$j+1 ))
            done
            rm lala.list
            awk -F "," '{print \$2}' tmp.info.csv >>tmporder.list
            for x in \$(cat tmporder.list);do
             	grep -w "\$x" tmpaminotype.info.csv | awk -F "," '{print \$2}' >>tmpder.list
            done
            paste -d "," tmpclusters.list tmplen.list tmphead.list tmpder.list >${params.projtag}_aminoClus_sum.csv
            rm tmp*
            rename_seq.py ${params.projtag}_aTypes.fasta ${params.projtag}_aminoheaders.list ${params.projtag}_AminoTypes_noTax.fasta
            stats.sh in=${params.projtag}_AminoTypes_noTax.fasta gc=${params.projtag}_clustered.gc gcformat=4
            awk '{print \$2}' ${params.projtag}_clustered.gc | sort | uniq -c | sort -bgr | awk '{print \$2}' >>tmp.lenny
            awk 'BEGIN{RS=">";ORS=""}length(\$2)<${params.minAA}{print ">"\$0}' ${fasta} >${params.projtag}_problematic_translations.fasta
            grep ">" ${params.projtag}_problematic_translations.fasta | awk -F ">" '{print \$2}' >problem_tmp.list
            seqtk subseq ${fasta} problem_tmp.list >${params.projtag}_problematic_nucleotides.fasta
            """
    }

    process clustalo_prot_VAP {

        label 'norm_cpus'

        publishDir "${params.mypwd}/${params.outdir}/AminoTypes/matrix", mode: "copy", overwrite: true

        input:
            file(prot) from aminotypesClustal

        output:
            file("*.matrix") into proclustmatrices

        script:
            """
            name=\$( echo ${prot} | awk -F ".fasta" '{print \$1}')
            clustalo -i ${prot} --distmat-out=\${name}_dist.matrix --full --force --threads=${task.cpus}
            clustalo -i ${prot} --distmat-out=\${name}_id.matrix --percent-id --full --force --threads=${task.cpus}
            """
    }

    process diamondaminotypes_VAP {

        label 'norm_cpus'

        publishDir "${params.mypwd}/${params.outdir}/Taxonomy/aminotypes", mode: "copy", overwrite: true, pattern: '*.{csv,tsv}'
        publishDir "${params.mypwd}/${params.outdir}/Taxonomy/aminotypes/diamond_out", mode: "copy", overwrite: true, pattern: '*dmd.out'

        input:
            file(reads) from aminotypesBlast

        output:
            file("*.fasta") into ( diamondAA_ch, diamondAA_clustal_ch, ntmafftAA )
            tuple file("*_sumPHYtab.csv"), file("*_summarytable.tsv"), file("*dmd.out") into summary_AAdiamond

        script:
            """
            virdb=${params.mypwd}/DBs/diamonddb_custom/${params.dbname}
            name=\$(ls ${reads} | awk -F ".fasta" '{print \$1}')
            diamond blastx -q ${reads} -d \${virdb} -p ${task.cpus} --min-score 50 --more-sensitive -o "\$name"_dmd.out -f 6 qseqid qlen sseqid qstart qend qseq sseq length qframe evalue bitscore pident btop --max-target-seqs 1
            echo "Preparing lists to generate summary .csv's"
            echo "[Best hit accession number]" >access.list
            echo "[pOTU sequence length]" >length.list
            echo "[e-value]" >evalue.list
            echo "[Bitscore]" >bit.list
            echo "[Percent ID (aa)]" >pid.list
            echo "[pOTU#]" >otu.list
            echo "[Virus ID]" >"\$name"_virus.list
            echo "[Gene]" >"\$name"_genes.list
            grep ">" ${reads} | awk -F ">" '{print \$2}' > seqids.lst
            echo "extracting genes and names"
            touch new_"\$name"_asvnames.txt
            j=1
            for s in \$(cat seqids.lst);do
                echo "Checking for \$s hit in diamond output"
                    if [[ ${params.refseq} == "T" ]];then
                        echo "RefSeq headers specified"
                        if [[ "\$(grep -wc "\$s" "\$name"_dmd.out)" -eq 1 ]];then
                            echo "Yep, there was a hit for \$s"
                            echo "Extracting the information now:"
                            acc=\$(grep -w "\$s" "\$name"_dmd.out | awk '{print \$3}')
                            echo "\$s" >> otu.list
                            echo "\$acc" >> access.list
                            line="\$(grep -w "\$s" "\$name"_dmd.out)"
                            echo "\$line" | awk '{print \$10}' >>evalue.list
                            echo "\$line" | awk '{print \$11}' >>bit.list
                            echo "\$line" | awk '{print \$12}' >>pid.list
                            echo "\$line" | awk '{print \$2}' >>length.list
                            echo "Extracting virus and gene ID for \$s now"
                            gene=\$(grep -w "\$acc" "\$headers" | awk -F "." '{ print \$2 }' | awk -F "[" '{ print \$1 }' | awk -F " " print substr(\$0, index(\$0,\$2)) | sed 's/ /_/g') &&
                            echo "\$gene" | sed 's/_/ /g' >> "\$name"_genes.list
                            virus=\$(grep -w "\$acc" "\$headers" | awk -F "[" '{ print \$2 }' | awk -F "]" '{ print \$1 }'| sed 's/ /_/g')
                            echo "\$virus" | sed 's/_/ /g' >> "\$name"_virus.list
                            echo ">pOTU\${j}_"\$virus"_"\$gene"" >> new_"\$name"_asvnames.txt
                            j=\$((\$j+1))
                            echo "\$s done."
                        else
                            echo "Ugh, there was no hit for \$s .."
                            echo "We still love \$s though and we will add it to the final fasta file"
                            echo "\$s" >> otu.list
                            echo "NO_HIT" >>access.list
                            echo "NO_HIT" >>"\$name"_genes.list
                            echo "NO_HIT" >>"\$name"_virus.list
                            echo "NO_HIT" >>evalue.list
                            echo "NO_HIT" >>bit.list
                            echo "NO_HIT" >>pid.list
                            echo "NO_HIT" >>length.list
                            virus="NO"
                            gene="HIT"
                            echo ">pOTU\${j}_"\$virus"_"\$gene"" >> new_"\$name"_asvnames.txt
                            j=\$((\$j+1))
                            echo "\$s done."
                    fi
                    else
                        echo "Using RVDB headers."
                        if [[ "\$(grep -wc "\$s" "\$name"_dmd.out)" -eq 1 ]];then
                            echo "Yep, there was a hit for \$s"
                            echo "Extracting the information now:"
                            acc=\$(grep -w "\$s" "\$name"_dmd.out | awk '{print \$3}' | awk -F "|" '{print \$3}')
                            echo "\$s" >>otu.list
                            echo "\$acc" >>access.list
                            line="\$(grep -w "\$s" "\$name"_dmd.out)"
                            echo "\$line" | awk '{print \$10}' >>evalue.list
                            echo "\$line" | awk '{print \$11}' >>bit.list
                            echo "\$line" | awk '{print \$12}' >>pid.list
                            echo "\$line" | awk '{print \$2}' >>length.list
                            echo "Extracting virus and gene ID for \$s now"
                            gene=\$(grep -w "\$acc" "\$headers" | awk -F "|" '{ print \$6 }' | awk -F "[" '{ print \$1 }' | sed 's/ /_/g') &&
                            echo "\$gene" | sed 's/_/ /g' >>"\$name"_genes.list
                            virus=\$(grep -w "\$acc" "\$headers" | awk -F "|" '{ print \$6 }' | awk -F "[" '{ print \$2 }' | awk -F "]" '{print \$1}' | sed 's/ /_/g') &&
                            echo "\$virus" | sed 's/_/ /g' >>"\$name"_virus.list
                            echo ">pOTU\${j}_"\$virus"_"\$gene"" >>new_"\$name"_asvnames.txt
                            j=\$((\$j+1))
                            echo "\$s done."
                        else
                            echo "Ugh, there was no hit for \$s .."
                            echo "We still love \$s though and we will add it to the final fasta file"
                            echo "\$s" >>otu.list
                            echo "NO_HIT" >>access.list
                            echo "NO_HIT" >>"\$name"_genes.list
                            echo "NO_HIT" >>"\$name"_virus.list
                            echo "NO_HIT" >>evalue.list
                            echo "NO_HIT" >>bit.list
                            echo "NO_HIT" >>pid.list
                            echo "NO_HIT" >>length.list
                            virus="NO"
                            gene="HIT"
                            echo ">pOTU\${j}_"\$virus"_"\$gene"" >>new_"\$name"_asvnames.txt
                            j=\$((\$j+1))
                            echo "\$s done."
                        fi
                    fi
            echo "Done with \$s"
            done
            echo "Now editing "\$name" fasta headers"
            ###### rename_seq.py
            rename_seq.py ${reads} new_"\$name"_asvnames.txt "\$name"_fin.fasta
            awk 'BEGIN {RS=">";FS="\\n";OFS=""} NR>1 {print ">"\$1; \$1=""; print}' "\$name"_fin.fasta > "\$name"_tmpssasv.fasta
            echo "[Sequence header]" > newnames.list
            cat new_"\$name"_asvnames.txt >> newnames.list
            touch sequence.list
            echo "     " > sequence.list
            grep -E '^[ACGT]+\$' "\$name"_tmpssasv.fasta >> sequence.list
            rm "\$name"_tmpssasv.fasta
            paste -d "," sequence.list "\$name"_virus.list "\$name"_genes.list otu.list newnames.list length.list bit.list evalue.list pid.list access.list >> "\$name"_sumPHYtab.csv
            paste -d"\t" otu.list access.list "\$name"_virus.list "\$name"_genes.list sequence.list length.list bit.list evalue.list pid.list >> "\$name"_summarytable.tsv
            rm evalue.list ; rm sequence.list ; rm bit.list ; rm pid.list ; rm length.list ;
            rm new_"\$name"_asvnames.txt
            rm "\$name"_virus.list
            rm "\$name"_genes.list
            rm newnames.list
            rm access.list
            echo "Done Assigning Taxonomy To : ${reads} "
            """
    }

    process mafft_prot_VAP {

        label 'norm_cpus'

        publishDir "${params.mypwd}/${params.outdir}/AminoTypes/alignment", mode: "copy", overwrite: true

        input:
            file(prot) from aminotypesMafft

        output:
            file("*_aln.fasta") into ( protmodeltest, protrax_align )
            file("${name}_aln.html") into prot_aln

        script:
            """
            name=\$( echo ${prot}  | awk -F ".fasta" '{print \$1}' )
            mafft --maxiterate 5000 --auto ${prot} >\${name}_ALN.fasta
            trimal -in \${name}_ALN.fasta -out \${name}_aln.fasta -keepheader -fasta -automated1 -htmlout \${name}_aln.html
            """
    }

    process modeltest_AminoTypes_VAP {

        label 'norm_cpus'

        publishDir "${params.mypwd}/${params.outdir}/AminoTypes/model_test", mode: "copy", overwrite: true, pattern: '*.{tree,log}'
        publishDir "${params.mypwd}/${params.outdir}/AminoTypes/model_test/summaryfiles", mode: "copy", overwrite: true, pattern: '*model.summary'

        input:
            file(prot) from protmodeltest

        output:
            file("*.tree") into protree_rax
            file("*.log") into prolog_rax
            file("*model.summary") into aamodelsum

        script:
            """
            name=\$( echo ${prot} | awk -F "_aln" '{print \$1}')
            modeltest-ng -i ${prot} -p ${task.cpus} -d aa -s 203 --disable-checkpoint
            echo "Modeltest analysis complete :)"
            echo "Modeltest-ng results summary:" > \${name}_model.summary
            tail -7 \${name}_aln.fasta.log >> \${name}_model.summary
            """
    }

    process raxml_prot_VAP {

        label 'rax_cpus'

        publishDir "${params.mypwd}/${params.outdir}/AminoTypes/trees", mode: "copy", overwrite: true

        input:
            file(palign) from protrax_align
            file(ptree) from protree_rax
            file(plog) from prolog_rax

        output:
            file("*raxml*") into raxml_prot_summary

        script:
            if (params.ptraxcust) {
                """
                pre=\$(echo ${palign} | awk -F "_aln.fasta" '{print \$1}' )
                raxml-ng --all --msa ${palign} --prefix \${pre} --threads ${task.cpus} ${params.ptraxcust}
                """
            } else if (params.ptmodeltrax) {
                """
                pre=\$(echo ${palign} | awk -F "_aln.fasta" '{print \$1}' )
                mod=\$(tail -14 ${plog} | head -1 | awk -F "--model " '{print \$2}')
                raxml-ng --all --msa ${palign} --model \${mod} --seed 2 --blopt nr_safe --threads ${task.cpus} --tree ${ptree} --bs-trees autoMRE --prefix \${pre} --redo
                """
            } else {
                """
                pre=\$(echo ${palign} | awk -F "_aln.fasta" '{print \$1}' )
                raxml-ng --all --msa ${palign} --model protGTR --seed 2 --blopt nr_safe --threads ${task.cpus} --tree ${ptree} --bs-trees autoMRE --prefix \${pre} --redo
                """
            }
    }

    process aminoTypes_counts_VAP {

        label 'norm_cpus'

        publishDir "${params.mypwd}/${params.outdir}/AminoTypes/counts", mode: "copy", overwrite: true

        input:
            file(fasta) from aminotypesCounts
            file(merged) from mergeforprotcounts
            file(samplist) from samplelist

        output:
            tuple file("*_protcounts.csv"), file("*dmd.out") into counts_summary

        script:
            """
            diamond makedb --in ${fasta} --db ${fasta}
            diamond blastx -q ${merged} -d ${fasta} -p ${task.cpus} --min-score ${params.ProtCountsBit} --id ${params.ProtCountID} -l ${params.ProtCountsLength} --more-sensitive -o ${params.projtag}_protCounts_dmd.out -f 6 qseqid qlen sseqid qstart qend qseq sseq length qframe evalue bitscore pident btop --max-target-seqs 1 --max-hsps 1
            echo "[Sequence]" >tmp.col1.txt
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
           paste -d "," tmp.col1.txt *col.txt > ${params.projtag}_protcounts.csv
           rm tmp*
           rm *col.txt
           """
    }

    if (params.clusteraaASV) {        // ASV_nucl -> ASV_aa -> clusteraa by %id with ch-hit -> extract representative nucl sequences to generate new OTU file

        process translateforclustering_VAP2 {

            label 'norm_cpus'

            conda 'python=2.7'

            publishDir "${params.mypwd}/${params.outdir}/Clustering/translation4clustering", mode: "copy", overwrite: true

            input:
                file(fasta) from nucl2aa

            output:
                file("*ASVprotforclust.fasta") into cluster_aa
                file("*_vr_report") into reportaa_VR
                file("*_ASVs_all.fasta") into asvfastaforaaclust

            script:
                """
                name=\$( echo ${fasta} | awk -F "_ASV" '{print \$1}')
                ${tools}/virtualribosomev2/dna2pep.py ${fasta} -r all -x -o none --fasta ${params.projtag}_allprot.fasta --report ${params.projtag}_translation_report
                """
        }

        process cluster_AA {

            label 'norm_cpus'

            publishDir "${params.mypwd}/${params.outdir}/Clustering/", mode: "copy", overwrite: true,

            input:
                file(fasta) from cluster_aa
                file(asvs) from asvfastaforaaclust

            output:
                file("${params.projtag}_nucleotides_pOTU\${id}.fasta") into pOTU_diamond_ch
                file("*_ASVs_all.fasta") into asvfastados
                file("${params.projtag}_aminos_pOTU\${id}_noTax.fasta") into ( pOTUaaforanalysis, pOTUaa_diamondbp, pOTUaaformafft, pOTUaaforcounts )

            script:
            // add awk script to count seqs
            if (params.aalist) {
                """
                for id in `echo ${params.aaidlist} | tr "," "\n"`;done
                    awk 'BEGIN{RS=">";ORS=""}length(\$2)>${params.minAA}{print ">"\$0}' ${fasta} > ${params.projtag}_filtered_proteins.fasta
                    cd-hit -i filtered.fasta -c \${id} --fasta ${params.projtag}_pOTU\${id}.fasta
                    sed 's/>Cluster />Cluster_/g' ${params.projtag}_pOTU\${id}.fasta.clstr >${params.projtag}_pOTU\${id}.clstr
                    grep ">Cluster_" ${params.projtag}_pOTU\${id}.clstr >tmpclusters.list
                        for x in \$(cat tmpclusters.list);
                        do 	    echo "Extracting \$x"
                                name=\$( echo \$x | awk -F ">" '{print \$2}')
	                            sed -n "/\$x/,/Cluster_/p" ${params.projtag}_pOTU\${id}.clstr > "\$name"_tmp.list
                        done
                        j=1
                        for x in *_tmp.list
                        do  	   name=\$(echo \$x | awk -F "_tmp" '{print \$1}')
	                               cluster="pOTU"\${j}""
	                               grep "ASV" \$x | awk -F ", " '{print \$2}' | awk -F "_" '{print \$1}' | awk -F ">" '{print \$2}' > seqs_tmps.list
	                               seqtk subseq ${asvs} seqs_tmps.list > \${cluster}_nucleotide_sequence.fasta
                                   vsearch --cluster_fast \${cluster}_nucleotide_sequences.fasta --id 0.2--centroids \${name}_centroids.fasta
                                   grep ">" \${name}_centroids.fasta >> tmp_cemtroids.list
                                   u="1"
                                   for y in \$( cat tmp_centroids.list)
                                   do   echo ">\${cluster}_type"\$u"" >> tmp_centroid.newheaders
                                        u=\$(( \$u+1 ))
                                   done
                                   rename_seq.py \${name}_centroids.fasta tmp_centroid.newheaders \${cluster}_types_labeled.fasta
                                   rm *tmp*
                        done
                        cat *_types_labeled.fasta >> ${params.projtag}_nucleotides_pOTU\${id}.fasta
                    grep -w "*" ${params.projtag}_pOTU\${id}.clstr | awk '{print \$3}' | awk -F "." '{print \$1}' >tmphead.list
                    grep -w "*" ${params.projtag}_pOTU\${id}.clstr | awk '{print \$2}' | awk -F "," '{print \$1}' >tmplen.list
                    paste -d"," tmpclusters.list tmphead.list >tmp.info.csv
                    grep ">" ${params.projtag}_pOTU\${id}.fasta >lala.list
                    j=1
                    for x in \$(cat lala.list);do
                        echo ">${params.projtag}_pOTU\${j}" >>${params.projtag}_aminoheaders.list
                        echo "\${x},>${params.projtag}_pOTU\${j}" >>tmpaminotype.info.csv
                        j=\$(( \$j+1 ))
                    done
                    rm lala.list
                    awk -F "," '{print \$2}' tmp.info.csv >>tmporder.list
                    for x in \$(cat tmporder.list);do
                     	grep -w "\$x" tmpaminotype.info.csv | awk -F "," '{print \$2}' >>tmpder.list
                    done
                    paste -d "," tmpclusters.list tmplen.list tmphead.list tmpder.list >${params.projtag}_pOTUCluster_summary.csv
                    rm tmp*
                    rename_seq.py ${params.projtag}_pOTU\${id}.fasta ${params.projtag}_aminoheaders.list ${params.projtag}_aminos_pOTU\${id}_noTax.fasta
                        stats.sh in=${params.projtag}_AminoTypes_noTax.fasta gc=${params.projtag}_clustered.gc gcformat=4
                    awk 'BEGIN{RS=">";ORS=""}length(\$2)<${params.minAA}{print ">"\$0}' ${fasta} >${params.projtag}_problematic_translations.fasta
                    grep ">" ${params.projtag}_problematic_translations.fasta | awk -F ">" '{print \$2}' > problem_tmp.list
                    seqtk subseq ${asvs} problem_tmp.list > ${params.projtag}_problematic_nucleotides.fasta
                    done
                done
                """
            } else if (params.aaid) {
                """
                id=${params.aaid}
                awk 'BEGIN{RS=">";ORS=""}length(\$2)>${params.minAA}{print ">"\$0}' ${fasta} > ${params.projtag}_filtered_proteins.fasta
                cd-hit -i filtered.fasta -c \${id} --fasta ${params.projtag}_pOTU\${id}.fasta
                sed 's/>Cluster />Cluster_/g' ${params.projtag}_pOTU\${id}.fasta.clstr >${params.projtag}_pOTU\${id}.clstr
                grep ">Cluster_" ${params.projtag}_pOTU\${id}.clstr >tmpclusters.list
                for x in \$(cat tmpclusters.list);
                do 	    echo "Extracting \$x"
                        name=\$( echo \$x | awk -F ">" '{print \$2}')
                        sed -n "/\$x/,/Cluster_/p" ${params.projtag}_pOTU\${id}.clstr > "\$name"_tmp.list
                done
                j=1
                for x in *_tmp.list
                do  	   name=\$(echo \$x | awk -F "_tmp" '{print \$1}')
                           cluster="pOTU"\${j}""
                           grep "ASV" \$x | awk -F ", " '{print \$2}' | awk -F "_" '{print \$1}' | awk -F ">" '{print \$2}' > seqs_tmps.list
                           seqtk subseq ${asvs} seqs_tmps.list > \${cluster}_nucleotide_sequence.fasta
                           vsearch --cluster_fast \${cluster}_nucleotide_sequences.fasta --id 0.2--centroids \${name}_centroids.fasta
                           grep ">" \${name}_centroids.fasta >> tmp_cemtroids.list
                           u="1"
                           for y in \$( cat tmp_centroids.list)
                           do   echo ">\${cluster}_type"\$u"" >> tmp_centroid.newheaders
                                u=\$(( \$u+1 ))
                           done
                           rename_seq.py \${name}_centroids.fasta tmp_centroid.newheaders \${cluster}_types_labeled.fasta
                           rm *tmp*
                done
                cat *_types_labeled.fasta >> ${params.projtag}_nucleotides_pOTU\${id}.fasta
                grep -w "*" ${params.projtag}_pOTU\${id}.clstr | awk '{print \$3}' | awk -F "." '{print \$1}' >tmphead.list
                grep -w "*" ${params.projtag}_pOTU\${id}.clstr | awk '{print \$2}' | awk -F "," '{print \$1}' >tmplen.list
                paste -d"," tmpclusters.list tmphead.list >tmp.info.csv
                grep ">" ${params.projtag}_pOTU\${id}.fasta >lala.list
                j=1
                for x in \$(cat lala.list);do
                    echo ">${params.projtag}_pOTU\${j}" >>${params.projtag}_aminoheaders.list
                    echo "\${x},>${params.projtag}_pOTU\${j}" >>tmpaminotype.info.csv
                    j=\$(( \$j+1 ))
                done
                rm lala.list
                awk -F "," '{print \$2}' tmp.info.csv >>tmporder.list
                for x in \$(cat tmporder.list);do
                    grep -w "\$x" tmpaminotype.info.csv | awk -F "," '{print \$2}' >>tmpder.list
                done
                paste -d "," tmpclusters.list tmplen.list tmphead.list tmpder.list >${params.projtag}_pOTUCluster_summary.csv
                rm tmp*
                rename_seq.py ${params.projtag}_pOTU\${id}.fasta ${params.projtag}_aminoheaders.list ${params.projtag}_aminos_pOTU\${id}_noTax.fasta
                stats.sh in=${params.projtag}_AminoTypes_noTax.fasta gc=${params.projtag}_clustered.gc gcformat=4
                awk 'BEGIN{RS=">";ORS=""}length(\$2)<${params.minAA}{print ">"\$0}' ${fasta} >${params.projtag}_problematic_translations.fasta
                grep ">" ${params.projtag}_problematic_translations.fasta | awk -F ">" '{print \$2}' > problem_tmp.list
                seqtk subseq ${asvs} problem_tmp.list > ${params.projtag}_problematic_nucleotides.fasta
                """
            }
        }

        process diamondpOTU_VAP2 {

            label 'norm_cpus'

            publishDir "${params.mypwd}/${params.outdir}/Taxonomy/", mode: "copy", overwrite: true, pattern: '*.{csv,tsv}'
            publishDir "${params.mypwd}/${params.outdir}/Taxonomy/diamond_out", mode: "copy", overwrite: true, pattern: '*dmd.out'

            input:
                file(reads) from pOTU_diamond_ch

            output:
                file("*.fasta") into ( diamondAA_ch, diamondAA_clustal_ch, ntmafftAA )
                tuple file("*_sumPHYtab.csv"), file("*_summarytable.tsv"), file("*dmd.out") into summary_AAdiamond

            script:
                """
                virdb=${params.mypwd}/DBs/diamonddb_custom/${params.dbname}
                name=\$(ls ${reads} | awk -F ".fasta" '{print \$1}')
                diamond blastx -q ${reads} -d \${virdb} -p ${task.cpus} --min-score 50 --more-sensitive -o "\$name"_dmd.out -f 6 qseqid qlen sseqid qstart qend qseq sseq length qframe evalue bitscore pident btop --max-target-seqs 1
                echo "Preparing lists to generate summary .csv's"
                echo "[Best hit accession number]" >access.list
                echo "[pOTU sequence length]" >length.list
                echo "[e-value]" >evalue.list
                echo "[Bitscore]" >bit.list
                echo "[Percent ID (aa)]" >pid.list
                echo "[pOTU#]" >otu.list
                echo "[Virus ID]" >"\$name"_virus.list
                echo "[Gene]" >"\$name"_genes.list
                grep ">" ${reads} | awk -F ">" '{print \$2}' > seqids.lst
                echo "extracting genes and names"
                touch new_"\$name"_asvnames.txt
                j=1
                for s in \$(cat seqids.lst);do
                    echo "Checking for \$s hit in diamond output"
                        if [[ ${params.refseq} == "T" ]];then
                            echo "RefSeq headers specified"
                            if [[ "\$(grep -wc "\$s" "\$name"_dmd.out)" -eq 1 ]];then
                                echo "Yep, there was a hit for \$s"
                                echo "Extracting the information now:"
                                acc=\$(grep -w "\$s" "\$name"_dmd.out | awk '{print \$3}')
                                echo "\$s" >> otu.list
                                echo "\$acc" >> access.list
                                line="\$(grep -w "\$s" "\$name"_dmd.out)"
                                echo "\$line" | awk '{print \$10}' >>evalue.list
                                echo "\$line" | awk '{print \$11}' >>bit.list
                                echo "\$line" | awk '{print \$12}' >>pid.list
                                echo "\$line" | awk '{print \$2}' >>length.list
                                echo "Extracting virus and gene ID for \$s now"
                                gene=\$(grep -w "\$acc" "\$headers" | awk -F "." '{ print \$2 }' | awk -F "[" '{ print \$1 }' | awk -F " " print substr(\$0, index(\$0,\$2)) | sed 's/ /_/g') &&
                                echo "\$gene" | sed 's/_/ /g' >> "\$name"_genes.list
                                virus=\$(grep -w "\$acc" "\$headers" | awk -F "[" '{ print \$2 }' | awk -F "]" '{ print \$1 }'| sed 's/ /_/g')
                                echo "\$virus" | sed 's/_/ /g' >> "\$name"_virus.list
                                echo ">pOTU\${j}_"\$virus"_"\$gene"" >> new_"\$name"_asvnames.txt
                                j=\$((\$j+1))
                                echo "\$s done."
                            else
                                echo "Ugh, there was no hit for \$s .."
                                echo "We still love \$s though and we will add it to the final fasta file"
                                echo "\$s" >> otu.list
                                echo "NO_HIT" >>access.list
                                echo "NO_HIT" >>"\$name"_genes.list
                                echo "NO_HIT" >>"\$name"_virus.list
                                echo "NO_HIT" >>evalue.list
                                echo "NO_HIT" >>bit.list
                                echo "NO_HIT" >>pid.list
                                echo "NO_HIT" >>length.list
                                virus="NO"
                                gene="HIT"
                                echo ">pOTU\${j}_"\$virus"_"\$gene"" >> new_"\$name"_asvnames.txt
                                j=\$((\$j+1))
                                echo "\$s done."
                        fi
                        else
                            echo "Using RVDB headers."
                            if [[ "\$(grep -wc "\$s" "\$name"_dmd.out)" -eq 1 ]];then
                                echo "Yep, there was a hit for \$s"
                                echo "Extracting the information now:"
                                acc=\$(grep -w "\$s" "\$name"_dmd.out | awk '{print \$3}' | awk -F "|" '{print \$3}')
                                echo "\$s" >>otu.list
                                echo "\$acc" >>access.list
                                line="\$(grep -w "\$s" "\$name"_dmd.out)"
                                echo "\$line" | awk '{print \$10}' >>evalue.list
                                echo "\$line" | awk '{print \$11}' >>bit.list
                                echo "\$line" | awk '{print \$12}' >>pid.list
                                echo "\$line" | awk '{print \$2}' >>length.list
                                echo "Extracting virus and gene ID for \$s now"
                                gene=\$(grep -w "\$acc" "\$headers" | awk -F "|" '{ print \$6 }' | awk -F "[" '{ print \$1 }' | sed 's/ /_/g') &&
                                echo "\$gene" | sed 's/_/ /g' >>"\$name"_genes.list
                                virus=\$(grep -w "\$acc" "\$headers" | awk -F "|" '{ print \$6 }' | awk -F "[" '{ print \$2 }' | awk -F "]" '{print \$1}' | sed 's/ /_/g') &&
                                echo "\$virus" | sed 's/_/ /g' >>"\$name"_virus.list
                                echo ">pOTU\${j}_"\$virus"_"\$gene"" >>new_"\$name"_asvnames.txt
                                j=\$((\$j+1))
                                echo "\$s done."
                            else
                                echo "Ugh, there was no hit for \$s .."
                                echo "We still love \$s though and we will add it to the final fasta file"
                                echo "\$s" >>otu.list
                                echo "NO_HIT" >>access.list
                                echo "NO_HIT" >>"\$name"_genes.list
                                echo "NO_HIT" >>"\$name"_virus.list
                                echo "NO_HIT" >>evalue.list
                                echo "NO_HIT" >>bit.list
                                echo "NO_HIT" >>pid.list
                                echo "NO_HIT" >>length.list
                                virus="NO"
                                gene="HIT"
                                echo ">pOTU\${j}_"\$virus"_"\$gene"" >>new_"\$name"_asvnames.txt
                                j=\$((\$j+1))
                                echo "\$s done."
                            fi
                        fi
                echo "Done with \$s"
                done
                echo "Now editing "\$name" fasta headers"
                ###### rename_seq.py
                rename_seq.py ${reads} new_"\$name"_asvnames.txt "\$name"_fin.fasta
                awk 'BEGIN {RS=">";FS="\\n";OFS=""} NR>1 {print ">"\$1; \$1=""; print}' "\$name"_fin.fasta > "\$name"_tmpssasv.fasta
                echo "[Sequence header]" > newnames.list
                cat new_"\$name"_asvnames.txt >> newnames.list
                touch sequence.list
                echo "     " > sequence.list
                grep -E '^[ACGT]+\$' "\$name"_tmpssasv.fasta >> sequence.list
                rm "\$name"_tmpssasv.fasta
                paste -d "," sequence.list "\$name"_virus.list "\$name"_genes.list otu.list newnames.list length.list bit.list evalue.list pid.list access.list >> "\$name"_sumPHYtab.csv
                paste -d"\t" otu.list access.list "\$name"_virus.list "\$name"_genes.list sequence.list length.list bit.list evalue.list pid.list >> "\$name"_summarytable.tsv
                rm evalue.list ; rm sequence.list ; rm bit.list ; rm pid.list ; rm length.list ;
                rm new_"\$name"_asvnames.txt
                rm "\$name"_virus.list
                rm "\$name"_genes.list
                rm newnames.list
                rm access.list
                echo "Done Assigning Taxonomy To : ${reads} "
                """
        }

        process vsearch_globalpOTU_VAP2 {

            label 'norm_cpus'

            publishDir "${params.mypwd}/${params.outdir}/Analyses/counts", mode: "copy", overwrite: true, pattern: '*.{biome,txt}'

            input:
                file(reads) from diamondAA_ch
                file("all_merged_clean.fasta") from mergeforpOTUcounts

            output:
                tuple file("*_counts.txt"), file("*_counts.biome") into pOTUcounts_vsearch

            script:
                """
                 if [ `echo ${reads} | grep -c "fin"` -eq 1 ];then
                    name=\$( echo ${reads} | awk -F "_fin" '{print \$1}')
                    vsearch --usearch_global all_merged_clean.fasta --db ${reads} --id ${params.asvID} --threads ${task.cpus} --otutabout "\$name"_counts.txt --biomout "\$name"_counts.biome
                else
                    ident=\$( echo ${reads} | awk -F "otu" '{print \$2}' | awk -F ".fasta" '{print \$1}')
                    name=\$( echo ${reads} | awk -F ".fasta" '{print \$1}')
                    vsearch --usearch_global all_merged_clean.fasta --db ${reads} --id \$ident --threads ${task.cpus} --otutabout "\$name"_counts.txt --biomout "\$name"_counts.biome
                fi
                """
        }

        process pOTUclustalo_VAP2 {

            label 'norm_cpus'

            publishDir "${params.mypwd}/${params.outdir}/Clustering/matrix", mode: "copy", overwrite: true

            input:
                file(reads) from diamondAA_clustal_ch

            output:
                file("*.matrix") into pOTUclustmatrices

            script:
                """
                if [ `echo ${reads} | grep -c "fin"` -eq 1 ];then
                    dent=\$( echo ${reads} | awk -F "otu" '{print \$2}' | awk -F "_fin" '{print \$1}')
                    name=\$( echo ${reads} | awk -F "_fin" '{print \$1}')
                    clustalo -i ${reads} --distmat-out=\${name}_dist.matrix --full --force --threads=${task.cpus}
                    clustalo -i ${reads} --distmat-out=\${name}_id.matrix --percent-id --full --force --threads=${task.cpus}
                else
                    ident=\$( echo ${reads} | awk -F "otu" '{print \$2}' | awk -F ".fasta" '{print \$1}')
                    name=\$( echo ${reads} | awk -F ".fasta" '{print \$1}')
                    clustalo -i ${reads} --distmat-out=\${name}_dist.matrix --full --force --threads=${task.cpus}
                    clustalo -i ${reads} --distmat-out=\${name}_id.matrix --percent-id --full --force --threads=${task.cpus}
                fi
                """
        }

        process pOTUmafft_VAP {

            label 'norm_cpus'

            publishDir "${params.mypwd}/${params.outdir}/Analyses/alignment", mode: "copy", overwrite: true

            input:
                file(reads) from ntmafftAA

            output:
                file("*_aln.fasta") into ( pOTUntmodeltest, pOTUntrax_align )
                file("${name}_aln.html") into align_html

            script:
                """
                if [ `echo ${reads} | grep -c "fin"` -eq 1 ];then
                    name=\$( echo ${reads} | awk -F "_fin" '{print \$1}' )
                    mafft --maxiterate 5000 --auto ${reads} >\${name}_ALN.fasta
                    trimal -in \${name}_ALN.fasta -out \${name}_aln.fasta -keepheader -fasta -automated1 -htmlout \${name}_aln.html
                fi
                """
        }

        process pOTUmodel_test_VAP2 {

            label 'norm_cpus'

            publishDir "${params.mypwd}/${params.outdir}/Analyses/modeltest", mode: "copy", overwrite: true, pattern: '*.{tree,log}'
            publishDir "${params.mypwd}/${params.outdir}/Analyses/modeltest/summaryfiles", mode: "copy", overwrite: true, pattern: '*model.summary'

            input:
                file(reads) from pOTUntmodeltest

            output:
                file("*.tree") into tree_raxpOTU
                file("*.log") into log_raxpOTU
                file("*model.summary") into ntmodelsumpOTU

            script:
                """
                name=\$( echo ${reads} | awk -F "_aln" '{print \$1}')
                modeltest-ng -i ${reads} -p ${task.cpus} -d nt -s 203 --disable-checkpoint
                echo "Modeltest analysis complete :)"
                echo "Modeltest-ng results summary:" > \${name}_model.summary
                tail -7 \${name}_aln.fasta.log >> \${name}_model.summary
                """
        }

        process pOTUraxml_VAP2 {

            label 'rax_cpus'

            publishDir "${params.mypwd}/${params.outdir}/Analyses/trees", mode: "copy", overwrite: true

            input:
                file(align) from ntrax_alignpOTU
                file(tree) from tree_raxpOTU
                file(log) from log_raxpOTU

            output:
                file("*raxml*") into raxml_summary

            script:
            if (params.ntraxcust) {
                """
                pre=\$(echo ${align} | awk -F "_aln.fasta" '{print \$1}' )
                raxml-ng --all --msa ${align} --prefix \${pre} --threads ${task.cpus} ${params.ntraxcust}
                """
            } else if (params.ntmodelt) {
                """
                pre=\$(echo ${align} | awk -F "_aln.fasta" '{print \$1}' )
                mod=\$(tail -14 ${log} | head -1 | awk -F "--model " '{print \$2}')
                raxml-ng --all --msa ${align} --model \${mod} --seed 2 --blopt nr_safe --threads ${task.cpus} --tree ${tree} --bs-trees autoMRE --prefix \${pre} --redo
                """
            } else {
                """
                pre=\$(echo ${align} | awk -F "_aln.fasta" '{print \$1}' )
                raxml-ng --all --msa ${align} --model GTR --seed 2 --blopt nr_safe --threads ${task.cpus} --tree ${tree} --bs-trees autoMRE --prefix \${pre} --redo
                """
            }
        }

        process pOTUclustalo_prot_VAP2 {

            label 'norm_cpus'

            publishDir "${params.mypwd}/${params.outdir}/Analyses/protein/matrix", mode: "copy", overwrite: true

            input:
                file(prot) from pOTUaaforanalysis

            output:
                file("*.matrix") into pOTUaaMatrix

            script:
                """
                name=\$( echo ${prot} | awk -F ".fasta" '{print \$1}')
                clustalo -i ${prot} --distmat-out=\${name}_dist.matrix --full --force --threads=${task.cpus}
                clustalo -i ${prot} --distmat-out=\${name}_id.matrix --percent-id --full --force --threads=${task.cpus}
                """
        }

        process pOTUaa_diamondblastp {

            label 'norm_cpus'

            publishDir "${params.mypwd}/${params.outdir}/Analyses/protein/taxonomy", mode: "copy", overwrite: true

            input:
                file(reads) from pOTUaa_diamondbp
            output:
                file("_dmd.out") into pOTUdiamblastp
                file("*newheaders.list") into pOTUnewheaders

            script:
                """
                virdb=${params.mypwd}/DBs/diamonddb_custom/${params.dbname}
                name=\$(ls ${reads} | awk -F ".fasta" '{print \$1}')
                diamond blastp -q ${reads} -d \${virdb} -p ${task.cpus} --min-score 50 --more-sensitive -o "\$name"_dmd.out -f 6 qseqid qlen sseqid qstart qend qseq sseq length qframe evalue bitscore pident btop --max-target-seqs 1
                echo "Preparing lists to generate summary .csv's"
                echo "[Best hit accession number]" >access.list
                echo "[pOTU sequence length]" >length.list
                echo "[e-value]" >evalue.list
                echo "[Bitscore]" >bit.list
                echo "[Percent ID (aa)]" >pid.list
                echo "[pOTUaa#]" >otu.list
                echo "[Virus ID]" >"\$name"_virus.list
                echo "[Gene]" >"\$name"_genes.list
                grep ">" ${reads} | awk -F ">" '{print \$2}' > seqids.lst
                echo "extracting genes and names"
                touch new_"\$name"_asvnames.txt
                j=1
                for s in \$(cat seqids.lst);do
                    echo "Checking for \$s hit in diamond output"
                        if [[ ${params.refseq} == "T" ]];then
                            echo "RefSeq headers specified"
                            if [[ "\$(grep -wc "\$s" "\$name"_dmd.out)" -eq 1 ]];then
                                echo "Yep, there was a hit for \$s"
                                echo "Extracting the information now:"
                                acc=\$(grep -w "\$s" "\$name"_dmd.out | awk '{print \$3}')
                                echo "\$s" >> otu.list
                                echo "\$acc" >> access.list
                                line="\$(grep -w "\$s" "\$name"_dmd.out)"
                                echo "\$line" | awk '{print \$10}' >>evalue.list
                                echo "\$line" | awk '{print \$11}' >>bit.list
                                echo "\$line" | awk '{print \$12}' >>pid.list
                                echo "\$line" | awk '{print \$2}' >>length.list
                                echo "Extracting virus and gene ID for \$s now"
                                gene=\$(grep -w "\$acc" "\$headers" | awk -F "." '{ print \$2 }' | awk -F "[" '{ print \$1 }' | awk -F " " print substr(\$0, index(\$0,\$2)) | sed 's/ /_/g') &&
                                echo "\$gene" | sed 's/_/ /g' >> "\$name"_genes.list
                                virus=\$(grep -w "\$acc" "\$headers" | awk -F "[" '{ print \$2 }' | awk -F "]" '{ print \$1 }'| sed 's/ /_/g')
                                echo "\$virus" | sed 's/_/ /g' >> "\$name"_virus.list
                                echo ">pOTUaa\${j}_"\$virus"_"\$gene"" >> new_"\$name"_asvnames.txt
                                j=\$((\$j+1))
                                echo "\$s done."
                            else
                                echo "Ugh, there was no hit for \$s .."
                                echo "We still love \$s though and we will add it to the final fasta file"
                                echo "\$s" >> otu.list
                                echo "NO_HIT" >>access.list
                                echo "NO_HIT" >>"\$name"_genes.list
                                echo "NO_HIT" >>"\$name"_virus.list
                                echo "NO_HIT" >>evalue.list
                                echo "NO_HIT" >>bit.list
                                echo "NO_HIT" >>pid.list
                                echo "NO_HIT" >>length.list
                                virus="NO"
                                gene="HIT"
                                echo ">pOTUaa\${j}_"\$virus"_"\$gene"" >> new_"\$name"_asvnames.txt
                                j=\$((\$j+1))
                                echo "\$s done."
                        fi
                        else
                            echo "Using RVDB headers."
                            if [[ "\$(grep -wc "\$s" "\$name"_dmd.out)" -eq 1 ]];then
                                echo "Yep, there was a hit for \$s"
                                echo "Extracting the information now:"
                                acc=\$(grep -w "\$s" "\$name"_dmd.out | awk '{print \$3}' | awk -F "|" '{print \$3}')
                                echo "\$s" >>otu.list
                                echo "\$acc" >>access.list
                                line="\$(grep -w "\$s" "\$name"_dmd.out)"
                                echo "\$line" | awk '{print \$10}' >>evalue.list
                                echo "\$line" | awk '{print \$11}' >>bit.list
                                echo "\$line" | awk '{print \$12}' >>pid.list
                                echo "\$line" | awk '{print \$2}' >>length.list
                                echo "Extracting virus and gene ID for \$s now"
                                gene=\$(grep -w "\$acc" "\$headers" | awk -F "|" '{ print \$6 }' | awk -F "[" '{ print \$1 }' | sed 's/ /_/g') &&
                                echo "\$gene" | sed 's/_/ /g' >>"\$name"_genes.list
                                virus=\$(grep -w "\$acc" "\$headers" | awk -F "|" '{ print \$6 }' | awk -F "[" '{ print \$2 }' | awk -F "]" '{print \$1}' | sed 's/ /_/g') &&
                                echo "\$virus" | sed 's/_/ /g' >>"\$name"_virus.list
                                echo ">pOTUaa\${j}_"\$virus"_"\$gene"" >>new_"\$name"_asvnames.txt
                                j=\$((\$j+1))
                                echo "\$s done."
                            else
                                echo "Ugh, there was no hit for \$s .."
                                echo "We still love \$s though and we will add it to the final fasta file"
                                echo "\$s" >>otu.list
                                echo "NO_HIT" >>access.list
                                echo "NO_HIT" >>"\$name"_genes.list
                                echo "NO_HIT" >>"\$name"_virus.list
                                echo "NO_HIT" >>evalue.list
                                echo "NO_HIT" >>bit.list
                                echo "NO_HIT" >>pid.list
                                echo "NO_HIT" >>length.list
                                virus="NO"
                                gene="HIT"
                                echo ">pOTUaa\${j}_"\$virus"_"\$gene"" >>new_"\$name"_asvnames.txt
                                j=\$((\$j+1))
                                echo "\$s done."
                            fi
                        fi
                echo "Done with \$s"
                done
                echo "Now editing "\$name" fasta headers"
                ###### rename_seq.py
                rename_seq.py ${reads} new_"\$name"_asvnames.txt "\$name"_fin.fasta
                awk 'BEGIN {RS=">";FS="\\n";OFS=""} NR>1 {print ">"\$1; \$1=""; print}' "\$name"_fin.fasta > "\$name"_tmpssasv.fasta
                echo "[Sequence header]" > newnames.list
                cat new_"\$name"_asvnames.txt >> newnames.list
                touch sequence.list
                awk 'BEGIN{RS=">";ORS=""}{print \$2"\n"}' \$name"_tmpssasv.fasta >>sequence.list
                rm "\$name"_tmpssasv.fasta
                paste -d "," sequence.list "\$name"_virus.list "\$name"_genes.list otu.list newnames.list length.list bit.list evalue.list pid.list access.list >> "\$name"_sumPHYtab.csv
                paste -d"\t" otu.list access.list "\$name"_virus.list "\$name"_genes.list sequence.list length.list bit.list evalue.list pid.list >> "\$name"_summarytable.tsv
                rm evalue.list ; rm sequence.list ; rm bit.list ; rm pid.list ; rm length.list ;
                rm new_"\$name"_asvnames.txt
                rm "\$name"_virus.list
                rm "\$name"_genes.list
                rm newnames.list
                rm access.list
                echo "Done Assigning Taxonomy To : ${reads} "
                """
            }

        process pOTUaamafft_prot_VAP2 {

            label 'norm_cpus'

            publishDir "${params.mypwd}/${params.outdir}/Analyses/protein/alignment", mode: "copy", overwrite: true

            input:
                file(prot) from pOTUaaformafft

            output:
                file("*_aln.fasta") into ( pOTUaamodeltest, pOTUaatrax_align )
                file("${name}_aln.html") into pOTUaa_aln

            script:
                """
                name=\$( echo ${prot}  | awk -F ".fasta" '{print \$1}' )
                mafft --maxiterate 5000 --auto ${prot} >\${name}_ALN.fasta
                trimal -in \${name}_ALN.fasta -out \${name}_aln.fasta -keepheader -fasta -automated1 -htmlout \${name}_aln.html
                """
        }

        process modeltest_prot_VAP {

            label 'norm_cpus'

            publishDir "${params.mypwd}/${params.outdir}/Analyses/protein/modeltest", mode: "copy", overwrite: true, pattern: '*.{tree,log}'
            publishDir "${params.mypwd}/${params.outdir}/Analyses/protein/modeltest/summaryfiles", mode: "copy", overwrite: true, pattern: '*model.summary'

            input:
                file(prot) from pOTUaamodeltest

            output:
                file("*.tree") into pOTUtree_rax
                file("*.log") into pOTUlog_rax
                file("*model.summary") into pOTUaamodelsum

            script:
                """
                name=\$( echo ${prot} | awk -F "_aln" '{print \$1}')
                modeltest-ng -i ${prot} -p ${task.cpus} -d aa -s 203 --disable-checkpoint
                echo "Modeltest analysis complete :)"
                echo "Modeltest-ng results summary:" > \${name}_model.summary
                tail -7 \${name}_aln.fasta.log >> \${name}_model.summary
                """
        }

        process raxml_pOTUaa_VAP2 {

            label 'rax_cpus'

            publishDir "${params.mypwd}/${params.outdir}/Analyses/protein/trees", mode: "copy", overwrite: true

            input:
                file(palign) from pOTUaatrax_align
                file(ptree) from pOTUaatree_rax
                file(plog) from pOTUaalog_rax

            output:
                file("*raxml*") into raxml_pOTUaa_summary

            script:
                if (params.ptraxcust) {
                    """
                    pre=\$(echo ${palign} | awk -F "_aln.fasta" '{print \$1}' )
                    raxml-ng --all --msa ${palign} --prefix \${pre} --threads ${task.cpus} ${params.ptraxcust}
                    """
                } else if (params.ptmodeltrax) {
                    """
                    pre=\$(echo ${palign} | awk -F "_aln.fasta" '{print \$1}' )
                    mod=\$(tail -14 ${plog} | head -1 | awk -F "--model " '{print \$2}')
                    raxml-ng --all --msa ${palign} --model \${mod} --seed 2 --blopt nr_safe --threads ${task.cpus} --tree ${ptree} --bs-trees autoMRE --prefix \${pre} --redo
                    """
                } else {
                    """
                    pre=\$(echo ${palign} | awk -F "_aln.fasta" '{print \$1}' )
                    raxml-ng --all --msa ${palign} --model protGTR --seed 2 --blopt nr_safe --threads ${task.cpus} --tree ${ptree} --bs-trees autoMRE --prefix \${pre} --redo
                    """
                }
        }

        process prot_counts_VAP {

            label 'norm_cpus'

            publishDir "${params.mypwd}/${params.outdir}/Analyses/protein/counts", mode: "copy", overwrite: true

            input:
                file(fasta) from pOTUaaforcounts
                file(merged) from mergeforpOTUaacounts

            output:
                tuple file("*_protcounts.csv"), file("*dmd.out") into potuaacounts_summary

            script:
                """
                potu="\$( echo ${fasta} | awk -F "_" 'print \$3')"
                diamond makedb --in ${fasta} --db ${fasta}
                diamond blastx -q ${merged} -d ${fasta} -p ${task.cpus} --min-score ${params.ProtCountsBit} --id ${params.ProtCountID} -l ${params.ProtCountsLength} --more-sensitive -o ${params.projtag}_\${potu}_Counts_dmd.out -f 6 qseqid qlen sseqid qstart qend qseq sseq length qframe evalue bitscore pident btop --max-target-seqs 1 --max-hsps 1
                echo "[Sequence]" >tmp.col1.txt
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
               paste -d "," tmp.col1.txt *col.txt > ${params.projtag}_\${potu}_counts.csv
               rm tmp*
               rm *col.txt
                """
        }
    }
} else {
    println("\n\t\033[0;31mMandatory argument not specified. For more info use `nextflow run vAMPirus.nf --help`\n\033[0m")
    exit 0
}

workflow.onComplete {
    log.info ( workflow.success ? \
        "---------------------------------------------------------------------------------" \
        + "\n\033[0;32mDone! Open the following report in your browser --> ${params.mypwd}/${params.outdir}/${params.tracedir}/vampirus_report.html\033[0m" : \
        "---------------------------------------------------------------------------------" \
        + "\n\033[0;31mSomething went wrong. Check error message below and/or log files.\033[0m" )
}
