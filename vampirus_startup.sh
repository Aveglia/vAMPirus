#!/usr/bin/env bash

usage () {
echo "

Run this script in the vAMPirus directory to set up Conda and the vAMPirus Conda environment

You can also use this script to download reference databases for taxonomy assignment by useing the -d option described below.

General execution:

vampirus_startup.sh -h -d [1|2|3|4]

    Command line options:

        [ -h ]                       	Print help information

        [ -d 1|2|3|4 ]                  Set this option to create a database directiory within the current working directory and download the following databases for taxonomy assignment:

                                                    1 - Download the proteic version of the Reference Viral DataBase (See the paper for more information on this database: https://f1000research.com/articles/8-530)
                                                    2 - Download only NCBIs Viral protein RefSeq database
                                                    3 - Download only the complete NCBI NR protein database
                                                    4 - Download all three databases

"

}

while getopts "hd:" OPTION; do
     case $OPTION in
         h) usage; exit;;
         d) DATABASE=${OPTARG};;
     esac
done
shift $((OPTIND-1))      # required, to "eat" the options that have been processed

export mypwd="$(pwd)"

sed "s|VAMPDIR|${mypwd}|g" ${mypwd}/vampirus.config > tmp.config
cat tmp.config > vampirus.config
rm tmp.config

os_c() {
    if [ -f /etc/os-release ];then
        echo -e "\n\t -- Downloading Linux Miniconda3 installation for Linux OS-- \n"
        curl -o Miniconda3-latest-Linux-x86_64.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    elif [ `which sw_vers | wc -l` -eq 1 ];then
        echo -e "\n\t -- Downloading Linux Miniconda3 installation for MacOS -- \n"
        curl -o Miniconda3-latest-MacOSX-x86_64.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
    else
        echo "Error in Miniconda installation"
        exit 0
    fi
}
source_c() {
    if [ -f ~/.bashrc ];then
        source ~/.bashrc
    else
        source ~/.bash_profile
    fi
}
conda_c() {
    source_c
    cd $mypwd
    #Check conda and environment
    check_conda=$( command -v conda )
    if [ "$check_conda" != "" ];then #&& [ "$ver" -gt "45" ];then
        echo -e "\n\t -- Conda seems to be installed in your system --\n"
        ver=$( conda -V | awk '{print $2}' | cut -f 1,2 -d "." )
        vern=4.8
        if [ $( echo "$ver >= $vern" | bc -l ) -eq 1 ];then
            echo -e "\n\t -- Conda is installed (v4.8 or higher). Checking environment... --\n"
            #Check environment
            check_env=$( conda info -e | grep -c "vAMPirus" )
	        if [ "$check_env" -eq 0 ];then
                echo -e "\n\t -- vAMPirus environment has not been created. Checking environment file... --\n"
                if [ -f vampirus_env.yml ];then
                    echo -e "\n\t -- vAMPirus environment file found. Creating environment... --\n"
                    conda env create -f vampirus_env.yml
                else
                    echo -e "\n\t\e[31m -- ERROR: vAMPirus environment file not found \(vAMPirus_env.yml\). Make sure you are running this from the vAMPirus program directory. --\e[39m\n"
                    exit 0
                fi
            elif [ "$check_env" -eq 1 ];then
                echo -e "\n\t -- vAMPirus environment is installed and ready to be used --\n"
            fi
        fi
    else
        echo -e "\n\t -- Conda is not intalled."
        echo -e -n "\n\t    Do you want to install Miniconda3? (y,n,exit): "
        read ans
        case $ans in
            [yY] | [yY][eE][sS])
                os_c
                echo -e "\n\t -- Starting Anaconda installation -- \n"
                bash Miniconda3*.sh
                echo -e "\n\t -- Installation done -- \n"
                rm Miniconda3*.sh
                source_c
                if [ -f vampirus_env.yml ];then
                    echo -e "\n\t -- vAMPirus environment file found. Creating environment... --\n"
                    conda env create -f vampirus_env.yml
                else
                    echo -e "\n\t\e[31m -- ERROR: vAMPirus environment file not found (vAMPirus_env.yml). Please check requirements and rerun the pre-check --\e[39m\n"
                    exit 0
                fi
            ;;
            [nN] | [nN][oO])
                echo -e "\n\t\e[31m -- ERROR: Download and Install Miniconda or Anaconda. Then rerun the pre-check  --\e[39m\n"
                exit 0
            ;;
            exit)
	           echo -e "\n\t -- Exiting -- \n"
               exit 0
            ;;
            *)
                echo -e "\n\n\t\e[31m -- Yes or No answer not specified. Try again --\e[39m\n"
	            conda_c
            ;;
        esac
    fi
}

nextflow_c() {
    source_c
    cd $mypwd
    echo "Checking if Nextflow installed system-wide.."
    check_nextflow=$( command -v nextflow )
    if [ "$check_nextflow" = "" ];then #&& [ "$ver" -gt "45" ];then
        echo "Nextflow is not system-wide, checking if its in the current working directory.."
        check_nextflow=$( ls nextflow )
    fi
    if [ "$check_nextflow" != "" ];then #&& [ "$ver" -gt "45" ];then
        echo "Nextflow seems to be installed in your system and in your current \$PATH, great!"
    else
        echo "Nextflow does not to seem to be downloaded or specified in your \$PATH.."
        echo "If you know you have it downloaded, answer n/N to the following question, add nextflow to your \$PATH variable and re-run this script to test. "
        echo "Otherwise, would you like me to install Nextflow for you now? (y,n,exit): "
        read ans
        case $ans in
            [yY] | [yY][eE][sS])
            echo "Awesome,starting Nextflow installation now ..."
            curl -s https://get.nextflow.io | bash
            echo "Nextflow installation finished, execultable in "$mypwd""
        ;;
        [nN] | [nN][oO])
            echo "ERROR: It seems like you answered with something other than y... Download and install Nextflow on your own and rerun re-run the script to make sure path is sepcified correctly."
            echo "If it was a mistake you regret your original answer, just re-run the script and all will be forgiven ^__^."
            exit 0
        ;;
        exit)
           echo "Alright, if you say so, exiting now."
           exit 0
        ;;
        *)
            echo "Yes or No answer not specified. Lets try this again..."
            nextflow_c
        ;;
        esac
    fi
}

echo "Alright, lets check your system for Conda..."
conda_c
echo "Editing path to conda directory in vampirus.config"
environment="$(conda env list  | sed 's/*//g' | grep "vAMPirus" | head -1 | awk '{print $2}')"
sed "s|CONDADIR|${environment}|g" "$mypwd"/vampirus.config > tmp1.config
cat tmp1.config > "$mypwd"/vampirus.config
rm tmp1.config

echo "-------------------------------------------------------------------------------- Conda loop done"

echo "Now lets check the status of Nextflow on your system..."
nextflow_c

echo "-------------------------------------------------------------------------------- nextflow loop done"

if [[ $DATABASE -eq 1 ]]
then    mkdir "$mypwd"/Databases
        cd "$mypwd"/Databases
        dir="$(pwd)"
        echo "Database installation: RVDB version 19.0 (latest as of 2020-06)"
        curl -o U-RVDBv19.0-prot.fasta.bz2  https://rvdb-prot.pasteur.fr/files/U-RVDBv19.0-prot.fasta.bz2
        bunzip2 U-RVDBv19.0-prot.fasta.bz2
        echo "Editing confiration file for you now..."
        sed 's/DATABASENAME/U-RVDBv19.0-prot.fasta/g' "$mypwd"/vampirus.config > tmp1.config
        sed "s|DATABASEDIR|${dir}|g" tmp1.config > tmp2.config
        rm tmp1.config
        cat tmp2.config > "$mypwd"/vampirus.config
        rm tmp2.config
        echo "Database downloaded, and configuration file edited, you still need to bunzip2 this database then just confirm the path was set correctly in the config file."
elif [[ $DATABASE -eq 2 ]]
then    mkdir "$mypwd"/Databases
        cd "$mypwd"/Databases
        dir="$(pwd)"
        echo "Database installation: Viral RefSeq database version 2.0 (latest as of 2020-07)"
        curl -o viral.2.protein.faa.gz https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.protein.faa.gz
        gunzip viral.2.protein.faa.gz
        echo "Editing confiration file for you now..."
        sed 's/DATABASENAME/viral.2.protein.faa/g' "$mypwd"/vampirus.config > tmp1.config
        sed "s|DATABASEDIR|${dir}|g" tmp1.config > tmp2.config
        rm tmp1.config
        cat tmp2.config > "$mypwd"/vampirus.config
        rm tmp2.config
        echo "Database downloaded, and configuration file edited, you still need to gunzip this database then just confirm the path was set correctly in the config file."
elif [[ $DATABASE -eq 3 ]]
then    mkdir "$mypwd"/Databases
        cd "$mypwd"/Databases
        dir="$(pwd)"
        echo "Database installation: NCBI NR protein database (should be the most up to date at time of running this script)"
        curl -o NCBI_nr_proteindb.faa.gz https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
        gunzip NCBI_nr_proteindb.faa.gz
        echo "Editing confiration file for you now..."
        sed 's/DATABASENAME/NCBI_nr_proteindb.faa/g' "$mypwd"/vampirus.config > tmp1.config
        sed "s|DATABASEDIR|${dir}|g" tmp1.config > tmp2.config
        rm tmp1.config
        cat tmp2.config > "$mypwd"/vampirus.config
        rm tmp2.config
        echo "Database downloaded, make sure you update the config file before running!"
elif [[ $DATABASE -eq 4 ]]
then    mkdir "$mypwd"/Databases
        cd "$mypwd"/Databases
        dir="$(pwd)"
        echo "Database installation: We want 'em all! Might take a little while....'"
        curl -o NCBI_nr_proteindb.faa.gz https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
        curl -o viral.2.protein.faa.gz https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.protein.faa.gz
        curl -o U-RVDBv19.0-prot.fasta.bz2 https://rvdb-prot.pasteur.fr/files/U-RVDBv19.0-prot.fasta.bz2
        sed "s|DATABASEDIR|${dir}|g" "$mypwd"/vampirus.config > tmp1.config
        cat tmp1.config > "$mypwd"/vampirus.config
        rm tmp1.config
        echo "Databases downloaded, make sure you update the config file with the one you would like to use and decompress the database before running."
elif [[ $DATABASE != "" ]]
then    echo "Error: Database download signaled but not given a value between 1-4"
        exit 1
fi
echo "-------------------------------------------------------------------------------- Database loop done"

cd "$mypwd"
echo "Ok, everything downloaded. To test installation, run the following commands and check for errors:"
echo "   "
echo "Checking DataCheck mode:"
echo ""$mypwd"/nextflow run  "$mypwd"/vAMPirusv0.1.0.nf -c  "$mypwd"/vampirus.config -profile conda,test --DataCheck"
echo "    "
echo "Next, test the analysis pipeline:"
echo ""$mypwd"/nextflow run  "$mypwd"/vAMPirusv0.1.0.nf -c  "$mypwd"/vampirus.config -profile conda,test --Analyze"
echo "   "
echo "If everything looks good, here are a example lanch commands to submit after testing installation and editing the paths to your data and other parameters for the run in the vampirus.config file."
echo "First, run the DataCheck:"
echo ""$mypwd"/nextflow run  "$mypwd"/vAMPirusv0.1.0.nf -c  "$mypwd"/vampirus.config -with-conda "$environment" --DataCheck"
echo "   "
echo "OR..."
echo ""$mypwd"/nextflow run  "$mypwd"/vAMPirusv0.1.0.nf -c  "$mypwd"/vampirus.config -profile conda --DataCheck"
echo "   "
echo "Then you can run the analysis:"
echo ""$mypwd"/nextflow run  "$mypwd"/vAMPirusv0.1.0.nf -c  "$mypwd"/vampirus.config -with-conda "$environment" --Analyze ..."
echo "   "
echo "OR..."
echo ""$mypwd"/nextflow run  "$mypwd"/vAMPirusv0.1.0.nf -c  "$mypwd"/vampirus.config -profile conda --Analyze ..."

touch STARTUP_HELP.txt
echo "-------------------------------------------------------------------------------------------------------------------------------- TESTING YOUR INSTALLATION" >> STARTUP_HELP.txt
echo "Ok, everything downloaded. To test installation, run the following commands and check for errors:" >> STARTUP_HELP.txt
echo "   " >> STARTUP_HELP.txt
echo "Checking DataCheck mode:" >> STARTUP_HELP.txt
echo ""$mypwd"/nextflow run  "$mypwd"/vAMPirusv0.1.0.nf -c  "$mypwd"/vampirus.config -profile conda,test --DataCheck" >> STARTUP_HELP.txt
echo "    " >> STARTUP_HELP.txt
echo "Next, test the analysis pipeline:" >> STARTUP_HELP.txt
echo ""$mypwd"/nextflow run  "$mypwd"/vAMPirusv0.1.0.nf -c  "$mypwd"/vampirus.config -profile conda,test --Analyze" >> STARTUP_HELP.txt
echo "   " >> STARTUP_HELP.txt
echo "--------------------------------------------------------------------------------------------------------------------------------" >> STARTUP_HELP.txt
echo "   " >> STARTUP_HELP.txt
echo "Ok, if everything went well (green text was spit out by Nextflow), now you can move on to the fun. First, you should review the help docs and the vampirus.config in the vAMPirus directory." >> STARTUP_HELP.txt
echo "   " >> STARTUP_HELP.txt
echo "-------------------------------------------------------------------------------------------------------------------------------- RUNNING DataCheck PIPELINE WITH YOUR DATA" >> STARTUP_HELP.txt
echo "If everything looks good, here are a example lanch commands to submit after testing installation and editing the paths to your data and other parameters for the run in the vampirus.config file:"
echo "   " >> STARTUP_HELP.txt
echo "First, run the DataCheck part of the pipeline using the -with-conda Nextflow option:" >> STARTUP_HELP.txt
echo "   " >> STARTUP_HELP.txt
echo ""$mypwd"/nextflow run  "$mypwd"/vAMPirusv0.1.0.nf -c  "$mypwd"/vampirus.config -with-conda "$environment" --DataCheck" >> STARTUP_HELP.txt
echo "   " >> STARTUP_HELP.txt
echo "OR using -profile option of Nextflow ..." >> STARTUP_HELP.txt
echo ""$mypwd"/nextflow run  "$mypwd"/vAMPirusv0.1.0.nf -c  "$mypwd"/vampirus.config -profile conda --DataCheck" >> STARTUP_HELP.txt
echo "   " >> STARTUP_HELP.txt
echo "--------------------------------------------------------------------------------------------------------------------------------" >> STARTUP_HELP.txt
echo "   " >> STARTUP_HELP.txt
echo "-------------------------------------------------------------------------------------------------------------------------------- RUNNING Analyze PIPELINE WITH YOUR DATA" >> STARTUP_HELP.txt
echo "Then you can run the analysis using the -with-conda Nextflow option, here is a launch command to run the complete analysis:" >> STARTUP_HELP.txt
echo "   " >> STARTUP_HELP.txt
echo ""$mypwd"/nextflow run  "$mypwd"/vAMPirusv0.1.0.nf -c  "$mypwd"/vampirus.config -with-conda "$environment" --Analyze --ncASV --pcASV" >> STARTUP_HELP.txt
echo "   " >> STARTUP_HELP.txt
echo "OR same command using -profile option of Nextflow ..." >> STARTUP_HELP.txt
echo "   " >> STARTUP_HELP.txt
echo ""$mypwd"/nextflow run  "$mypwd"/vAMPirusv0.1.0.nf -c  "$mypwd"/vampirus.config -profile conda --Analyze --ncASV --pcASV" >> STARTUP_HELP.txt
echo "   " >> STARTUP_HELP.txt
echo "--------------------------------------------------------------------------------------------------------------------------------" >> STARTUP_HELP.txt
echo "    "
echo "Setup script is complete!"
echo "Check out the STARTUP_HELP.txt file for more information on how to move forward with the analysis."
