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

os_c() {
    if [ -f /etc/os-release ];then
        echo -e "\n\t -- Downloading Linux Miniconda3 installation -- \n"
        curl -o Miniconda3-latest-Linux-x86_64.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    else
        echo -e "\n\t\e[31m -- ERROR: Are you in a Unix system? Please check requirements for vAMPirus --\e[39m\n"
        exit 0
    fi
}
source_c() {
    if [ -f ~/.bashrc ];then
        source ~/.bashrc
    fi
}
conda_c() {
    source ~/.bashrc
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
                if [ -f vAMPirus_env.yml ];then
                    echo -e "\n\t -- vAMPirus environment file found. Creating environment... --\n"
                    conda env create -f vAMPirus_env.yml
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
                if [ -f vAMPirus_env.yml ];then
                    echo -e "\n\t -- vAMPirus environment file found. Creating environment... --\n"
                    conda env create -f vAMPirus_env.yml
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
    source ~/.bashrc
    cd $mypwd
    echo "Checking if Nextflow installed system-wide.."
    check_nextflow=$( command -v nextflow )
    if [ "$check_nextflow" = "" ];then #&& [ "$ver" -gt "45" ];then
        echo -e "\n\t -- Nextflow is not system-wide, checking if its in the current working directory.. --\n"
        check_nextflow=$( ls nextflow)
    fi
    if [ "$check_nextflow" != "" ];then #&& [ "$ver" -gt "45" ];then
        echo -e "\n\t -- Nextflow seems to be installed in your system and in your current \$PATH, great! --\n"
    else
        echo -e "\n\t -- Nextflow does not to seem to be downloaded or it specified in your \$PATH --\n"
        echo -e "\n\t -- If you know you have it downloaded, answer n/N to the following question, add nextflow to your \$PATH variable and re-run this script to test --\n"
        echo -e -n "\n\t If not, would you like me to install Nextflow for you? (y,n,exit): "
        read ans
        case $ans in
            [yY] | [yY][eE][sS])
            echo -e "\n\t Starting Nextflow installation \n"
            curl -s https://get.nextflow.io | bash
            echo -e "\n\t Nextflow installation finished, execultable in "$mypwd" \n"
        ;;
        [nN] | [nN][oO])
            echo -e "\n\t\e[31m -- ERROR: Download and Install Nextflow. Then re-run the script. --\e[39m\n"
            exit 0
        ;;
        exit)
           echo -e "\n\t -- Exiting -- \n"
           exit 0
        ;;
        *)
            echo -e "\n\n\t\e[31m -- Yes or No answer not specified. Try again --\e[39m\n"
            nextflow_c
        ;;
        esac
    fi
}

echo "Alright, lets check your system for Conda..."
conda_c
echo "Now lets check the status of Nextflow on your system..."
nextflow_c

if [[ $DATABASE -eq 1 ]]
then    mkdir "$mypwd"/DATABASES
        cd "$mypwd"/DATABASES
        echo "Database installation: RVDB version 19.0 (latest as of 2020-06)"
        curl -o U-RVDBv19.0-prot.fasta.bz2  https://rvdb-prot.pasteur.fr/files/U-RVDBv19.0-prot.fasta.bz2
        echo "Database downloaded, make sure you update the config file before running!"
elif [[ $DATABASE -eq 2 ]]
then    mkdir "$mypwd"/DATABASES
        cd "$mypwd"/DATABASES
        echo "Database installation: Viral RefSeq database version 2.0 (latest as of 2020-07)"
        curl -o viral.2.protein.faa.gz https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.protein.faa.gz
        echo "Database downloaded, make sure you update the config file before running!"
elif [[ $DATABASE -eq 3 ]]
then    mkdir "$mypwd"/DATABASES
        cd "$mypwd"/DATABASES
        echo "Database installation: NCBI NR protein database (should be the most up to date at time of running this script)"
        curl -o NCBI_nr_proteindb.faa.gz https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
        echo "Database downloaded, make sure you update the config file before running!"
elif [[ $DATABASE -eq 4 ]]
then    mkdir "$mypwd"/DATABASES
        cd "$mypwd"/DATABASES
        echo "Database installation: We want 'em all! Might take a little while....'"
        curl -o NCBI_nr_proteindb.faa.gz https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
        curl -o viral.2.protein.faa.gz https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.protein.faa.gz
        curl -o U-RVDBv19.0-prot.fasta.bz2 https://rvdb-prot.pasteur.fr/files/U-RVDBv19.0-prot.fasta.bz2
        echo "Databases downloaded, make sure you update the config file with the one you would like to use before running!"
elif [[ $DATABASE != "" ]]
then    echo "Error: Database download signaled but not given a value between 1-4"
        exit 1
fi

echo "Setup script is complete!"
