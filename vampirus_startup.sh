#!/usr/bin/env bash

usage () {
echo "

Run this script in the vAMPirus

General exicution:



"

}

while getopts "hrbp:laguf" OPTION; do
     case $OPTION in
         h) usage; exit;;
         r)
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
                    echo -e "\n\t\e[31m -- ERROR: vAMPirus environment file not found \(vAMPirus_env.yml\). Please check requirements and rerun the pre-check --\e[39m\n"
                    exit 0
                fi
            elif [ "$check_env" -eq 1 ];then
                echo -e "\n\t -- vAMPirus environment is installed and ready to be used --\n"
            fi
        fi
    else
        echo -e "\n\t -- Conda is not intalled. Please install Anaconda (https://www.anaconda.com) and rerun this script --\n"
        echo -e -n "\n\t    Do you want to install Anaconda? (y,n,exit): "
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
                echo -e "\n\t\e[31m -- ERROR: Download and Install Anaconda. Then rerun the pre-check  --\e[39m\n"
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
