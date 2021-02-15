FROM continuumio/miniconda3

LABEL authors="Alex J. Veglia and Ramon Rivera-Vicens" \
      description="Docker image containing all requirements for vAMPirus pipeline" \
      version="1.0dev"

RUN apt update; apt install -y gcc bc procps

COPY vampirus_env.yml /
RUN conda env create -f /vampirus_env.yml  && conda clean -a

RUN conda create -n virtualribosome -y python=2.7

ENV PATH /opt/conda/envs/vAMPirus/bin:$PATH
