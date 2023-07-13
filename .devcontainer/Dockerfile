FROM --platform=linux/amd64 mcr.microsoft.com/vscode/devcontainers/python:3.10
RUN python3.10 -m pip install ipykernel 
RUN python3.10 -m pip install --no-warn-conflicts 'colabfold[alphafold-minus-jax] @ git+https://github.com/sokrypton/ColabFold'
RUN python3.10 -m pip install -q  streamlit fair-esm langchain openai Bio torch 


# download MAFFT image

# For easy upgrade in the future. ARG variables only persist during build time.
ARG MAFFT_VER="7.475"

# metadata
LABEL base.image="ubuntu:bionic"
LABEL dockerfile.version="1"
LABEL software="MAFFT"
LABEL software.version="7.475"
LABEL description="Multiple alignment program for amino acid or nucleotide sequences"
LABEL website="https://mafft.cbrc.jp/alignment/software/"
LABEL license="https://mafft.cbrc.jp/alignment/software/license.txt"
LABEL maintainer="Erin Young"
LABEL maintainer.email="eriny@utah.gov"

# install wget, remove apt garbage
RUN apt-get update && apt-get install -y wget curl ncbi-blast+ && \
 rm -rf /var/lib/apt/lists/* && apt-get autoclean

# install mafft and BLAST
RUN wget https://mafft.cbrc.jp/alignment/software/mafft_${MAFFT_VER}-1_amd64.deb && \
    dpkg -i mafft_${MAFFT_VER}-1_amd64.deb && \
    rm mafft_${MAFFT_VER}-1_amd64.deb && \
    mkdir /data

# clone MPNN 
RUN git clone https://github.com/dauparas/ProteinMPNN.git


