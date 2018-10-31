Bootstrap: docker
From: continuumio/miniconda3

%labels
    MAINTAINER Carlos Guzman <cag104@eng.ucsd.edu>
    DESCRIPTION Container image containing all requirements for the adapted heinzlab/chip-seq-pipeline
    VERSION 0.1dev

%files
    environment.yml /

%environment
	PATH=/opt/conda/envs/chipseq-0.1dev/bin:$PATH
	export PATH

%post
    apt-get -y update
    apt-get -y install build-essential libboost-all-dev libgsl-dev libz-dev

    /opt/conda/bin/conda env create -f /environment.yml
    /opt/conda/bin/conda clean -a

    git clone https://github.com/rnakato/SSP.git
    cd SSP
    make
    mv bin/ssp /opt/conda/envs/chipseq-0.1dev/bin
