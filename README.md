[![singularity](https://img.shields.io/badge/singularity-%3E%3D%202.4.2-blue.svg)](http://singularity.lbl.gov/)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.30.2-brightgreen.svg)](https://www.nextflow.io/)

**heinzlab/chip-seq-pipeline** is a bioinformatics best-practice analysis pipeline used for ChIP sequencing data in the Heinz Lab at UCSD. The pipeline is adapted from the nf-core set of Nextflow pipelines.

The pipeline uses [Nextflow](https://www.nextflow.io), a bioinformatics workflow tool. It pre-processes raw data from FastQ inputs, aligns the reads and performs extensive quality-control on the results.

This pipeline can be run on any Unix based environment that has Nextflow installed.

## What is here?
* **README.md:** is what you are reading, which has a complete walkthrough of building and running the container.
* **Singularity:** includes the build recipe for the main Singularity container. Can more or less be copied over for implementation of other genomics pipelines.
* **nextflow.config:** the main nextflow config file with profile and default parameter information.

## Installation
### NextFlow installation
See https://github.com/SciLifeLab/NGI-NextflowDocs for instructions on how to install and configure Nextflow. The easiest method involves using Anaconda (https://www.anaconda.com/).

```
conda install -c bioconda nextflow
```

### Singularity installation
Singularity is a container technology similar to Docker that has become popular for use on HPCs to handle multiple software dependencies. Admin assistance might be required for installation if your cluster does not already support Singularity containerization.

See https://singularity.lbl.gov/ for instructions on how to install Singularity.

See https://singularity.lbl.gov/install-request for information on how to request for Singularity to be installed on your university's HPC if not already available.

### Pipeline installation
This pipeline itself needs no installation - NextFlow will automatically fetch it from GitHub when `heinzlab/chip-seq-pipeline` is specified as the pipeline name.

## Running the pipeline
The typical command for running the pipeline is as follows:

```
nextflow run heinzlab/chip-seq-pipeline --reads '*.fastq.gz' --genome hg38
```

Information regarding the mandatory and optional parameters available can be found as follows:

```
nextflow run heinzlab/chip-seq-pipeline --help
```

The following files will be created in the directory the pipeline is run:

```
work            # Directory containing the nextflow working files
results         # Finished results for each sample, one directory per pipeline step
.nextflow.log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

## Reference genomes available
The following genomes are available:

* UCSC hg38 (--genome hg38)
* Ensembl GRCh37 (--genome GRCh37)
* Ensembl GRCm38 (--genome GRCm38)
* Ensembl Rnor_6.0 (--genome Rnor_6.0)

## Setting up and running the pipeline tutorial
Step by step guide to setting up and running the pipeline.

1. Install Anaconda (https://www.anaconda.com/download/).

2. Install Nextflow:

```
conda install --yes -c bioconda nextflow
```

3. Download and unzip references from igenomes (https://support.illumina.com/sequencing/sequencing_software/igenome.html)

4. Copy and paste the `conf/base.config` file at https://github.com/heinzlab/chip-seq-pipeline/blob/master/conf/base.config into a text editor and edit the `igenomes_base` parameter at the bottom to the path that holds your igenome files. Name it whatever you want.

5. Set up cluster environment config file if running on a cluster (see section below).

6. Run the pipeline:

```
nextflow run heinzlab/chip-seq-pipeline --reads '/path/to/reads/*.fastq.gz' --genome your_genome --singleEnd -c your_config_file.txt
```
