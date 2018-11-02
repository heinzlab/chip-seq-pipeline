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
nextflow run heinzlab/chip-seq-pipeline --reads '*.fastq.gz' --genome hg38 --singleEnd -ansi
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
Specifying long paths every time you run the pipeline is a pain. To make this easier, the pipeline comes configured to understand reference genome keywords which correspond to preconfigured paths, meaning that you can just specify `--genome ID` when running the pipeline.

Note that this genome key can also be specified in a config file if you always use the same genome.

To use this system, add paths to your config file using the following template:

```groovy
params {
  genomes {
    'YOUR-ID' {
      bed12  = '<PATH TO BED FILE>/genes.bed'
      fasta  = '<PATH TO FASTA FILE>/genome.fa'
      gtf    = '<PATH TO GTF FILE>/genes.gtf'
      star    = '<PATH TO STAR INDEX>/STARIndex/'
    }
    'OTHER-GENOME' {
      // [..]
    }
  }
  // Optional - default genome. Ignored if --genome 'OTHER-GENOME' specified on command line
  genome = 'YOUR-ID'
}
```

You can add as many genomes as you like as long as they have unique IDs.

For those in the Heinz lab using the Marvin cluster, you can see the currently available genomes by doing the following:

```
nextflow run heinzlab/chip-seq-pipeline --list_genomes
```

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
nextflow run heinzlab/chip-seq-pipeline --reads '/path/to/reads/*.fastq.gz' --genome your_genome --singleEnd -c your_config_file.txt --ansi
```

## Cluster Environment
By default, pipeline uses the `local` Nextflow executor - in other words, all jobs are run in the login session. If you're using a simple server, this may be fine. If you're using a compute cluster, this is bad as all jobs will run on the head node.

To specify your cluster environment, add the following line to your config file:

```groovy
process {
  executor = 'YOUR_SYSTEM_TYPE'
}
```

Many different cluster types are supported by Nextflow. For more information, please see the [Nextflow documentation](https://www.nextflow.io/docs/latest/executor.html).

Note that you may need to specify cluster options, such as a project or queue. To do so, use the `clusterOptions` config option:

```groovy
process {
  executor = 'SLURM'
  clusterOptions = '-A myproject'
}
```
