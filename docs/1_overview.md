---
sort: 1
permalink: /overview
---

# Overview

## What is it

A reads to report workflow for metavirome analysis:

* Raw reads cleanup and profiling with [Kraken2](https://ccb.jhu.edu/software/kraken2/)
* Assembly with [SPAdes](https://github.com/ablab/spades#readme) or [Megahit](https://github.com/voutcn/megahit#readme) (overview with [MetaQuast](http://quast.sourceforge.net/metaquast))
* Phage mining with:
  * [VirFinder](https://github.com/jessieren/VirFinder#readme)
  * [VirSorter](https://github.com/simroux/VirSorter#readme)
  * [Phigaro](https://github.com/bobeobibo/phigaro#readme)
  * [Vibrant](https://github.com/AnantharamanLab/VIBRANT#readme)
* Consensus (cd-hit) and backmapping ([bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)), and quantification ([bamtocov](https://github.com/telatin/bamtocov))
* Genefinding and viral taxonomy ([vContact2](https://bitbucket.org/MAVERICLab/vcontact2/wiki/Home))
* Comprehensive report

## Workflow manager

The workflow is written in [Nextflow](https://nextflow.io/), a DSL and task
orchestrator that allows the reproducible execution and scale up from:
* local execution (e.g. virtual machine)
* HPC (Slurm, PBS...)
* Cloud (Amazon, Google,...)

The dependencies can be easily installed, via:
* Docker
* Singularity
* Conda environment

## To get started

A complex pipeline requires a lot of dependencies and the appropriate resources 
to parallelize the execution of concurrent tasks.

[Nextflow](https://nextflow.io/) allows the user to have a custom configuration specifying how to execute the
tasks (locally, in a cluster, on the cloud) and where to find the packages (Singularity, Docker).

Finally, MetaPhage needs a set of databases that can be easily downloaded with a script provided 
in the repository.

## The workflow


![MetaPhage Schematics]({{ site.baseurl }}{% link imgs/workflow.svg %})


