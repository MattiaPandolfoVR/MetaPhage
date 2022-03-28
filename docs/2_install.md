---
sort: 2
permalink: /installation
---

# Installation

[Nextflow](https://nextflow.io/) allows to execute pipelines locally, in a cluster with a scheduler ([Slurm](https://slurm.schedmd.com/documentation.html), [PBS](https://www.openpbs.org/), ...) or on the cloud (AWS, Azure...).

The dependencies can be installed in a [Conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html), or using a [Docker](https://www.docker.com/) or [Singularity](https://apptainer.org/) container.

## Prerequisites

The pipeline requires:

* A set of dependencies (provided via Conda or Singularity or Docker)
* A Linux system capable of performing *de novo* assemblies (a single local machine or an HPC cluster)
* A set of databases (installable with a script provided in the repository)

## Installation

:bulb: Ensure you have Nextflow installed otherwise
[install it first](https://www.nextflow.io/docs/latest/getstarted.html)

1. Clone the repository (and enter it):
   
```bash
# Download a copy of the repository
git clone https://github.com/MattiaPandolfoVR/MetaPhage
# Enter the repository
cd MetaPhage
```

2. Download the database (the default location can be `db/` inside the repository)

```bash
# Install a required module
pip install wget
# Download the databases in the "db" subdirectory
python bin/python/db_manager.py -o ./db
```

## Dependencies

If using Docker or Singularity, the pipeline is capable of retrieving the appropriate
image automatically. For local executions or tests, it is possible to generate a single
conda environment (see _Miniconda_).

### Miniconda

:bulb: Ensure that Miniconda is [already installed](https://telatin.github.io/microbiome-bioinformatics/Install-Miniconda/) in the system.

```bash
cd deps
conda env create -n MetaPhage --file env.yaml
```

To execute the pipeline, remember to activate the environment first:
```bash
conda activate MetaPhage
```

### Docker

A Docker image is available from `andreatelatin/metaphage:1.0`.
Nextflow can fetch it before running the pipeline but if you want to download
it before:

```bash
docker pull andreatelatin/metaphage:1.0
```

### Singularity

A Docker image is available from `docker://andreatelatin/metaphage:1.0`.
Nextflow can fetch it before running the pipeline but if you want to download
it before:

```bash
wget "https://s3.climb.ac.uk/ifrqmra-metaphage/v1.0/metaphage.simg"
```
