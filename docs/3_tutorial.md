---
sort: 3
permalink: /tutorial
---

# Tutorial

## Before we start

This tutorial will show how to run the pipeline using either Conda or Singularity as package providers. **Linux** is required.

* Ensure you have Nextflow installed (try `nextflow -version`), or install it [as described here](https://www.nextflow.io/docs/latest/getstarted.html#installation)
* Ensure you have Miniconda or Singularity installed (try `conda info` or `singularity --version`)

## Download MetaPhage

```bash
# Download MetaPhage and expand it
wget https://github.com/MattiaPandolfoVR/MetaPhage/archive/refs/tags/v0.2.0.tar.gz
tax xvf v0.2.0.tar.gz

# Keep note of the installation directory:
# If you move away you can return with `cd "$METAPHAGE_DIR"`
export METAPHAGE_DIR=$(readlink -f MetaPhage-0.2.0)
```

## Download the databases

We will download the databases before starting the pipeline, as
[described here]({{ site.baseurl }}{% link notes/databases.md %}):

```bash
# This will download the database in "$METAPHAGE_DIR/db/"
cd $METEAPHAGE_DIR
./bin/python/db_manager.py -o ./db/ -m 6
```

:bulb: The `-m INT` parameter specifies the number of concurrent downloads.

## Get the raw data

We will analyse the 20 samples from
"*The stepwise assembly of the neonatal virome is modulated by breastfeeding*"
([Liang et al 2020](https://www.nature.com/articles/s41586-020-2192-1)).

A downloader script is provided, which (unless otherwise specified,
use `--help` for details), will download the reads into the
`demo/` subdirectory of the installation directory.

:bulb: The script will download the samples specified in the *demo* subdirectory
(see [list](https://github.com/MattiaPandolfoVR/MetaPhage/blob/main/demo/infant-metadata.csv)).

```bash
# Always from the $METAPHAGE_DIR directory
./bin/getExample.py --verbose -t 8
```

Again, you can specify the number of concurrent downloads via `-t INT`.
The whole process takes approximately 10-15 minutes, depending on the available bandwitdth.

## Dependencies

### Miniconda

If you plan to use miniconda, create an environment as specified in the repository and activate it:

```bash
conda env create -n metaphage --file deps/env.yaml 
conda activate metaphage
```

### Singularity

You can pre-download the singularity image with the following command:

```bash
mkdir -p $METAPHAGE_DIR/containers/
wget -O $METAPHAGE_DIR/containers/metaphage.simg "https://s3.climb.ac.uk/ifrqmra-metaphage/v1.0/metaphage.simg"
```

## Create the project configuration file

We will use the [newProject.py]({{ site.baseurl }}{% link 5_newproject.md %})
script to generate the configuration file,
as we would do in real-life usage.

```bash
python ./bin/newProject.py -i demo \
    -m demo/infant-metadata.csv \
    -v Infant_delivery_type \
    -s demo.conf
```

Where:

* `-i` is the input directory containing the reads
* `-m` is the metadata file in CSV format (see [details]({{ site.baseurl }}{% link notes/metadata.md %}))
* `-v` is the main metadata variable. By default is *Treatment*, and you can run the program without specifying it as you will get an error and a full list of valid variables (column names) to pick the correct one.
* `-d` is the database directory (by default will be ./db in the installation directory, so it can be omitted in our case)


## Run the pipeline

To run the pipeline locally you will need at least 64 Gb of RAM and 16 CPUs,
provided that the appropriate *conda environment is active*:


:bulb: If you have less, see 
[here]({{ site.baseurl }}{% link notes/resources.md %}))
```bash
cd $METAPHAGE_DIR
nextflow run main.nf -c demo.conf
```

## With singularity

```bash
cd $METAPHAGE_DIR
nextflow run main.nf -c demo.conf -with-singularity ./containers/metaphage.simg
```

## With a scheduler

:bulb: If using a scheduler (in a cluster, for example), we recommend using a Singularity
rather than Miniconda.

You can add the following to your *demo.conf* file, to specify a specific scheduler (example, Slurm)
and to drive the choice of the queue (can be a fixed value, or a conditional value depending on the time as in the
example).

:book: See [executor documentation](https://www.nextflow.io/docs/latest/executor.html) from Nextflow.

```text
process {
    executor     = 'slurm'
    queue        = { task.time <= 2.h ? 'qib-short' : task.time <= 48.h ? 'qib-medium' : 'qib-long' }
    clusterOptions = ' --constraint=intel '
}
```
