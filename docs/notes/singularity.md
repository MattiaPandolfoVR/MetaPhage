---
sort: 4
---

# Singularity

[Singularity](https://apptainer.org/) is the most widely used container system for HPC.
The configuration file of the pipeline allows the automatic download of the image
from [Docker hub](https://hub.docker.com/r/andreatelatin/metaphage).

## Offline use

Several HPC are kept off-line, so it's convenient to pre-download the Singularity image.

```bash
wget "https://s3.climb.ac.uk/ifrqmra-metaphage/v1.0/metaphage.simg"
```

## Running with Singularity (manually)

Nextflow can be run using a Singularity image adding
`-with-singularity $PATH_TO_IMG` to the command,
for example:

```bash
nextflow run MetaPhage/main.nf -c project.config -with-singularity $PATH_TO_SIMG
```

## Running with Singularity (configuration)

A simpler way to specify the container is with a configuration file, for example:

```nextflow
singularity.enabled    = true
process.container      = 'metaphage.simg'
singularity.autoMounts = true
```

:bulb: If the image is saved in the directory specified in `$NXF_SINGULARITY_CACHEDIR` it's not
necessary to add the path to the container, but just the container filename.

## Nextflow documentation

* [Nextflow: Singularity containers](https://www.nextflow.io/docs/latest/singularity.html)