---
sort: 2
---

# Configuration file

Nextflow allow to separate the configuration for the pipeline execution
from the workflow logic.

MetaPhage comes with some pre-configured profiles, but it can be convenient
to simply use a configuration file for the setup and then add `-c mycluster.conf` to
the pipeline execution.

## What to configure

* Executor (local, slurm, PBS...)
* Dependencies (Singularity, Docker,...)

## Executor

### Using a local machine

When running the pipeline in a local machine, you can use the local executor,
supplying the details about your available CPUs and memory.
MetaPhage will require at least 16 CPUs and 64 Gb of RAM.

```groovy
params {
    max_cpus = 64
    max_memory = 128.GB
    max_time = 72.h
}
process {
  executor='local'
}
```

### Using slurm

A predefined template is available for slurm, but should be customized with
the available queues in your HPC.

```groovy
process {
  executor='slurm'
  withLabel:big_res {
        cpus = 32
        memory = 128.GB
        time = 72.h
        queue = "your_queue_name"
    }
  withLabel: med_res {
        cpus = 16
        memory = 32.GB
        time = 48.h
        queue = "your_queue_name"
    }
  withLabel: low_res {
        cpus = 8
        memory = 16.GB
        time = 24.h
        queue = "your_queue_name"
    }
}
```

## Dependencies

### Miniconda

If you created a Miniconda environment with the requested dependencies, just
actvitate it before running the workflow.


### Singularity

If you want to use Singularity, you have to
[download the image]({{ site.baseurl }}{% link notes/singularity.md %}) first and then
either add `-with-singularity PATH_TO_IMAGE` to the pipeline execution, or add

```groovy
// Change the value to the path of your singulariity image
process {
    container      = '$projectDir/containers/metaphage.simg'
}
```

### Docker

A Docker image is also available, as [described here]({{ site.baseurl }}{% link notes/docker.md %}).

## Nextflow documentation

* [Nextflow: Configuration](https://www.nextflow.io/docs/latest/config.html)