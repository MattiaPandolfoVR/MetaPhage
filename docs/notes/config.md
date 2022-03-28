---
sort: 4
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

## Example

```nextflow
```

## Nextflow documentation

* [Nextflow: Configuration](https://www.nextflow.io/docs/latest/config.html)