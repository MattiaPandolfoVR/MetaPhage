---
sort: 3
---

# Docker

A Docker image with all the dependencies of MetaPhage is available from [Docker hub](https://hub.docker.com/r/andreatelatin/metaphage), and can be pre-downloaded via:

```bash
docker pull andreatelatin/metaphage:1.0
```

## Running with Docker (manually)

Add `-with-docker andreatelatin/metaphage:1.0` to the pipeline to manually specify that
the dependencies will be found in the MetaPhage docker image.

## Running with Docker (configuration)


A simpler way to specify the container is with a configuration file, for example:

```nextflow
docker.enabled     = true
process.container  = 'andreatelatin/metaphage:1.0'
```

## Nextflow documentation

* [Nextflow: Docker containers](https://www.nextflow.io/docs/latest/docker.html)