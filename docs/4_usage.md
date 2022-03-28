---
sort: 4
permalink: /usage
---

# Usage

To run the pipeline you need:

* A [configuration file]({{ site.baseurl }}{% link notes/config.md %}) specific for your execution environment (optional, but recommended). We will assume it's called _cluster.conf_.
* The [databases]({{ site.baseurl }}{% link notes/databases.md %})
* A project configuration file, created with [newProject.py]({{ site.baseurl }}{% link 5_newproject.md %}), that includes information on the reads and metadata (we will assume the file is called _project.conf_)


```bash
nextflow run MetaPhage/main.nf -c project.conf -c cluster.conf
```
  