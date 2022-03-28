---
sort: 1
---

# MetaPhage databases

MetaPhage requires a set of databases:
* Inphared
* Kraken2/minikraken
* Phigaro
* PhiX
* Vibrant
* Virsorter

## Download the databases

To download all the databases, use the **db_manager.py** script, available in the
_bin/python/_ subdirectory in this repository.

```text
usage: db_manager.py [-h] -o OUTDIR [-m MAXTHREADS]

This script is designed to manage all the databases used in MetaPhage. This script is placed in MetaPhage/bin.

optional arguments:
  -h, --help            show this help message and exit

Options:
  -o OUTDIR, --outdir OUTDIR
                        Output directory
  -m MAXTHREADS, --maxthreads MAXTHREADS
                        Maximum number of concurrent downloads
```

:bulb: The default location for the databases is in the `./db` subdirectory
of the repository.

```bash
# $MFDIR is the directory of the MetaPhage cloned repository
pip install wget
python ${MFDIR}/bin/python/db_manager.py -o ${MFDIR}/db/ -m 4
```