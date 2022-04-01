---
sort: 5
permalink: /output
---

# Pipeline output

The pipeline will produce a directory with the following subdirectories:

* **fastp_qc** (with the QC metrics)
* **assembly** (with the contigs produced by Megahit)
* **mining** (with the output of phigaro, vibrant, virfinder, virsorter)
* **prodigal** (with the gene prediction)
* **cd-hit** (with the dereplicated vOTUs)
* **taxonomy** (with the output of Kraken2 and Krona)
* **report** (with the interactive HTML report)

:mag: an [**example output**](https://telatin.github.io/microbiome-bioinformatics/attachments/metaphage/demo/report/index.html)
is available for a preview, although for space reasons
not all attachments are available.

The report not only summarises graphically the output of the pipeline, but also contains hyperlinks
to relevant files:

* The count table
* The taxonomy table
* PhyloSeq object
* PhyloSeq object, with normalized counts (css)
* The multifasta set of viral OTUs
* The multifasta set of proteins predicted in the viral OTUs
* The annotation (GFF) of the viral OTUs
