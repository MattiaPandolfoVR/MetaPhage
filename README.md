This is MetaPhage, a nextflow pipeline for automatic phage discovery. MetaPhage can be run on Linux or MacOS.

# Structure

This pipeline consists of several modules. The image below summarizes them all.

<p align="center">
  <img src="./figures/metaphage.drawio.svg">
</p>


# Installation

## With Conda

1. Conda is a handy package manager. First of all, install Conda following these instructions https://docs.conda.io/projects/conda/en/latest/user-guide/install/.

2. Install mamba
```
conda install mamba -c conda-forge -y
```

3. Download and extract the repo into the previously specified folder:
```
git clone --single-branch --branch pre_release https://github.com/MattiaPandolfoVR/MetaPhage.git

user: MattiaPandolfoVR
pwd: 011811037f4123b6929d11d0ea1d5363cdb1c394
```

4. Using `cd` command, set your working directory where you want to download the pipeline:
```
cd your/path
```

5. Create the MetaPhage core environment using mamba and activate it
```
mamba env create -n metaphage --file metaphage_core_env.yml
conda activate metaphage
```

## With Docker

Not implemented yet.

## With Singularity (recommended)

Not implemented yet.

# Usage

## With Conda

1. First of all, either copy your paired datasets in the `datasets/base` folder or specify the path where they are stored with the `--readPath` option (see below).

2. Open your terminal and activate the environment previously created:
```
conda activate metaphage
```

3. Using `cd` command, set your working directory to the previously downloaded MetaPhage folder:
```
cd /your/path/MetaPhage
```

4. Start the pipeline with the following command, selecting the right combination of configs (one for datasets and one for system resources):
```
nextflow run main.nf -profile dataset,load
```
At this level you can specify all your **custom options** (read the dedicated section), for example the `--readPath` option previously mentioned:
```
nextflow run main.nf -profile base --readPath /your/path/
```

# Options

## General options

### `--readPath`

Specify the folder where your datasets are stored. Default is `$projectDir/dataset`. Note that $projectDir correspond to the MetaPhage folder.

### `--metaPath`

Specify the folder where the dataset's metadata are stored. Default is `$readPath/metadata/`.

### `--dbPath`

Specify the folder where databases are stored. Default is `$projectDir/db`.

### `--virome_dataset`

(Boolean) Specify if the dataset provided is a virome. This will affect different tools usage. Default is `false`.

### `--singleEnd`

(Boolean) Specify if your datasets are in single-end mode. Default is `false`. Please note that single-end mode is not supported yet.

### `--outdir`

Specify the folder where your results are stored. Default is `$projectDir/output`.

### `--temp_dir`

Specify the folder where your temporary files are stored. Default is `$projectDir/temp`.

### `--workDir`

Specify the directory where tasks temporary files are created. Default is `$projectDir/work`.

## Quality check and trimming

### `--skip_qtrimming`

(Boolean) Specify whether to perform the quality trimming of your reads or not. Default is `false`.

### `--adapter_forward` and `--adapter_reverse`

Specify the adapter sequences. Deafault are `AGATCGGAAGAGCACACGTCTGAACTCCAGTCA` for forward and `AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT` for reverse. 

### `--mean_quality`

Given a read, every base having quality less than `--mean_quality` (default is `15`) is marked as "unqualified". If the percentage of unqualified in the read is more that 40%, that that read if excluded from the analysis. See <https://github.com/OpenGene/fastp> for more informations.

### `--trimming_quality`

Sliding window trimming is enabled in 5'→3' and in 3'→5' with a window 4bp large. Bases inside the window are trimmed if their mean quality is less then `--trimming_quality` (default is 15). See <https://github.com/OpenGene/fastp> for more informations.

### `--keep_phix`

(Boolean) Specify whether to remove the phix or not. Default is `false`. 

### `--mod_phix`

Specify the modality of phix removal. There are 3 possibilities:

- `phiX174` (default) search and remove the complete genome of Coliphage phiX174 isolate S1 (GenBank: AF176027.1, <https://www.ncbi.nlm.nih.gov/nuccore/AF176027>). Genome is automatically downloaded if not already present in `./db/phix/`.
- `WA11` search and remove the complete genome of Coliphage WA11 (GenBank: DQ079895.1, <https://www.ncbi.nlm.nih.gov/nuccore/DQ079895>). Genome is automatically downloaded if not already present in `./db/phix/`.
- `custom` search and remove the sequence specified with `--file_phix_alone` (path to the .fasta file; the path is relative to the pipeline's root directory, for example `--file_phix_alone ./db/phix/genome.fasta`).

## Microbial taxonomy

### `--skip_bacterial_taxo`

(Boolean) Specify whether to skip microbial taxonomy classification step (kraken2 and krona) or not. Default is `false`.

### `--skip_kraken2`

(Boolean) Specify whether to skip the microbial taxonomy classification with Kraken2 or not. Default is `false`.

### `--mod_kraken2`

Specify the modality of the short read alignment with Kraken2. There are 3 possibilities:

- `miniBAV` (default) align against RefSeq bacteria, archaea, and viral libraries. Pre-built database taken from <https://ccb.jhu.edu/software/kraken2/downloads.shtml>.
- `miniBAVH` align against RefSeq bacteria, archaea, and viral libraries, and against the GRCh38 human genome. Pre-built database taken from <https://ccb.jhu.edu/software/kraken2/downloads.shtml>.
- `custom` align using your custom database. With this modality you have to specify also `--file_kraken2_db`, which is the path to the folder containing the `hash.k2d`, `opts.k2d` and `taxo.k2d` files. The path is relative to the pipeline's root directory, for example `--file_kraken2_db ./db/kraken2/folder/`.

### `--skip_krona`

(Boolean) Specify whether to generate the krona-compatible TEXT file using KrakenTools/kreport2krona.py or not. Default is `false`.

## Assembly

### `--skip_megahit`

(Boolean) Specify whether to skip the assembly with MEGAHIT or not. Default is `true`.

### `--skip_metaquast`

(Boolean) Specify whether to skip the assembly evaluation with metaQUAST or not. Default is `false`.

## Phage mining

### `--skip_mining`

(Boolean) Specify whether to skip the phage mining entire step (VIBRANT, phigaro, VirSorter, VirFinder) or not. Default is `false`. If you want to exclude a single or multiple tools from this step, use the specific skip parameter instead.

### `--skip_vibrant`

(Boolean) Specify whether to skip the phage mining with VIBRANT or not. Default is `false`.

### `--mod_vibrant`

Specify the modality of the phage mining with VIBRANT. There are 2 possibilities:

- `legacy` (default) use VIBRANT 1.0.1.

- `standard` use VIBRANT 1.2.1. **Not working yet** (does not produce output).

### `--skip_phigaro`

(Boolean) Specify whether to skip the phage mining with Phigaro or not. Default is `false`.

### `--mod_phigaro`

Specify the modality of the phage mining with phigaro. There are X possibilities:

- `standard` (default) use phigaro 2.3.0.

- `custom` use phigaro using your custom config file. With this modality you have to specify also `--file_figaro_config`, which is the path to the .yml config file containing yout custom parameters to run the miner.

### `--skip_virsorter`

(Boolean) Specify whether to skip the phage mining with VirSorter or not. Default is `false`.

### `--mod_virsorter`

Specify the modality of the phage mining with VirSorter. There are 2 possibilities:

- `legacy` (default) use VirSorter 1.0.6.

- `standard` use VirSorter 2.0.beta. **Not working yet** (does not produce output).

-  `custom` use VirSorter with your custom database. With this modality you have to specify also `--file_virsorter_db`, wich is the path to the folder containing the database file. Please verify that your database files match the required files requested by the VirSorter version (1.0.6).

### `--skip_virfinder`

(Boolean) Specify whether to skip the phage mining with VirFinder or not. Default is `false`.

## Dereplication and reads mapping

### `--skip_dereplication`

(Boolean) Specify whether to skip the dereplication of viral scaffolds or not. Default is `false`.

### `--minlen`

Specify the minimal length in bp for a viral _consensus_ scaffold. Default is `1000`.

## Viral taxonomy

### `--skip_viral_taxo`

(Boolean) Specify whether to skip the viral taxonomy classification step (vcontact2, graphanalyzer) or not. Default is `false`.

### `--skip_vcontact2`

(Boolean) Specify whether to skip the phage taxonomy classification with vcontact2. Default is `false`.

### `--mod_vcontact2`

Specify the modality of the phage taxonomy analysis with vConTACT2. Currently there are 2 possibilities (new modalities will be added periodically):

- `Jan2022` (default) use vConTACT2 with the outputs generated by <https://github.com/RyanCook94/inphared.pl> script runned on January 2022. 
- `custom` use vConTACT2 using your custom reference genomes and taxonomy. With this modality you have to specify also `--file_vcontact2_db`, which is the path to the folder containing the `vConTACT2_proteins.faa`, `vConTACT2_gene_to_genome.csv` and `data_excluding_refseq.tsv` files. The path is relative to the pipeline's root directory, for example `--file_vcontact2_db ./db/inphared/custom/`.

### `--vcontact2_file_head`

Specify the INphared files prefix (vcontact2 db). Default is `20Jan2022_vConTACT2_` If you use a different version of inphared, specify the file prefix string (usually `dd/mm/yyyy_vConTACT2_`). 

### `--skip_graphanalyzer`

(Boolean) Specify wheter to skip the automatic phage taxonomy assignment with graphanalyzer and taxonomy table csv file. Default is `false`.

## Plots and report

### `--skip_miner_comparison`

(Boolean) Specify whether to skip the miner comparison plot (upSet plot) or not. Default is `false`.

### `--skip_summary`

(Boolean) Specify whether to skip the summary table generation and single (for each sample) violin plots creation or not. Default is `false`.

### `--skip_taxonomy_table`

(Boolean) Specify whether to skip the taxonomy table generation or not. Default is `false`.

### `--skip_heatmap`

(Boolean) Specify whether to skip the heatmap plot creation or not. Default is `false`.

### `--heatmap_var`

Specify the name of the variable (metadata column) to use for the heatmap top dendrogram subdivision. (Requested) in order to generate the heatmap.

### `--skip_alpha_diversity`

(Boolean) Specify whether to skip the alpha-diversity plots creation or not. Default is `false`.

### `--alpha_var1`

Specify the name of the variable (metadata column) to use for the alpha-diversity sample clustering (on the x-axis). (Requested) in order to generate the alpha-diversity plots

### `--alpha_var2`

Specify the name of the variable (metadata column) to use for the alpha-diversity color mapping. (Requested) in order to generate the alpha-diversity plots.

### `--skip_beta_diversity`

(Boolean) Specify whether to skip the beta-diversity plots creation or not. Default is `false`.

### `--beta_var`

Specify the name of the variable (metadata column) to use for the beta-diversity color mapping. (Requested) in order to generate the beta-diversity plots.

### `--skip_violin_plots`

(Boolean) Specify whether to skip the violin plot (samples clustered for a variable) or not. Default is `false`.

### `--violin_var`

Specify the name of the variable (metadata column) to use for the violin plot clustering. (Requested) in order to generate the violin plot.

### `--skip_report`

(Boolean) Specify whether to skip the report step (Multiqc report) or not. Default is `false`.
