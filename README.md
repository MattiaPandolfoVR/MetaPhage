This is MetaPhage, a nextflow pipeline for phage discovery. MetaPhage can be run on Linux or MacOS.

# Installation

## With Conda

1. Conda is a handy package manager. First of all, install Conda following these instructions https://docs.conda.io/projects/conda/en/latest/user-guide/install/.

2. Open your terminal, create a new Conda environment named `nf` and activate it:
```
conda create --name nf
conda activate nf
```

3. Install nextflow:
```
conda install nextflow==20.10.0 -c bioconda
```

4. Install **wget** and **unzip**:
```
conda install wget==1.20.1 unzip==6.0 -c anaconda
```

5. Using `cd` command, set your working directory where you want to download the pipeline:
```
cd your/path
```

6. Download and extract the repo into the previously specified folder:
```
wget -O MetaPhage.zip https://github.com/MattiaPandolfoVR/MetaPhage/archive/main.zip && unzip MetaPhage.zip
```

## With Docker

Not implemented yet.

## With Singularity (recommended)

Not implemented yet.

## Dependencies

The execution of the pipeline depends on some dependencies that are automatically resolved by Conda. However, if you prefer, can resolve them manually with these commands:
```
conda install python==3.7.8 -c conda-forge
conda install pandas==1.1.4 -c conda-forge
conda install wget==1.20.1 -c anaconda
conda install graphviz==2.42.3 -c conda-forge
conda install fastp==0.20.1 -c bioconda
conda install htstream==1.0.0 -c bioconda
```

# Usage

## With Conda

1. First of all, either copy your paired datasets in the `datasets/base` folder or specify the path where they are stored with the `--readPaths` option (see below).

2. Open your terminal and activate the environment previously created:
```
conda activate nf
```

3. Using `cd` command, set your working directory to the previously downloaded MetaPhage folder:
```
cd /your/path/MetaPhage
```

4. Start the pipeline with the following command:
```
nextflow run main.nf -profile base
```
At this level you can specify all your **custom options** (read the dedicated section), for example the `--readPaths` option previously mentioned:
```
nextflow run main.nf -profile base --readPath /your/path/
```

## Custom options

### `--readPaths`

Specify the folder where your datasets are stored. Default is `./datasets/base/`.

### `--singleEnd`

Specify if your datasets are in single-end mode. Default is `false`. Please note that single-end mode is not supported yet.

### `--adapter_forward` and `--adapter_reverse`

Specify the adapter sequences. Deafault are `AGATCGGAAGAGCACACGTCTGAACTCCAGTCA` for forward and `AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT` for reverse. 

### `--mean_quality`

Given a read, every base having quality less than `--mean_quality` (default is `15`) is marked as "unqualified". If the percentage of unqualified in the read is more that 40%, that that read if excluded from the analysis. See <https://github.com/OpenGene/fastp> for more informations.

### `--trimming_quality`

Sliding window trimming is enabled in 5'→3' and in 3'→5' with a window 4bp large. Bases inside the window are trimmed if their mean quality is less then `--trimming_quality` (default is 15). See <https://github.com/OpenGene/fastp> for more informations.

### `--keep_phix`

Specify whether to remove the phix or not. Default is `false`. 

### `--mod_phix`

Specify the modality of phix removal. There are 3 possibilities:

- `phiX174` (default) search and remove the complete genome of Coliphage phiX174 isolate S1 (GenBank: AF176027.1, <https://www.ncbi.nlm.nih.gov/nuccore/AF176027>). Genome is automatically downloaded if not already present in `./db/phix/`.
- `WA11` search and remove the complete genome of Coliphage WA11 (GenBank: DQ079895.1, <https://www.ncbi.nlm.nih.gov/nuccore/DQ079895>). Genome is automatically downloaded if not already present in `./db/phix/`.
- `custom` search and remove the sequence specified with `--file_phix_alone` (path to the .fasta file).

### `--skip_kraken2`

Specify whether to perform the short read alignment with Kraken2 or not. Default is `false`.

### `--mod_kraken2`

Specify the modality of the short read alignment with Kraken2. There are 3 possibilities:

- `miniBAV` (default) align against RefSeq bacteria, archaea, and viral libraries. Pre-built database taken from <https://ccb.jhu.edu/software/kraken2/downloads.shtml>.
- `miniBAVH` align against RefSeq bacteria, archaea, and viral libraries, and against the GRCh38 human genome. Pre-built database taken from <https://ccb.jhu.edu/software/kraken2/downloads.shtml>.
- `custom` align using your custom database. With this modality you have to specify also `--file_kraken2_db`, which is the path to the `hash.k2d` file. `opts.k2d` and `taxo.k2d` files MUST be placed in the same folder of `hash.k2d`.

### `--skip_bracken`

Specify whether to perform the quantification with Bracken or not. Default is `false`.

### `--bracken_read_length`

Specify the read length to be used in Bracken (default is `100`). The databases provided with `--mod_kraken2 miniBAV` and `--mod_kraken2 miniBAVH` include files for read lengths 100, 150, or 200.

### `--bracken_abundance_level`

Specifies the taxonomic rank to analyze (default is `S`). Options are `D`, `P`, `C`, `O`, `F`, `G`, and `S`. Each classification at this specified rank will receive an estimated number of reads belonging to that rank after abundance estimation.

# Structure

This pipeline consists of several modules. The image below summarizes them all.

<p align="center">
  <img src="./slides/pipeline_2020_11_18.drawio.svg">
</p>
