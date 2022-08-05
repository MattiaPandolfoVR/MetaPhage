# MetaPhage tutorial (v2 beta)

## Requirements

* An Intel 64-bit Linux machine with:
  * Internet access
  * Conda (see [how to install](https://telatin.github.io/microbiome-bioinformatics/Install-Miniconda/))
  * Nextflow (see [how to install](https://www.nextflow.io/docs/latest/getstarted.html#installation))
  * Compiling tools (g++ and make), in Ubuntu `sudo apt install -y build-essential`
* Resources will vary depending on the project size, we recommend at least 128 Gb RAM and 64 CPU cores

## Installation

1. Ensure to have conda and nextflow available, install mamba

```bash
# Check required tools
nextflow -version   # Otherwise install Nextflow: https://www.nextflow.io/
conda --version     # Otherwise install Miniconda: https://www.nextflow.io/
g++ --version       # Otherwise install the compiler, "sudo apt install build-essential" from ubuntu

# Install mamba
conda install -y -c conda-forge mamba
```

2. Clone this repository (*dev* branch):

```bash
# Download the repository
git clone --branch dev https://github.com/MattiaPandolfoVR/MetaPhage.git

# IMPORTANT: Set the installation directory now or the following stesp won't work
export METAPHAGE_DIR="$PWD"/MetaPhage
```

3. Create an environment for MetaPhage2:

```bash
# Create the environment (needed once)
mamba env create -n metaphage2 --file "$METAPHAGE_DIR"/deps/env-v2.yaml

# ⚠️ Activate the environment to use MetaPhage
conda activate metaphage2

# Update/download krona and quast datasets
ktUpdateTaxonomy.sh
quast-download-silva
quast-download-busco
```

4. Clone cd-hit and compile it to support more sequences:

```bash 
# Clone and compile
git clone https://github.com/weizhongli/cdhit.git
cd cdhit
make MAX_SEQ=2000000 

# Move binaries to the directory with MetaPhage binaries
mv cd-hit cd-hit-2d  cd-hit-454 cd-hit-div cd-hit-est "$METAPHAGE_DIR"/bin/

cd ..
rm -rf cdhit
```

5. Download databases and setup VirSorter2

```bash
# Download up to 6 files simultaneously to "$METAPHAGE_DIR"/DB/ (default location)
# It's important to specify the database release 2022.1 as by default you will get the bundle for v1
"$METAPHAGE_DIR"/bin/python/db_manager.py -o "$METAPHAGE_DIR"/DB/ -m 6 -r 2022.1

# Download and setup virsorter 2
virsorter setup --db-dir "$METAPHAGE_DIR"/DB/virsorter/virsorter2 --jobs 4
```

## Usage

At this point you can use MetaPhage v2 in a similar fashion thant MetaPhage v1:
* Generate a Project file starting with a directory with your reads and a metadata file ([Generate a new project](https://mattiapandolfovr.github.io/MetaPhage/new))
* then execute the pipeline ensuring the conda environment is active (see [run the pipeline](https://mattiapandolfovr.github.io/MetaPhage/tutorial#create-the-project-configuration-file))

