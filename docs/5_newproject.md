---
sort: 5
permalink: /new
---

# Start new project

MetaPhage comes with an utility script to initialize a new project.
The utility, **newProject.py**, is located in the _bin_ subdirectory.

*newProject.py* creates a project configuration file to be fed to
nextflow. The configuration file acts as an analysis protocol, and makes
it easier to generate reproducible results.

## Help screen

```text
usage: newProject.py [-h] [-s SAVE] -i READSDIR [-o OUTPUT_DIR]
                     [-d DATABASE_DIR] [-m METADATA] [-p PROJECT]
                     [-v MAIN_VARIABLE] [-a ALPHA_DIV_1] [-A ALPHA_DIV_2]
                     [-b BETA_DIV] [-V VIOLIN] [-S SUM_VIOL_VAR] [-H HEATMAP]
                     [--single-end] [-l] [--img IMG] [--tmp TEMPDIR]
                     [--work WORKINGDIR] [--verbose]

Generate a Metaphage configuration file for a new project.

optional arguments:
  -h, --help            show this help message and exit

Main arguments:
  -s SAVE, --save SAVE  Configuration file output [default: None]
  -i READSDIR, --reads-dir READSDIR
                        Directory containing the reads.
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Output directory [default: ./MetaPhage]
  -d DATABASE_DIR, --database-dir DATABASE_DIR
                        Database directory
  -m METADATA, --metadata-file METADATA
                        Metadata file.
  -p PROJECT, --project PROJECT
                        Project name

Metadata arguments:
  -v MAIN_VARIABLE, --main-variable MAIN_VARIABLE
                        Default variable of the metadata table for comparisons
  -a ALPHA_DIV_1, --alpha-div-1 ALPHA_DIV_1
                        Variable for alpha diversity (otherwise -v)
  -A ALPHA_DIV_2, --alpha-div-2 ALPHA_DIV_2
                        Secondary variable for alpha diversity (otherwise -v)
  -b BETA_DIV, --beta-div BETA_DIV
                        Variable for alpha diversity (otherwise -v)
  -V VIOLIN, --violin VIOLIN
                        Variable for violin plots (otherwise -v)
  -S SUM_VIOL_VAR, --sum-viol-var SUM_VIOL_VAR
                        Variable for total violin plots (otherwise -v)
  -H HEATMAP, --heatmap HEATMAP
                        Variable for heatmap (otherwise -v)

Metadata arguments:
  --single-end          Single end reads (by default is inferred)
  -l, --local-run       Configure for local execution
  --img IMG             Singularity image [default: None]
  --tmp TEMPDIR         Temporary directory [default: /tmp]
  --work WORKINGDIR     Nextflow work directory [default: /tmp]
  --verbose             Enable verbose output
```

## Main parameters

### Input files

* **`-i`**, **`--reads-dir`** _DIRECTORY_: path to the input reads. 
* **`-m`**, **`--metadata-file`** _FILE_: CSV file with the [metadata]({{ site.baseurl }}{% link notes/metadata.md %}). Will contain the columns used for diversity analyses.

### Other paths
* **`-d`**, **`--database-dir`** _DIRECTORY_: path to the downloaded [databases]({{ site.baseurl }}{% link notes/databases.md %}). By default will select the `./db` subdirectory of the MetaPhage installation directory.
* **`-o`**, **`--output-dir`** _DIRECTORY_: path to the output directory (default: `./MetaPhage`).
* **`-s`**, **`--save`** _FILE_: configuration file created by the script (default: stdout)
* **`--tmp`** _TEMPDIR_: Temporary directory [default: /tmp]
* **`--work`** _WORKINGDIR_: Nextflow work directory [default: /tmp]

### Metadata variables

To produce plots and diversity analyses we can specify the variables (column headers in the metadata file). To perform a primary analysis, it is possible to simply specify the main variable that will 
be used for all the plots:

* **`-v`**, **`--main-variable`** _NAME_: metadata variable to be used in all the plots. By default will be _Treatment_ but an error will be thrown if no "Treatment" column is provided.

Other variables:

* **`-a`**, **`--alpha-div-1`**, _ALPHA_DIV_1_: Variable for alpha diversity (otherwise -v)
* **`-A`**, **`--alpha-div-2`**, _ALPHA_DIV_2_: Secondary variable for alpha diversity (otherwise -v)
* **`-b`**, **`--beta-div`**, _BETA_DIV_:  Variable for alpha diversity (otherwise -v)
* **`-V`**, **`--violin`**, _VIOLIN_: Variable for violin plots (otherwise -v)
* **`-S`**, **`--sum-viol-var`**, _SUM_VIOL_VAR_: Variable for total violin plots (otherwise -v)
* **`-H`**, **`--heatmap`**, _HEATMAP_: Variable for heatmap (otherwise -v)


## Assumptions

The program will scan the reads directory to check the number of samples and if they are
single end or paired ends (based on the presence of "_1"/"_R1" and "_2"/"_R2" tags).
The list of samples is then compared with the "Sample" colum in the metadata file.

The "Metadata" variables are checked and if not present in the metadata file will throw and error
listing the detected columns to simplify the troubleshooting.

If it is planned to run the pipeline in a single node / virtual machine, the script can be used 
to launch it with `-l` (`--local-run`), and in this case the total RAM and cores available will be
also written in the configuration file.

## Example output

```text
params {    
    config_profile_name = 'MetaPhage project'    
    config_profile_description = 'MetaPhage analysis configuration'

    // INPUT PATHS    
    readPath = "/home/ubuntu/volume/metaphage-test/input-change"    
    fqpattern = "_{1,2}.fastq.gz"    
    metaPath = "/home/ubuntu/volume/metaphage-test/input-change/metadata"    
    dbPath = "/qib/platforms/Informatics/telatin/git/MetaPhage/db"

    // OUTPUT/WORKING PATHS    
    outdir = "/home/ubuntu/volume/metaphage-test/MetaPhage"    
    temp_dir = "/tmp"

    // METADATA     
    metadata = true    
    virome_dataset = true    
    singleEnd = false    
    sum_viol_var = "Infant_delivery_type"     
    heatmap_var = "Infant_delivery_type"     
    alpha_var1 = "Infant_delivery_type"     
    alpha_var2 = "Infant_formula_type"     
    beta_var = "Infant_delivery_type"     
    violin_var = "Infant_delivery_type" 
}
```

## How to use it

```bash
nextflow run main.nf -c project.conf [other options]
```

The **project.conf`** file is the output of _newProject.py_, and contains the
information required to locate the files and generate the plots.
It is missing the information on how to run the pipeline (using a scheduler...)
or where to locate the required programs.

If running locally, you can run the pipeline inside the appropriate **conda environment**
and this will remove the need to specify otherwise.
Alternatively, you can use Singularity or Docker. If you pre-downloaded the singularity
image (for example for off-line execution) you can add `-with-singularity PATH_TO_IMAGE`.

Similarily, you can pre-download the Docker image and then add 
`-with-docker andreatelatin/metaphage:1.0`.

See also:

* [MetaPhage and Docker]({{ site.baseurl }}{% link notes/docker.md %})
* [MetaPhage and Singularity]({{ site.baseurl }}{% link notes/singularity.md %})
* [MetaPhage Configuration]({{ site.baseurl }}{% link notes/config.md %})