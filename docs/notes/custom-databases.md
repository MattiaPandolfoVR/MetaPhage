---
sort: 9
---

# Custom Databases

MetaPhage requires a set of databases:
* PhiX
* Kraken2/minikraken
* Phigaro
* Vibrant
* Virsorter2
* CheckV
* AllVSAll
* Inphared

## Use a custom database
By default, **db_manager.py** will download all the databases the first time MetaPhage is launched.
Each database-requiring tool has its own folder in db folder:
db/
    checkv/
            checkv-db-v1.2/
    diamond/
            June2022/
    inphared/
            1Jun2022_data_excluding_refseq.tsv
            1Jun2022_vConTACT2_proteins.faa
            1Jun2022_vConTACT2_gene_to_genome.cs
    kraken2/
            miniBAV/
    phigaro/
            standard/
    phix/
            phiX174.fasta
    vibrant/
            legacy/
    virsorter/
            virsorter2/

If you want to use a different database for any of the tool:
1. Move your custom database folder in the relative tool db/ folder.
2. Modify your specific project.conf located in conf/datasets (by default), setting the correct tool database variable with the name of your custon database folder. Tool database variables with their default values are reported below:
mod_phix = "phiX174"
mod_kraken2 = "miniBAV"
mod_phigaro = "standard"
mod_vibrant = "legacy"
mod_virsorter2 = "virsorter2"
mod_checkv = "checkv-db-v1.2"
mod_vcontact2 = "Jun2022"
vcontact2_file_head = "1Jun2022_vConTACT2_"

For example, if you want to change CheckV default db (checkv-db-v1.2) with your custom one (checkv-custom), move the checkv-custom folder to db/checkv/, and set mod_checkv = "checkv-custom" in your config file.

## Phix and vcontact2
To change the phix genome to use, simply copy the .fasta file in db/phix folder.
Similarly, copy your custom vConTACT2 files in db/inphared.
By default, vConTACT2 make use of two INPHARED files: *_vConTACT2_proteins.faa and *_vConTACT2_gene_to_genome.csv, while graphanalyzer make use of *_data_excluding_refseq.tsv.
If you want to use a custom database for vConTACT2 of graphanalyzer, in addition to copy the files in the db/inphared folder, be sure to match the prefix of your custom files with the parameter vcontact2_file_head (e.g. vcontact2_file_head = "my_custom_" in order to use my_custom_proteins.faa etc). The prefix must be the same for all the file in the folder.