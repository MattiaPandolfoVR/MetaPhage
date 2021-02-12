---
marp: true
theme: gaia
paginate: true
backgroundImage: url('https://encrypted-tbn0.gstatic.com/images?q=tbn%3AANd9GcTIqXpTkNoyuasJ53q6CWyEssP1dqp3u23pGQ&usqp=CAU')
---

<!-- _class: lead -->
<!-- _footer: https://github.com/MattiaPandolfoVR/MetaPhage -->

# **MetaPhage**

A nextflow pipeline 
for phage discovery :smile:

---

### Update **12 october**

![bg left:40% 95%](pipeline_2020_10_12.drawio.svg)

:white_check_mark: minimal pipeline working

running at Quadram:
:x: databases to be manually downloaded 
:x: absoluth paths everywhere
:x: conda packages to be manually installed

---

### Update **18 november**

![bg left:40% 80%](pipeline_2020_11_18.drawio.svg)


:white_check_mark: only relative paths
:white_check_mark: auto-installed conda packages
:white_check_mark: wiki manual under development
:white_check_mark: **auto-handled dbs**
:x: no Docker/Singularity support

:exclamation: nf-core template removed

---

### **TO-DO on december-january**

:x: complete miners implementation (virFinder & virSorter)
:x: clusterizzation -> discuss with Evelyin!
:x: count matrixes (seq e prot)
:x: output magaement ?
:x: start the paper draft 
:x: complete the wiki
:x: MultiQC implementation!

---

### Update **december**

![bg left:40% 95%](null)

:cry: lost!

---

### Update **15 january**

![bg left:60% 100%](pipeline_2020_12_01.drawio.svg)

:white_check_mark: binning implemented (maxbin2/metabat2/das_tool)
:white_check_mark: all conded miners are implemented
:white_check_mark: paper draft writing started
:x: still no Docker/Singularity support

---

### Update **15 january**

![bg left:60% 100%](dereplic_2020_12_01.drawio.svg)

:white_check_mark: "proto" dereplication with CD-HIT (only Vibrant)
:construction: vConTACT2 has to run only once, with dereplication results!

---

### Update **15 gennaio**

![bg left:70% 100%](images_2020_12_01/screen2.png)

:white_check_mark: bowtie2 vcs (viral consensus scaffolding)
:white_check_mark: MultiQC count matrix creation!

---

### Update **15 gennaio**

![bg left:70% 100%](images_2020_12_01/screen1.png)

:white_check_mark: MultiQC counts barplot creation!

---

### Update **21 January**

![bg left:60% 100%](pipeline_2020_12_21.drawio.svg)

:white_check_mark: dereplication now includes Vibrant, Phigaro, VirSorter1, VirFinder and Marvel

---

### Update **21 January**

![bg left:60% 100%](dereplic_2020_12_21.drawio.svg)

:white_check_mark: rough dereplication with CD-HIT