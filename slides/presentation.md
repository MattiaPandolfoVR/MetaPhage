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

### Update **12 ottobre**

![bg left:40% 95%](pipeline_2020_10_12.drawio.svg)

:white_check_mark: pipeline minima funzionante

Gira solo al Quadram:
:x: database da scaricare a mano
:x: path assoluti ovunque
:x: pkg conda da installare a mano

---

### Update **18 novembre**

![bg left:40% 80%](pipeline_2020_11_18.drawio.svg)


:white_check_mark: esclusivamente path relativi
:white_check_mark: pkg conda si installano da soli
:white_check_mark: iniziata la scrittura del manuale
:white_check_mark: **dbs gestiti automaticamente**
:x: ancora nessun supporto Docker

:exclamation: rimozione delle strutture nf-core
<!-- spiegarne tutti i benefici -->

---

### **TO-DO**

:x: finire implementazione (virFinder & virSorter)
:x: clusterizzazione -> discutere con Evelyin
:x: matrici conte (seq e prot)
:x: ragionamento output
:x: iniziare scrittura manoscritto
:x: finire implementazione manuale
:x: MultiQC!

---

### Update **dicembre**

![bg left:40% 95%](null)

:cry: perso!

---

### Update **15 gennaio**

![bg left:60% 100%](pipeline_2020_12_01.drawio.svg)

:white_check_mark: implementata la parte di binning
:white_check_mark: implementati tutti i miner disponibili in conda
:x: iniziata la bozza del paper
:x: ancora nessun supporto Docker/Singularity

---

### Update **15 gennaio**

![bg left:60% 100%](dereplic_2020_12_01.drawio.svg)

:white_check_mark: dereplicazione embrionale con CD-HIT (solo Vibrant)
:construction: vConTACT2 deve girare una volta sola, prendendo il risultato della dereplicazione

---

### Update **15 gennaio**

![bg left:70% 100%](images_2020_12_01/screen2.png)

:white_check_mark: conteggio con bowtie2
:white_check_mark: creazione tabella di conte in MultiQC

---

### Update **15 gennaio**

![bg left:70% 100%](images_2020_12_01/screen1.png)

:white_check_mark: creazione barplot di conte in MultiQC

---

### Update **21 January**

![bg left:60% 100%](pipeline_2020_12_21.drawio.svg)

:white_check_mark: dereplication includes Vibrant and Phigaro