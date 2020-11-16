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

![bg left:40% 95%](pipeline_2020_11_18.drawio.svg)


:white_check_mark: esclusivamente path relativi
:white_check_mark: pkg conda si installano da soli
:white_check_mark: iniziata la scrittura del manuale
:white_check_mark: **dbs gestiti automaticamente**
:x: ancora nessun supporto Docker

:exclamation: rimozione delle strutture nf-core
<!-- spiegarne tutti i benefici -->