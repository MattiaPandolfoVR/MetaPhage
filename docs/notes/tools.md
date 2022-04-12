---
sort: 7
---

# Tools implemented


## Miners 

**VirSorter** (Roux et al., 2015) was one of the first tools to focus on homology-based bacteriophage identification rather than prophage regions,
detecting viral hallmarks genes homologous to sequences in an archaea and bacteria reference database.
These genes are then used to build probabilistic models exploiting different metrics to measure a confidence score, consisting of PFAM genes,
uncharacterised genes, short genes, strand switching and viral-like genes.

**VirFinder** (Ren et al., 2017) was the first machine learning approach tool, based on k-mer signatures rather than
reference-based viral identification, on which a model is trained using viral and host sequences in equal numbers.
The sequence itself is predicted using a logistic regression model trained on previously known sequences,
which may introduce biases leading to its variable performance when applied to samples collected from different environments. 

More recent tools, such as **VIBRANT** (Kieft et al., 2020), use an alternative, hybrid approach;
employing both methods discussed above, the tool is able to recover a wide array of phages infecting bacteria and archaea, including prophages,
additionally characterizing the metabolic genes and pathways of the identified phages. The prophage detection approach,
implemented by **Phigaro** (Starikova et al., 2020) makes use of a gene predictor to retrieve a list of genes, proteins and coordinates,
together with related metrics. Protein sequences are then annotated using phage-specific profile hidden Markov models HMMs) from a prokaryotic
Virus Orthologous Group (pVOGs).
A smoothing window algorithm then determines phageal genes high density regions and prophage regions and boundaries,
considering the GC content and pVOGs annotations. For non-bioinformatician scientists, the setting up of a pipeline for phage identification
and classification may be quite demanding.
