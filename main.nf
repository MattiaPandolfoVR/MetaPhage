#!/usr/bin/env nextflow


/* CONFIGURATION VARIABLES */

// quality check 
params.skip_qc = false

// Trimming 
params.adapter_forward = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
params.adapter_reverse = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
params.mean_quality = 15
params.trimming_quality = 15
params.keep_phix = false
params.phix_reference = "./db/phix/GCA_002596845.1_ASM259684v1_genomic.fna.gz" 

// microbial taxonomy 
params.skip_kraken = false 
params.kraken2_db = false                                                                                                                                          

// Assembly 
params.skip_spades = false
params.skip_quast = false

// Phage-hunting - Vibrant 
params.skip_vibrant = false
params.skip_vibrant_db = false  
params.vibrant_db = './db/vibrant'                                                          

// Viral Taxonomy - vContact2 
params.skip_vcontact2 = false                                                                               

  
/* FILE INPUT */

if (params.readPaths) {                 // declared in profile config
    if (params.singleEnd == true) {     // declared in profile config
        error "singleEnd mode not supported yet"
    } else {
        ch_reads_fastqc = Channel.fromFilePairs("${params.readPaths}/*_{R1,R2}.fastq.gz")
    }
}
else
    error "No input files supplied, please check your config file"


/* PROCESSES */

/* STEP 1 - quality check */
process fastqc {
    tag "${reads.baseName}" 
    publishDir "${params.outdir}/fastqc/", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") == -1 ? "$filename" : null}

    when:
    !params.skip_qc

    input:
    tuple val(name), file(reads) from ch_reads_fastqc

    output:
    file("*_fastqc.{zip,html}") into ch_fastqc_results_raw

    script:
    """
    fastqc \
    -t ${task.cpus} \
    -q $reads
    """
}
