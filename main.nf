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
        Channel.fromFilePairs("${params.readPaths}/*_{R1,R2}.fastq.gz")
            .into {ch_reads_fastp}
    }
}
else
    error "No input files supplied, please check your config file"


/* PROCESSES */

/* STEP 1 - quality check  and trimming */
process fastp {
    tag "$seqID"
    publishDir "${params.outdir}/fastp/qc/", mode: 'copy',
        saveAs: {filename -> filename.endsWith(".html") ? "$filename" : null}

    input:
    tuple val(seqID), file(reads) from ch_reads_fastp

    output:
    tuple val(seqID), file("*_trimmed.fastq.gz") into ch_fastp_phix
    file("${seqID}_qc_report.html")

    script:
    """
    fastp \
    -w ${task.cpus} \
    -q ${params.mean_quality} \
    -5 \
    -3 \
    --cut_mean_quality ${params.trimming_quality} \
    --adapter_sequence=${params.adapter_forward} \
    --adapter_sequence_r2=${params.adapter_reverse} \
    -i ${reads[0]} \
    -I ${reads[1]} \
    -o "${seqID}_R1_trimmed.fastq.gz" \
    -O "${seqID}_R2_trimmed.fastq.gz" \
    -h "${seqID}_qc_report.html"
    """
}
