#!/usr/bin/env nextflow


/* CONFIGURATION VARIABLES */

// Trimming 
params.adapter_forward = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
params.adapter_reverse = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
params.mean_quality = 15
params.trimming_quality = 15
params.keep_phix = false
params.mod_phix = "WA11" 
params.file_phix_alone = "-" 

// microbial taxonomy 
params.skip_kraken = false 
params.mod_kraken2 = "-"                                                                                                                                          

// Assembly 
params.skip_spades = false
params.skip_quast = false

// Phage-hunting - Vibrant 
params.skip_vibrant = false 
params.mod_vibrant = "-"                                                          

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

/* STEP 0 - check presence and download required files */
process db_manager {
    conda "anaconda::pandas==1.1.3 anaconda::wget==1.20.1"
    
    publishDir "${params.outdir}/", mode: 'copy',
        saveAs: {filename -> filename.endsWith(".log") ? "$filename" : null}

    output:
    file("db_manager.log")
    file("file_phix_alone") into ch_file_phix_alone

    script:
    """
    python $workflow.projectDir/bin/db_manager.py \
    --mod_phix ${params.mod_phix} \
    --file_phix_alone ${params.file_phix_alone}
    """
   
}

/* STEP 1 - quality check  and trimming */
process fastp {
    conda "bioconda::fastp==0.20.1"
    
    tag "$seqID"
    publishDir "${params.outdir}/fastp/qc/", mode: 'copy',
        saveAs: {filename -> filename.endsWith(".html") ? "$filename" : null}

    input:
    tuple val(seqID), file(reads) from ch_reads_fastp

    output:
    tuple val(seqID), file("*_trimmed.fastq") into ch_fastp_phix
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
    -o ${seqID}_R1_trimmed.fastq \
    -O ${seqID}_R2_trimmed.fastq \
    -h ${seqID}_qc_report.html
    """
}

process remove_phix {
    conda "bioconda::htstream==1.0.0"

    tag "$seqID"
    publishDir "${params.outdir}/hts_SeqScreener/", mode: 'copy',
        saveAs: {filename -> filename.endsWith(".fastq.gz") ? "$filename" : null}

    when:
    !params.keep_phix

    input:
    file file_phix_alone from ch_file_phix_alone
    tuple val(seqID), file(reads) from ch_fastp_phix

    output:
    tuple val(seqID), file("dephixed*.fastq.gz") into ch_phix_kraken2

    script:
    path_file_phix_alone = file("$workflow.projectDir/db/groovy_vars/${file_phix_alone}").text
    """
    hts_SeqScreener \
    -1 ${reads[0]} \
    -2 ${reads[1]} \
    --seq $workflow.projectDir/${path_file_phix_alone} \
    --check-read-2 \
    --gzip-output \
    --prefix dephixed_${seqID} \
    --force
    """
}