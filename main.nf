#!/usr/bin/env nextflow


/* CONFIGURATION VARIABLES */

// Trimming 
params.adapter_forward = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
params.adapter_reverse = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
params.mean_quality = 15
params.trimming_quality = 15
params.keep_phix = false
params.mod_phix = "phiX174" 
params.file_phix_alone = "-" 

// microbial taxonomy 
params.skip_kraken2 = false 
params.mod_kraken2 = "miniBAV"
params.file_kraken2_db = "-"                                                                                                                                     

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
            .set {ch_reads_fastp}
    }
}
else
    error "No input files supplied, please check your config file"



/* PROCESSES */

/* STEP 0 - check presence and download required files */
process db_manager {
    echo true
    conda "anaconda::pandas==1.1.3 anaconda::wget==1.20.1 conda-forge::tar==1.29"
    
    publishDir "${params.outdir}/", mode: 'copy',
        saveAs: {filename -> filename.endsWith(".log") ? "$filename" : null}

    output:
    file("db_manager.log")
    file("file_phix_alone") into ch_file_phix_alone
    file("file_kraken2_db") into ch_file_kraken2_db

    script:
    println "Checking presence of required databases. Downloading missing databases… (a detailed log will be created at ./output/db_manager.log)"
    """
    python $workflow.projectDir/bin/db_manager.py \
    --mod_phix ${params.mod_phix} \
    --file_phix_alone ${params.file_phix_alone} \
    --mod_kraken2 ${params.mod_kraken2} \
    --file_kraken2_db ${params.file_kraken2_db}  
    """
   
}

/* STEP 1a - quality check  and trimming */
process fastp {
    conda "bioconda::fastp==0.20.1"
    
    tag "$seqID"
    publishDir "${params.outdir}/fastp/qc/", mode: 'copy',
        saveAs: {filename -> filename.endsWith(".html") ? "$filename" : null}

    input:
    tuple val(seqID), file(reads) from ch_reads_fastp

    output:
    tuple val(seqID), file("*_trimmed.fastq*") into ch_fastp_phix
    file("${seqID}_qc_report.html")

    script:
    { ext = (params.keep_phix == true) ? ".gz" : "" } // hts_SeqScreener accepts only .fastq (NOT .fastq.gz)
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
    -o ${seqID}_R1_trimmed.fastq${ext} \
    -O ${seqID}_R2_trimmed.fastq${ext} \
    -h ${seqID}_qc_report.html
    """
}

/* STEP 1b - phix removal */
if(!params.keep_phix) {
    process remove_phix {
        conda "bioconda::htstream==1.0.0"

        tag "$seqID"
        publishDir "${params.outdir}/hts_SeqScreener/", mode: 'copy',
            saveAs: {filename -> filename.endsWith(".fastq.gz") ? "$filename" : null}

        input:
        file file_phix_alone from ch_file_phix_alone
        tuple val(seqID), file(reads) from ch_fastp_phix

        output:
        tuple val(seqID), file("dephixed*.fastq.gz") into (ch_trimm_kraken2, ch_trimm_spades)

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
}
else {
    ch_fastp_phix.into {ch_trimm_spades; ch_trimm_kraken2}
}

/* STEP 2 - short reads alignment */

process kraken2 {
    conda "bioconda::kraken2==2.1.0 conda-forge::llvm-openmp==11.0.0"
    
    tag "$seqID"
    publishDir "${params.outdir}/taxonomy/kraken2/", mode: 'copy',
        saveAs: {filename -> (filename.endsWith(".kraken") || filename.endsWith(".txt")) ? "$filename" : null}

    when:
    !params.skip_kraken2

    input:
    file file_kraken2_db from ch_file_kraken2_db
    tuple val(seqID), file(reads) from ch_trimm_kraken2

    output:
    file("*.txt")

    script:
    path_file_kraken2_db = file("$workflow.projectDir/db/groovy_vars/${file_kraken2_db}").text.replace("hash.k2d", "")
    """
    kraken2 \
    --report-zero-counts \
    --threads ${task.cpus} \
    --db $workflow.projectDir/${path_file_kraken2_db} \
    --report ${seqID}_report.txt \
    --paired ${reads[0]} ${reads[1]} 
    """
}
