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
params.skip_bracken = false 
params.bracken_read_length = 100
params.bracken_abundance_level = "S"                                                                                                                                    

// Assembly 
params.skip_metaspades = false
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
    file("file_phix_alone") into (ch_file_phix_alone)
    file("file_kraken2_db") into (ch_file_kraken2_db, ch_file_bracken_db)

    script:
    println "\n\nChecking presence of required databases. Downloading missing databases… (a detailed log will be created at ./output/db_manager.log). Several GB may to be downloaded: this could take long time!\nWait please...\n\n"
    """
    python -u $workflow.projectDir/bin/db_manager.py \
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
    tuple val(seqID), file("*_trimmed.fastq*") into (ch_fastp_phix)
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
        conda "bioconda::htstream==1.3.3 conda-forge::boost==1.70.0"

        tag "$seqID"
        publishDir "${params.outdir}/hts_SeqScreener/", mode: 'copy',
            saveAs: {filename -> filename.endsWith(".fastq.gz") ? "$filename" : null}

        input:
        file file_phix_alone from ch_file_phix_alone
        tuple val(seqID), file(reads) from ch_fastp_phix

        output:
        tuple val(seqID), file("dephixed*.fastq.gz") into (ch_trimm_kraken2, ch_trimm_metaspades)

        script:
        path_file_phix_alone = file("$workflow.projectDir/db/groovy_vars/${file_phix_alone}").text
        """
        hts_SeqScreener \
        -1 ${reads[0]} \
        -2 ${reads[1]} \
        --seq $workflow.projectDir/${path_file_phix_alone} \
        --check-read-2 \
        --fastq-output dephixed_${seqID} \
        --force
        """
    }
}
else {
    ch_fastp_phix.into {ch_trimm_metaspades; ch_trimm_kraken2}
}

/* STEP 2 - short reads alignment */
process kraken2 {
    conda "bioconda::kraken2==2.1.0 conda-forge::llvm-openmp==11.0.0"
    
    tag "$seqID"
    publishDir "${params.outdir}/taxonomy/kraken2/", mode: 'copy',
        saveAs: {filename -> (filename.endsWith(".kraken2") || filename.endsWith(".txt")) ? "$filename" : null}

    when:
    !params.skip_kraken2

    input:
    file file_kraken2_db from ch_file_kraken2_db
    tuple val(seqID), file(reads) from ch_trimm_kraken2

    output:
    file("${seqID}_output.kraken2") 
    file("${seqID}_report.txt") into (ch_kraken2_bracken)
    val(seqID) into (ch_seqID_bracken)

    script:
    path_file_kraken2_db = file("$workflow.projectDir/db/groovy_vars/${file_kraken2_db}").text.replace("hash.k2d", "")
    """
    kraken2 \
    --report-zero-counts \
    --threads ${task.cpus} \
    --db $workflow.projectDir/${path_file_kraken2_db} \
    --output ${seqID}_output.kraken2 \
    --report ${seqID}_report.txt \
    --paired ${reads[0]} ${reads[1]} 
    """
}

process bracken {
    conda "bioconda::bracken==2.5.3 conda-forge::libcxx==9.0.1 conda-forge::llvm-openmp==10.0.1 conda-forge::python=3.7 conda-forge::python_abi==3.7=1_cp37m"

    tag "$seqID"
    publishDir "${params.outdir}/taxonomy/bracken/", mode: 'copy',
        saveAs: {filename -> filename.endsWith(".txt") ? "$filename" : null}

    when:
    !params.skip_bracken

    input:
    file file_bracken_db from ch_file_bracken_db
    file report from ch_kraken2_bracken
    val seqID from ch_seqID_bracken

    output:
    file("${seqID}_abundancies.txt")
    file("${seqID}_report_bracken_species.txt") into (ch_bracken_krona)
    val(seqID) into (ch_seqID_krona)

    script:
    path_file_bracken_db = file("$workflow.projectDir/db/groovy_vars/${file_bracken_db}").text.replace("hash.k2d", "")
    """
    bracken \
    -d $workflow.projectDir/${path_file_bracken_db} \
    -i ${report} \
    -o ${seqID}_abundancies.txt \
    -r ${params.bracken_read_length} \
    -l ${params.bracken_abundance_level} 
    """
}

process krona {
    conda "bioconda::krona==2.7.1 anaconda::python==3.7"

    tag "$seqID"
    publishDir "${params.outdir}/taxonomy/krona/", mode: 'copy',
        saveAs: {filename -> filename.endsWith(".html") ? "$filename" : null}

    when:
    !params.skip_kraken2 && !params.skip_bracken

    input:
    file report_bracken from ch_bracken_krona
    val seqID from ch_seqID_krona

    output:
    file("${seqID}_krona_abundancies.html")

    script:
    """
    python $workflow.projectDir/bin/kreport2krona.py \
    --report-file ${report_bracken} \
    --output to_krona.txt 

    ktImportText \
    to_krona.txt \
    -o ${seqID}_krona_abundancies.html
    """
}

/* STEP 3 - assembly */
process metaSPAdes {
    conda "bioconda::spades==3.14.1 conda-forge::llvm-openmp==8.0.0"

    tag "$seqID"
    publishDir "${params.outdir}/assembly/$seqID", mode: 'copy',
        saveAs: {filename -> filename.endsWith(".fasta") ? "$filename" : null}
    
    when:
    !params.skip_metaspades

    input:
    tuple val(seqID), file(reads) from ch_trimm_metaspades

    output:
    tuple val(seqID), file("${seqID}_scaffolds.fasta") into ch_metaspades_quast
    tuple val(seqID), file("${seqID}_contigs.fasta")

    script:
    """    
    spades.py \
    --meta \
    --threads ${task.cpus} \
    --memory ${task.memory.toGiga()} \
    --pe1-1 ${reads[0]} \
    --pe1-2 ${reads[1]} \
    -o ./

    mv scaffolds.fasta ${seqID}_scaffolds.fasta
    mv contigs.fasta ${seqID}_contigs.fasta
    """
}

process quast {
    tag "$seqID"
    publishDir "${params.outdir}/assembly/$seqID/quast", mode: 'copy',
        saveAs: { filename -> filename.endsWith(".tsv") ? "$filename" : null }

    when:
    !params.skip_spades || !params.skip_quast 

    input:
    tuple val(seqID), file(scaffold) from ch_metaspades_quast

    output:
    file("$seqID/quast/*") into ch_quast_results

    script:
    """
    metaquast.py \
    -t ${task.cpus} \
    --rna-finding \
    --max-ref-number 0 \
    -l ${seqID} ${scaffold} \
    -o "${seqID}/quast/"
    """
} 
