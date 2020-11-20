#!/usr/bin/env nextflow


/* CONFIGURATION VARIABLES */
params.db_manager_reports = false

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
params.skip_megahit = false
params.skip_quast = false

// Phage-hunting 
params.skip_vibrant = false 
params.mod_vibrant = "legacy"    
params.file_vibrant_db = "-" 
params.skip_phigaro = false 
params.mod_phigaro = "standard"    
params.file_phigaro_config = "-" 
params.skip_virsorter2 = false 
params.mod_virsorter2 = "legacy"    
params.file_virsorter2_db = "-"  
params.virsorter2_viromes = false  
params.skip_virfinder = false 
params.mod_virfinder = "standard"    
params.file_virfinder_db = "-"                                                  

// Viral Taxonomy - vContact2 
params.skip_vcontact2 = true // true while debugging the pipeline                                                                               


  
/* FILE INPUT */

if (params.readPaths) {         // declared in profile config
    if (params.singleEnd) {     // declared in profile config
        Channel.fromPath("${params.readPaths}/*_single.fastq.gz", checkIfExists: true)
            .map { file -> tuple(file.simpleName, file) }
            .ifEmpty { exit 1, "No input files supplied! Please check params.readPaths in your config file!" }
            .set { ch_reads_fastp } 
    } else {
        Channel.fromFilePairs("${params.readPaths}/*_{R1,R2}.fastq.gz", checkIfExists: true)
            .ifEmpty { exit 1, "No input files supplied! Please check params.readPaths in your config file!" }
            .set { ch_reads_fastp }
    }
} else
    error "No input files supplied! Please check params.readPaths in your config file!"



/* PROCESSES */

def welcomeScreen() {
    println """
        ====================================================
         __  __      _        _____  _                      
        |  \\/  |    | |      |  __ \\| |                     
        | \\  / | ___| |_ __ _| |__) | |__   __ _  __ _  ___ 
        | |\\/| |/ _ \\ __/ _` |  ___/| '_ \\ / _` |/ _` |/ _ \\
        | |  | |  __/ || (_| | |    | | | | (_| | (_| |  __/
        |_|  |_|\\___|\\__\\__,_|_|    |_| |_|\\__,_|\\__, |\\___|
                                                  __/ |     
                                                 |___/      
        ====================================================
        """.stripIndent()
}
def cursystem = System.properties['os.name']
welcomeScreen()

/* STEP 0 - check presence and download required files */
process db_manager {
    if (params.db_manager_reports) { echo true }
    conda "bioconda::vibrant==1.2.1 conda-forge::tar==1.29 bioconda::phigaro==2.3.0"

    tag "Downloading missing databases..."

    output:
    file("file_phix_alone") into (ch_file_phix_alone)
    file("file_kraken2_db") into (ch_file_kraken2_db, ch_file_bracken_db)
    file("file_vibrant_db") into (ch_file_vibrant_db)
    file("file_phigaro_config") into (ch_file_phigaro_config)
    file("file_virsorter2_db") into (ch_file_virsorter2_db)

    script:
    """
    python $workflow.projectDir/bin/db_manager.py \
    --mod_phix ${params.mod_phix} \
    --file_phix_alone ${params.file_phix_alone} \
    --mod_kraken2 ${params.mod_kraken2} \
    --file_kraken2_db ${params.file_kraken2_db} \
    --mod_vibrant ${params.mod_vibrant} \
    --file_vibrant_db ${params.file_vibrant_db} \
    --mod_phigaro ${params.mod_phigaro} \
    --file_phigaro_config ${params.file_phigaro_config} \
    --mod_virsorter2 ${params.mod_virsorter2} \
    --file_virsorter2_db ${params.file_virsorter2_db}
    """
}

/* STEP 1a - quality check and trimming */
process fastp {
    conda "bioconda::fastp==0.20.1"

    tag "$seqID"
    publishDir "${params.outdir}/fastp_qc/", mode: 'copy',
        saveAs: {filename -> filename.endsWith(".html") ? "$filename" : null}

    input:
    tuple val(seqID), file(reads) from ch_reads_fastp

    output:
    tuple val(seqID), file("*_trimmed.fastq*") into (ch_fastp_phix)
    file("${seqID}_qc_report.html")

    script:
    def ext = params.keep_phix ? ".gz" : "" // hts_SeqScreener accepts only .fastq (NOT .fastq.gz)
    def inp = params.singleEnd ? "-i ${reads[0]}" : "-i ${reads[0]} -I ${reads[1]}"
    def out = params.singleEnd ? "-o ${seqID}_trimmed.fastq${ext}" : "-o ${seqID}_R1_trimmed.fastq${ext} -O ${seqID}_R2_trimmed.fastq${ext}"
    """
    fastp \
    -w ${task.cpus} \
    -q ${params.mean_quality} \
    -5 \
    -3 \
    --cut_mean_quality ${params.trimming_quality} \
    --adapter_sequence=${params.adapter_forward} \
    --adapter_sequence_r2=${params.adapter_reverse} \
    $inp \
    $out \
    -h ${seqID}_qc_report.html
    """
}

/* STEP 1b - phix removal */
if(!params.keep_phix) {
    process remove_phix {
        conda "bioconda::htstream==1.3.3 conda-forge::boost==1.70.0"

        tag "$seqID"

        input:
        file file_phix_alone from ch_file_phix_alone
        tuple val(seqID), file(reads) from ch_fastp_phix

        output:
        tuple val(seqID), file("dephixed*.fastq.gz") into (ch_trimm_kraken2, ch_trimm_metaspades, ch_trimm_megahit)

        script:
        path_file_phix_alone = file("$workflow.projectDir/bin/groovy_vars/${file_phix_alone}").text
        def inp = params.singleEnd ? "-U ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
        def check = params.singleEnd ? "" : "--check-read-2"
        """
        hts_SeqScreener \
        $inp \
        --seq $workflow.projectDir/${path_file_phix_alone} \
        $check \
        --fastq-output dephixed_${seqID} \
        --force
        """
    }
}
else {
    ch_fastp_phix.into {ch_trimm_kraken2; ch_trimm_metaspades; ch_trimm_megahit}
}

/* STEP 2 - short reads alignment */
process kraken2 {
    conda "bioconda::kraken2==2.1.0 conda-forge::llvm-openmp==11.0.0"
    
    tag "$seqID"
    publishDir "${params.outdir}/taxonomy/kraken2/", mode: 'copy'

    when:
    !params.skip_kraken2

    input:
    file file_kraken2_db from ch_file_kraken2_db
    tuple val(seqID), file(reads) from ch_trimm_kraken2

    output:
    file("${seqID}_output.txt") 
    file("${seqID}_report.txt") into (ch_kraken2_bracken)
    val(seqID) into (ch_seqID_bracken)

    script:
    path_file_kraken2_db = file("$workflow.projectDir/bin/groovy_vars/${file_kraken2_db}").text
    def input = params.singleEnd ? "${reads}" :  "--paired ${reads[0]} ${reads[1]}"
    """
    kraken2 \
    --report-zero-counts \
    --threads ${task.cpus} \
    --db $workflow.projectDir/${path_file_kraken2_db} \
    --output ${seqID}_output.txt \
    --report ${seqID}_report.txt \
    $input
    """
}

process bracken {
    conda "bioconda::bracken==2.5.3 conda-forge::libcxx==9.0.1 conda-forge::llvm-openmp==10.0.1 conda-forge::python=3.7 conda-forge::python_abi==3.7=1_cp37m"

    tag "$seqID"
    publishDir "${params.outdir}/taxonomy/bracken/", mode: 'copy'

    when:
    !params.skip_kraken2 && !params.skip_bracken

    input:
    file file_bracken_db from ch_file_bracken_db
    file report from ch_kraken2_bracken
    val seqID from ch_seqID_bracken

    output:
    file("${seqID}_abundancies.txt")
    file("${seqID}_report_bracken_species.txt") into (ch_bracken_krona)
    val(seqID) into (ch_seqID_krona)

    script:
    path_file_bracken_db = file("$workflow.projectDir/bin/groovy_vars/${file_bracken_db}").text
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
    publishDir "${params.outdir}/taxonomy/krona/", mode: 'copy'

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
    if (cursystem.contains('Mac')) {
        conda "bioconda::spades==3.14.1 conda-forge::llvm-openmp==10.0.1"
    }
    else { // spades on Linux has slightly different dependencies
        conda "bioconda::spades==3.14.1 conda-forge::llvm-openmp==8.0.0"
    }

    tag "$seqID"
    publishDir "${params.outdir}/assembly/metaspades/${seqID}", mode: 'copy'
    
    when:
    !params.singleEnd && !params.skip_metaspades

    input:
    tuple val(seqID), file(reads) from ch_trimm_metaspades

    output:
    tuple val("metaspades"), val(seqID), file("${seqID}_scaffolds.fasta") into (ch_metaspades_quast, ch_metaspades_vibrant, ch_metaspades_phigaro, ch_metaspades_virsorter2, ch_metaspades_virfinder)
    tuple val(seqID), file("${seqID}_contigs.fasta")

    script:
    def input = params.singleEnd ? "" : "--pe1-1 ${reads[0]} --pe1-2 ${reads[1]}"  //check when metaspades accepts single-end read libraries!
    """    
    spades.py \
    --meta \
    --threads ${task.cpus} \
    --memory ${task.memory.toGiga()} \
    $input \
    -o ./

    mv scaffolds.fasta ${seqID}_scaffolds.fasta
    mv contigs.fasta ${seqID}_contigs.fasta
    """
}

process megahit {
    conda "bioconda::megahit==1.2.9"

    tag "$seqID"
    publishDir "${params.outdir}/assembly/megahit/${seqID}", mode: 'copy'

    when:
    !params.skip_megahit

    input:
    tuple val(seqID), file(reads) from ch_trimm_megahit

    output:
    tuple val("megahit"), val(seqID), file("${seqID}_contigs.fasta") into (ch_megahit_quast, ch_megahit_vibrant, ch_megahit_phigaro, ch_megahit_virsorter2, ch_megahit_virfinder)

    script:
    def input = params.singleEnd ? "--read ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    megahit \
    --num-cpu-threads ${task.cpus} \
    --memory ${task.memory.toBytes()} \
    $input \
    --out-dir result \
    --out-prefix ${seqID}

    mv result/${seqID}.contigs.fa ${seqID}_contigs.fasta
    """
}

process quast {
    conda "bioconda::quast==5.0.2"

    tag "$assembler-$seqID"
    publishDir "${params.outdir}/assembly/${assembler}/quast/${seqID}", mode: 'copy'

    when:
    !params.skip_metaspades && !params.skip_megahit && !params.skip_quast 

    input:
    tuple val(assembler), val(seqID), file(scaffold) from Channel.empty().mix(ch_metaspades_quast, ch_megahit_quast)

    output:
    file("report.html")
    file("report.pdf")
    file("report.tsv")

    script:
    """
    metaquast.py \
    --threads ${task.cpus} \
    --rna-finding \
    --max-ref-number 0 \
    --labels ${seqID} \
    -o ./ \
    ${scaffold}
    """
} 

/* STEP 4 - phage mining */
process vibrant {
    if (params.mod_vibrant == "legacy") {
        conda "bioconda::vibrant==1.0.1"
    }
    else {
       conda "bioconda::vibrant==1.2.1"
    }
    

    tag "$assembler-$seqID"
    publishDir "${params.outdir}/mining/vibrant/${assembler}", mode: 'copy'

    when:
    !params.skip_vibrant

    input:
    file file_vibrant_db from ch_file_vibrant_db
    tuple val(assembler), val(seqID), file(scaffold) from Channel.empty().mix(ch_metaspades_vibrant, ch_megahit_vibrant)

    output:
    file("*")
    tuple val(assembler), val(seqID), file("**/*.phages_combined.faa") into (ch_vibrant_vcontact2)

    script:
    path_file_vibrant_db = file("$workflow.projectDir/bin/groovy_vars/${file_vibrant_db}").text
    if (params.mod_vibrant == "legacy")
        """
        VIBRANT_run.py \
        -t ${task.cpus} \
        -i ${scaffold} \
        -k $workflow.projectDir/${path_file_vibrant_db}KEGG_profiles_prokaryotes.HMM \
        -p $workflow.projectDir/${path_file_vibrant_db}Pfam-A_v32.HMM \
        -v $workflow.projectDir/${path_file_vibrant_db}VOGDB94_phage.HMM \
        -e $workflow.projectDir/${path_file_vibrant_db}Pfam-A_plasmid_v32.HMM \
        -a $workflow.projectDir/${path_file_vibrant_db}Pfam-A_phage_v32.HMM \
        -c $workflow.projectDir/${path_file_vibrant_db}VIBRANT_categories.tsv \
        -n $workflow.projectDir/${path_file_vibrant_db}VIBRANT_names.tsv \
        -s $workflow.projectDir/${path_file_vibrant_db}VIBRANT_KEGG_pathways_summary.tsv \
        -m $workflow.projectDir/${path_file_vibrant_db}VIBRANT_machine_model.sav \
        -g $workflow.projectDir/${path_file_vibrant_db}VIBRANT_AMGs.tsv
        """
    else 
        """
        VIBRANT_run.py \
        -t ${task.cpus} \
        -i ${scaffold} \
        -folder ./ \
        -d $workflow.projectDir/${path_file_vibrant_db}databases/ \
        -m $workflow.projectDir/${path_file_vibrant_db}files/ 
        """
}

process phigaro {
    conda "bioconda::phigaro==2.3.0"

    tag "$assembler-$seqID"
    publishDir "${params.outdir}/mining/phigaro/${assembler}", mode: 'copy'

    when:
    !params.skip_phigaro

    input:
    file file_phigaro_config from ch_file_phigaro_config // this act just like a timer
    tuple val(assembler), val(seqID), file(scaffold) from Channel.empty().mix(ch_metaspades_phigaro, ch_megahit_phigaro)

    output:
    file("*")

    script:
    path_file_phigaro_config = file("$workflow.projectDir/bin/groovy_vars/${file_phigaro_config}").text
    """
    python $workflow.projectDir/bin/phigaro_config_creator.py

    printf 'Y\n' | phigaro \
    --threads ${task.cpus} \
    --fasta-file ${scaffold} \
    --config config.yml \
    --print-vogs \
    --extension html \
    --output phigaro_${seqID} \
    --not-open \
    --save-fasta \
    --mode basic
    """
}

process virsorter2 {
    if (params.mod_virsorter2 == "legacy") {
        conda "bioconda::virsorter==1.0.6"
    }
    else {
       conda "bioconda::virsorter==2.0.beta"
    }

    tag "$assembler-$seqID"
    publishDir "${params.outdir}/mining/virsorter2/${assembler}/${seqID}", mode: 'copy'

    when:
    !params.skip_virsorter2

    input:
    file file_virsorter2_db from ch_file_virsorter2_db 
    tuple val(assembler), val(seqID), file(scaffold) from Channel.empty().mix(ch_metaspades_virsorter2, ch_megahit_virsorter2)

    output:
    file("*")

    script:
    path_file_virsorter2_db = file("$workflow.projectDir/bin/groovy_vars/${file_virsorter2_db}").text
    def viromes = params.virsorter2_viromes ? "2" : "1"
    if (params.mod_virsorter2 == "legacy")
        """
        wrapper_phage_contigs_sorter_iPlant.pl \
        -f ${scaffold} \
        --db $viromes \
        --wdir ./ \
        --ncpu ${task.cpus} \
        --data-dir $workflow.projectDir/${path_file_virsorter2_db}
        """
    else 
        """
        echo $workflow.projectDir/${path_file_virsorter2_db}
        """
}

process virfinder {
    conda "bioconda::r-virfinder==1.1"

    tag "$assembler-$seqID"
    publishDir "${params.outdir}/mining/virfinder/${assembler}", mode: 'copy'

    when:
    !params.skip_virfinder

    input:
    tuple val(assembler), val(seqID), file(scaffold) from Channel.empty().mix(ch_metaspades_virfinder, ch_megahit_virfinder)

    output:
    file("*")

    script:
    """
    Rscript $workflow.projectDir/bin/virfinder_execute.R ${scaffold}
    mv results.txt ${seqID}_results.txt
    """
}

/* STEP 5 - viral taxonomy */
process vcontact2 {
    conda "bioconda::vcontact2==0.9.19"

    tag "$assembler-$seqID"
    publishDir "${params.outdir}/taxonomy/vcontact2/${assembler}/${seqID}", mode: 'copy'

    when:
    !params.skip_vcontact2

    input:
    tuple val(assembler), val(seqID), file(phages_combined) from ch_vibrant_vcontact2

    output:
    file("*")

    script:
    """
    $workflow.projectDir/bin/simplify_faa-ffn.py ${phages_combined}

    $workflow.projectDir/bin/vcontact2_gene2genome.py \
    -p ${phages_combined}.simple.faa \
    -o viral_genomes_g2g.csv \
    -s 'Prodigal-FAA'

    vcontact2 \
    -t ${task.cpus} \
    --raw-proteins ${phages_combined}.simple.faa \
    --proteins-fp viral_genomes_g2g.csv \
    --db 'ProkaryoticViralRefSeq94-Merged' \
    --pcs-mode MCL \
    --vcs-mode ClusterONE \
    --c1-bin $workflow.projectDir/bin/cluster_one-1.0.jar \
    --output-dir ./
    """
}