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
params.skip_virsorter = false 
params.mod_virsorter = "legacy"    
params.file_virsorter_db = "-"  
params.virsorter_viromes = false  
params.skip_virfinder = false                                                 

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
            .into { ch_reads_fastp }
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
    file("file_virsorter_db") into (ch_file_virsorter_db)

    script:
    """
    python $workflow.projectDir/bin/db_manager.py \
    --mod_phix ${params.mod_phix} --keep_phix ${params.keep_phix} --file_phix_alone ${params.file_phix_alone} \
    --mod_kraken2 ${params.mod_kraken2} --skip_kraken2 ${params.skip_kraken2} --file_kraken2_db ${params.file_kraken2_db} \
    --mod_vibrant ${params.mod_vibrant} --skip_vibrant ${params.skip_vibrant} --file_vibrant_db ${params.file_vibrant_db} \
    --mod_phigaro ${params.mod_phigaro} --skip_phigaro ${params.skip_phigaro} --file_phigaro_config ${params.file_phigaro_config} \
    --mod_virsorter ${params.mod_virsorter} --skip_virsorter ${params.skip_virsorter} --file_virsorter_db ${params.file_virsorter_db}
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
        tuple val(seqID), file("dephixed*.fastq.gz") into (ch_trimm_kraken2, ch_trimm_metaspades, ch_trimm_megahit, ch_trimm_mapping)

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
    ch_fastp_phix.into {ch_trimm_kraken2; ch_trimm_metaspades; ch_trimm_megahit; ch_trimm_mapping}
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
    tuple val("metaspades"), val(seqID), file("${seqID}_scaffolds.fasta") into (ch_metaspades_quast, ch_metaspades_vibrant, ch_metaspades_phigaro, ch_metaspades_virsorter, ch_metaspades_virfinder)
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
    tuple val("megahit"), val(seqID), file("${seqID}_contigs.fasta") into (ch_megahit_quast, ch_megahit_vibrant, ch_megahit_phigaro, ch_megahit_virsorter, ch_megahit_virfinder)

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
    file "**/*.phages_combined.fna" into (ch_vibrant_collect)

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

process virsorter {
    if (params.mod_virsorter == "legacy") {
        conda "bioconda::virsorter==1.0.6 bioconda::perl-bioperl==1.7.2"
    }
    else {
       conda "bioconda::virsorter==2.0.beta"
    }

    tag "$assembler-$seqID"
    publishDir "${params.outdir}/mining/virsorter/${assembler}", mode: 'copy'

    when:
    !params.skip_virsorter

    input:
    file file_virsorter_db from ch_file_virsorter_db 
    tuple val(assembler), val(seqID), file(scaffold) from Channel.empty().mix(ch_metaspades_virsorter, ch_megahit_virsorter)

    output:
    file("*")

    script:
    path_file_virsorter_db = file("$workflow.projectDir/bin/groovy_vars/${file_virsorter_db}").text
    def viromes = params.virsorter_viromes ? "2" : "1"
    if (params.mod_virsorter == "legacy")
        """
        wrapper_phage_contigs_sorter_iPlant.pl \
        -f ${scaffold} \
        --db $viromes \
        --wdir ${seqID}_virsorter \
        --ncpu ${task.cpus} \
        --data-dir $workflow.projectDir/${path_file_virsorter_db}
        """
    else 
        """
        echo $workflow.projectDir/${path_file_virsorter_db}
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

process cd_hit {
    conda "bioconda::cd-hit==4.8.1 bioconda::seqkit==0.14.0"
    tag "all"
    publishDir "${params.outdir}/CD-HIT", mode: 'copy'

    input:
    file '*.phages_combined.fna' from ch_vibrant_collect.collect()

    output:
    file("*")
    file "splitted83/*.fasta" into (ch_collect_bowtie2)

    script:
    """
    cat *.phages_combined.fna > concat.fasta

    cd-hit-est \
    -T ${task.cpus} \
    -M ${task.memory.toMega()} \
    -i concat.fasta \
    -o derep83.fasta \
    -c 0.83 \
    -n 5 

    seqkit split derep83.fasta \
    --by-id \
    --force \
    --out-dir splitted83
    """
}
/*
process bowtie2 {
    conda "bioconda::bowtie2==2.4.1 bioconda::samtools==1.9 bioconda::qualimap==2.2.2a=1"

    tag "${seqID}-all"
    publishDir "${params.outdir}/bowtie2", mode: 'copy'

    input:
    tuple val(seqID), file(reads) from ch_trimm_mapping
    file consensus_scaffolds from Channel.value(ch_collect_bowtie2.collect())

    output:
    file("*")
    //file "${seqID}.sorted.bam" into (ch_bowtie2_collect2)

    script:
    """
    echo "${seqID}" > ${seqID}.txt

    counter=0
    for scaffold in ${consensus_scaffolds}
    do
        counter=\$((counter+1))
        mkdir \$counter

        bowtie2-build --threads ${task.cpus} \$scaffold \${scaffold}_index

        bowtie2 -p ${task.cpus} -x \${scaffold}_index -1 ${reads[0]} -2 ${reads[1]} -S ${seqID}_\$scaffold.sam 

        samtools view -@ ${task.cpus} -S -b ${seqID}_\$scaffold.sam > ${seqID}_\$scaffold.bam

        samtools sort -@ ${task.cpus} ${seqID}_\$scaffold.bam -o ${seqID}_\$scaffold.sorted.bam

        #samtools index -@ ${task.cpus} ${seqID}_\$scaffold.sorted.bam

        #samtools flagstat -@ ${task.cpus} ${seqID}_\$scaffold.sorted.bam > ${seqID}_mappingstats.txt

        #qualimap bamqc -nt ${task.cpus} -outdir qualimap_bamqc_${seqID} -bam ${seqID}_\$scaffold.sorted.bam

        mv \${scaffold}* \$counter
        mv ${seqID}* \$counter
    done
    """
}
*/

/*
process collect2 {
    conda "bioconda::samtools==1.9 bioconda::qualimap==2.2.2a=1"

    tag "all"
    publishDir "${params.outdir}/collect2", mode: 'copy'

    input:
    file '*.sorted.bam' from ch_bowtie2_collect2.collect()

    output:
    file("*")
    file("collect.sorted.bam") into (ch_collect2_metabat2)


    script:
    """
    samtools merge -@ ${task.cpus} collect.sorted.bam *.sorted.bam

    samtools index -@ ${task.cpus} collect.sorted.bam

    samtools flagstat -@ ${task.cpus} collect.sorted.bam > mappingstats.txt

    qualimap bamqc -nt ${task.cpus} -outdir qualimap_bamqc -bam collect.sorted.bam 
    """
}
*/
/*
process drep {
    //conda "bioconda::drep==2.6.2 bioconda::checkm-genome bioconda::mash bioconda::mummer bioconda::prodigal"
    conda "bioconda::drep==2.6.2"

    tag "all"
    publishDir "${params.outdir}/dereplication", mode: 'copy'

    input:
    file "concat.fasta" from ch_collect_drep

    output:
    file("*")

    script:
    """
    which dRep > drep_v.txt
    which checkm > checkm_v.txt

    dRep dereplicate ./ \
    -p ${task.cpus} \
    -g concat.fasta 
    """
}
*/
/*
process metabat2 {
    conda "bioconda::metabat2==2.14"

    tag "$seqID"
    publishDir "${params.outdir}/binning/metabat2", mode: 'copy'

    input:
    file(concat) from ch_bowtie_metabat
    file(alignment) from ch_collect2_metabat2
    
    output:
    file("*")

    script:
    """

    jgi_summarize_bam_contig_depths \
    --outputDepth depth.txt \
    --percentIdentity 95 \
    ${alignment}
    
    metabat2 \
    --numThreads ${task.cpus} \
    -i ${concat} \
    -a depth.txt \
    -o metabat \
    --minContig 2000 \
    --seed 1 \
    --unbinned \
    --verbose > log.txt
    """
}
*/

/* STEP X - binning */
/*
process bowtie2_samtools {
    conda "bioconda::bowtie2==2.4.1 bioconda::samtools==1.9 bioconda::qualimap==2.2.2a=1"

    tag "$seqID"
    publishDir "${params.outdir}/mapping/bowtie2", mode: 'copy'

    input:
    tuple val(seqID), file(reads) from ch_trimm_bowtie2
    tuple val(assembler), val(seqID), file(scaffold) from ch_metaspades_bowtie2

    output:
    file("*")
    tuple val(assembler), val(seqID), file("${seqID}.sorted.bam"), file("${seqID}.sorted.bam.bai") into (ch_bowtie2_metabat2)

    script:
    """
    echo "${seqID} : ${scaffold} : ${reads[0]} : ${reads[1]}" > echo_bowtie2_${seqID}.txt

    bowtie2-build --threads ${task.cpus} ${scaffold} ${seqID}_index

    bowtie2 -p ${task.cpus} -x ${seqID}_index -1 ${reads[0]} -2 ${reads[1]} -S ${seqID}.sam 

    samtools view -S -b ${seqID}.sam > ${seqID}.bam

    samtools sort -@ ${task.cpus} ${seqID}.bam -o ${seqID}.sorted.bam

    samtools index -@ ${task.cpus} ${seqID}.sorted.bam

    samtools flagstat -@ ${task.cpus} ${seqID}.sorted.bam > mappingstats_${seqID}.txt

    qualimap bamqc -nt ${task.cpus} -outdir qualimap_bamqc_${seqID} -bam ${seqID}.sorted.bam 
    """
}

*/
/*
process metabat2 {
    conda "bioconda::metabat2==2.14"

    tag "$seqID"
    publishDir "${params.outdir}/binning/metabat2", mode: 'copy'

    input:
    tuple val(assembler), val(seqID), file(scaffold) from ch_metaspades_metabat2
    tuple val(assembler), val(seqID), file(alignment), file(index) from ch_bowtie2_metabat2
    
    output:
    file("*")

    script:
    """
    echo "${seqID} : ${scaffold} : ${alignment} : ${index}" > echo_metabat2_${seqID}.txt

    jgi_summarize_bam_contig_depths \
    --outputDepth ${seqID}_depth.txt \
    --percentIdentity 95 \
    ${alignment}
    
    metabat2 \
    --numThreads ${task.cpus} \
    -i ${scaffold} \
    -a ${seqID}_depth.txt \
    -o ${seqID}.metabat \
    --minContig 2000 \
    --seed 1 \
    --unbinned \
    --verbose > ${seqID}_log.txt
    """
}
*/

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