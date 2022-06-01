// Default command to launch is:
//                              nextflow run main.nf -profile [load_profile],[dataset_profile]
//                              where [load_profile] is the profile for your machine and [dataset_profile] the one for your datasets paths and variables!
params.fqpattern = "_{1,2}.fastq.gz"

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
    
    println " Input:      " + params.readPath
    println " Metadata:   " + params.metaPath
    println " Databases:  " + params.dbPath
    println " Output:     " + params.outdir
    println " Resources:  " + params.max_cpus + " CPUs, " + params.max_memory + " RAM"
}
def cursystem = System.properties['os.name']
welcomeScreen()

/* INPUT FILES */

// Sequences
if (params.readPath) {
    if (params.singleEnd) {
        Channel.fromPath("${params.readPath}/*.fastq.gz", checkIfExists: true)
            .map { file -> tuple(file.simpleName, file) }
            .ifEmpty { exit 1, "No input files supplied! Please check params.readPath in your config file!" }
            .set { ch_reads_fastp } 
    } else {
        Channel.fromFilePairs("${params.readPath}/*" + params.fqpattern, checkIfExists: true)
            .ifEmpty { exit 1, "No input files supplied! Please check params.readPath in your config file!" }
            .set { ch_reads_fastp }
    }
} else
    error "No input files supplied! Please check the env variable readPath in your config file!"

// Metadata
if (params.metadata == true) {
    Channel.fromPath("${params.metaPath}/*.csv", checkIfExists: true)
        .ifEmpty { exit 1, "No metadata file found! Please check params.metadata in your config file!"}
        .set { ch_metadata_checker }
    //println " Metadata loaded."
} else {
    println "No metadata supplied! Metadata are mandatory for beta diversity and heatmaps plots!"
    println "If you are interested in these plots, simply provide a metadata .csv file using params.metadata in your config file!"
    println "For more informations, check the metadata paragraph on the wiki at https://github.com/MattiaPandolfoVR/MetaPhage#input_files"
    ch_metadata_kraken_files = Channel.empty()
    params.skip_kraken_files = true
    ch_metadata_summary = Channel.empty()
    params.skip_taxonomy_table = true
    ch_metadata_taxonomy = Channel.empty()
    params.skip_alpha_diversity = true
    ch_metadata_alphadiversity = Channel.empty()
    params.skip_beta_diversity = true
    ch_metadata_betadiversity = Channel.empty()
    params.skip_heatmap = true
    ch_metadata_heatmap = Channel.empty()
    params.skip_violin_plots = true
    ch_metadata_violin = Channel.empty()
}

/* STEP 0 - check metadata and db presence and download required files */

process csv_validator {
    tag "$code"
    label 'min_res'

    input:
    val code from 'Checking metadata format...'
    path(metadata) from ch_metadata_checker

    output:
    path("metadata.csv") into (ch_metadata_kraken_files, ch_metadata_summary, ch_metadata_phylo,
    ch_metadata_taxonomy, ch_metadata_alphadiversity, ch_metadata_betadiversity, ch_metadata_heatmap,
    ch_metadata_violin)

    script:
    """
    mv ${metadata} metadata_to_check.csv
    python3 $projectDir/bin/python/csv_checker.py
    """
}

process db_manager {
    tag "$code"
    label 'low_res'

    input:
    val code from 'Downloading missing databases...'

    output:
    path("db_path.csv") into (ch_dbm_phix, ch_dbm_kraken2, ch_dbm_vibrant, ch_dbm_virsorter, ch_dbm_vcontact2)

    script:
    """
    python3 $projectDir/bin/python/db_manager.py \
    -o ${params.dbPath} \
    -m ${task.cpus}
    """
}

/* STEP 1 - quality check and trimming */

process fastp {
    tag "$seqID"
    label 'med_res'
    publishDir "${params.outdir}/fastp_qc/${seqID}", mode: 'copy', pattern: "*.html"

    when:
    !params.skip_qtrimming

    input:
    tuple val(seqID), file(reads) from ch_reads_fastp

    output:
    tuple val(seqID), file("*_trimmed.fastq*") into (ch_fastp_phix)
    file("${seqID}_qc_report.html")
    file("${seqID}_fastp.json") into ch_fastp_multiqc

    script: 
    def ext = params.keep_phix ? ".gz" : ""
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
    -h ${seqID}_qc_report.html \
    -j ${seqID}_fastp.json
    """
}

if(!params.keep_phix) {

    ch_phix = ch_dbm_phix.splitCsv().flatMap { it -> "${it[0]}" + "/${params.mod_phix}.fasta" }
    
    process remove_phix {
        tag "$seqID"
        label 'low_res'

        when:
        !params.skip_qtrimming

        input:
        each file_phix_alone from ch_phix
        tuple val(seqID), file(reads) from ch_fastp_phix
        
        output:
        tuple val(seqID), file("*.gz") into (ch_trimm_kraken2, ch_trimm_derep)
        tuple val(seqID), file("*.fastq") into ch_trimm_megahit

        script:
        def inp = params.singleEnd ? "-U ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
        def check = params.singleEnd ? "" : "--check-read-2"
        """
        hts_SeqScreener \
        $inp \
        --seq ${file_phix_alone} \
        $check \
        --fastq-output ${seqID}_dephixed \
        --force
	    cp ${seqID}_dephixed_R1.fastq.gz ${seqID}_2_dephixed_R1.fastq.gz
	    cp ${seqID}_dephixed_R2.fastq.gz ${seqID}_2_dephixed_R2.fastq.gz
	    gzip -d ${seqID}_2_dephixed_R1.fastq.gz > ${seqID}_dephixed_R1.fastq
	    gzip -d ${seqID}_2_dephixed_R2.fastq.gz > ${seqID}_dephixed_R2.fastq
        """
    }
}
else {
    ch_fastp_phix.into {ch_trimm_megahit; ch_trimm_derep}
}

/* STEP 2 - taxonomy classification */

ch_kraken2 = ch_dbm_kraken2.splitCsv().flatMap { it -> "${it[1]}" + "/${params.mod_kraken2}/" }

process kraken2 {
    tag "$seqID"
    label 'big_res'
    publishDir "${params.outdir}/taxonomy/kraken2/${seqID}", mode: 'copy'

    when:
    !params.skip_bacterial_taxo && !params.skip_kraken2

    input:
    each file_kraken2_db from ch_kraken2
    tuple val(seqID), file(reads) from ch_trimm_kraken2

    output:
    file("${seqID}_output.txt") into ch_kraken2_output
    tuple val(seqID), file("${seqID}_report.txt") into (ch_kraken2_krona, ch_kraken2_multiqc)

    script:
    def inp = params.singleEnd ? "${reads}" :  "--paired ${reads[0]} ${reads[1]}"
    //println "running"
    """
    kraken2 \
    --report-zero-counts \
    --threads ${task.cpus} \
    --db ${file_kraken2_db} \
    --output ${seqID}_output.txt \
    --report ${seqID}_report.txt \
    $inp
    """
}

process krona {
    cache 'lenient'
    tag "$seqID"
    label 'big_res'
    publishDir "${params.outdir}/taxonomy/krona/${seqID}", mode: 'copy'

    when:
    !params.skip_bacterial_taxo && !params.skip_kraken2 && !params.skip_krona

    input:
    tuple val(seqID), file(report) from ch_kraken2_krona

    output:
    file("${seqID}_*_krona_abundancies.html") into ch_krona_output

    script:
    """
    echo "bracken was not run, using kracken2 output to generate krona report!"
    python $projectDir/bin/python/kreport2krona.py \
    --report-file ${report} \
    --output to_krona.txt 
    ktImportText \
    to_krona.txt \
    -o ${seqID}_krak_krona_abundancies.html
    """
    }

/* STEP 3 - assembly */

process megahit {
    cache 'lenient'
    tag "$seqID"
    label 'big_res'
    publishDir "${params.outdir}/assembly/megahit/${seqID}", mode: 'copy'

    when:
    !params.skip_megahit

    input:
    tuple val(seqID), file(reads) from ch_trimm_megahit

    output:
    tuple val(seqID), val("megahit"), file("${seqID}_megahit_contigs.fasta") into (ch_megahit_quast, ch_megahit_vibrant, ch_megahit_phigaro,
     ch_megahit_virsorter, ch_megahit_virfinder)
    
    script:
    def inp = params.singleEnd ? "--read ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    megahit \
    --num-cpu-threads ${task.cpus} \
    --memory ${task.memory.toBytes()} \
    $inp \
    --out-dir result \
    --out-prefix ${seqID}
    mv result/${seqID}.contigs.fa ${seqID}_megahit_contigs.fasta
    """
}

process metaquast {
    cache 'lenient'
    tag "$assembler-$seqID"
    label 'med_res'
    publishDir "${params.outdir}/assembly/${assembler}/metaquast", mode: 'copy'

    when:
    !params.skip_megahit || !params.skip_metaquast 

    input:
    tuple val(seqID), val(assembler), file(scaffold) from ch_megahit_quast

    output:
    file("*")
    file("${seqID}") into ch_metaquast_multiqc

    script:
    """
    metaquast.py \
    --threads ${task.cpus} \
    --rna-finding \
    --max-ref-number 0 \
    --labels ${seqID} \
    -o ${seqID} \
    ${scaffold}
    """
}

/* STEP 4 - phage mining */

ch_vibrant = ch_dbm_vibrant.splitCsv().flatMap { it -> "${it[2]}" + "/${params.mod_vibrant}/" }

process vibrant {
    cache 'lenient'
    tag "$assembler-$seqID"
    label 'med_res'
    publishDir "${params.outdir}/mining/vibrant/${assembler}", mode: 'copy'

    when:
    !params.skip_mining && !params.skip_vibrant

    input:
    each file_vibrant_db from ch_vibrant
    tuple val(seqID), val(assembler), file(scaffold) from ch_megahit_vibrant

    output:
    file("*")
    tuple val(seqID), val(assembler), val("vibrant"), file("${seqID}_vibrant/VIBRANT_phages_${seqID}_megahit_contigs/${seqID}_correct_vibrant.fna") into ch_vibrant_cdhit optional true

    script:
    if(!params.virome_dataset)
        """
        VIBRANT_run.py \
        -t ${task.cpus} \
        -i ${scaffold} \
        -d ${file_vibrant_db} \
        -m ${file_vibrant_db}
        mv VIBRANT_*/ ${seqID}_vibrant/
        sed 's/^>/>vibrant_/' ${seqID}_vibrant/VIBRANT_phages_${seqID}_megahit_contigs/${seqID}_megahit_contigs.phages_combined.fna >\
         ${seqID}_vibrant/VIBRANT_phages_${seqID}_megahit_contigs/${seqID}_correct_vibrant.fna
        rm -f temp*
        """
    else
        """
        VIBRANT_run.py \
        -t ${task.cpus} \
        -i ${scaffold} \
        -d ${file_vibrant_db} \
        -m ${file_vibrant_db} \
        -virome
        mv VIBRANT_*/ ${seqID}_vibrant/
        sed 's/^>/>vibrant_/' ${seqID}_vibrant/VIBRANT_phages_${seqID}_megahit_contigs/${seqID}_megahit_contigs.phages_combined.fna >\
         ${seqID}_vibrant/VIBRANT_phages_${seqID}_megahit_contigs/${seqID}_correct_vibrant.fna
        rm -f temp*
        """
}

process phigaro {
    errorStrategy 'ignore'
    tag "$assembler-$seqID"
    label 'med_res'
    publishDir "${params.outdir}/mining/phigaro/${assembler}", mode: 'copy', pattern: "${seqID}_phigaro"

    when:
    !params.skip_mining && !params.skip_phigaro

    input:
    tuple val(seqID), val(assembler), file(scaffold) from ch_megahit_phigaro

    output:
    file("*")
    tuple val(seqID), val(assembler), val("phigaro"), file("${seqID}_phigaro/${seqID}_correct_phigaro.fasta") into ch_phigaro_cdhit optional true

    script:
    """
    python $projectDir/bin/python/phigaro_config_creator.py ${params.dbPath}
    printf 'Y\n' | phigaro \
    --threads ${task.cpus} \
    --fasta-file ${scaffold} \
    --config config.yml \
    --print-vogs \
    --output ${seqID}_phigaro \
    --not-open \
    --save-fasta \
    --mode basic
    # Add a control for unexistent phigaro files
    if [ ! -e ${seqID}_phigaro/*.phigaro.fasta ]; then
        mkdir -p ${seqID}_phigaro && touch ${seqID}_phigaro/${seqID}_megahit_contigs.phigaro.fasta
    fi
    # add miner flag at each fasta header
    sed 's/^>/>phigaro_/' ${seqID}_phigaro/${seqID}_megahit_contigs.phigaro.fasta > ${seqID}_phigaro/${seqID}_correct_tmp_phigaro.fasta
    sed 's/_prophage/ prophage/' ${seqID}_phigaro/${seqID}_correct_tmp_phigaro.fasta > ${seqID}_phigaro/${seqID}_correct_phigaro.fasta
    rm ${seqID}_phigaro/${seqID}_correct_tmp_phigaro.fasta
    """
}

ch_virsorter = ch_dbm_virsorter.splitCsv().flatMap { it -> "${it[3]}" + "/${params.mod_virsorter}/" }

process virsorter {
    cache 'lenient'
    tag "$assembler-$seqID"
    label 'med_res'
    publishDir "${params.outdir}/mining/virsorter/${assembler}", mode: 'copy'

    when:
    !params.skip_mining && !params.skip_virsorter

    input:
    each file_virsorter_db from ch_virsorter
    tuple val(seqID), val(assembler), file(scaffold) from ch_megahit_virsorter

    output:
    file("*")
    tuple val(seqID), val(assembler), val("1virsorter"), file("${seqID}_virsorter/Predicted_viral_sequences/${seqID}_correct_virsorter_cat-1.fasta") into ch_1virsorter_cdhit optional true
    tuple val(seqID), val(assembler), val("4virsorter"), file("${seqID}_virsorter/Predicted_viral_sequences/${seqID}_correct_virsorter_prophages_cat-4.fasta") into ch_4virsorter_cdhit optional true

    script:
    def viromes = params.virome_dataset ? "2" : "1"
    if (params.mod_virsorter == "legacy")
        """
        wrapper_phage_contigs_sorter_iPlant.pl \
        -f ${scaffold} \
        --db $viromes \
        --wdir ${seqID}_virsorter \
        --ncpu ${task.cpus} \
        --data-dir ${file_virsorter_db}
        cd ${seqID}_virsorter/Predicted_viral_sequences/
        for FILENAME in *;
        do
            mv \$FILENAME ${seqID}_\$FILENAME;
            if [ ! -e ${seqID}_VIRSorter_cat-1.fasta ]; then
                touch ${seqID}_VIRSorter_cat-1.fasta
            elif [ ! -e  ${seqID}_VIRSorter_prophages_cat-4.fasta ]; then
                touch ${seqID}_VIRSorter_prophages_cat-4.fasta
            else
                break
            fi
        done
        # add a space after sequenceID in header (before flag)
        sed 's/_flag=/ flag=/' ${seqID}_VIRSorter_cat-1.fasta > ${seqID}_correct_virsorter_cat-1.fasta;
        sed 's/_flag=/ flag=/' ${seqID}_VIRSorter_prophages_cat-4.fasta > ${seqID}_correct_virsorter_prophages_cat-4.fasta
        """
    else 
        """
        echo ${file_virsorter_db}
        """
}

process virfinder {
    cache 'lenient'
    tag "$assembler-$seqID"
    label 'big_res'
    publishDir "${params.outdir}/mining/virfinder/${assembler}/${seqID}_virfinder", mode: 'copy'

    when:
    !params.skip_mining && !params.skip_virfinder

    input:
    tuple val(seqID), val(assembler), file(scaffold) from ch_megahit_virfinder

    output:
    file("*")
    tuple val(seqID), val(assembler), val("virfinder"), file("${seqID}_viral_sequences.fasta") into ch_virfinder_cdhit optional true

    script:
    """
    Rscript $projectDir/bin/Rscript/virfinder_execute.R ${scaffold} ${task.cpus} $projectDir
    mv results.txt ${seqID}_results.txt
    mv viral_sequences.fasta ${seqID}_viral_sequences.fasta
    """
}

/* STEP 5 - dereplication and reads mapping */

process cdhit {
    cache 'lenient'
    tag "$assembler"
    label 'med_res'
    publishDir "${params.outdir}/cd-hit/", mode: 'copy'

    when:
    !params.skip_dereplication

    input:
    tuple val(seqID), val(assembler), val(miner), file(scaffolds) from Channel.empty().mix(ch_vibrant_cdhit, ch_phigaro_cdhit, ch_1virsorter_cdhit, ch_4virsorter_cdhit, ch_virfinder_cdhit).groupTuple(by: 1)

    output:
    file("*")
    tuple val(assembler), file("${assembler}_vOTUs_consensus.fasta") into (ch_cdhit_prodigal, ch_cdhit_bowtie2, ch_cdhit_chopper_c, ch_cdhit_taxonomytable_c, ch_cdhit_summary)
    tuple val(assembler), file("${assembler}_min_comp.tsv") into (ch_cdhit_mincomp, ch_cdhit_taxonomytable_min)

    script:
    """
    # append all the mined viral scaffolds
    cat ${scaffolds} > concat_${assembler}.fasta
    # set viral scaffolds reads header to be unique
    seqkit rename concat_${assembler}.fasta > concat_unique_${assembler}.fasta
    # CD-HIT-EST
    cd-hit-est \
    -T ${task.cpus} \
    -M ${task.memory.toMega()} \
    -i concat_unique_${assembler}.fasta \
    -o derep95_${assembler}.fasta \
    -c 0.95 \
    -aS 0.85 \
    -n 9 \
    -d 0 \
    -p 1 \
    -t 4 \
    -g 1
    # filter dereplicated sequences by length
    seqkit seq \
    --min-len ${params.minlen} \
    --out-file filtered_derep95_${assembler}.fasta \
    derep95_${assembler}.fasta
    # obtain which mined viral scaffold collapse in a vOTU
    python $projectDir/bin/python/miner_comparison.py ${assembler}
    # sort vOTU file by name
    seqkit sort --id-regexp ">vOTU_([0-9]+)" vOTU_${assembler}.fasta > ${assembler}_vOTUs_consensus.fasta
    """
}

process prodigal {
    cache 'lenient'
    tag "$assembler"
    label 'med_res'
    publishDir "${params.outdir}/prodigal", mode: 'copy'

    when:
    !params.skip_dereplication

    input:
    tuple val(assembler), file(consensus) from ch_cdhit_prodigal

    output:
    file("*")
    tuple val(assembler), file("${assembler}_vOTUs_proteins.faa") into (ch_prodigal_diamond, ch_prodigal_chopper_p, ch_prodigal_taxonomytable_p) optional true
    tuple val(assembler), file("${assembler}_vOTUs_coords.gff") into (ch_prodigal_chopper_g, ch_prodigal_taxonomytable_g) optional true

    script:
    """
    prodigal \
    -i ${consensus} \
    -o ${assembler}_vOTUs_coords.gff \
    -a ${assembler}_vOTUs_proteins.faa \
    -f gff \
    -p meta
    """
}

process bowtie2_derep {
    cache 'lenient'
    tag "$assembler"
    label 'med_res'

    when:
    !params.skip_dereplication

    input:
    tuple val(assembler), file(consensus) from ch_cdhit_bowtie2.collect()
    tuple val(seqID), file(reads) from ch_trimm_derep

    output:
    tuple val(assembler), file("${seqID}_${consensus}.sorted.bam") into ch_bowtie2bam_covtocounts2 
    tuple val(assembler), file("${seqID}_${consensus}.sorted.bam.bai") into ch_bowtie2bai_covtocounts2

    script:
    """
    bowtie2-build \
    --threads ${task.cpus} ${consensus} ${consensus}_index
    bowtie2 \
    -p ${task.cpus} \
    -x ${consensus}_index \
    -1 ${reads[0]} \
    -2 ${reads[1]} \
    -S ${seqID}_${consensus}.sam 
    samtools view \
    -@ ${task.cpus} \
    -S -b ${seqID}_${consensus}.sam > ${seqID}_${consensus}.bam
    samtools sort \
    -@ ${task.cpus} ${seqID}_${consensus}.bam \
    -o ${seqID}_${consensus}.sorted.bam
    samtools index \
    -@ ${task.cpus} ${seqID}_${consensus}.sorted.bam
    samtools flagstat \
    -@ ${task.cpus} ${seqID}_${consensus}.sorted.bam > mappingstats_${seqID}_${consensus}.txt
    qualimap bamqc \
    -nt ${task.cpus} \
    -outdir qualimap_bamqc_${seqID}_${consensus}.folder \
    -bam ${seqID}_${consensus}.sorted.bam
    """
}

process covtocounts2 {
    cache 'lenient'
    tag "$assembler"
    label 'low_res'
    publishDir "${params.outdir}/report/tables", mode: 'copy'

    when:
    !params.skip_dereplication

    input:
    tuple val(assembler), file(sortedbam) from ch_bowtie2bam_covtocounts2.groupTuple()
    tuple val(assembler), file(sortedbambai) from ch_bowtie2bai_covtocounts2.groupTuple()

    output:
    file("count_table.csv")
    tuple val(assembler), file("count_table.csv") into (ch_covtocounts2_summary, ch_covtocounts2_phylo, ch_covtocounts2_heatmap,
    ch_covtocounts2_alphadiversity, ch_covtocounts2_betadiversity, ch_covtocounts2_violin)

    script:
    """
    bamcountrefs \
    --threads ${task.cpus} \
    --multiqc \
    ${sortedbam} > count_table_mqc.txt
    python $projectDir/bin/multiqc/multiqc_model_editor.py ${assembler}
    cp count_table_mqc.txt count_table.csv 
    sed -i 1,4d count_table.csv
    """
}

/* STEP 6 - viral taxonomy */

ch_diamond = Channel.fromPath("${params.dbPath}/diamond/${params.mod_vcontact2}/allVSall_${params.mod_vcontact2}.csv")
ch_vcontact2_db = ch_dbm_vcontact2.splitCsv().flatMap { it -> "${it[4]}"  }

process diamond_vcontact2 {
    cache 'lenient'
    tag "$assembler"
    label 'big_res'
    publishDir "${params.outdir}/taxonomy/vcontact2/${assembler}", mode: 'copy',
        saveAs: {filename -> filename.startsWith("allVSall_") ? "$filename" : null}

    when:
    !params.skip_dereplication && !params.skip_viral_taxo && !params.skip_vcontact2
    
    input:
    path allVSall from ch_diamond
    val file_vcontact2_db from ch_vcontact2_db
    tuple val(assembler), file(viral_orfs) from ch_prodigal_diamond

    output:
    tuple val(assembler), file("allVSall_${params.mod_vcontact2}_app.csv"), file("viral_genomes_combined.csv") into ch_diamond_vcontact2 optional true

    script:  
    reference_db_faa = file("${file_vcontact2_db}/${params.vcontact2_file_head}proteins.faa")
    gene_to_genome = file("${file_vcontact2_db}/${params.vcontact2_file_head}gene_to_genome.csv")
    reference_db = file("${file_vcontact2_db}/${params.mod_vcontact2}_reference_db.dmnd")

    if(allVSall.exists() == false){
        """
        # Vcontact db: ${file_vcontact2_db}
        # Vcontact file head: ${params.vcontact2_file_head}
        # Vcontact mod: ${params.mod_vcontact2}
        # Reference DB faa: ${reference_db_faa}
        # Gene to genome: ${gene_to_genome}
        # Reference DB: ${reference_db}

        # Create the db/vcontact2 directory to store diamond processed files
        mkdir -p ${params.dbPath}/diamond/${params.mod_vcontact2}/
        # Format reference_db
        diamond makedb -p ${task.cpus} --in $reference_db_faa -d $reference_db
        # Align the reference_db against itself to create allVSall
        diamond blastp -p ${task.cpus} --sensitive -d $reference_db -q $reference_db_faa -o $allVSall
        # Format the viral_orf in the correct way
        $projectDir/bin/python/simplify_faa-ffn_derep.py $viral_orfs
        # Create gene2genome file
        $projectDir/bin/python/vcontact2_gene2genome.py -p ${viral_orfs}.simple.faa -o viral_genomes_g2g.csv -s 'Prodigal-FAA'
        # Combine the mined viral_orfs.faa with inphared protein.faa
        cat $reference_db_faa ${viral_orfs}.simple.faa > ${viral_orfs}.merged.faa 
        # Replace "None_provided" to "none" in viral_genomes_g2g.csv column
        sed -i 's/,None_provided/,none/g' viral_genomes_g2g.csv
        # Combine viral_genomes_g2g.csv with inphared gene_to_genome.csv, removing the header line of the second
        sed -i 1d $gene_to_genome
        cat viral_genomes_g2g.csv $gene_to_genome > viral_genomes_combined.csv
        # Align the viral orfs against the reference_db
        diamond blastp -p ${task.cpus} --sensitive -d $reference_db -q ${viral_orfs}.merged.faa -o viral_orfs_alignment.csv
        # Append the viral orfs alignment to a copy of the allVSall file 
        cat viral_orfs_alignment.csv $allVSall > allVSall_${params.mod_vcontact2}_app.csv
        """
    }
    else {
        """
        # Format the viral_orf in the correct way
        $projectDir/bin/python/simplify_faa-ffn_derep.py $viral_orfs
        # create gene2genome file
        $projectDir/bin/python/vcontact2_gene2genome.py -p ${viral_orfs}.simple.faa -o viral_genomes_g2g.csv -s 'Prodigal-FAA'
        # Combine the mined viral_orfs.faa with inphared protein.faa
        cat $reference_db_faa ${viral_orfs}.simple.faa > ${viral_orfs}.merged.faa 
        # Replace "None_provided" to "none" in viral_genomes_g2g.csv column
        sed -i 's/,None_provided/,none/g' viral_genomes_g2g.csv
        # Combine viral_genomes_g2g.csv with inphared gene_to_genome.csv, removing the header line of the second
        sed -i 1d $gene_to_genome
        cat viral_genomes_g2g.csv $gene_to_genome > viral_genomes_combined.csv
        # Align the viral orfs against the reference_db
        diamond blastp -p ${task.cpus} --sensitive -d $reference_db -q ${viral_orfs}.merged.faa -o viral_orfs_alignment.csv
        # Append the viral orfs alignment to a copy of the allVSall file
        cat viral_orfs_alignment.csv $allVSall > allVSall_${params.mod_vcontact2}_app.csv
        """
    }
}

process vcontact2 {
    cache 'lenient'
    tag "$assembler"
    label 'big_res'
    publishDir "${params.outdir}/taxonomy/vcontact2/${assembler}", mode: 'copy'

    when:
    !params.skip_dereplication && !params.skip_viral_taxo && !params.skip_vcontact2

    input:
    each clusterONE from Channel.fromPath("$projectDir/bin/cluster_one-1.0.jar")
    tuple val(assembler), file(alignment), file(gene2genome) from ch_diamond_vcontact2

    output:
    tuple val(assembler), file("c1.ntw"), file("genome_by_genome_overview.csv") into ch_vcontact2_extender optional true

    script:
    """
    # run vConTACT2
    vcontact2_local \
    -t ${task.cpus} \
    --blast-fp ${alignment} \
    --rel-mode 'Diamond' \
    --proteins-fp ${gene2genome} \
    --db 'None' \
    --pcs-mode MCL \
    --vcs-mode ClusterONE \
    --c1-bin ${clusterONE} \
    --output-dir ./
    """
    stub:
    """
    if [[ -e ${projectDir}/tests/assets/vcontact2/c1/ntw.gz ]]; then
        cp ${projectDir}/tests/assets/vcontact2/*gz .
        gunzip *gz
    fi
    """
}

process graphanalyzer {
    cache 'lenient'
    tag "$assembler"
    label 'med_res'
    publishDir "${params.outdir}/taxonomy/vcontact2", mode: 'copy'
    publishDir "${params.outdir}/report/tables/", mode: 'copy', pattern: '*.csv'

    when:
    !params.skip_dereplication && !params.skip_viral_taxo && !params.skip_vcontact2 && !params.skip_graphanalyzer

    input:
    each file_vcontact2_db from Channel.fromPath("${params.dbPath}/inphared/")
    tuple val(assembler), file(netfile), file(csvfile) from ch_vcontact2_extender


    output:
    file("*") optional true
    tuple val(assembler), file("taxonomy_table.csv") into (ch_vcontact2_phylo, ch_vcontact2_taxonomytable,
    ch_vcontact2_heatmap, ch_vcontact2_alphadiversity, ch_vcontact2_betadiversity, ch_vcontact2_violin) optional true

    script:
    """
    python $projectDir/bin/python/graph_analyzer.py \
    --graph ${netfile} \
    --csv ${csvfile} \
    --metas ${file_vcontact2_db}/*data_excluding_refseq.tsv \
    --output ./ \
    --suffix ${assembler}
    cp custom_taxonomy_table_${assembler}_mqc.txt custom_taxonomy_table_mqc.txt
    python $projectDir/bin/python/taxonomy_table_namer.py
    cp custom_taxonomy_table_mqc.txt custom_taxonomy_table_${assembler}_mqc.csv
    sed -i 1,5d custom_taxonomy_table_${assembler}_mqc.csv
    mv custom_taxonomy_table_${assembler}_mqc.csv taxonomy_table.csv
    rm custom_taxonomy_table_megahit_mqc.txt
    rm custom_taxonomy_table_mqc.txt
    """
}

/* STEP 7 - statistics and report section generation */

process kraken_file {
    cache 'lenient'
    tag "Creating the table..."
    label 'low_res'
    publishDir "${params.outdir}/report/tables", mode: 'copy'
    
    when:
    !params.skip_dereplication && !skip_bacterial_taxo && !params.skip_kraken2 && !params.skip_kraken_files

    input:
    path metadata from ch_metadata_kraken_files
    file(kraken) from ch_kraken2_output.collect()
    file(krona) from ch_krona_output.collect().ifEmpty([])

    output:
    file("kraken_files.html") into ch_krakenfiles_multiqc

    script:
    def path = "${params.outdir}"
    """
    Rscript $projectDir/bin/Rscript/kraken_files.R $path ${metadata}
    """
}

process miner_comparison {
    cache 'lenient'
    tag "$assembler"
    label 'low_res'
    publishDir "${params.outdir}/report/plots", mode: 'copy', pattern: "*comparison.png"
    
    when:
    !params.skip_miner_comparison
    
    input:
    tuple val(assembler), file(min_comp) from ch_cdhit_mincomp.collect().ifEmpty([])
    
    output:
    file("miner_comparison_mqc.png") into ch_mincomp_multiqc
    
    script:
    """
    Rscript $projectDir/bin/Rscript/miner_comparison.R ${min_comp}
    cp miner_comparison.png miner_comparison_mqc.png
    """
}

process summary {
    cache 'lenient'
    tag "$assembler"
    label 'low_res'
    publishDir "${params.outdir}/report/tables", mode: 'copy', pattern: "*_plot.html"

    when:
    !params.skip_violin_plots

    input:
    tuple val(assembler), file(count_table) from ch_covtocounts2_summary.collect().ifEmpty([])
    tuple val(assembler), file(vOTUs_consensus) from ch_cdhit_summary.collect().ifEmpty([])
    path metadata from ch_metadata_summary

    output:
    file("vOTUs_summary.html") into ch_summary_multiqc

    script:
    def path = "${params.outdir}"
    """
    Rscript $projectDir/bin/Rscript/summary_report.R ${count_table} ${vOTUs_consensus} ${metadata} ${params.sum_viol_var} $path
    """
}

process phylo_obj {
    cache 'lenient'
    tag "$assembler"
    label 'min_res'
    publishDir "${params.outdir}/report/tables", mode: 'copy', pattern: "*.rds"

    when:
    params.metadata && !params.skip_dereplication && !params.skip_viral_taxo && !params.skip_vcontact2 && !params.skip_graphanalyzer

    input:
    tuple val(assembler), file(count_table) from ch_covtocounts2_phylo.collect().ifEmpty([])
    tuple val(assembler), file(taxonomy_table) from ch_vcontact2_phylo.collect().ifEmpty([])
    path metadata from ch_metadata_phylo

    output:
    file("phyloseq.rds")
    file("phyloseq_filt_css_norm.rds")

    script:
    """
    Rscript $projectDir/bin/Rscript/phyloseq.R ${count_table} ${taxonomy_table} ${metadata}
    """
}

process file_chopper {
    cache 'lenient'
    tag "$assembler"
    label 'low_res'
    publishDir "${params.outdir}/cd-hit/vOTUs_consensus", mode: 'copy', pattern: "vOTU_*.fasta"
    publishDir "${params.outdir}/cd-hit/vOTUs_proteins", mode: 'copy', pattern: "vOTU_*.faa"
    publishDir "${params.outdir}/cd-hit/vOTUs_coords", mode: 'copy', pattern: "vOTU_*.gff"

    when:
    !params.skip_dereplication

    input:
    tuple val(assembler), path(vOTUs_consensus) from ch_cdhit_chopper_c.collect().ifEmpty([])
    tuple val(assembler), path(vOTUs_proteins) from ch_prodigal_chopper_p.collect().ifEmpty([])
    tuple val(assembler), path(vOTUs_coords) from ch_prodigal_chopper_g.collect().ifEmpty([])

    output:
    file("*")

    script:
    """
    Rscript $projectDir/bin/Rscript/file_chopper.R ${vOTUs_consensus} ${vOTUs_proteins} ${vOTUs_coords}
    """
}

process taxonomy_table {
    cache 'lenient'
    tag "$assembler"
    label 'low_res'
    publishDir "${params.outdir}/report/tables", mode: 'copy'
    
    when:
    !params.skip_dereplication && !skip_viral_taxo && !params.skip_taxonomy_table

    input:
    tuple val(assembler), path(taxonomy_table) from ch_vcontact2_taxonomytable.collect().ifEmpty([])
    tuple val(assembler), path(vOTUs_consensus) from ch_cdhit_taxonomytable_c.collect().ifEmpty([])
    tuple val(assembler), path(vOTUs_proteins) from ch_prodigal_taxonomytable_p.collect().ifEmpty([])
    tuple val(assembler), path(vOTUs_coords) from ch_prodigal_taxonomytable_g.collect().ifEmpty([])
    tuple val(assembler), path(min_comp) from ch_cdhit_taxonomytable_min.collect().ifEmpty([])
    path metadata from ch_metadata_taxonomy

    output:
    file("taxonomy_table.html") into ch_taxonomytable_multiqc

    script:
    def path = "${params.outdir}"
    """
    Rscript $projectDir/bin/Rscript/taxonomy_table.R ${taxonomy_table} ${vOTUs_consensus} ${vOTUs_proteins} ${vOTUs_coords} $path ${min_comp} ${metadata}
    """
}

process heatmap {
    cache 'lenient'
    tag "$assembler"
    label 'med_res'
    publishDir "${params.outdir}/report/plots", mode: 'copy'

    when:
    params.metadata && !params.skip_heatmap

    input:
    tuple val(assembler), file(count_table) from ch_covtocounts2_heatmap.collect().ifEmpty([])
    tuple val(assembler), file(taxonomy_table) from ch_vcontact2_heatmap.collect().ifEmpty([])
    path metadata from ch_metadata_heatmap

    output:
    file("*.html")

    script:
    """
    Rscript $projectDir/bin/Rscript/heatmap.R ${count_table} ${taxonomy_table} ${metadata} ${params.heatmap_var}
    """
}

process alpha_diversity {
    cache 'lenient'
    tag "$assembler"
    label 'med_res'
    publishDir "${params.outdir}/report/plots/alpha_div", mode: 'copy', pattern: "*.html"

    when:
    params.metadata && !params.skip_alpha_diversity

    input:
    tuple val(assembler), file(count_table) from ch_covtocounts2_alphadiversity.collect().ifEmpty([])
    tuple val(assembler), file(taxonomy_table) from ch_vcontact2_alphadiversity.collect().ifEmpty([])
    path metadata from ch_metadata_alphadiversity

    output:
    file("*.html")

    script:
    """
    Rscript $projectDir/bin/Rscript/alpha_diversity.R ${count_table} ${taxonomy_table} ${metadata} ${params.alpha_var1} ${params.alpha_var2}
    """
}

process betadiversity {
    cache 'lenient'
    tag "$assembler"
    label 'med_res'
    publishDir "${params.outdir}/report/plots/beta_div", mode: 'copy', pattern: "*.html"
    when:
    params.metadata && !params.skip_beta_diversity

    input:
    tuple val(assembler), file(count_table) from ch_covtocounts2_betadiversity.collect().ifEmpty([])
    tuple val(assembler), file(taxonomy_table) from ch_vcontact2_betadiversity.collect().ifEmpty([])
    path metadata from ch_metadata_betadiversity

    output:
    file("*.html")

    script:
    """
    Rscript $projectDir/bin/Rscript/beta_diversity.R ${count_table} ${taxonomy_table} ${metadata} ${params.beta_var}
    """
}

process violin_plots {
    cache 'lenient'
    tag "$assembler"
    label 'med_res'
    publishDir "${params.outdir}/report/plots", mode: 'copy', pattern: "*.html"

    when:
    !params.skip_violin_plots

    input:
    tuple val(assembler), path(count_table) from ch_covtocounts2_violin.collect().ifEmpty([])
    tuple val(assembler), file(taxonomy_table) from ch_vcontact2_violin.collect().ifEmpty([])
    path metadata from ch_metadata_violin

    output: 
    file("*.html")

    script:
    """
    Rscript $projectDir/bin/Rscript/violin_plot.R ${count_table} ${taxonomy_table} ${metadata} ${params.violin_var}
    """
}

/* STEP 8 summary report */

process multiqc {
    cache 'lenient'
    tag "Generating the report..."
    label 'low_res'
    publishDir "${params.outdir}/report/", mode: 'copy'

    when:
    !params.skip_report

    input:
    // id:miner_comparison
    file(min_comp) from ch_mincomp_multiqc.collect().ifEmpty([])
    // id:vOTUs_summary
    file(summary_table) from ch_summary_multiqc.collect().ifEmpty([])
    // id:taxonomy_table
    file(taxo_table) from ch_taxonomytable_multiqc.collect().ifEmpty([])
    // id:fastp
    file(fastp) from ch_fastp_multiqc.collect().ifEmpty([])
    // id:kraken2
    file(kraken2) from ch_kraken2_multiqc.collect().ifEmpty([])
    // id:kraken_files
    file(krakenfiles) from ch_krakenfiles_multiqc.collect().ifEmpty([])
    // id: quast
    file(metaquast) from ch_metaquast_multiqc.collect().ifEmpty([])

    output:
    file("MetaPhage_report.html")

    script:
    """
    python $projectDir/bin/python/spacer_taxo.py
    multiqc \
    --config $projectDir/bin/multiqc/multiqc_config.yaml \
    --filename "MetaPhage_report.html" \
    --exclude general_stats $projectDir/bin/multiqc/. -f \
    ${min_comp} vOTUs_summary_mqc.html taxonomy_table_mqc.html ${fastp} ${kraken2} kraken_files_mqc.html ${metaquast}
    """
}

