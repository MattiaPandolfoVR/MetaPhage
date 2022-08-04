// Coded by Mattia Pandolfo (mattia.pandolfo@univr.it)
// Default command to launch is:
//
// nextflow run metaphage.nf -profile [load_profile],[dataset_profile]
//
// where [load_profile] is the profile for your machine and [dataset_profile] the one for your datasets paths and variables!

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
        Channel.fromFilePairs("${params.readPath}/*" + "${params.fqpattern}", checkIfExists: true)
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
    cache 'lenient'
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
    python $projectDir/bin/python/csv_checker.py
    """
}

process db_manager {
    cache 'lenient'
    tag "$code"
    label 'low_res'

    input:
    val code from 'Downloading missing databases...'

    output:
    path("db_path.csv") into (ch_dbm_phix, ch_dbm_kraken2, ch_dbm_vibrant, ch_dbm_virsorter2, ch_dbm_vcontact2) //ch_dbm_checkv ch_set_virsorter2 ch_dbm_allvsall

    script:
    """
    python $projectDir/bin/python/db_manager.py \
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
    publishDir "${params.outdir}/taxonomy/kraken2/${seqID}", mode: 'copy', pattern: "*.txt"
    publishDir "${params.outdir}/report/taxonomy/kraken2/${seqID}", mode: 'copy', pattern: "*_report.txt"

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
    publishDir "${params.outdir}/report/taxonomy/krona/${seqID}", mode: 'copy', pattern: "*_krak_krona_abundancies.html"

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
    tuple val(seqID), val("megahit"), file("${seqID}_megahit_contigs.fasta") into (ch_megahit_quast, ch_megahit_deepvirfinder, ch_megahit_phigaro, ch_megahit_vibrant,
     ch_megahit_virsorter2, ch_megahit_virfinder)
    
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
    !params.skip_megahit && !params.skip_metaquast 

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
    --max-ref-number ${params.metaquast_max_ref_num} \
    --labels ${seqID} \
    -o ${seqID} \
    ${scaffold}
    """
}

/* STEP 4 - phage mining */

process deepvirfinder {
    cache 'lenient'
    tag "$assembler-$seqID"
    label 'med_res'
    publishDir "${params.outdir}/mining/deepvirfinder/${assembler}/${seqID}_deepvirfinder", mode: 'copy',
        saveAs: {filename -> filename.endsWith("_results.txt") | filename.endsWith("_viral_sequences.fasta") ? "$filename" : null}

    when:
    !params.skip_mining && !params.skip_deepvirfinder

    input:
    tuple val(seqID), val(assembler), file(scaffold) from ch_megahit_deepvirfinder

    output:
    file("*")
    tuple val(seqID), val(assembler), val("deepvirfinder"), file("${seqID}_correct_viral_sequences.fasta") into ch_deepvirfinder_cdhit optional true

    script:
    """
    # copy all files from deepvirfinder master and create a temporary output dir
    cp -r $projectDir/bin/DeepVirFinder-master/* .
    mkdir ./output
    # launch deepvirfinder
    python dvf.py \
    -c ${task.cpus} \
    -m ./models \
    -i ./${scaffold} \
    -o ./output \
    -l 0
    # rename result file
    mv ./output/*dvfpred.txt ./${seqID}_results.txt
    # since deepvirfinder outputs a list of viral scaffolds headers, we need to collect these and extract the related viral sequence only if they respect the pvalue threshold
	python << END
	import pandas as pnd
	df = pnd.read_csv(
	    './${seqID}_results.txt', sep='\\t', header=0,
	    names=['header', 'length', 'score', 'pvalue']
	)
	f = open('./${seqID}_filtered_headers.txt', 'w')
	for index, row in df.iterrows():
	    if row.pvalue <= ${params.dvf_maxpvalue}:
	        f.write(row.header + '\\n')
	f.close()
	END
    # join the viral scaffold header with sequence
    seqtk subseq ./${scaffold} ./${seqID}_filtered_headers.txt > ${seqID}_viral_sequences.fasta
    # add miner flag at each fasta header
    sed 's/^>/>deepvirfinder_/' ${seqID}_viral_sequences.fasta > ${seqID}_correct_viral_sequences.fasta
    """
}

process phigaro {
    cache 'lenient'
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

process virfinder {
    cache 'lenient'
    tag "$assembler-$seqID"
    label 'med_res'
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
    # since virfinder outputs a list of viral scaffolds headers, we need to collect these and extract the related viral sequence only if they respect the pvalue threshold
    python $projectDir/bin/python/pvalue_virfinder.py ${seqID}
    # join the viral scaffold header with sequence
    seqtk subseq ./${scaffold} ./${seqID}_filtered_headers.txt > ${seqID}_viral_tmp_sequences.fasta
    # add miner flag at each fasta header
    sed 's/^>/>virfinder_/' ${seqID}_viral_tmp_sequences.fasta > ${seqID}_viral_sequences.fasta
    rm ${seqID}_viral_tmp_sequences.fasta
    """
}

ch_virsorter2 = ch_dbm_virsorter2.splitCsv().flatMap { it -> "${it[3]}" + "/${params.mod_virsorter2}/" }

process virsorter2 {
    cache 'lenient'
    tag "$assembler-$seqID"
    label 'med_res'
    publishDir "${params.outdir}/mining/virsorter2/${assembler}", mode: 'copy', pattern: "${seqID}_virsorter2"

    when:
    !params.skip_mining && !params.skip_virsorter2

    input:
    each file_virsorter2_db from ch_virsorter2
    tuple val(seqID), val(assembler), file(scaffold) from ch_megahit_virsorter2

    output:
    file("*")
    tuple val(seqID), val(assembler), val("virsorter2"), file("${seqID}_virsorter2/${seqID}-correct-final-viral-combined.fa") into ch_virsorter2_cdhit optional true

    script:
    def viromes = params.virome_dataset ? "2" : "1"
    """
    # run virsorter2
    virsorter run \
    -j ${task.cpus} \
    -d ${file_virsorter2_db} \
    --min-length 0 \
    -i ${scaffold} \
    --include-groups ${params.virsorter2_include_groups} \
    -w ${seqID}_virsorter2 \
    -l ${seqID} \
    --rm-tmpdir
    # add miner flag at each fasta header
    sed 's/^>/>virsorter2_/' ${seqID}_virsorter2/${seqID}-final-viral-combined.fa > ${seqID}_virsorter2/${seqID}-tmp-final-viral-combined.fa
    sed 's/||/ /' ${seqID}_virsorter2/${seqID}-tmp-final-viral-combined.fa > ${seqID}_virsorter2/${seqID}-correct-final-viral-combined.fa
    rm ${seqID}_virsorter2/${seqID}-tmp-final-viral-combined.fa
    """
}//${params.virsorter2_min-length}

/* STEP 5 - dereplication and reads mapping */

process cdhit {
    cache 'lenient'
    tag "$assembler"
    label 'big_res'
    publishDir "${params.outdir}/cd-hit/", mode: 'copy'
    publishDir "${params.outdir}/report/files/", mode: 'copy', pattern: "*_vOTUs_consensus.fasta"

    when:
    !params.skip_dereplication

    input:
    tuple val(seqID), val(assembler), val(miner), file(scaffolds) from Channel.empty().mix(ch_deepvirfinder_cdhit, ch_phigaro_cdhit, ch_vibrant_cdhit, ch_virsorter2_cdhit, ch_virfinder_cdhit).groupTuple(by: 1)

    output:
    file("*")
    tuple val(assembler), file("${assembler}_vOTUs_consensus.fasta") into (ch_cdhit_prodigal, ch_cdhit_checkV, ch_cdhit_bowtie2, ch_cdhit_chopper_c, ch_cdhit_taxonomytable_c, ch_cdhit_summary)
    tuple val(assembler), file("${assembler}_min_comp.tsv") into (ch_cdhit_votus_mincomp, ch_cdhit_taxonomytable_min)
    file("${assembler}_mined_viral_contigs_comp.tsv") into ch_cdhit_vc_mincomp

    script:
    """
    # append all the mined viral scaffolds
    cat ${scaffolds} > concat_${assembler}.fasta
    # set viral scaffolds reads header to be unique
    seqkit rename concat_${assembler}.fasta > concat_unique_${assembler}.fasta
    # CD-HIT-EST
    # copy all files from cdhit in bin. This workaround is required for max sequence length > 653560
    cp -r $projectDir/bin/cdhit/* .
    # launch cd-hit-est
    ./cd-hit-est \
    -T ${task.cpus} \
    -M 0 \
    -i concat_unique_${assembler}.fasta \
    -o derep95_${assembler}.fasta \
    -c ${params.cdhit_sequence_identity} \
    -aS ${params.cdhit_alignment_coverage} \
    -d 0 \
    -p 1 \
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

//ch_checkv = ch_dbm_checkv.splitCsv().flatMap { it -> "${it[4]}" + "/${params.mod_checkv}/" } to be added whith new db_manager.py

process checkV {
    cache 'lenient'
    tag "$assembler"
    label 'low_res'
    publishDir "${params.outdir}/cd-hit/checkV", mode: 'copy'

    when:
    !params.skip_dereplication && !params.skip_checkv

    input:
    //each file_checkv_db from ch_checkv to be added whith new db_manager.py
    tuple val(assembler), file(consensus) from ch_cdhit_checkV


    output:
    file("*")
    file("quality_summary.tsv") into ch_checkv_checkvtable

    script:
    """
    checkv end_to_end \
    -t ${task.cpus} \
    ${consensus} \
    -d ${params.dbPath}/checkv/${params.mod_checkv} \
    .
    """ // -d ${file_checkv_db} to be added with new db_manager.
    // If you don't want to set the environmet variable, you can still use the database through the -d parameter.
}

process prodigal {
    cache 'lenient'
    tag "$assembler"
    label 'med_res'
    publishDir "${params.outdir}/prodigal", mode: 'copy'
    publishDir "${params.outdir}/report/files/", mode: 'copy', pattern: "*_vOTUs_proteins.faa"
    publishDir "${params.outdir}/report/files/", mode: 'copy', pattern: "*_vOTUs_coords.gff"

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

//ch_diamond = ch_dbm_diamond.splitCsv().flatMap { it -> "${it[5]}" + "/${params.mod_vcontact2}" } to be added whith new db_manager.py
//ch_vcontact2_db = ch_dbm_vcontact2.splitCsv().flatMap { it -> "${it[6]}" + "/${params.mod_vcontact2}/" } to be added whith new db_manager.py
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
    label 'vcontact2'
    publishDir "${params.outdir}/taxonomy/vcontact2/${assembler}", mode: 'copy'
    publishDir "${params.outdir}/report/taxonomy/vcontact2", mode: 'copy', pattern: 'c1.ntw'
    publishDir "${params.outdir}/report/taxonomy/vcontact2", mode: 'copy', pattern: 'genome_by_genome_overview.csv'

    when:
    !params.skip_dereplication && !params.skip_viral_taxo && !params.skip_vcontact2

    input:
    each clusterONE from Channel.fromPath("$projectDir/bin/cluster_one-1.0.jar")
    tuple val(assembler), file(alignment), file(gene2genome) from ch_diamond_vcontact2

    output:
    file("*")
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
    label 'graphanalyzer'
    publishDir "${params.outdir}/taxonomy/vcontact2", mode: 'copy', pattern: 'results_vcontact2_*.csv'
    publishDir "${params.outdir}/report/taxonomy/vcontact2", mode: 'copy', pattern: 'results_vcontact2_*.csv'
    publishDir "${params.outdir}/taxonomy/vcontact2", mode: 'copy', pattern: '*.{xlxs,log}'
    publishDir "${params.outdir}/report/tables/", mode: 'copy', pattern: 'taxonomy_table.csv'
    publishDir "${params.outdir}/report/taxonomy/vcontact2/", mode: 'copy', pattern: 'single-views_*/*.html'

    when:
    !params.skip_dereplication && !params.skip_viral_taxo && !params.skip_vcontact2 && !params.skip_graphanalyzer

    input:
    each file_vcontact2_db from Channel.fromPath("${params.dbPath}/inphared/")
    tuple val(assembler), file(netfile), file(csvfile) from ch_vcontact2_extender

    output:
    file("*")
    file("single-views_*/*.html")
    tuple val(assembler), file("taxonomy_table.csv") into (ch_vcontact2_phylo, ch_vcontact2_taxonomytable,
    ch_vcontact2_heatmap, ch_vcontact2_alphadiversity, ch_vcontact2_betadiversity, ch_vcontact2_violin) optional true

    script:
    """
    python $projectDir/bin/python/graphanalyzer.py \
    --threads ${task.cpus} \
    --graph ${netfile} \
    --csv ${csvfile} \
    --metas ${file_vcontact2_db}/*data_excluding_refseq.tsv \
    --prefix vOTU \
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
    tag "Creating kraken2 table..."
    label 'low_res'
    publishDir "${params.outdir}/report/tables", mode: 'copy'
    
    when:
    !params.skip_dereplication && !params.skip_bacterial_taxo && !params.skip_kraken2 && !params.skip_kraken_files

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
    publishDir "${params.outdir}/report/plots", mode: 'copy', pattern: "*_mqc.png"
    publishDir "${params.outdir}/report/plots", mode: 'copy', pattern: "*_mqc.html"
    publishDir "${params.outdir}/report/tables", mode: 'copy', pattern: "*table.html"
    
    when:
    !params.skip_miner_comparison
    
    input:
    tuple val(assembler), file(vOTUs_min_comp) from ch_cdhit_votus_mincomp.collect().ifEmpty([])
    file(vc_min_comp) from ch_cdhit_vc_mincomp.collect().ifEmpty([])

    output:
    file("votus_miner_comparison_mqc.png") into ch_mincomp_multiqc_mincomp
    file("vc_distribution_plot_mqc.html")
    file("vc_miner_comparison_table.html") into ch_mincomp_multiqc_vctable
    
    script:
    """
    Rscript $projectDir/bin/Rscript/miner_comparison.R ${vc_min_comp} ${vOTUs_min_comp} ${params.minlen}
    """
}

process checkv_table {
    cache 'lenient'
    tag "Creating checkV table..."
    label 'low_res'
    publishDir "${params.outdir}/report/tables", mode: 'copy'

    when:
    !params.skip_dereplication && !params.skip_checkv && !params.skip_checkv_table

    input:
    file(quality_summary) from ch_checkv_checkvtable

    output:
    file("checkv_table.html") into ch_checkVtable_multiqc

    script:
    """
    Rscript $projectDir/bin/Rscript/checkv_table.R ${quality_summary}
    """
}

process summary {
    cache 'lenient'
    tag "$assembler"
    label 'low_res'
    publishDir "${params.outdir}/report/tables", mode: 'copy', pattern: "*.html"

    when:
    !params.skip_violin_plots

    input:
    tuple val(assembler), file(count_table) from ch_covtocounts2_summary.collect().ifEmpty([])
    tuple val(assembler), file(vOTUs_consensus) from ch_cdhit_summary.collect().ifEmpty([])
    path metadata from ch_metadata_summary

    output:
    file("vOTUs_summary.html")
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
    publishDir "${params.outdir}/report/files/phyloseq", mode: 'copy', pattern: "*.rds"

    when:
    params.metadata && !params.skip_dereplication && !params.skip_viral_taxo && !params.skip_vcontact2 && !params.skip_graphanalyzer && !params.skip_phylo_obj

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
    publishDir "${params.outdir}/report/files/vOTUs_consensus", mode: 'copy', pattern: "vOTU_*.fasta"
    publishDir "${params.outdir}/report/files/vOTUs_proteins", mode: 'copy', pattern: "vOTU_*.faa"
    publishDir "${params.outdir}/report/files/vOTUs_coords", mode: 'copy', pattern: "vOTU_*.gff"

    when:
    !params.skip_dereplication && !params.skip_file_chopper

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
    params.metadata && !params.skip_dereplication && !params.skip_viral_taxo && !params.skip_taxonomy_table

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

process alpha_diversity {
    cache 'lenient'
    tag "$assembler"
    label 'med_res'
    publishDir "${params.outdir}/report/plots/alpha_diversity", mode: 'copy', pattern: "*.html"

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

process beta_diversity {
    cache 'lenient'
    tag "$assembler"
    label 'med_res'
    publishDir "${params.outdir}/report/plots/beta_diversity", mode: 'copy', pattern: "*.html"

    when:
    params.metadata && !params.skip_beta_diversity

    input:
    tuple val(assembler), file(count_table) from ch_covtocounts2_betadiversity.collect()
    tuple val(assembler), file(taxonomy_table) from ch_vcontact2_betadiversity.collect()
    path metadata from ch_metadata_betadiversity

    output:
    file("beta*.html") optional true

    script:
    """
    Rscript $projectDir/bin/Rscript/beta_diversity.R ${count_table} ${taxonomy_table} ${metadata} ${params.beta_var}
    """
}

process heatmap {
    cache 'lenient'
    tag "$assembler"
    label 'med_res'
    publishDir "${params.outdir}/report/plots/heatmap", mode: 'copy', pattern: "*.html"

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

process violin_plots {
    cache 'lenient'
    tag "$assembler"
    label 'med_res'
    publishDir "${params.outdir}/report/plots/total_violin", mode: 'copy', pattern: "*.html"

    when:
    params.metadata && !params.skip_violin_plots

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
    file(votus_min_comp) from ch_mincomp_multiqc_mincomp.collect().ifEmpty([])
    file(vc_min_comp_table) from ch_mincomp_multiqc_vctable.collect().ifEmpty([])
    file(summary_table) from ch_summary_multiqc.collect().ifEmpty([])
    file(taxo_table) from ch_taxonomytable_multiqc.collect().ifEmpty([])
    file(fastp) from ch_fastp_multiqc.collect().ifEmpty([])
    file(kraken2) from ch_kraken2_multiqc.collect().ifEmpty([])
    file(krakenfiles) from ch_krakenfiles_multiqc.collect().ifEmpty([])
    file(checkV) from ch_checkVtable_multiqc.collect().ifEmpty([])
    file(metaquast) from ch_metaquast_multiqc.collect().ifEmpty([])

    output:
    file("MetaPhage_${params.run_name}_report.html")

    script:
    """
    python $projectDir/bin/python/spacer_taxo.py
    multiqc \
    --config $projectDir/bin/multiqc/multiqc_config.yaml \
    --filename "MetaPhage_report.html" \
    --exclude general_stats $projectDir/bin/multiqc/. -f \
    ${votus_min_comp} vc_miner_comparison_table_mqc.html vOTUs_summary_mqc.html taxonomy_table_mqc.html checkv_table_mqc.html ${fastp} ${kraken2} kraken_files_mqc.html ${metaquast}
    python $projectDir/bin/python/disable_toolbox.py MetaPhage_report.html ${params.run_name}
    """
}
