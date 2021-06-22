/* CONFIGURATION VARIABLES */
params.db_manager_reports = false

// Trimming 
params.skip_qtrimming = false
params.adapter_forward = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
params.adapter_reverse = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
params.mean_quality = 15
params.trimming_quality = 15
params.keep_phix = false
params.mod_phix = "phiX174" 
params.file_phix_alone = "-" 

// microbial taxonomy
params.skip_bacterial_taxo = false
params.skip_kraken2 = false 
params.mod_kraken2 = "miniBAV"
params.file_kraken2_db = "-"
params.skip_bracken = false 
params.bracken_read_length = 100
params.bracken_abundance_level = "S"                                                                                                                                    

// Assembly (default is metaspades)
params.multiple_alliners = false
params.skip_metaspades = false
params.skip_megahit = true // add a control to chose which assembler and deactivate the other (two_flag = params.one && params.two ? false : params.two)
params.skip_metaquast = false

// Phage-mining
params.skip_mining = false
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

// Dereplication
params.skip_dereplication = false
params.minlen = 1000

// Viral Taxonomy - vContact2
params.skip_viral_taxo = false
params.skip_vcontact2 = false // true while debugging the pipeline
params.mod_vcontact2 = "Jan2021"
params.file_vcontact2_db = "-"

// Report
params.skip_error_barplots = false
//params.beta_diversity = true
//params.miner_comparison = true
params.skip_taxonomy_table = false
//params.heatmap = true
params.skip_report = false

  
/* FILE INPUT */

if (params.readPaths) {         // declared in profile config
    if (params.singleEnd) {     // declared in profile config
        Channel.fromPath("${params.readPaths}/*_single.fastq.gz", checkIfExists: true)
            .map { file -> tuple(file.simpleName, file) }
            .ifEmpty { exit 1, "No input files supplied! Please check params.readPaths in your config file!" }
            .set { ch_reads_fastp } 
    } else {
        Channel.fromFilePairs("${params.readPaths}/*{R1,R2}*.fastq.gz", checkIfExists: true)
            .ifEmpty { exit 1, "No input files supplied! Please check params.readPaths in your config file!" }
            .set { ch_reads_fastp }
    }
} else
    error "No input files supplied! Please check params.readPaths in your config file!"

/* METADATA */

if (params.metadata) {
        Channel.fromPath("${params.metadata}/*_metadata.csv", checkIfExists: true)
            println "found metadata!"
            //.set { ch_metadata_betadiverity, ch_metadata_heatmap }
    } else {
        println "No metadata supplied! Metadata are mandatory for beta diversity and heatmaps plots!"
        println "If you are interested in these plots, simply provide a metadata .csv file using params.metadata in your config file!"
        println "For more informations, check the metadata paragraph on the wiki at _______________"     //provide wiki url
        params.betadiversity = false
        params.heatmap = false
}

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

    tag "Downloading missing databases..."

    output:
    file("file_phix_alone") into (ch_file_phix_alone)
    file("file_kraken2_db") into (ch_file_kraken2_db, ch_file_bracken_db)
    file("file_vibrant_db") into (ch_file_vibrant_db)
    file("file_phigaro_config") into (ch_file_phigaro_config)
    file("file_virsorter_db") into (ch_file_virsorter_db)
    file("file_vcontact2_db") into (ch_file_vcontact2_db, ch_file_extender_db)

    script:
    """
    python $workflow.projectDir/bin/db_manager.py \
    --mod_phix ${params.mod_phix} --keep_phix ${params.keep_phix} --file_phix_alone ${params.file_phix_alone} \
    --mod_kraken2 ${params.mod_kraken2} --skip_kraken2 ${params.skip_kraken2} --file_kraken2_db ${params.file_kraken2_db} \
    --mod_vibrant ${params.mod_vibrant} --skip_vibrant ${params.skip_vibrant} --file_vibrant_db ${params.file_vibrant_db} \
    --mod_phigaro ${params.mod_phigaro} --skip_phigaro ${params.skip_phigaro} --file_phigaro_config ${params.file_phigaro_config} \
    --mod_virsorter ${params.mod_virsorter} --skip_virsorter ${params.skip_virsorter} --file_virsorter_db ${params.file_virsorter_db} \
    --mod_vcontact2 ${params.mod_vcontact2} --skip_vcontact2 ${params.skip_vcontact2} --file_vcontact2_db ${params.file_vcontact2_db}
    """
}

/* STEP 1 - quality check and trimming */
process fastp {
    tag "$seqID"
    publishDir "${params.outdir}/fastp_qc/", mode: 'copy',
        saveAs: {filename -> filename.endsWith(".html") ? "$filename" : null}

    when:
    !params.skip_qtrimming

    input:
    tuple val(seqID), file(reads) from ch_reads_fastp

    output:
    tuple val(seqID), file("*_trimmed.fastq*") into (ch_fastp_phix)
    file("${seqID}_qc_report.html")
    file("${seqID}_fastp.json") into (ch_fastp_multiqc)

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
    -h ${seqID}_qc_report.html \
    -j ${seqID}_fastp.json
    """
}

if(!params.keep_phix) {
    process remove_phix {
        tag "$seqID"

        when:
        !params.skip_qtrimming

        input:
        file file_phix_alone from ch_file_phix_alone
        tuple val(seqID), file(reads) from ch_fastp_phix

        output:
        tuple val(seqID), file("*.fastq.gz") into (ch_trimm_kraken2, ch_trimm_metaspades, ch_trimm_megahit, ch_trimm_mapping, ch_trimm_derep)

        script:
        path_file_phix_alone = file("$workflow.projectDir/bin/.groovy_vars/${file_phix_alone}").text
        def inp = params.singleEnd ? "-U ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
        def check = params.singleEnd ? "" : "--check-read-2"
        """
        hts_SeqScreener \
        $inp \
        --seq $workflow.projectDir/${path_file_phix_alone} \
        $check \
        --fastq-output ${seqID}_dephixed \
        --force
        """
    }
}
else {
    ch_fastp_phix.into {ch_trimm_kraken2; ch_trimm_metaspades; ch_trimm_megahit; ch_trimm_mapping; ch_trimm_derep}
}

/* STEP 2 - taxonomy classification */
process kraken2 {
    tag "$seqID"
    publishDir "${params.outdir}/taxonomy/kraken2/", mode: 'copy'

    when:
    !params.skip_bacterial_taxo && !params.skip_kraken2

    input:
    file file_kraken2_db from ch_file_kraken2_db
    tuple val(seqID), file(reads) from ch_trimm_kraken2

    output:
    file("${seqID}_output.txt")
    tuple val(seqID), file("${seqID}_report.txt") into (ch_kraken2_bracken)
    tuple val(seqID), file("${seqID}_report.txt") into (ch_kraken2_krona)

    script:
    path_file_kraken2_db = file("$workflow.projectDir/bin/.groovy_vars/${file_kraken2_db}").text
    def inp = params.singleEnd ? "${reads}" :  "--paired ${reads[0]} ${reads[1]}"
    println "running"
    """
    kraken2 \
    --report-zero-counts \
    --threads ${task.cpus} \
    --db $workflow.projectDir/${path_file_kraken2_db} \
    --output ${seqID}_output.txt \
    --report ${seqID}_report.txt \
    $inp
    """
}

process bracken {
    tag "$seqID"
    publishDir "${params.outdir}/taxonomy/bracken/", mode: 'copy'

    when:
    !params.skip_bacterial_taxo && (!params.skip_kraken2 && !params.skip_bracken)

    input:
    file file_bracken_db from ch_file_bracken_db
    tuple val(seqID), file(report_kraken2) from ch_kraken2_bracken

    output:
    file("${seqID}_abundancies.txt")
    tuple val(seqID), file("${seqID}_report_bracken_species.txt") into (ch_bracken_krona)

    script:
    path_file_bracken_db = file("$workflow.projectDir/bin/.groovy_vars/${file_bracken_db}").text
    """
    bracken \
    -d $workflow.projectDir/${path_file_bracken_db} \
    -i ${report_kraken2} \
    -o ${seqID}_abundancies.txt \
    -r ${params.bracken_read_length} \
    -l ${params.bracken_abundance_level} 
    """
}

process krona {
    tag "$seqID"
    publishDir "${params.outdir}/taxonomy/krona/", mode: 'copy'

    when:
    !params.skip_bacterial_taxo && !params.skip_kraken2

    input:
    tuple val(seqID), file(report) from Channel.empty().mix(ch_kraken2_krona, ch_bracken_krona)

    output:
    file("${seqID}_*_krona_abundancies.html")

    script:
    if (params.skip_bracken == true && params.skip_kraken2 == false)
        """
        echo "bracken was not run, using kracken2 output to generate krona report!"
        python $workflow.projectDir/bin/kreport2krona.py \
        --report-file ${report} \
        --output to_krona.txt 
        ktImportText \
        to_krona.txt \
        -o ${seqID}_krak_krona_abundancies.html
        """
    else
        """
        python $workflow.projectDir/bin/kreport2krona.py \
        --report-file ${report} \
        --output to_krona.txt 
        ktImportText \
        to_krona.txt \
        -o ${seqID}_brack_krona_abundancies.html
        """
}

/* STEP 3 - assembly */
process metaspades {
    tag "$seqID"
    publishDir "${params.outdir}/assembly/metaspades/${seqID}", mode: 'copy'
    
    when:
    !params.singleEnd && !params.skip_metaspades

    input:
    tuple val(seqID), file(reads) from ch_trimm_metaspades

    output:
    tuple val(seqID), val("metaspades"), file("${seqID}_metaspades_scaffolds.fasta") into (ch_metaspades_quast, ch_metaspades_mapping, ch_metaspades_vibrant, ch_metaspades_phigaro, ch_metaspades_virsorter, ch_metaspades_virfinder, ch_metaspades_virfinderproc)
    tuple val(seqID), val("metaspades"), file("${seqID}_metaspades_contigs.fasta")

    script:
    def inp = params.singleEnd ? "" : "--pe1-1 ${reads[0]} --pe1-2 ${reads[1]}"
    """    
    spades.py \
    --meta \
    --threads ${task.cpus} \
    --memory ${task.memory.toGiga()} \
    $inp \
    -o ./
    mv scaffolds.fasta ${seqID}_metaspades_scaffolds.fasta
    mv contigs.fasta ${seqID}_metaspades_contigs.fasta
    """
}

process megahit {
    tag "$seqID"
    publishDir "${params.outdir}/assembly/megahit/${seqID}", mode: 'copy'

    when:
    !params.skip_megahit

    input:
    tuple val(seqID), file(reads) from ch_trimm_megahit

    output:
    tuple val(seqID), val("megahit"), file("${seqID}_megahit_contigs.fasta") into (ch_megahit_quast, ch_megahit_vibrant, ch_megahit_phigaro, ch_megahit_virsorter, ch_megahit_virfinder, ch_megahit_virfinderproc)
    
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
    tag "$assembler-$seqID"
    publishDir "${params.outdir}/assembly/${assembler}/metaquast", mode: 'copy'

    when:
    (!params.skip_metaspades && !params.skip_megahit) || !params.skip_metaquast 

    input:
    tuple val(seqID), val(assembler), file(scaffold) from Channel.empty().mix(ch_metaspades_quast, ch_megahit_quast)

    output:
    file("${seqID}") into (ch_metaquast_multiqc) // here ${seqID} is a folder 

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

/* STEP 6 - phage mining */

process vibrant_legacy {
    tag "$assembler-$seqID"
    publishDir "${params.outdir}/mining/vibrant/${assembler}", mode: 'copy', overwrite: true

    when:
    !params.skip_mining && !params.skip_vibrant

    input:
    file file_vibrant_db from ch_file_vibrant_db
    tuple val(seqID), val(assembler), file(scaffold) from Channel.empty().mix(ch_metaspades_vibrant, ch_megahit_vibrant)

    output:
    file("*")
    tuple val(seqID), val(assembler), val("vibrant"), file("**/*.phages_combined.fna") into (ch_vibrant_cdhit)

    script:
    path_file_vibrant_db = file("$workflow.projectDir/bin/.groovy_vars/${file_vibrant_db}").text
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
        mv VIBRANT_*/ ${seqID}_vibrant/
        rm -f temp*
        """
    else 
        """
        VIBRANT_run.py \
        -t ${task.cpus} \
        -i ${scaffold} \
        -d $workflow.projectDir/${path_file_vibrant_db}databases/ \
        -m $workflow.projectDir/${path_file_vibrant_db}files/
        mv VIBRANT_*/ ${seqID}_vibrant/
        rm -f temp*
        """
}

process phigaro {
    tag "$assembler-$seqID"
    publishDir "${params.outdir}/mining/phigaro/${assembler}", mode: 'copy', overwrite: true

    when:
    !params.skip_mining && !params.skip_phigaro

    input:
    file file_phigaro_config from ch_file_phigaro_config // this acts just like a timer
    tuple val(seqID), val(assembler), file(scaffold) from Channel.empty().mix(ch_metaspades_phigaro, ch_megahit_phigaro)

    output:
    file("*")
    tuple val(seqID), val(assembler), val("phigaro"), file("**/*.phigaro.fasta") into (ch_phigaro_cdhit)

    script:
    path_file_phigaro_config = file("$workflow.projectDir/bin/.groovy_vars/${file_phigaro_config}").text
    """
    python $workflow.projectDir/bin/phigaro_config_creator.py
    printf 'Y\n' | phigaro \
    --threads ${task.cpus} \
    --fasta-file ${scaffold} \
    --config config.yml \
    --print-vogs \
    --output ${seqID}_phigaro \
    --not-open \
    --save-fasta \
    --mode basic
    """
}    

process virsorter_legacy {
    tag "$assembler-$seqID"
    publishDir "${params.outdir}/mining/virsorter/${assembler}", mode: 'copy', overwrite: true

    when:
    !params.skip_mining && !params.skip_virsorter

    input:
    file file_virsorter_db from ch_file_virsorter_db 
    tuple val(seqID), val(assembler), file(scaffold) from Channel.empty().mix(ch_metaspades_virsorter, ch_megahit_virsorter)

    output:
    file("*")
    tuple val(seqID), val(assembler), val("1virsorter"), file("**/${seqID}_VIRSorter_cat-1.fasta") into (ch_1virsorter_cdhit)
    tuple val(seqID), val(assembler), val("4virsorter"), file("**/${seqID}_VIRSorter_prophages_cat-4.fasta") into (ch_4virsorter_cdhit)

    script:
    path_file_virsorter_db = file("$workflow.projectDir/bin/.groovy_vars/${file_virsorter_db}").text
    def viromes = params.virsorter_viromes ? "2" : "1"
    if (params.mod_virsorter == "legacy")
        """
        wrapper_phage_contigs_sorter_iPlant.pl \
        -f ${scaffold} \
        --db $viromes \
        --wdir ${seqID}_virsorter \
        --ncpu ${task.cpus} \
        --data-dir $workflow.projectDir/${path_file_virsorter_db}
        cd ${seqID}_virsorter/Predicted_viral_sequences/
        for FILENAME in *; do mv \$FILENAME ${seqID}_\$FILENAME; done
        mv * ../
        cd ../
        rmdir Predicted_viral_sequences/
        """
    else 
        """
        echo $workflow.projectDir/${path_file_virsorter_db}
        """
}

process virfinder {
    tag "$assembler-$seqID"
    publishDir "${params.outdir}/mining/virfinder/${assembler}/${seqID}_virfinder", mode: 'copy', overwrite: true

    when:
    !params.skip_mining && !params.skip_virfinder

    input:
    tuple val(seqID), val(assembler), file(scaffold) from Channel.empty().mix(ch_metaspades_virfinder, ch_megahit_virfinder)

    output:
    file("*")
    tuple val(seqID), val(assembler), val("virfinder"), file("*_viral_sequences.fasta") into (ch_virfinder_cdhit)

    script:
    """
    Rscript $workflow.projectDir/bin/Rscript/virfinder_execute.R ${scaffold} ${task.cpus} $workflow.projectDir
    mv results.txt ${seqID}_results.txt
    mv viral_sequences.fasta ${seqID}_viral_sequences.fasta
    """
}

/* STEP 7 - dereplication and reads mapping */

process cdhit {
    tag "$assembler"
    publishDir "${params.outdir}/CD-HIT/", mode: 'copy', overwrite: true

    when:
    !params.skip_dereplication

    input:
    tuple val(seqID), val(assembler), val(miner), file(scaffolds) from Channel.empty().mix(ch_vibrant_cdhit, ch_phigaro_cdhit, ch_1virsorter_cdhit, ch_4virsorter_cdhit, ch_virfinder_cdhit).groupTuple(by: 1)

    output:
    file("*")
    file("renamed_filtered_derep95_${assembler}.fasta") into (ch_cdhit_bowtie2)
    tuple val(assembler), file("renamed_filtered_derep95_${assembler}.fasta") into (ch_cdhit_prodigal)

    script:
    if (params.skip_metaspades == false && params.skip_megahit == false)
        error "Dereplication works with one assembler at a time!"
    else
        """
        cat ${scaffolds} > concat_${assembler}.fasta
        cd-hit-est \
        -T ${task.cpus} \
        -M ${task.memory.toMega()} \
        -i concat_${assembler}.fasta \
        -o derep95_${assembler}.fasta \
        -c 0.95 \
        -aS 0.85 \
        -n 9 \
        -d 0 \
        -p 1 \
        -g 1
        seqkit seq \
        --min-len ${params.minlen} \
        --out-file filtered_derep95_${assembler}.fasta \
        derep95_${assembler}.fasta
        python << END
        from Bio import SeqIO
        from Bio.SeqRecord import SeqRecord
        corresp = open('corresp.txt', 'w')
        renamed = []
        counter = 1
        for seqrec in SeqIO.parse('filtered_derep95_${assembler}.fasta', 'fasta'):
            newid = 'VCS_' + str(counter) + '_length_' + str(len(seqrec))
            corresp.write(seqrec.id + ' '+seqrec.description +'\\n' + newid + '\\n\\n')
            renamed.append(SeqRecord(seqrec.seq, id=newid, description=''))
            counter += 1
        SeqIO.write(renamed, 'renamed_filtered_derep95_${assembler}.fasta', 'fasta')
        corresp.close()
        END
        """
}

process prodigal { // prokka is better! or even balrog
    tag "$assembler"
    publishDir "${params.outdir}/CD-HIT/", mode: 'copy', overwrite: true

    when:
    !params.skip_dereplication

    input:
    tuple val(assembler), file(vcs) from ch_cdhit_prodigal

    output:
    file("*")
    tuple val(assembler), file("derep_prots_${assembler}.faa") into ch_prodigal_diamond

    script:
    """
    prodigal \
    -i ${vcs} \
    -o derep_coords_${assembler}.gff \
    -a derep_prots_${assembler}.faa \
    -f gff
    """
}

process bowtie2_derep {
    tag "${seqID}"
    publishDir "${params.outdir}/bowtie2", mode: 'copy', overwrite: true

    when:
    !params.skip_dereplication

    input:
    file consensus from ch_cdhit_bowtie2.collect()
    tuple val(seqID), file(reads) from ch_trimm_derep

    output:
    file("*")
    file("*.sorted.bam") into (ch_bowtie2bam_covtocounts2)
    file("*.sorted.bam.bai") into (ch_bowtie2bai_covtocounts2)

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
    publishDir "${params.outdir}/bowtie2", mode: 'copy', overwrite: true

    when:
    !params.skip_report

    input:
    file(sortedbam) from ch_bowtie2bam_covtocounts2.collect()
    file(sortedbambai) from ch_bowtie2bai_covtocounts2.collect()

    output:
    file("*")
    tuple file("custom_count_table_mqc.txt"), file("custom_count_plot_mqc.txt") into (ch_covtocounts2_multiqc)
    file("variance_table.txt") into ch_covtocounts2_barplots

    script:
    """
    bamcountrefs \
    --threads ${task.cpus} \
    --multiqc \
    ${sortedbam} > multiqc_model.txt
    cp multiqc_model.txt variance_table.txt
    sed -i 1,3d variance_table.txt
    python $workflow.projectDir/bin/multiqc_model_editor.py
    """
}

params.allVSall = "$workflow.projectDir/db/vcontact2/allVSall.csv"
allVSall = file(params.allVSall)

/* STEP 8 - viral taxonomy */

process diamond_vcontact2 {
    tag "$assembler"
    publishDir "${params.outdir}/taxonomy/vcontact2/${assembler}", mode: 'copy', overwrite: true

    when:
    !params.skip_dereplication && !params.skip_viral_taxo && !params.skip_vcontact2
    
    input:
    file file_vcontact2_db from ch_file_vcontact2_db
    tuple val(assembler), file(viral_orfs) from ch_prodigal_diamond

    output:
    file("*")
    tuple val(assembler), file("allVSall.csv"), file("viral_genomes_g2g.csv") into ch_diamond_vcontact2

    script:
    path_file_vcontact2_db = file("$workflow.projectDir/bin/.groovy_vars/${file_vcontact2_db}").text
    reference_db_faa = file("$workflow.projectDir/${path_file_vcontact2_db}vConTACT2_proteins.faa")
    reference_db = file("$workflow.projectDir/db/vcontact2/reference_db.dmnd")
    if(allVSall.exists() == false){
        """
        # Format reference_db
        diamond makedb \
        -p ${task.cpus} \
        --in ${reference_db_faa} \
        -d ${reference_db}
        # Align the reference_db against itself to create allVSall
        diamond blastp \
        -p ${task.cpus} \
        --sensitive \
        -d ${reference_db} \
        -q ${reference_db_faa} \
        -o ${allVSall}
        # Format the viral_orf in the correct way
        $workflow.projectDir/bin/simplify_faa-ffn_derep.py ${viral_orfs}
        # create gene2genome file
        $workflow.projectDir/bin/vcontact2_gene2genome.py \
        -p ${viral_orfs}.simple.faa \
        -o viral_genomes_g2g.csv \
        -s 'Prodigal-FAA'
        # Align the viral orfs against the reference_db
        diamond blastp \
        -p ${task.cpus} \
        --sensitive \
        -d ${reference_db} \
        -q ${viral_orfs}.simple.faa \
        -o viral_orfs_alignment.csv
        # Append the viral orfs alignment to the allVSall file
        cat viral_orfs_alignment.csv >> ${allVSall}
        cp ${allVSall} ./allVSall.csv
        """
    }
    else {
        """
        # Format reference_db
        diamond makedb \
        -p ${task.cpus} \
        --in ${reference_db_faa} \
        -d ${reference_db}
        # Format the viral_orf in the correct way
        $workflow.projectDir/bin/simplify_faa-ffn_derep.py ${viral_orfs}
        # create gene2genome file
        $workflow.projectDir/bin/vcontact2_gene2genome.py \
        -p ${viral_orfs}.simple.faa \
        -o viral_genomes_g2g.csv \
        -s 'Prodigal-FAA'
        # Align the viral orfs against the reference_db
        diamond blastp \
        -p ${task.cpus} \
        --sensitive \
        -d ${reference_db} \
        -q ${viral_orfs}.simple.faa \
        -o viral_orfs_alignment.csv
        # Append the viral orfs alignment to the allVSall file
        cat viral_orfs_alignment.csv >> ${allVSall}
        cp ${allVSall} ./allVSall.csv
        """
    }
}

params.clusterONE = "$workflow.projectDir/bin/cluster_one-1.0.jar"
clusterONE = file(params.clusterONE)

process vcontact2 {
    tag "$assembler"
    publishDir "${params.outdir}/taxonomy/vcontact2/${assembler}", mode: 'copy', overwrite: true

    when:
    !params.skip_dereplication && !params.skip_viral_taxo && !params.skip_vcontact2

    input:
    tuple val(assembler), file(total_alignment), file(gene2genome) from ch_diamond_vcontact2

    output:
    file("*")
    tuple val(assembler), file("c1.ntw"), file("genome_by_genome_overview.csv") into ch_vcontact2_extender

    script:
    """
    # run vConTACT2
    vcontact \
    -t ${task.cpus} \
    --blast-fp ${total_alignment} \
    --rel-mode 'Diamond' \
    --proteins-fp ${gene2genome} \
    --db None \
    --pcs-mode MCL \
    --vcs-mode ClusterONE \
    --c1-bin ${clusterONE} \
    --output-dir ./
    """
}

process vcontact2_extender {
    tag "$assembler"
    publishDir "${params.outdir}/report", mode: 'copy', overwrite: true

    when:
    !params.skip_dereplication && !params.skip_viral_taxo && !params.skip_vcontact2

    input:
    file file_vcontact2_db from ch_file_extender_db 
    tuple val(assembler), file(netfile), file(csvfile) from ch_vcontact2_extender
    
    output:
    file("*")
    file("custom_taxonomy_table_mqc.txt") into ch_vcontact2_taxonomytable
    file("custom_graph_plot_mqc.html") into ch_vcontact2_multiqc

    script:
    path_file_vcontact2_db = file("$workflow.projectDir/bin/.groovy_vars/${file_vcontact2_db}").text
    """
    python $workflow.projectDir/bin/graph_analyzer.py \
    --input-graph ${netfile} \
    --input-csv ${csvfile} \
    --input-metas $workflow.projectDir/${path_file_vcontact2_db}data_excluding_refseq.tsv \
    --output ./ \
    --suffix ${assembler}
    """
}

/* STEP 9 - statistics and report generation */

process error_barplots {
    tag "all"
    publishDir "${params.outdir}/report/plots", mode: 'copy', overwrite: true

    when:
    !params.skip_report && !params.skip_error_barplots

    input:
    file(variance_table) from ch_covtocounts2_barplots

    output: 
    file("*")
    file("error_barplots.html") into ch_errorbarplots_multiqc

    script:
    """
    Rscript $workflow.projectDir/bin/Rscript/error_barplots.R ${variance_table}
    """
}

process taxonomy_table {
    tag "all"
    publishDir "${params.outdir}/report/plots", mode: 'copy', overwrite: true

    when:
    !params.skip_report && !params.skip_taxonomy_table

    input:
    file(taxonomy_table_vcontact2) from ch_vcontact2_taxonomytable

    output: 
    file("*")
    file("taxonomy_table_mqc.html") into ch_taxonomytable_multiqc

    script:
    """
    Rscript $workflow.projectDir/bin/Rscript/taxonomy_table.R ${taxonomy_table_vcontact2}
    """
}

// process taxonomy_table_debug {
//     tag "all"
//     publishDir "${params.outdir}/report/plots", mode: 'copy', overwrite: true

//     when:
//     !params.skip_report && !params.skip_taxonomy_table

//     input:
//     file(taxonomy_table_vcontact2) from Channel.fromPath("${params.taxoTable}", checkIfExists: true)

//      output: 
//      file("*")
//      file("taxonomy_table_mqc.html") into ch_taxonomytable_multiqc

//     script:
//     """
//     Rscript $workflow.projectDir/bin/Rscript/taxonomy_table.R ${taxonomy_table_vcontact2}
//     """
// }

/* STEP 10 summary report */

process multiqc {
    tag "all"
    publishDir "${params.outdir}/report/", mode: 'copy', overwrite: true 

    when:
    !params.skip_report

    input:
    file(fastp) from ch_fastp_multiqc.collect().ifEmpty([])
    file(metaquast) from ch_metaquast_multiqc.collect().ifEmpty([])
    tuple file(custom_count_table), file(custom_count_plot) from ch_covtocounts2_multiqc.collect().ifEmpty([])
    file(error_barplot) from ch_errorbarplots_multiqc.collect().ifEmpty([])
    file(taxonomy_table) from ch_taxonomytable_multiqc.collect().ifEmpty([])
    file(graph_plot) from ch_vcontact2_multiqc.collect().ifEmpty([])

    output:
    file("*")

    script:
    """
    multiqc \
    --config $workflow.projectDir/bin/multiqc/multiqc_config.yaml \
    --filename "MetaPhage_report.html" \
    --exclude general_stats $workflow.projectDir/bin/multiqc/. -f \
    . ${error_barplot} ${custom_count_plot} ${custom_count_table} ${taxonomy_table} #${graph_plot}
    $workflow.projectDir/bin/multiqc/remove_footer.py \
    MetaPhage_report.html MetaPhage_report.html
    """
}
