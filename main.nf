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
params.skip_megahit = false
params.skip_quast = false

// Mapping
params.skip_mapping = false

// Binning
params.skip_binning = false
params.skip_metabat2 = false
params.skip_maxbin2 = false
params.min_contig_size = 2000
params.skip_das_tool = false

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
params.skip_virfinderproc = false
params.skip_marvel = false                                            

// Dereplication
params.skip_dereplication = false

// Viral Taxonomy - vContact2
params.skip_viral_taxo = true
params.skip_vcontact2 = false // true while debugging the pipeline                                                                               

// Report
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

/* STEP 1 - quality check and trimming */
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
        conda "bioconda::htstream==1.3.3 conda-forge::boost==1.70.0"

        tag "$seqID"

        input:
        file file_phix_alone from ch_file_phix_alone
        tuple val(seqID), file(reads) from ch_fastp_phix

        output:
        tuple val(seqID), file("*.fastq.gz") into (ch_trimm_kraken2, ch_trimm_metaspades, ch_trimm_megahit, ch_trimm_mapping, ch_trimm_derep)

        script:
        path_file_phix_alone = file("$workflow.projectDir/bin/groovy_vars/${file_phix_alone}").text
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
    conda "bioconda::kraken2==2.1.0 conda-forge::llvm-openmp==11.0.0"
    
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

    script:
    path_file_kraken2_db = file("$workflow.projectDir/bin/groovy_vars/${file_kraken2_db}").text
    def inp = params.singleEnd ? "${reads}" :  "--paired ${reads[0]} ${reads[1]}"
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
    conda "bioconda::bracken==2.5.3 conda-forge::libcxx==9.0.1 conda-forge::llvm-openmp==10.0.1 conda-forge::python=3.7 conda-forge::python_abi==3.7=1_cp37m"

    tag "$seqID"
    publishDir "${params.outdir}/taxonomy/bracken/", mode: 'copy'

    when:
    !params.skip_bacterial_taxo && (!params.skip_kraken2 && !params.skip_bracken)

    input:
    file file_bracken_db from ch_file_bracken_db
    tuple val(seqID), file(report) from ch_kraken2_bracken

    output:
    file("${seqID}_abundancies.txt")
    tuple val(seqID), file("${seqID}_report_bracken_species.txt") into (ch_bracken_krona)

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
    !params.skip_bacterial_taxo && (!params.skip_kraken2 && !params.skip_bracken)

    input:
    tuple val(seqID), file(report_bracken) from ch_bracken_krona

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
    conda "bioconda::megahit==1.2.9"

    tag "$seqID"
    publishDir "${params.outdir}/assembly/megahit/${seqID}", mode: 'copy'

    when:
    params.multiple_alliners && !params.skip_megahit

    input:
    tuple val(seqID), file(reads) from ch_trimm_megahit

    output:
    tuple val(seqID), val("megahit"), file("${seqID}_megahit_contigs.fasta") into (ch_megahit_quast, ch_megahit_mapping, ch_megahit_vibrant, ch_megahit_phigaro, ch_megahit_virsorter, ch_megahit_virfinder, ch_megahit_virfinderproc)
    
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

process quast {
    conda "bioconda::quast==5.0.2"

    tag "$assembler-$seqID"
    publishDir "${params.outdir}/assembly/${assembler}/quast/${seqID}", mode: 'copy'

    when:
    (!params.skip_metaspades && !params.skip_megahit) && !params.skip_quast 

    input:
    tuple val(seqID), val(assembler), file(scaffold) from Channel.empty().mix(ch_metaspades_quast, ch_megahit_quast)

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

/* STEP 4 - mapping */ 

process bowtie2_mapping {
    conda "bioconda::samtools==1.11 bioconda::bowtie2==2.4.2"
    
    tag "$seqID-$assembler"
    publishDir "${params.outdir}/mapping/${assembler}", mode: 'copy',
        saveAs: {filename ->
                    if(filename.endsWith(".fasta")) null
                    else if(filename.endsWith(".fastq.gz")) null
                    else filename
                }

    input:
    tuple val(seqID), val(assembler), file(assembly), val(mapID), file(mapReads) from Channel.empty().mix(ch_metaspades_mapping, ch_megahit_mapping).combine(ch_trimm_mapping).flatMap{ tup -> if(tup[0] == tup[3]){ [tup] } }

    output:
    tuple val(seqID), val(assembler), file(assembly), file("*.bam") into ch_bowtie2_metabat2
    tuple val(seqID), val(assembler), file(assembly), file(mapReads) into ch_bowtie2_maxbin2
    file(assembly) into ch_bowtie2_das_tool
    tuple val(seqID), file("*.bam"), file("*.bam.bai")

    when:
    !params.skip_mapping

    script:
    def name = "${assembler}_${seqID}_${mapID}"
    def inp = params.singleEnd ? "-U ${mapReads}" : "-1 ${mapReads[0]} -2 ${mapReads[1]}"
    """
    bowtie2-build \
        --threads ${task.cpus} ${assembly} ref

    bowtie2 \
        -p ${task.cpus} \
        -x ref $inp | \
    samtools view -@ ${task.cpus} -bS | \
    samtools sort -@ ${task.cpus} -o "${name}.bam"

    samtools index "${name}.bam"
    """
}

/* Step 5 - binning */ 

process metabat2 {
    conda "bioconda::metabat2==2.15"

    tag "$seqID-$assembler"
    publishDir "${params.outdir}/binning/metabat2/${assembler}", mode: 'copy'

    when:
    !params.skip_binning && !params.skip_metabat2

    input:
    tuple val(seqID), val(assembler), file(assembly), file(bam) from ch_bowtie2_metabat2
    val(min_size) from params.min_contig_size

    output:
    tuple val("metabat2"), file("${seqID}_${assembler}_metabat2_scaffolds2bin.tsv") into ch_metabat2_das_tool
    file("*.fa")
    file("*.txt")

    script:
    def name = "${seqID}_${assembler}"
    """
    jgi_summarize_bam_contig_depths \
    --outputDepth ${name}.depth.txt ${bam}

    metabat2 \
    -t ${task.cpus} \
    -i ${assembly} \
    -o "${name}.metabat" \
    -m ${min_size} \
    -v --unbinned
    
    $workflow.projectDir/bin/Fasta_to_Scaffolds2Bin.sh -e fasta > ${name}_metabat2_scaffolds2bin.tsv
    """
}

process maxbin2 {
    conda "bioconda::maxbin2==2.2.7"

    tag "$seqID-$assembler"
    publishDir "${params.outdir}/binning/maxbin2/${assembler}", mode: 'copy'

    when:
    !params.skip_binning && !params.skip_maxbin2

    input:
    tuple val(seqID), val(assembler), file(assembly), file(mapReads) from ch_bowtie2_maxbin2
    val(min_size) from params.min_contig_size
    
    output:
    tuple val("maxbin2"), file("${seqID}_${assembler}_maxbin2_scaffolds2bin.tsv") into ch_maxbin2_das_tool
    file("*")

    script:
    def name = "${seqID}_${assembler}"
    """
    gunzip -dc ${seqID}_dephixed_R1.fastq.gz > ${seqID}_R1.fastq
    gunzip -dc ${seqID}_dephixed_R2.fastq.gz > ${seqID}_R2.fastq

    run_MaxBin.pl \
    -thread ${task.cpus} \
    -min_contig_length ${min_size} \
    -contig ${assembly} \
    -reads ${seqID}_R1.fastq \
    -reads2 ${seqID}_R2.fastq \
    -out ${name}.fasta

    $workflow.projectDir/bin/Fasta_to_Scaffolds2Bin.sh -e fasta > ${name}_maxbin2_scaffolds2bin.tsv
    """    
}

process das_tool {
    conda "bioconda::diamond==2.0.6"
    if (cursystem.contains('Mac')) {
        conda "bioconda::das_tool==1.1.2"
    }
    else { // das_tool on Linux has slightly different dependencies
        conda "bioconda::das_tool==1.1.1"
    }

    tag "$binner"
    publishDir "${params.outdir}/binning/das_tool/", mode: 'copy'

    when:
    !params.skip_binning && !params.skip_das_tool

    input:
    tuple val(binner), file(metabat2_bins) from ch_metabat2_das_tool
    tuple val(binner), file(maxbin2_bins) from ch_maxbin2_das_tool
    file(assembly) from ch_bowtie2_das_tool

    output:
    tuple val(binner), file("*,contigs.fa") into ch_das_tool_marvel
    file("*")

    script:
    """
    DAS_Tool \
    -i ${maxbin2_bins},${metabat2_bins} \
    -l MAXBIN2,METABAT2 \
    -c ${assembly} \
    -o ./ \
    -t ${task.cpus} \
    --search_engine blast \
    --write_bins 1
    """
}

/* STEP 6 - phage mining */

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
    !params.skip_mining && !params.skip_vibrant

    input:
    file file_vibrant_db from ch_file_vibrant_db
    tuple val(seqID), val(assembler), file(scaffold) from Channel.empty().mix(ch_metaspades_vibrant, ch_megahit_vibrant)

    output:
    file("*")
    tuple val(seqID), val(assembler), val("vibrant"), file("**/*.phages_combined.fna") into (ch_vibrant_cdhit)

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
    !params.skip_mining && !params.skip_phigaro

    input:
    file file_phigaro_config from ch_file_phigaro_config // this acts just like a timer
    tuple val(seqID), val(assembler), file(scaffold) from Channel.empty().mix(ch_metaspades_phigaro, ch_megahit_phigaro)

    output:
    file("*")
    tuple val(seqID), val(assembler), val("phigaro"), file("**/*.phigaro.fasta") into (ch_phigaro_cdhit)

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

    if [ ! -e phigaro_${seqID}/*.phigaro.fasta ]; then
        mkdir -p phigaro_${seqID}
        touch phigaro_${seqID}/empty.phigaro.fasta
    fi
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
    !params.skip_mining && !params.skip_virsorter

    input:
    file file_virsorter_db from ch_file_virsorter_db 
    tuple val(seqID), val(assembler), file(scaffold) from Channel.empty().mix(ch_metaspades_virsorter, ch_megahit_virsorter)

    output:
    file("*")
    tuple val(seqID), val(assembler), val("1virsorter"), file("**/${seqID}_VIRSorter_cat-1.fasta") into (ch_1virsorter_cdhit)
    tuple val(seqID), val(assembler), val("4virsorter"), file("**/${seqID}_VIRSorter_prophages_cat-4.fasta") into (ch_4virsorter_cdhit)

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

        cd ${seqID}_virsorter/Predicted_viral_sequences/
        for FILENAME in *; do mv \$FILENAME ${seqID}_\$FILENAME; done
        cd ../../
        """
    else 
        """
        echo $workflow.projectDir/${path_file_virsorter_db}
        """
}

process virfinder {
    conda "bioconda::r-virfinder==1.1"

    tag "$assembler-$seqID"
    publishDir "${params.outdir}/mining/virfinder/${assembler}/${seqID}", mode: 'copy'

    when:
    !params.skip_mining && !params.skip_virfinder

    input:
    tuple val(seqID), val(assembler), file(scaffold) from Channel.empty().mix(ch_metaspades_virfinder, ch_megahit_virfinder)

    output:
    file("*")
    tuple val(seqID), val(assembler), file("*_virus.csv") into (ch_virfinder_virfinderproc)

    script:
    """
    Rscript $workflow.projectDir/bin/virfinder_execute.R ${scaffold}

    mv results.txt ${seqID}_results.txt
    mv virus.csv ${seqID}_virus.csv
    """
}

process virfinderproc {
    conda "bioconda::r-seqinr==3.1_3"

    tag "$assembler-$seqID"
    publishDir "${params.outdir}/mining/virfinder/${assembler}/${seqID}", mode: 'copy'

    when:
    !params.skip_mining && !params.skip_virfinderproc

    input:
    tuple val(seqID), val(assembler), file(scaffold) from Channel.empty().mix(ch_metaspades_virfinderproc, ch_megahit_virfinderproc)
    tuple val(seqID), val(assembler), file(csvfile) from ch_virfinder_virfinderproc

    output:
    file("*")
    tuple val(seqID), val(assembler), val("virfinder"), file("*_viral_sequences.fasta") into (ch_virfinder_cdhit)

    script:
    """
    Rscript $workflow.projectDir/bin/virfinder_process.R ${scaffold} ${csvfile}

    mv viral_sequences.fasta ${seqID}_viral_sequences.fasta
    """
}

process marvel {
    conda "bioconda::marvel==0.2"

    tag "$binner"
    publishDir "${params.outdir}/mining/marvel/${binner}", mode: 'copy'

    when:
    !params.skip_mining && !params.skip_mapping && !params.skip_marvel

    input:
    tuple val(binner), file(bins) from ch_das_tool_marvel
    val(min_size) from params.min_contig_size

    output:
    file("*")
    
    script:
    """
    marvel_prokka  \
    -i $bins \
    -t ${task.cpus} \
    -o ./ \
    -m ${min_size} \
    """
}

/* STEP 7 - dereplication and reads mapping */

process cdhit {
    conda "bioconda::cd-hit==4.8.1 bioconda::seqkit==0.14.0"
    
    tag "$assembler"
    publishDir "${params.outdir}/CD-HIT/", mode: 'copy'

    when:
    !params.skip_dereplication

    input:
    tuple val(seqID), val(assembler), val(miner), file(scaffolds) from Channel.empty().mix(ch_vibrant_cdhit, ch_phigaro_cdhit, ch_1virsorter_cdhit, ch_4virsorter_cdhit, ch_virfinder_cdhit).groupTuple(by: 1)

    output:
    file("*")
    file("splitted83_${assembler}/*.fasta") into (ch_cdhit_bowtie2)
    tuple val(assembler), file("derep83_${assembler}.fasta") into (ch_cdhit_prodigal)

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
        -o derep83_${assembler}.fasta \
        -c 0.83 \
        -n 5 

        seqkit split derep83_${assembler}.fasta \
        --by-id \
        --force \
        --out-dir splitted83_${assembler}
        """
}

process prodigal { // prokka is better! or even balrog
    conda "bioconda::prodigal==2.6.3"
    
    tag "$assembler"
    publishDir "${params.outdir}/CD-HIT/", mode: 'copy'

    when:
    !params.skip_dereplication

    input:
    tuple val(assembler), file(vcs) from ch_cdhit_prodigal

    output:
    file("*")
    tuple val(assembler), file("derep_prots_${assembler}.faa") into (ch_prodigal_vcontact2)

    script:
    """
    prodigal -i ${vcs} -o derep_coords_${assembler}.gbk -a derep_prots_${assembler}.faa
    """
}

process bowtie2_derep {
    conda "bioconda::bowtie2==2.4.1 bioconda::samtools==1.9 bioconda::qualimap==2.2.2a=1"

    tag "${seqID}"
    publishDir "${params.outdir}/bowtie2", mode: 'copy'

    when:
    !params.skip_dereplication

    input:
    file consensus from ch_cdhit_bowtie2.collect()
    tuple val(seqID), file(reads) from ch_trimm_derep

    output:
    file("*")
    file("count___*") into ch_bowtie2_collector

    script:
    """
    for scaffold in ${consensus}
    do
        bowtie2-build --threads ${task.cpus} \$scaffold \${scaffold}_index

        bowtie2 -p ${task.cpus} -x \${scaffold}_index -1 ${reads[0]} -2 ${reads[1]} -S ${seqID}_\$scaffold.sam 

        samtools view -@ ${task.cpus} -S -b ${seqID}_\$scaffold.sam > ${seqID}_\$scaffold.bam

        samtools sort -@ ${task.cpus} ${seqID}_\$scaffold.bam -o ${seqID}_\$scaffold.sorted.bam

        samtools index -@ ${task.cpus} ${seqID}_\$scaffold.sorted.bam

        samtools flagstat -@ ${task.cpus} ${seqID}_\$scaffold.sorted.bam > mappingstats_${seqID}_\$scaffold.txt

        qualimap bamqc -nt ${task.cpus} -outdir qualimap_bamqc_${seqID}_\$scaffold.folder -bam ${seqID}_\$scaffold.sorted.bam
    
        samtools view -c -F 260 ${seqID}_\$scaffold.sorted.bam > count___${seqID}___\$scaffold.txt
    done
    """
}

process collector {
    conda "anaconda::pandas==1.1.3"

    tag "all"
    publishDir "${params.outdir}/report/", mode: 'copy'

    when:
    !params.skip_dereplication

    input:
    file values from ch_bowtie2_collector.collect()

    output:
    tuple file("custom_plot_mqc.yaml"), file("custom_table_mqc.txt") into ch_collector_multiqc

    script:
    """
    python $workflow.projectDir/bin/collector.py \
    --alignments ${values}
    """
}

/* STEP 8 - viral taxonomy */

process vcontact2 {
    conda "bioconda::vcontact2==0.9.19"

    tag "$assembler"
    publishDir "${params.outdir}/taxonomy/vcontact2/${assembler}", mode: 'copy'

    when:
    !params.skip_dereplication && !params.skip_viral_taxo && !params.skip_vcontact2

    input:
    tuple val(assembler), file(phages_combined) from ch_prodigal_vcontact2

    output:
    file("*")
    tuple val(assembler), file("c1.ntw"), file("genome_by_genome_overview.csv") into (ch_vcontact2_extender)

    script:
    """
    # $workflow.projectDir/bin/simplify_faa-ffn.py ${phages_combined}
    $workflow.projectDir/bin/simplify_faa-ffn_derep.py ${phages_combined}

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

process vcontact2_extender {
    conda "anaconda::python=3.7 conda-forge::pandas==1.1.4 anaconda::networkx==2.5 conda-forge::graphviz==2.42.3 conda-forge::pygraphviz==1.6 pyviz::hvplot==0.6.0 pyviz::panel==0.10.1 conda-forge::plotly==4.12.0"

    tag "$assembler"
    publishDir "${params.outdir}/report", mode: 'copy'

    input:
    tuple val(assembler), file(netfile), file(csvfile) from ch_vcontact2_extender
    //file(netfile) from Channel.fromPath( 'extra/c1.ntw' )
    //file(csvfile) from Channel.fromPath( 'extra/genome_by_genome_overview.csv' )

    output:
    file("*")

    script:
    """
    python $workflow.projectDir/bin/graph_analyzer.py \
    --input-graph ${netfile} \
    --input-csv ${csvfile} \
    --output ./ \
    --suffix ${assembler}
    """
}

/* STEP 9 - report generation */
process multiqc {
    conda "bioconda::multiqc==1.9=py_1 conda-forge::python==3.9.0"

    tag "all"
    publishDir "${params.outdir}/report/", mode: 'copy'

    when:
    !params.skip_report

    input:
    file("*_fastp.json") from ch_fastp_multiqc.collect().ifEmpty([])
    tuple file(custom_plot), file(custom_table) from ch_collector_multiqc

    output:
    file("*")

    script:
    """
    multiqc \
    --config $workflow.projectDir/bin/multiqc_config.yaml \
    --filename "MultiPhate_report.html" \
    . ${custom_plot} ${custom_table}
    """
}
