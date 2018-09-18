#!/usr/bin/env nextflow

import Helper

// Pipeline version
if (workflow.commitId){
    version = "0.1 $workflow.revision"
} else {
    version = "0.1 (local version)"
}

params.help = false
if (params.help){
    Help.print_help(params)
    exit 0
}

def infoMap = [:]
if (params.containsKey("fastq")){
    infoMap.put("fastq", file(params.fastq).size())
}
if (params.containsKey("fasta")){
    if (file(params.fasta) instanceof LinkedList){
        infoMap.put("fasta", file(params.fasta).size())
    } else {
        infoMap.put("fasta", 1) 
    }
}
if (params.containsKey("accessions")){
    // checks if params.accessions is different from null
    if (params.accessions) {
        BufferedReader reader = new BufferedReader(new FileReader(params.accessions));
        int lines = 0;
        while (reader.readLine() != null) lines++;
        reader.close();
        infoMap.put("accessions", lines)
    }
}

Help.start_info(infoMap, "$workflow.start", "$workflow.profile")
    

// Placeholder for main input channels
if (params.fastq instanceof Boolean){exit 1, "'fastq' must be a path pattern. Provide value:'$params.fastq'"}
if (!params.fastq){ exit 1, "'fastq' parameter missing"}
IN_fastq_raw = Channel.fromFilePairs(params.fastq).ifEmpty { exit 1, "No fastq files provided with pattern:'${params.fastq}'" }

// Placeholder for secondary input channels


// Placeholder for extra input channels


// Placeholder to fork the raw input channel

IN_fastq_raw.set{ integrity_coverage_in_1_0 }


IN_genome_size_1_1 = Channel.value(params.genomeSize_1_1)
    .map{it -> it.toString().isNumber() ? it : exit(1, "The genomeSize parameter must be a number or a float. Provided value: '${params.genomeSize__1_1}'")}

IN_min_coverage_1_1 = Channel.value(params.minCoverage_1_1)
    .map{it -> it.toString().isNumber() ? it : exit(1, "The minCoverage parameter must be a number or a float. Provided value: '${params.minCoverage__1_1}'")}

process integrity_coverage_1_1 {

    // Send POST request to platform
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_1 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_1 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_1 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId integrity_coverage_1_1 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { sample_id }
    // This process can only use a single CPU
    cpus 1

    input:
    set sample_id, file(fastq_pair) from integrity_coverage_in_1_0
    val gsize from IN_genome_size_1_1
    val cov from IN_min_coverage_1_1
    // This channel is for the custom options of the integrity_coverage.py
    // script. See the script's documentation for more information.
    val opts from Channel.value('')

    output:
    set sample_id,
        file(fastq_pair),
        file('*_encoding'),
        file('*_phred'),
        file('*_coverage'),
        file('*_max_len') into MAIN_integrity_1_1
    file('*_report') optional true into LOG_report_coverage1_1_1
    set sample_id, val("1_1_integrity_coverage"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_integrity_coverage_1_1
set sample_id, val("integrity_coverage_1_1"), val("1_1"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_integrity_coverage_1_1
file ".versions"

    script:
    template "integrity_coverage.py"

}

// TRIAGE OF CORRUPTED SAMPLES
LOG_corrupted_1_1 = Channel.create()
MAIN_PreCoverageCheck_1_1 = Channel.create()
// Corrupted samples have the 2nd value with 'corrupt'
MAIN_integrity_1_1.choice(LOG_corrupted_1_1, MAIN_PreCoverageCheck_1_1) {
    a -> a[2].text == "corrupt" ? 0 : 1
}

// TRIAGE OF LOW COVERAGE SAMPLES
integrity_coverage_out_1_0 = Channel.create()
SIDE_phred_1_1 = Channel.create()
SIDE_max_len_1_1 = Channel.create()

MAIN_PreCoverageCheck_1_1
// Low coverage samples have the 4th value of the Channel with 'fail'
    .filter{ it[4].text != "fail" }
// For the channel to proceed with FastQ in 'sample_good' and the
// Phred scores for each sample in 'SIDE_phred'
    .separate(integrity_coverage_out_1_0, SIDE_phred_1_1, SIDE_max_len_1_1){
        a -> [ [a[0], a[1]], [a[0], a[3].text], [a[0], a[5].text]  ]
    }

/** REPORT_COVERAGE - PLUG-IN
This process will report the expected coverage for each non-corrupted sample
and write the results to 'reports/coverage/estimated_coverage_initial.csv'
*/
process report_coverage_1_1 {

    // This process can only use a single CPU
    cpus 1
    publishDir 'reports/coverage_1_1/'

    input:
    file(report) from LOG_report_coverage1_1_1.filter{ it.text != "corrupt" }.collect()

    output:
    file 'estimated_coverage_initial.csv'

    """
    echo Sample,Estimated coverage,Status >> estimated_coverage_initial.csv
    cat $report >> estimated_coverage_initial.csv
    """
}

/** REPORT_CORRUPT - PLUG-IN
This process will report the corrupted samples and write the results to
'reports/corrupted/corrupted_samples.txt'
*/
process report_corrupt_1_1 {

    // This process can only use a single CPU
    cpus 1
    publishDir 'reports/corrupted_1_1/'

    input:
    val sample_id from LOG_corrupted_1_1.collect{it[0]}

    output:
    file 'corrupted_samples.txt'

    """
    echo ${sample_id.join(",")} | tr "," "\n" >> corrupted_samples.txt
    """

}


SIDE_phred_1_1.set{ SIDE_phred_1_2 }


SIDE_max_len_1_1.set{ SIDE_max_len_1_7 }


// Check sliding window parameter
if ( params.trimSlidingWindow_1_2.toString().split(":").size() != 2 ){
    exit 1, "'trimSlidingWindow_1_2' parameter must contain two values separated by a ':'. Provided value: '${params.trimSlidingWindow_1_2}'"
}
if ( !params.trimLeading_1_2.toString().isNumber() ){
    exit 1, "'trimLeading_1_2' parameter must be a number. Provide value: '${params.trimLeading_1_2}'"
}
if ( !params.trimTrailing_1_2.toString().isNumber() ){
    exit 1, "'trimTrailing_1_2' parameter must be a number. Provide value: '${params.trimTrailing_1_2}'"
}
if ( !params.trimMinLength_1_2.toString().isNumber() ){
    exit 1, "'trimMinLength_1_2' parameter must be a number. Provide value: '${params.trimMinLength_1_2}'"
}

IN_trimmomatic_opts_1_2 = Channel.value([params.trimSlidingWindow_1_2,params.trimLeading_1_2,params.trimTrailing_1_2,params.trimMinLength_1_2])
IN_adapters_1_2 = Channel.value(params.adapters_1_2)

clear = params.clearAtCheckpoint ? "true" : "false"
checkpointClear_1_2 = Channel.value(clear)

process fastqc_1_2 {

    // Send POST request to platform
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_2 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_2 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_2 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId fastqc_trimmomatic_1_2 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { sample_id }
    publishDir "reports/fastqc_1_2/", pattern: "*.html"

    input:
    set sample_id, file(fastq_pair) from integrity_coverage_out_1_0
    val ad from Channel.value('None')

    output:
    set sample_id, file(fastq_pair), file('pair_1*'), file('pair_2*') into MAIN_fastqc_out_1_2
    file "*html"
    set sample_id, val("1_2_fastqc"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_fastqc_1_2
set sample_id, val("fastqc_1_2"), val("1_2"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_fastqc_1_2
file ".versions"

    script:
    template "fastqc.py"
}

/** FASTQC_REPORT - MAIN
This process will parse the result files from a FastQC analyses and output
the optimal_trim information for Trimmomatic
*/
process fastqc_report_1_2 {

    // Send POST request to platform
    
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_2 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_2 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_2 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId fastqc_trimmomatic_1_2 \"$params.platformSpecies\" false"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }
    

    tag { sample_id }
    // This process can only use a single CPU
    cpus 1
    publishDir 'reports/fastqc_1_2/run_1/', pattern: '*summary.txt', mode: 'copy'

    input:
    set sample_id, file(fastq_pair), file(result_p1), file(result_p2) from MAIN_fastqc_out_1_2
    val opts from Channel.value("--ignore-tests")

    output:
    set sample_id, file(fastq_pair), 'optimal_trim', ".status" into _MAIN_fastqc_trim_1_2
    file '*_trim_report' into LOG_trim_1_2
    file "*_status_report" into LOG_fastqc_report_1_2
    file "${sample_id}_*_summary.txt" optional true
    set sample_id, val("1_2_fastqc_report"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_fastqc_report_1_2
set sample_id, val("fastqc_report_1_2"), val("1_2"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_fastqc_report_1_2
file ".versions"

    script:
    template "fastqc_report.py"

}

MAIN_fastqc_trim_1_2 = Channel.create()
_MAIN_fastqc_trim_1_2
        .filter{ it[3].text == "pass" }
        .map{ [it[0], it[1], file(it[2]).text] }
        .into(MAIN_fastqc_trim_1_2)


/** TRIM_REPORT - PLUG-IN
This will collect the optimal trim points assessed by the fastqc_report
process and write the results of all samples in a single csv file
*/
process trim_report_1_2 {

    publishDir 'reports/fastqc_1_2/', mode: 'copy'

    input:
    file trim from LOG_trim_1_2.collect()

    output:
    file "FastQC_trim_report.csv"

    """
    echo Sample,Trim begin, Trim end >> FastQC_trim_report.csv
    cat $trim >> FastQC_trim_report.csv
    """
}


process compile_fastqc_status_1_2 {

    publishDir 'reports/fastqc_1_2/', mode: 'copy'

    input:
    file rep from LOG_fastqc_report_1_2.collect()

    output:
    file 'FastQC_1run_report.csv'

    """
    echo Sample, Failed? >> FastQC_1run_report.csv
    cat $rep >> FastQC_1run_report.csv
    """

}


/** TRIMMOMATIC - MAIN
This process will execute trimmomatic. Currently, the main channel requires
information on the trim_range and phred score.
*/
process trimmomatic_1_2 {

    // Send POST request to platform
    
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_2 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_2 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_2 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId fastqc_trimmomatic_1_2 \"$params.platformSpecies\" false"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }
    

    tag { sample_id }
    publishDir "results/trimmomatic_1_2", pattern: "*.gz"

    input:
    set sample_id, file(fastq_pair), trim_range, phred from MAIN_fastqc_trim_1_2.join(SIDE_phred_1_2)
    val opts from IN_trimmomatic_opts_1_2
    val ad from IN_adapters_1_2
    val clear from checkpointClear_1_2

    output:
    set sample_id, "${sample_id}_*trim.fastq.gz" into fastqc_trimmomatic_out_1_1
    file 'trimmomatic_report.csv'
    set sample_id, val("1_2_trimmomatic"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_trimmomatic_1_2
set sample_id, val("trimmomatic_1_2"), val("1_2"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_trimmomatic_1_2
file ".versions"

    script:
    template "trimmomatic.py"

}



IN_adapter_1_3 = Channel.value(params.adapter_1_3)


process filter_poly_1_3 {

    // Send POST request to platform
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_3 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_3 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_3 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId filter_poly_1_3 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { sample_id }

    errorStrategy { task.exitStatus == 120 ? 'ignore' : 'retry' }

    input:
    set sample_id, file(fastq_pair) from fastqc_trimmomatic_out_1_1
    val adapter from IN_adapter_1_3

    output:
    set sample_id , file("${sample_id}_filtered_{1,2}.fastq.gz") into filter_poly_out_1_2
    set sample_id, val("1_3_filter_poly"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_filter_poly_1_3
set sample_id, val("filter_poly_1_3"), val("1_3"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_filter_poly_1_3
file ".versions"

    script:
    """
    gunzip -c ${fastq_pair[0]} >  ${sample_id}_1.fq
    gunzip -c ${fastq_pair[1]} >  ${sample_id}_2.fq

    for seqfile in *.fq;
    do if [ ! -s \$seqfile  ]
    then
        echo \$seqfile is empty && exit 120
    fi
    done

    prinseq-lite.pl --fastq ${sample_id}_1.fq  --fastq2 ${sample_id}_2.fq  --custom_params "${adapter}" -out_format 3 -out_good ${sample_id}_filtered

    gzip ${sample_id}_filtered_*.fastq

    rm *.fq *.fastq

    """
}



IN_index_files_1_4 = Channel.value(params.refIndex_1_4)


process remove_host_1_4 {

    // Send POST request to platform
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_4 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_4 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_4 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId remove_host_1_4 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { sample_id }
    publishDir 'results/mapping/remove_host_1_4/', pattern: '*_bowtie2.log', mode: 'copy'

    input:
    set sample_id, file(fastq_pair) from filter_poly_out_1_2
    val bowtie2Index from IN_index_files_1_4

    output:
    set sample_id , file("${sample_id}*.headersRenamed_*.fq.gz") into remove_host_out_1_3
    set sample_id, file("*_bowtie2.log") into into_json_1_4
    set sample_id, val("1_4_remove_host"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_remove_host_1_4
set sample_id, val("remove_host_1_4"), val("1_4"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_remove_host_1_4
file ".versions"

    script:
    """
    bowtie2 -x ${bowtie2Index} -1 ${fastq_pair[0]} -2 ${fastq_pair[1]} -p $task.cpus 1> ${sample_id}.bam 2> ${sample_id}_bowtie2.log

    samtools view -buh -f 12 -o ${sample_id}_samtools.bam -@ $task.cpus ${sample_id}.bam

    rm ${sample_id}.bam

    samtools fastq -1 ${sample_id}_unmapped_1.fq -2 ${sample_id}_unmapped_2.fq ${sample_id}_samtools.bam

    rm ${sample_id}_samtools.bam

    renamePE_samtoolsFASTQ.py -1 ${sample_id}_unmapped_1.fq -2 ${sample_id}_unmapped_2.fq

    gzip *.headersRenamed_*.fq

    rm *.fq
    """
}



process report_remove_host_1_4 {

        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_4 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_4 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_4 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId remove_host_1_4 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { sample_id }
    
    input:
    set sample_id, file(bowtie_log) from into_json_1_4

    output:
    set sample_id, val("1_4_report_remove_host"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_report_remove_host_1_4
set sample_id, val("report_remove_host_1_4"), val("1_4"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_report_remove_host_1_4
file ".versions"

    script:
    template "process_mapping.py"

}


// Check for the presence of absence of both index and fasta reference
if (params.index_1_5 == null && params.reference_1_5 == null){
    exit 1, "An index or a reference fasta file must be provided."
} else if (params.index_1_5 != null && params.reference_1_5 != null){
    exit 1, "Provide only an index OR a reference fasta file."
}


if (params.reference_1_5){

    reference_in_1_5 = Channel.fromPath(params.reference_1_5)
        .map{it -> file(it).exists() ? [it.toString().tokenize('/').last().tokenize('.')[0..-2].join('.') ,it] : null}
        .ifEmpty{ exit 1, "No fasta file was provided"}

    process bowtie_build_1_5 {

        // Send POST request to platform
            if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_5 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_5 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_5 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId bowtie_1_5 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

        tag { build_id }
        storeDir 'bowtie_index/'
        maxForks 1

        input:
        set build_id, file(fasta) from reference_in_1_5

        output:
        val build_id into bowtieIndexId_1_5
        file "${build_id}*.bt2" into bowtieIndex_1_5

        script:
        """
        bowtie2-build ${fasta} $build_id > ${build_id}_bowtie2_build.log
        """
    }
} else {
    bowtieIndexId_1_5 = Channel.value(params.index_1_5.split("/").last())
    bowtieIndex_1_5 = Channel.fromPath("${params.index_1_5}*.bt2").collect().toList()
}


process bowtie_1_5 {

    // Send POST request to platform
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_5 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_5 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_5 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId bowtie_1_5 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { sample_id }
    publishDir 'results/mapping/bowtie_1_5/'

    input:
    set sample_id, file(fastq_pair) from remove_host_out_1_3
    each index from bowtieIndexId_1_5
    each file(index_files) from bowtieIndex_1_5

    output:
    set sample_id , file("*.bam") into bowtie_out_1_4
    set sample_id, file("*_bowtie2.log") into into_json_1_5
    set sample_id, val("1_5_bowtie"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_bowtie_1_5
set sample_id, val("bowtie_1_5"), val("1_5"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_bowtie_1_5
file ".versions"

    script:
    """
    bowtie2 -x $index -1 ${fastq_pair[0]} -2 ${fastq_pair[1]} -p $task.cpus 1> ${sample_id}.bam 2> ${sample_id}_bowtie2.log
    """
}


process report_bowtie_1_5 {

        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_5 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_5 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_5 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId bowtie_1_5 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }
    tag { sample_id }
    
    input:
    set sample_id, file(bowtie_log) from into_json_1_5

    output:
    set sample_id, val("1_5_report_bowtie"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_report_bowtie_1_5
set sample_id, val("report_bowtie_1_5"), val("1_5"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_report_bowtie_1_5
file ".versions"

    script:
    template "process_mapping.py"

}


process retrieve_mapped_1_6 {

    // Send POST request to platform
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_6 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_6 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_6 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId retrieve_mapped_1_6 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { sample_id }
    publishDir 'results/mapping/retrieve_mapped_1_6/'

    input:
    set sample_id, file(bam) from bowtie_out_1_4

    output:
    set sample_id , file("*.headersRenamed_*.fq.gz") into _retrieve_mapped_out_1_5
    set sample_id, val("1_6_retrieve_mapped"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_retrieve_mapped_1_6
set sample_id, val("retrieve_mapped_1_6"), val("1_6"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_retrieve_mapped_1_6
file ".versions"

    script:
    """
    samtools view -buh -F 12 -o ${sample_id}_samtools.bam -@ $task.cpus ${bam}

    rm ${bam}

    samtools fastq -1 ${sample_id}_mapped_1.fq -2 ${sample_id}_mapped_2.fq ${sample_id}_samtools.bam

    rm ${sample_id}_samtools.bam

    renamePE_samtoolsFASTQ.py -1 ${sample_id}_mapped_1.fq -2 ${sample_id}_mapped_2.fq

    gzip *.headersRenamed_*.fq

    rm *.fq
    """
}


_retrieve_mapped_out_1_5.into{ retrieve_mapped_out_1_5;_LAST_fastq_1_8 }

//MAIN INPUT - FASTQ FILES
spades_in = Channel.create()
megahit_in = Channel.create()
retrieve_mapped_out_1_5.into{ spades_in; megahit_in }

//EXPECTED GENOME SIZE
if ( !params.minimumContigSize_1_7.toString().isNumber() ){
    exit 1, "'minimumContigSize_1_7' parameter must be a number. Provided value: '${params.minimumContigSize_1_7}'"
}

//SPADES OPTIONS
if ( !params.spadesMinCoverage_1_7.toString().isNumber() ){
    exit 1, "'spadesMinCoverage_1_7' parameter must be a number. Provided value: '${params.spadesMinCoverage_1_7}'"
}
if ( !params.spadesMinKmerCoverage_1_7.toString().isNumber()){
    exit 1, "'spadesMinKmerCoverage_1_7' parameter must be a number. Provided value: '${params.spadesMinKmerCoverage_1_7}'"
}

if ( params.spadesKmers_1_7.toString().split(" ").size() <= 1 ){
    if (params.spadesKmers_1_7.toString() != 'auto'){
        exit 1, "'spadesKmers_1_7' parameter must be a sequence of space separated numbers or 'auto'. Provided value: ${params.spadesKmers_1_7}"
    }
}

clear = params.clearAtCheckpoint ? "true" : "false"
checkpointClear_1_7 = Channel.value(clear)

//MEGAHIT OPTIONS
if ( params.megahitKmers_1_7.toString().split(" ").size() <= 1 ){
    if (params.megahitKmers_1_7.toString() != 'auto'){
        exit 1, "'megahitKmers_1_7' parameter must be a sequence of space separated numbers or 'auto'. Provided value: ${params.megahitKmers_1_7}"
    }
}

//SPADES INPUT CHANNELS
IN_spades_opts_1_7 = Channel.value([params.spadesMinCoverage_1_7,params.spadesMinKmerCoverage_1_7])
IN_spades_kmers_1_7 = Channel.value(params.spadesKmers_1_7)

//MEGAGIT INPUT CHANNELS
IN_megahit_kmers_1_7 = Channel.value(params.megahitKmers_1_7)

SIDE_max_len_spades = Channel.create()
SIDE_max_len_megahit = Channel.create()
SIDE_max_len_1_7.into{SIDE_max_len_spades ; SIDE_max_len_megahit}

process va_spades_1_7 {

    // Send POST request to platform
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_7 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_7 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_7 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId viral_assembly_1_7 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    validExitStatus 0,1

    tag { sample_id }
    publishDir 'results/assembly/spades_1_7/', pattern: '*_spades*.fasta', mode: 'copy'

    input:
    set sample_id, file(fastq_pair), max_len from spades_in.join(SIDE_max_len_spades)
    val opts from IN_spades_opts_1_7
    val kmers from IN_spades_kmers_1_7
    val clear from checkpointClear_1_7

    output:
    set sample_id, file({task.exitStatus == 1 ? ".exitcode" : '*_spades*.fasta'}) into assembly_spades
    set sample_id, val("1_7_va_spades"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_va_spades_1_7
set sample_id, val("va_spades_1_7"), val("1_7"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_va_spades_1_7
file ".versions"

    script:
    template "spades.py"

}

class VerifyCompletness {

    public static boolean contigs(String filename, int threshold){
        BufferedReader reader = new BufferedReader(new FileReader(filename));
        boolean result = processContigs(reader, threshold);
        reader.close()

        return result;
    }

    private static boolean processContigs(BufferedReader reader, int threshold){
        String line;
        int lineThreshold = 0;
        List splittedLine

        while ((line = reader.readLine()) != null) {
            if (line.startsWith('>')) {
                splittedLine = line.split('_')
                lineThreshold = splittedLine[3].toInteger()
                if(lineThreshold >= threshold) {
                    return true;
                }
             }
        }

        return false;
    }
}

megahit = Channel.create()
good_assembly = Channel.create()
assembly_spades.choice(good_assembly, megahit){a -> a[1].toString() == "null" ? false : VerifyCompletness.contigs(a[1].toString(), params.minimumContigSize_1_7.toInteger()) == true ? 0 : 1}


process va_megahit_1_7  {

    // Send POST request to platform
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_7 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_7 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_7 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId viral_assembly_1_7 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { sample_id }
    publishDir 'results/assembly/megahit_1_7/', pattern: '*_megahit*.fasta', mode: 'copy'

    input:
    set sample_id, file(fastq_pair), max_len from megahit_in.join(megahit).map{ ot -> [ot[0], ot[1]] }.join(SIDE_max_len_megahit)
    val kmers from IN_megahit_kmers_1_7

    output:
    set sample_id, file('*megahit*.fasta') into megahit_assembly
    set sample_id, val("1_7_va_megahit"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_va_megahit_1_7
set sample_id, val("va_megahit_1_7"), val("1_7"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_va_megahit_1_7
file ".versions"

    script:
    template "megahit.py"

}


good_assembly.mix(megahit_assembly).into{ to_report_1_7 ; viral_assembly_out_1_6 }
orf_size = Channel.value(params.minimumContigSize_1_7)


process report_viral_assembly_1_7 {

        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_7 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_7 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_7 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId viral_assembly_1_7 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }
    
    tag { sample_id }
    
    input:
    set sample_id, file(assembly) from to_report_1_7
    val min_size from orf_size

    output:
    set sample_id, val("1_7_report_viral_assembly"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_report_viral_assembly_1_7
set sample_id, val("report_viral_assembly_1_7"), val("1_7"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_report_viral_assembly_1_7
file ".versions"

    script:
    template "process_viral_assembly.py"

}



if ( !params.minAssemblyCoverage_1_8.toString().isNumber() ){
    if (params.minAssemblyCoverage_1_8.toString() != 'auto'){
        exit 1, "'minAssemblyCoverage_1_8' parameter must be a number or 'auto'. Provided value: ${params.minAssemblyCoverage_1_8}"
    }
}
if ( !params.AMaxContigs_1_8.toString().isNumber() ){
    exit 1, "'AMaxContigs_1_8' parameter must be a number. Provide value: '${params.AMaxContigs_1_8}'"
}

IN_assembly_mapping_opts_1_8 = Channel.value([params.minAssemblyCoverage_1_8,params.AMaxContigs_1_8])
IN_genome_size_1_8 = Channel.value(params.genomeSize_1_8)


process assembly_mapping_1_8 {

    // Send POST request to platform
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_8 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_8 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_8 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId assembly_mapping_1_8 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { sample_id }

    input:
    set sample_id, file(assembly), file(fastq) from viral_assembly_out_1_6.join(_LAST_fastq_1_8)

    output:
    set sample_id, file(assembly), 'coverages.tsv', 'coverage_per_bp.tsv', 'sorted.bam', 'sorted.bam.bai' into MAIN_am_out_1_8
    set sample_id, file("coverage_per_bp.tsv") optional true into SIDE_BpCoverage_1_8
    set sample_id, val("1_8_assembly_mapping"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_assembly_mapping_1_8
set sample_id, val("assembly_mapping_1_8"), val("1_8"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_assembly_mapping_1_8
file ".versions"

    script:
    """
    {
        echo [DEBUG] BUILDING BOWTIE INDEX FOR ASSEMBLY: $assembly >> .command.log 2>&1
        bowtie2-build --threads ${task.cpus} $assembly genome_index >> .command.log 2>&1
        echo [DEBUG] MAPPING READS FROM $fastq >> .command.log 2>&1
        bowtie2 -q --very-sensitive-local --threads ${task.cpus} -x genome_index -1 ${fastq[0]} -2 ${fastq[1]} -S mapping.sam >> .command.log 2>&1
        echo [DEBUG] CONVERTING AND SORTING SAM TO BAM >> .command.log 2>&1
        samtools sort -o sorted.bam -O bam -@ ${task.cpus} mapping.sam && rm *.sam  >> .command.log 2>&1
        echo [DEBUG] CREATING BAM INDEX >> .command.log 2>&1
        samtools index sorted.bam >> .command.log 2>&1
        echo [DEBUG] ESTIMATING READ DEPTH >> .command.log 2>&1
        parallel -j ${task.cpus} samtools depth -ar {} sorted.bam \\> {}.tab  ::: \$(grep ">" $assembly | cut -c 2- | tr " " "_")
        # Insert 0 coverage count in empty files. See Issue #2
        echo [DEBUG] REMOVING EMPTY FILES  >> .command.log 2>&1
        find . -size 0 -print0 | xargs -0 -I{} sh -c 'echo -e 0"\t"0"\t"0 > "{}"'
        echo [DEBUG] COMPILING COVERAGE REPORT  >> .command.log 2>&1
        parallel -j ${task.cpus} echo -n {.} '"\t"' '&&' cut -f3 {} '|' paste -sd+ '|' bc >> coverages.tsv  ::: *.tab
        cat *.tab > coverage_per_bp.tsv
        rm *.tab
        if [ -f "coverages.tsv" ]
        then
            echo pass > .status
        else
            echo fail > .status
        fi
        echo -n "" > .report.json
        echo -n "" > .versions
    } || {
        echo fail > .status
    }
    """
}


/** PROCESS_ASSEMBLY_MAPPING -  MAIN
Processes the results from the assembly_mapping process and filters the
assembly contigs based on coverage and length thresholds.
*/
process process_assembly_mapping_1_8 {

    // Send POST request to platform
    
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_8 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_8 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_8 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId assembly_mapping_1_8 \"$params.platformSpecies\" false"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }
    

    tag { sample_id }
    // This process can only use a single CPU
    cpus 1

    input:
    set sample_id, file(assembly), file(coverage), file(coverage_bp), file(bam_file), file(bam_index) from MAIN_am_out_1_8
    val opts from IN_assembly_mapping_opts_1_8
    val gsize from IN_genome_size_1_8

    output:
    set sample_id, '*_filt.fasta', 'filtered.bam', 'filtered.bam.bai' into assembly_mapping_out_1_7
    set sample_id, val("1_8_process_am"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_process_am_1_8
set sample_id, val("process_am_1_8"), val("1_8"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_process_am_1_8
file ".versions"

    script:
    template "process_assembly_mapping.py"

}


SIDE_BpCoverage_1_8.set{ SIDE_BpCoverage_1_9 }



process pilon_1_9 {

    // Send POST request to platform
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_9 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_9 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_9 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId pilon_1_9 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { sample_id }
    echo false
    publishDir 'results/assembly/pilon_1_9/', mode: 'copy', pattern: "*.fasta"

    input:
    set sample_id, file(assembly), file(bam_file), file(bam_index) from assembly_mapping_out_1_7

    output:
    set sample_id, '*_polished.fasta' into pilon_out_1_8, pilon_report_1_9
    set sample_id, val("1_9_pilon"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_pilon_1_9
set sample_id, val("pilon_1_9"), val("1_9"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_pilon_1_9
file ".versions"

    script:
    """
    {
        pilon_mem=${String.valueOf(task.memory).substring(0, String.valueOf(task.memory).length() - 1).replaceAll("\\s", "")}
        java -jar -Xms256m -Xmx\${pilon_mem} /NGStools/pilon-1.22.jar --genome $assembly --frags $bam_file --output ${assembly.name.replaceFirst(~/\.[^\.]+$/, '')}_polished --changes --threads $task.cpus >> .command.log 2>&1
        echo pass > .status
    } || {
        echo fail > .status
    }
    """

}

process pilon_report_1_9 {

    
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh"
        afterScript "report_POST.sh $params.projectId $params.pipelineId 1_9 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId pilon_1_9 \"$params.platformSpecies\" false"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh"
        }
    

    tag { sample_id }

    input:
    set sample_id, file(assembly), file(coverage_bp) from pilon_report_1_9.join(SIDE_BpCoverage_1_9)

    output:
    file "*_assembly_report.csv" into pilon_report_out_1_9
    set sample_id, val("1_9_pilon_report"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_pilon_report_1_9
set sample_id, val("pilon_report_1_9"), val("1_9"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_pilon_report_1_9
file ".versions"

    script:
    template "assembly_report.py"

}


process compile_pilon_report_1_9 {

    publishDir "reports/assembly/pilon_1_9/", mode: 'copy'

    input:
    file(report) from pilon_report_out_1_9.collect()

    output:
    file "pilon_assembly_report.csv"

    """
    echo Sample,Number of contigs,Average contig size,N50,Total assembly length,GC content,Missing data > pilon_assembly_report.csv
    cat $report >> pilon_assembly_report.csv
    """
}



// Check for the presence of absence of the minimum contig size parameter
if (params.size_1_10 == null){
    exit 1, "A minimum contig size must be provided."
}

IN_min_contig_size_1_10 = Channel.value(params.size_1_10)

process split_assembly_1_10 {

        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_10 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_10 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_10 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId split_assembly_1_10 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { sample_id }

    publishDir "results/assembly/split_assembly_1_10/${sample_id}/"

    input:
    set sample_id, file(assembly) from pilon_out_1_8
    val min_contig_size from IN_min_contig_size_1_10

    output:
    file '*split.fasta' into splitCh_1_10 optional true
    set sample_id, val("1_10_split_assembly"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_split_assembly_1_10
set sample_id, val("split_assembly_1_10"), val("1_10"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_split_assembly_1_10
file ".versions"

    script:
    template "split_fasta.py"


}
_split_assembly_out_1_9 = Channel.create()

splitCh_1_10.flatMap().map{ it -> [it.toString().tokenize('/').last().tokenize('.')[0..-2].join('.'), it] }
    .into(_split_assembly_out_1_9)


_split_assembly_out_1_9.into{ split_assembly_out_1_9;dengue_typing_in_1_10;mafft_in_1_11 }

file(params.BD_sequence_file_2_11) ? params.BD_sequence_file_2_11 : exit(1, "'BD_sequence_file_2_11' parameter missing")

process dengue_typing_2_11 {

    // Send POST request to platform
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 2_11 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 2_11 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 2_11 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId dengue_typing_2_11 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { sample_id }
    errorStrategy "ignore"
    publishDir "results/dengue_typing/${sample_id}/"

    input:
    set sample_id, file(assembly) from dengue_typing_in_1_10

    output:
    file "seq_typing*"
    set sample_id, val("2_11_dengue_typing"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_dengue_typing_2_11
set sample_id, val("dengue_typing_2_11"), val("2_11"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_dengue_typing_2_11
file ".versions"

    script:
    """
    {
        # Prevents read-only issues
        mkdir rematch_temp
        cp -r /NGStools/ReMatCh rematch_temp
        export PATH="\$(pwd)/rematch_temp/ReMatCh:\$PATH"

        seq_typing.py assembly -f ${assembly} -b ${ params.BD_sequence_file_2_11 } -o ./ -j $task.cpus -t nucl

        # Add information to dotfiles
        json_str="{'tableRow':[{'sample':'${sample_id}','data':[{'header':'seqtyping','value':'\$(cat seq_typing.report.txt)','table':'typing'}]}],'metadata':[{'sample':'${sample_id}','treeData':'\$(cat seq_typing.report.txt)'}]}"
        echo \$json_str > .report.json
        version_str="[{'program':'seq_typing.py','version':'0.1'}]"
        echo \$version_str > .versions

        rm -r rematch_temp

        if [ -s seq_typing.report.txt ];
        then
            echo pass > .status
        else
            echo fail > .status
        fi
    } || {
        echo fail > .status
        json_str="{'tableRow':[{'sample':'${sample_id}','data':[{'header':'seqtyping','value':'NA','table':'typing'}]}]}"
        echo \$json_str > .report.json
    }
    """

}

process mafft_3_12 {

        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 3_12 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 3_12 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 3_12 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId mafft_3_12 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { 'mafft' }

    publishDir "results/alignment/mafft_3_12/"

    input:
    file(assembly) from mafft_in_1_11.map{ it[1] }.collect()

    output:
    file ("*.align") into mafft_out_3_11
    set val('single'), val("3_12_mafft"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_mafft_3_12
set val('single'), val("mafft_3_12"), val("3_12"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_mafft_3_12
file ".versions"

    script:
    """

    cat ${assembly} > all_assemblies.fasta

    mafft --adjustdirection --thread $task.cpus --auto all_assemblies.fasta > ${workflow.scriptName}.align
    """

}


IN_substitution_model_3_13 = Channel.value(params.substitutionModel_3_13)
IN_seed_number_3_13 = Channel.value(params.seedNumber_3_13)
IN_bootstrap_number_3_13 = Channel.value(params.bootstrap_3_13)

process raxml_3_13 {

        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 3_13 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 3_13 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 3_13 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId raxml_3_13 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { 'raxml' }

    publishDir "results/phylogeny/raxml_3_13/"

    input:
    file(alignment) from mafft_out_3_11
    val substitution_model from IN_substitution_model_3_13
    val seednumber from IN_seed_number_3_13
    val bootstrapnumber from IN_bootstrap_number_3_13

    output:
    file ("RAxML_*") into raxml_out_3_12
    set val('single'), val("3_13_raxml"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_raxml_3_13
set val('single'), val("raxml_3_13"), val("3_13"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_raxml_3_13
file ".versions"

    script:
    """
    raxmlHPC -s ${alignment} -p 12345 -m ${substitution_model} -T $task.cpus -n $workflow.scriptName -f a -x ${seednumber} -N ${bootstrapnumber}

    # Add information to dotfiles
    json_str="{'treeData':[{'trees':['\$(cat RAxML_bipartitions.*.nf)', 'bootstrap': '${bootstrapnumber}']}]}"

    echo \$json_str > .report.json

    version_str="[{'program':'raxmlHPC','version':'8.2.11'}]"
    echo \$version_str > .versions
    """

}



/** STATUS
Reports the status of a sample in any given process.
*/
process status {

    tag { sample_id }
    publishDir "pipeline_status/$task_name"

    input:
    set sample_id, task_name, status, warning, fail, file(log) from STATUS_integrity_coverage_1_1.mix(STATUS_fastqc_1_2,STATUS_fastqc_report_1_2,STATUS_trimmomatic_1_2,STATUS_filter_poly_1_3,STATUS_remove_host_1_4,STATUS_report_remove_host_1_4,STATUS_bowtie_1_5,STATUS_report_bowtie_1_5,STATUS_retrieve_mapped_1_6,STATUS_va_spades_1_7,STATUS_va_megahit_1_7,STATUS_report_viral_assembly_1_7,STATUS_assembly_mapping_1_8,STATUS_process_am_1_8,STATUS_pilon_1_9,STATUS_pilon_report_1_9,STATUS_split_assembly_1_10,STATUS_dengue_typing_2_11,STATUS_mafft_3_12,STATUS_raxml_3_13)

    output:
    file '*.status' into master_status
    file '*.warning' into master_warning
    file '*.fail' into master_fail
    file '*.log'

    """
    echo $sample_id, $task_name, \$(cat $status) > ${sample_id}_${task_name}.status
    echo $sample_id, $task_name, \$(cat $warning) > ${sample_id}_${task_name}.warning
    echo $sample_id, $task_name, \$(cat $fail) > ${sample_id}_${task_name}.fail
    echo "\$(cat .command.log)" > ${sample_id}_${task_name}.log
    """
}

process compile_status_buffer {

    input:
    file status from master_status.buffer( size: 5000, remainder: true)
    file warning from master_warning.buffer( size: 5000, remainder: true)
    file fail from master_fail.buffer( size: 5000, remainder: true)

    output:
    file 'master_status_*.csv' into compile_status_buffer
    file 'master_warning_*.csv' into compile_warning_buffer
    file 'master_fail_*.csv' into compile_fail_buffer

    """
    cat $status >> master_status_${task.index}.csv
    cat $warning >> master_warning_${task.index}.csv
    cat $fail >> master_fail_${task.index}.csv
    """
}

process compile_status {

    publishDir 'reports/status'

    input:
    file status from compile_status_buffer.collect()
    file warning from compile_warning_buffer.collect()
    file fail from compile_fail_buffer.collect()

    output:
    file "*.csv"

    """
    cat $status >> master_status.csv
    cat $warning >> master_warning.csv
    cat $fail >> master_fail.csv
    """

}


/** Reports
Compiles the reports from every process
*/
process report {

    tag { sample_id }

    input:
    set sample_id,
            task_name,
            pid,
            report_json,
            version_json,
            trace from REPORT_integrity_coverage_1_1.mix(REPORT_fastqc_1_2,REPORT_fastqc_report_1_2,REPORT_trimmomatic_1_2,REPORT_filter_poly_1_3,REPORT_remove_host_1_4,REPORT_report_remove_host_1_4,REPORT_bowtie_1_5,REPORT_report_bowtie_1_5,REPORT_retrieve_mapped_1_6,REPORT_va_spades_1_7,REPORT_va_megahit_1_7,REPORT_report_viral_assembly_1_7,REPORT_assembly_mapping_1_8,REPORT_process_am_1_8,REPORT_pilon_1_9,REPORT_pilon_report_1_9,REPORT_split_assembly_1_10,REPORT_dengue_typing_2_11,REPORT_mafft_3_12,REPORT_raxml_3_13)

    output:
    file "*" optional true into master_report

    """
    prepare_reports.py $report_json $version_json $trace $sample_id $task_name 1 $pid $workflow.scriptId $workflow.runName
    """

}


process compile_reports {

    publishDir "pipeline_report/", mode: "copy"

    if ( params.reportHTTP != null ){
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH;"
        afterScript "metadata_POST.sh $params.projectId $params.pipelineId 0 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId 0 \"$params.platformSpecies\""
    }

    input:
    file report from master_report.collect()
    file forks from Channel.fromPath(".forkTree.json")
    file dag from Channel.fromPath(".treeDag.json")
    file js from Channel.fromPath("${workflow.projectDir}/resources/main.js.zip")

    output:
    file "pipeline_report.json"
    file "pipeline_report.html"
    file "src/main.js"

    script:
    template "compile_reports.py"
}

workflow.onComplete {
  // Display complete message
  log.info "Completed at: " + workflow.complete
  log.info "Duration    : " + workflow.duration
  log.info "Success     : " + workflow.success
  log.info "Exit status : " + workflow.exitStatus
}

workflow.onError {
  // Display error message
  log.info "Workflow execution stopped with the following message:"
  log.info "  " + workflow.errorMessage
}
