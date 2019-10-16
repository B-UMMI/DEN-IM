#!/usr/bin/env nextflow

import Helper
import CollectInitialMetadata

// Pipeline version
if (workflow.commitId){
    version = "2.2 $workflow.revision"
} else {
    version = "2.2 (local version)"
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

Help.start_info(infoMap, "$workflow.start", "$workflow.profile", version)
CollectInitialMetadata.print_metadata(workflow)
    

// Placeholder for main input channels
if (params.fastq instanceof Boolean){exit 1, "'fastq' must be a path pattern. Provide value:'$params.fastq'"}
if (!params.fastq){ exit 1, "'fastq' parameter missing"}
IN_fastq_raw = Channel.fromFilePairs(params.fastq, size: -1).ifEmpty { exit 1, "No fastq files provided with pattern:'${params.fastq}'" }

// Placeholder for secondary input channels


// Placeholder for extra input channels


// Placeholder to fork the raw input channel

IN_fastq_raw.set{ integrity_coverage_in_1_0 }


IN_genome_size_1_1 = Channel.value(params.genomeSize)
    .map{it -> it.toString().isNumber() ? it : exit(1, "The genomeSize parameter must be a number or a float. Provided value: '${params.genomeSize_}'")}

IN_min_coverage_1_1 = Channel.value(params.minCoverage)
    .map{it -> it.toString().isNumber() ? it : exit(1, "The minCoverage parameter must be a number or a float. Provided value: '${params.minCoverage_}'")}

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


// Check sliding window parameter
if ( params.trimSlidingWindow.toString().split(":").size() != 2 ){
    exit 1, "'trimSlidingWindow' parameter must contain two values separated by a ':'. Provided value: '${params.trimSlidingWindow}'"
}
if ( !params.trimLeading.toString().isNumber() ){
    exit 1, "'trimLeading' parameter must be a number. Provide value: '${params.trimLeading}'"
}
if ( !params.trimTrailing.toString().isNumber() ){
    exit 1, "'trimTrailing' parameter must be a number. Provide value: '${params.trimTrailing}'"
}
if ( !params.trimMinLength.toString().isNumber() ){
    exit 1, "'trimMinLength' parameter must be a number. Provide value: '${params.trimMinLength}'"
}

IN_trimmomatic_opts_1_2 = Channel.value([params.trimSlidingWindow,params.trimLeading,params.trimTrailing,params.trimMinLength])
IN_adapters_1_2 = Channel.value(params.adapters)

clear = params.clearInput ? "true" : "false"
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
    set sample_id, file(fastq_pair), file('pair_*') into MAIN_fastqc_out_1_2
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
    set sample_id, file(fastq_pair), file(results) from MAIN_fastqc_out_1_2
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
    set sample_id, "*trim.fastq.gz" into fastqc_trimmomatic_out_1_1
    file 'trimmomatic_report.csv'
    set sample_id, val("1_2_trimmomatic"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_trimmomatic_1_2
set sample_id, val("trimmomatic_1_2"), val("1_2"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_trimmomatic_1_2
file ".versions"

    script:
    template "trimmomatic.py"

}



IN_adapter_1_3 = Channel.value(params.pattern)

clear = params.clearInput ? "true" : "false"
checkpointClear_1_3 = Channel.value(clear)

process filter_poly_1_3 {

    // Send POST request to platform
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_3 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_3 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_3 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId filter_poly_1_3 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { sample_id }
    echo true

    errorStrategy { task.exitStatus == 120 ? 'ignore' : 'retry' }

    input:
    set sample_id, file(fastq_pair) from fastqc_trimmomatic_out_1_1
    val adapter from IN_adapter_1_3
    val clear from checkpointClear_1_3

    output:
    set sample_id , file("${sample_id}_filtered*") into filter_poly_out_1_2
    set sample_id, val("1_3_filter_poly"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_filter_poly_1_3
set sample_id, val("filter_poly_1_3"), val("1_3"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_filter_poly_1_3
file ".versions"

    script:
    """
    a=(${fastq_pair})

    if ((\${#a[@]} > 1));
    then
        gunzip -c ${fastq_pair[0]} >  ${sample_id}_1.fq
        gunzip -c ${fastq_pair[1]} >  ${sample_id}_2.fq
    else
        gunzip -c ${fastq_pair[0]} >  ${sample_id}.fq
    fi

    for seqfile in *.fq;
    do if [ ! -s \$seqfile  ]
    then
        echo \$seqfile is empty
        echo 'No data left after polymorphic sequence filtering' > .fail
        exit 120
    fi
    done

    if ((\${#a[@]} > 1));
    then
        prinseq-lite.pl --fastq ${sample_id}_1.fq  --fastq2 ${sample_id}_2.fq  --custom_params "${adapter}" -out_format 3 -out_good ${sample_id}_filtered
    else
        prinseq-lite.pl --fastq ${sample_id}.fq --custom_params "${adapter}" -out_format 3 -out_good ${sample_id}_filtered
    fi

    gzip ${sample_id}_filtered*

    if [ "$clear" = "true" ];
    then
        work_regex=".*/work/.{2}/.{30}/.*"
        file_source1=\$(readlink -f \$(pwd)/${fastq_pair[0]})
        file_source2=\$(readlink -f \$(pwd)/${fastq_pair[1]})
        if [[ "\$file_source1" =~ \$work_regex ]]; then
            rm \$file_source1 \$file_source2
        fi
    fi
    """
}



// Check for the presence of absence of both index and fasta reference
if (params.index == null && params.reference == null){
    exit 1, "An index or a reference fasta file must be provided."
} else if (params.index != null && params.reference != null){
    exit 1, "Provide only an index OR a reference fasta file."
}

clear = params.clearInput ? "true" : "false"
checkpointClear_1_4 = Channel.value(clear)

if (params.reference){

    reference_in_1_4 = Channel.fromPath(params.reference)
        .map{it -> file(it).exists() ? [it.toString().tokenize('/').last().tokenize('.')[0..-2].join('.') ,it] : null}

    process bowtie_build_1_4 {

        // Send POST request to platform
            if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_4 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_4 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_4 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId bowtie_1_4 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

        tag { build_id }
        storeDir 'bowtie_index/'
        maxForks 1

        input:
        set build_id, file(fasta) from reference_in_1_4

        output:
        val build_id into bowtieIndexId_1_4
        file "${build_id}*.bt2" into bowtieIndex_1_4

        script:
        """
        # checking if reference file is empty. Moved here due to allow reference file to be inside the container.
        if [ ! -f "$fasta" ]
        then
            echo "Error: ${fasta} file not found."
            exit 1
        fi

        bowtie2-build ${fasta} $build_id > ${build_id}_bowtie2_build.log
        """
    }
} else {
    bowtieIndexId_1_4 = Channel.value(params.index.split("/").last())
    bowtieIndex_1_4 = Channel.fromPath("${params.index}*.bt2").collect().toList()
}


process bowtie_1_4 {

    // Send POST request to platform
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_4 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_4 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_4 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId bowtie_1_4 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { sample_id }
    publishDir 'results/mapping/bowtie_1_4/'

    input:
    set sample_id, file(fastq_pair) from filter_poly_out_1_2
    each index from bowtieIndexId_1_4
    each file(index_files) from bowtieIndex_1_4

    output:
    set sample_id , file("pair_info.txt"), file("*.bam") into bowtie_out_1_3
    set sample_id, file("*_bowtie2.log") into into_json_1_4
    set sample_id, val("1_4_bowtie"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_bowtie_1_4
set sample_id, val("bowtie_1_4"), val("1_4"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_bowtie_1_4
file ".versions"

    script:
    """
    {
        a=(${fastq_pair})

        if ((\${#a[@]} > 1));
        then
            echo "True" > pair_info.txt
            bowtie2 -x $index -1 ${fastq_pair[0]} -2 ${fastq_pair[1]} -p $task.cpus 1> ${sample_id}.bam 2> ${sample_id}_bowtie2.log

            if [ "$clear" = "true" ];
            then
                work_regex=".*/work/.{2}/.{30}/.*"
                file_source1=\$(readlink -f \$(pwd)/${fastq_pair[0]})
                file_source2=\$(readlink -f \$(pwd)/${fastq_pair[1]})
                if [[ "\$file_source1" =~ \$work_regex ]]; then
                    rm \$file_source1 \$file_source2
                fi
            fi

            echo pass > .status
        else
            echo "False" > pair_info.txt
            bowtie2 -x $index -U ${fastq_pair[0]} -p $task.cpus 1> ${sample_id}.bam 2> ${sample_id}_bowtie2.log

            if [ "$clear" = "true" ];
            then
                work_regex=".*/work/.{2}/.{30}/.*"
                file_source1=\$(readlink -f \$(pwd)/${fastq_pair[0]})
                if [[ "\$file_source1" =~ \$work_regex ]]; then
                    rm \$file_source1
                fi
            fi

            echo pass > .status
        fi
    } || {
        echo fail > .status
    }
    """
}


process report_bowtie_1_4 {

        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_4 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_4 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_4 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId bowtie_1_4 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { sample_id }

    input:
    set sample_id, file(bowtie_log) from into_json_1_4

    output:
    set sample_id, val("1_4_report_bowtie"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_report_bowtie_1_4
set sample_id, val("report_bowtie_1_4"), val("1_4"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_report_bowtie_1_4
file ".versions"

    script:
    template "process_mapping.py"

}


process retrieve_mapped_1_5 {

    // Send POST request to platform
    if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_5 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_5 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_5 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId retrieve_mapped_1_5 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { sample_id }
    publishDir 'results/mapping/retrieve_mapped_1_5/'

    input:
    set sample_id, file(is_pair), file(bam) from bowtie_out_1_3

    output:
    set sample_id , file("*_mapped*") into OUT_retrieve_mapped_1_4
    set sample_id, val("1_5_retrieve_mapped"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_retrieve_mapped_1_5
set sample_id, val("retrieve_mapped_1_5"), val("1_5"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_retrieve_mapped_1_5
    file ".versions"

    script:
    """
    if [[ \$(cat ${is_pair}) == "True" ]];
    then
        samtools view -buh -F 12 -o ${sample_id}_samtools.bam -@ $task.cpus ${bam}
        rm ${bam}

        samtools fastq -1 ${sample_id}_mapped_1.fq -2 ${sample_id}_mapped_2.fq ${sample_id}_samtools.bam
        rm ${sample_id}_samtools.bam

    else
        samtools view -buh -F 4 -o ${sample_id}_samtools.bam -@ $task.cpus ${bam}
        rm ${bam}

        samtools fastq ${sample_id}_samtools.bam > ${sample_id}_mapped.fq
        rm ${sample_id}_samtools.bam

    fi
    """
}

process renamePE_1_4 {

    tag { sample_id }
    publishDir 'results/mapping/retrieve_mapped_{{ pid }}/'

    if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_5 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_5 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_5 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId renamePE_1_5 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
    }

    input:
    set sample_id, file(fastq_pair) from OUT_retrieve_mapped_1_4

    output:
    set sample_id , file("*.headersRenamed*") into retrieve_mapped_out_1_4
    set sample_id, val("1_5_renamePE"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_renamePE_1_5
    set sample_id, val("renamePE_1_5"), val("1_5"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_renamePE_1_5
    file ".versions"

    script:
    template "renamePE_samtoolsFASTQ.py"

}

IN_genome_size_1_6 = Channel.value(params.genomeSize)
    .map{it -> it.toString().isNumber() ? it : exit (1, "The genomeSize parameter must be a number or a float. Provided value: '${params.genomeSize}'")}
IN_min_coverage_1_6 = Channel.value(params.minCoverage)
    .map{it -> it.toString().isNumber() ? it : exit (1, "The minCoverage parameter must be a number or a float. Provided value: '${params.minCoverage}'")}

process integrity_coverage2_1_6 {

    // Send POST request to platform
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_6 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_6 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_6 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId check_coverage_1_6 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { sample_id }
    cpus 1

    input:
    set sample_id, file(fastq_pair) from retrieve_mapped_out_1_4
    val gsize from IN_genome_size_1_6
    val cov from IN_min_coverage_1_6
    // Use -e option for skipping encoding guess
    val opts from Channel.value('-e')

    output:
    set sample_id,
        file(fastq_pair),
        file('*_coverage'),
        file('*_max_len') optional true into MAIN_integrity_1_6
    file('*_report') into LOG_report_coverage_1_6
    set sample_id, val("1_6_check_coverage"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_check_coverage_1_6
set sample_id, val("check_coverage_1_6"), val("1_6"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_check_coverage_1_6
file ".versions"

    script:
    template "integrity_coverage.py"
}

_check_coverage_out_1_5 = Channel.create()
SIDE_max_len_1_6 = Channel.create()

MAIN_integrity_1_6
    .filter{ it[2].text != "fail" }
    .separate(_check_coverage_out_1_5, SIDE_max_len_1_6){
        a -> [ [a[0], a[1]], [a[0], a[3].text]]
    }


process report_coverage2_1_6 {

    // This process can only use a single CPU
    cpus 1
    publishDir 'reports/coverage_1_6/'

    input:
    file(report) from LOG_report_coverage_1_6.filter{ it.text != "corrupt" }.collect()

    output:
    file 'estimated_coverage_second.csv'

    """
    echo Sample,Estimated coverage,Status >> estimated_coverage_second.csv
    cat $report >> estimated_coverage_second.csv
    """
}


_check_coverage_out_1_5.into{ check_coverage_out_1_5;_LAST_fastq_1_8;_LAST_fastq_1_11 }


SIDE_max_len_1_6.set{ SIDE_max_len_1_7 }


//MAIN INPUT - FASTQ FILES
spades_in = Channel.create()
megahit_in = Channel.create()

check_coverage_out_1_5.into{ spades_in; megahit_in }

//EXPECTED GENOME SIZE
if ( !params.minimumContigSize.toString().isNumber() ){
    exit 1, "'minimumContigSize' parameter must be a number. Provided value: '${params.minimumContigSize}'"
}

//SPADES OPTIONS
if ( !params.spadesMinCoverage.toString().isNumber() ){
    exit 1, "'spadesMinCoverage' parameter must be a number. Provided value: '${params.spadesMinCoverage}'"
}
if ( !params.spadesMinKmerCoverage.toString().isNumber()){
    exit 1, "'spadesMinKmerCoverage' parameter must be a number. Provided value: '${params.spadesMinKmerCoverage}'"
}

if ( params.spadesKmers.toString().split(" ").size() <= 1 ){
    if (params.spadesKmers.toString() != 'auto'){
        exit 1, "'spadesKmers' parameter must be a sequence of space separated numbers or 'auto'. Provided value: ${params.spadesKmers}"
    }
}

clear = params.clearInput ? "true" : "false"
checkpointClearSpades_1_7 = Channel.value(clear)
checkpointClearMegahit_1_7 = Channel.value(clear)

//MEGAHIT OPTIONS
if ( params.megahitKmers.toString().split(" ").size() <= 1 ){
    if (params.megahitKmers.toString() != 'auto'){
        exit 1, "'megahitKmers' parameter must be a sequence of space separated numbers or 'auto'. Provided value: ${params.megahitKmers}"
    }
}

//SPADES INPUT CHANNELS
IN_spades_opts_1_7 = Channel.value([params.spadesMinCoverage,params.spadesMinKmerCoverage])
IN_spades_kmers_1_7 = Channel.value(params.spadesKmers)

//MEGAGIT INPUT CHANNELS
IN_megahit_kmers_1_7 = Channel.value(params.megahitKmers)

SIDE_max_len_spades = Channel.create()
SIDE_max_len_megahit = Channel.create()
SIDE_max_len_1_7.into{SIDE_max_len_spades ; SIDE_max_len_megahit}

disableRR_1_7 = "false"

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
    val clear from checkpointClearSpades_1_7
    val disable_rr from disableRR_1_7

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
assembly_spades.choice(good_assembly, megahit){a -> a[1].toString() == "null" ? false : VerifyCompletness.contigs(a[1].toString(), params.minimumContigSize.toInteger()) == true ? 0 : 1}


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
    val clear from checkpointClearSpades_1_7

    output:
    set sample_id, file('*megahit*.fasta') into megahit_assembly
    set sample_id, val("1_7_va_megahit"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_va_megahit_1_7
set sample_id, val("va_megahit_1_7"), val("1_7"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_va_megahit_1_7
file ".versions"

    script:
    template "megahit.py"

}

good_assembly.mix(megahit_assembly).into{ to_report_1_7 ; viral_assembly_out_1_6 }
orf_size = Channel.value(params.minimumContigSize)


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



if ( !params.minAssemblyCoverage.toString().isNumber() ){
    if (params.minAssemblyCoverage.toString() != 'auto'){
        exit 1, "'minAssemblyCoverage' parameter must be a number or 'auto'. Provided value: ${params.minAssemblyCoverage}"
    }
}
if ( !params.AMaxContigs.toString().isNumber() ){
    exit 1, "'AMaxContigs' parameter must be a number. Provide value: '${params.AMaxContigs}'"
}

IN_assembly_mapping_opts_1_8 = Channel.value([params.minAssemblyCoverage,params.AMaxContigs])
IN_genome_size_1_8 = Channel.value(params.genomeSize)


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
        a=(${fastq})

        echo [DEBUG] BUILDING BOWTIE INDEX FOR ASSEMBLY: $assembly >> .command.log 2>&1
        bowtie2-build --threads ${task.cpus} $assembly genome_index >> .command.log 2>&1

        if ((\${#a[@]} > 1));
        then
            echo [DEBUG] MAPPING READS FROM $fastq >> .command.log 2>&1
            bowtie2 -q --very-sensitive-local --threads ${task.cpus} -x genome_index -1 ${fastq[0]} -2 ${fastq[1]} -S mapping.sam >> .command.log 2>&1

        else
            echo [DEBUG] MAPPING READS FROM $fastq >> .command.log 2>&1
            bowtie2 -q --very-sensitive-local --threads ${task.cpus} -x genome_index -U ${fastq[0]} -S mapping.sam >> .command.log 2>&1
        fi

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



clear = params.clearInput ? "true" : "false"
checkpointClear_1_9 = Channel.value(clear)

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
    val clear from checkpointClear_1_9

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

        if [ "$clear" = "true" ];
        then
            work_regex=".*/work/.{2}/.{30}/.*"
            assembly_file=\$(readlink -f \$(pwd)/${assembly})
            bam_file=\$(readlink -f \$(pwd)/${bam_file})
            if [[ "\$assembly_file" =~ \$work_regex ]]; then
                rm \$assembly_file \$bam_file
            fi
        fi

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
if (params.splitSize == null){
    exit 1, "A minimum contig size must be provided."
}

IN_min_contig_size_1_10 = Channel.value(params.splitSize)

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
    file('*.fasta') into splitCh_1_10
    set sample_id, val("1_10_split_assembly"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_split_assembly_1_10
set sample_id, val("split_assembly_1_10"), val("1_10"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_split_assembly_1_10
file ".versions"

    script:
    template "split_fasta.py"


}

split_assembly_out_1_9 = Channel.create()

splitCh_1_10.flatMap().map{ it -> [it.toString().tokenize('/').last().tokenize('.')[0..-2].join('.'), it]}.into( split_assembly_out_1_9 )


// Check for the presence of absence of fasta reference
if (params.typingReference == null) {
    exit 1, "Dengue_typing: A reference fasta file must be provided."
}

getRef_1_11 = params.getGenome ? "true" : "false"
checkpointReferenceGenome_1_11 = Channel.value(getRef_1_11)
checkpointReferenceGenome_1_11.into{ reference_reads_1_11 ; reference_assembly_1_11 }

reference_1_11 = Channel.fromPath(params.typingReference)

class VerifyCompletnessTyping {

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
                lineThreshold = 0
            } else {
                lineThreshold += line.length()
                if(lineThreshold >= threshold) {
                    return true;
                }
             }
        }

        return false;
    }
}


type_reads_1_11 = Channel.create()
type_assembly_1_11 = Channel.create()
split_assembly_out_1_9.choice(type_assembly_1_11, type_reads_1_11){a -> a[1].toString() == "null" ? false : VerifyCompletnessTyping.contigs(a[1].toString(), 10000) == true ? 0 : 1}

process dengue_typing_assembly_1_11 {

    // Send POST request to platform
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_11 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_11 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_11 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId dengue_typing_1_11 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { sample_id }

    publishDir "results/dengue_typing/${sample_id}/"


    input:
    set sample_id, file(assembly) from type_assembly_1_11
    val get_reference from reference_assembly_1_11
    each file(reference) from Channel.fromPath("${params.typingReference}")

    output:
    file "seq_typing*"
    set sample_id, file(assembly) into out_typing_assembly_1_11
    file("*.fa") optional true into _ref_seqTyping_assembly_1_11
    set sample_id, val("1_11_dengue_typing_assembly"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_dengue_typing_assembly_1_11
set sample_id, val("dengue_typing_assembly_1_11"), val("1_11"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_dengue_typing_assembly_1_11
file ".versions"

    script:
    template "dengue_typing_assembly.py"

}


process dengue_typing_reads_1_11 {

// Send POST request to platform
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_11 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_11 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_11 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId dengue_typing_1_11 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { sample_id }

    publishDir "results/dengue_typing/${sample_id}/"

    errorStrategy { task.exitStatus == 120 ? 'ignore' : 'retry' }

    input:
    set sample_id, file(assembly), file(fastq_pair) from type_reads_1_11.join(_LAST_fastq_1_11)
    val get_reference from reference_reads_1_11
    each file(reference) from Channel.fromPath("${params.typingReference}")

    output:
    file "seq_typing*"
    set sample_id, file("*consensus.fasta") into out_typing_reads_1_11
    file("*.fa") optional true into _ref_seqTyping_reads_1_11
    set sample_id, val("1_11_dengue_typing_reads"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_dengue_typing_reads_1_11
set sample_id, val("dengue_typing_reads_1_11"), val("1_11"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_dengue_typing_reads_1_11
file ".versions"

    script:
    template "dengue_typing_reads.py"

}

out_typing_assembly_1_11.mix(out_typing_reads_1_11).set{ dengue_typing_out_1_10 }

_ref_seqTyping_assembly_1_11.mix(_ref_seqTyping_reads_1_11).set{ _ref_seqTyping_1_11 }


_ref_seqTyping_1_11.set{ _ref_seqTyping_1_12 }


// True when a dengue_typing secondary channel is connected
has_ref_1_12 = binding.hasVariable('_ref_seqTyping_1_12')

if ( has_ref_1_12 ){
    dengue_typing_out_1_10.map{ it[1] }.collect().mix(_ref_seqTyping_1_12.unique{it.name}).set{mafft_input}
} else {
    dengue_typing_out_1_10.map{ it[1] }.collect().set{mafft_input}
}

//dengue_typing_out_1_10.map{ it[1] }.mix(_ref_seqTyping_1_12.unique()).set{mafft_input}

include_ncbi = params.includeNCBI ? "true" : "false"
IN_ncbi_1_12 = Channel.value(include_ncbi)

process mafft_1_12 {

        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_12 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_12 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_12 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId mafft_1_12 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { 'mafft' }

    publishDir "results/alignment/mafft_1_12/"

    input:
    file(assembly) from mafft_input.collect()
    val ncbi from IN_ncbi_1_12

    output:
    file ("*.align") into mafft_out_1_11
    set val('single'), val("1_12_mafft"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_mafft_1_12
set val('single'), val("mafft_1_12"), val("1_12"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_mafft_1_12
file ".versions"

    script:
    """
    cat ${assembly} > all_assemblies.fasta

    samples=`cat all_assemblies.fasta | grep ">" | wc -l`;

    if (( \$samples < 4 )) || [ $ncbi = "true" ]
    then
        cat ${workflow.projectDir}/ref/NCBI.fasta >> all_assemblies.fasta
        echo '{"metadata":[{"sample":"NCBI-DENV-1","treeData":"1-IV","column":"typing"},{"sample":"NCBI-DENV-2","treeData":"2-V(AsianI)","column":"typing"},{"sample":"NCBI-DENV-3","treeData":"3-III","column":"typing"},{"sample":"NCBI-DENV-4","treeData":"4-II","column":"typing"}]}' > .report.json
    fi


    mafft --adjustdirection --thread $task.cpus --auto all_assemblies.fasta > ${workflow.scriptName}.align
    """

}


IN_substitution_model_1_13 = Channel.value(params.substitutionModel)
IN_seed_number_1_13 = Channel.value(params.seedNumber)
IN_bootstrap_number_1_13 = Channel.value(params.bootstrap)
IN_simple_label_1_13 = Channel.value(params.simpleLabel)

process raxml_1_13 {

        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_13 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_13 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_13 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId raxml_1_13 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { 'raxml' }

    publishDir "results/phylogeny/raxml_1_13/"

    errorStrategy { task.exitStatus == 120 ? 'ignore' : 'retry' }

    input:
    file(alignment) from mafft_out_1_11
    val substitution_model from IN_substitution_model_1_13
    val seednumber from IN_seed_number_1_13
    val bootstrapnumber from IN_bootstrap_number_1_13

    output:
    file ("RAxML_*") into raxml_out_1_12
    file ("RAxML_bipartitions.*.nf") into into_json_1_13
    set val('single'), val("1_13_raxml"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_raxml_1_13
set val('single'), val("raxml_1_13"), val("1_13"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_raxml_1_13
file ".versions"

    script:
    """
    samples=`cat ${alignment} | grep ">" | wc -l`;
    if (( \$samples < 4 ))
    then
        echo ERROR: Too few species! RAxML is very unhappy!
        exit 120
    fi

    raxmlHPC -s ${alignment} -p 12345 -m ${substitution_model} -T $task.cpus -n $workflow.scriptName -f a -x ${seednumber} -N ${bootstrapnumber}

    # Add information to dotfiles
    version_str="[{'program':'raxmlHPC','version':'8.2.11'}]"
    echo \$version_str > .versions
    """

}

process report_raxml_1_13 {

        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_13 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_13 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_13 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId raxml_1_13 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { 'raxml' }

    input:
    file(newick) from into_json_1_13
    val label from IN_simple_label_1_13

    output:
    set val('single'), val("1_13_report_raxml"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_report_raxml_1_13
set val('single'), val("report_raxml_1_13"), val("1_13"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_report_raxml_1_13
file ".versions"

    script:
    template "process_newick.py"

}




/** STATUS
Reports the status of a sample in any given process.
*/
process status {

    tag { sample_id }
    publishDir "pipeline_status/$task_name"

    input:
    set sample_id, task_name, status, warning, fail, file(log) from STATUS_integrity_coverage_1_1.mix(STATUS_fastqc_1_2,STATUS_fastqc_report_1_2,STATUS_trimmomatic_1_2,STATUS_filter_poly_1_3,STATUS_bowtie_1_4,STATUS_report_bowtie_1_4,STATUS_retrieve_mapped_1_5,STATUS_renamePE_1_5,STATUS_check_coverage_1_6,STATUS_va_spades_1_7,STATUS_va_megahit_1_7,STATUS_report_viral_assembly_1_7,STATUS_assembly_mapping_1_8,STATUS_process_am_1_8,STATUS_pilon_1_9,STATUS_pilon_report_1_9,STATUS_split_assembly_1_10,STATUS_dengue_typing_assembly_1_11,STATUS_dengue_typing_reads_1_11,STATUS_mafft_1_12,STATUS_raxml_1_13,STATUS_report_raxml_1_13)

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
            trace from REPORT_integrity_coverage_1_1.mix(REPORT_fastqc_1_2,REPORT_fastqc_report_1_2,REPORT_trimmomatic_1_2,REPORT_filter_poly_1_3,REPORT_bowtie_1_4,REPORT_report_bowtie_1_4,REPORT_retrieve_mapped_1_5,REPORT_renamePE_1_5,REPORT_check_coverage_1_6,REPORT_va_spades_1_7,REPORT_va_megahit_1_7,REPORT_report_viral_assembly_1_7,REPORT_assembly_mapping_1_8,REPORT_process_am_1_8,REPORT_pilon_1_9,REPORT_pilon_report_1_9,REPORT_split_assembly_1_10,REPORT_dengue_typing_assembly_1_11,REPORT_dengue_typing_reads_1_11,REPORT_mafft_1_12,REPORT_raxml_1_13,REPORT_report_raxml_1_13)

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
    file forks from Channel.fromPath("${workflow.projectDir}/.forkTree.json")
    file dag from Channel.fromPath("${workflow.projectDir}/.treeDag.json")
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