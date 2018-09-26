class Help {

    static def start_info(Map info, String time, String profile) {

        println ""
        println "============================================================"
        println "                F L O W C R A F T"
        println "============================================================"
        println "Built using flowcraft v1.3.0"
        println ""
        if (info.containsKey("fastq")){
        int nsamples = info.fastq / 2
        println " Input FastQ                 : $info.fastq"
        println " Input samples               : $nsamples"
        }
        if (info.containsKey("fasta")){
        println " Input Fasta                 : $info.fasta"
        }
        if (info.containsKey("accessions")){
        println " Input accessions            : $info.accessions"
        }
        println " Reports are found in        : ./reports"
        println " Results are found in        : ./results"
        println " Profile                     : $profile"
        println ""
        println "Starting pipeline at $time"
        println ""

    }

    static void complete_info(nextflow.script.WorkflowMetadata wf) {

        println ""
        println "Pipeline execution summary"
        println "=========================="
        println "Completed at                 : $wf.complete"
        println "Duration                     : $wf.duration"
        println "Success                      : $wf.success"
        println "Work directory               : $wf.workDir"
        println "Exit status                  : $wf.exitStatus"
        println ""

    }

    static def print_help(Map params) {

        println ""
        println "============================================================"
        println "                F L O W C R A F T"
        println "============================================================"
        println "Built using flowcraft v1.3.0"
        println ""
        println ""
        println "Usage: "
        println "    nextflow run DEN-IM.nf"
        println ""
        println "       --fastq                     Path expression to paired-end fastq files. (default: $params.fastq) (default: 'fastq/*_{1,2}.*')"
        println "       "
        println "       Component 'INTEGRITY_COVERAGE_1_1'"
        println "       ----------------------------------"
        println "       --genomeSize_1_1            Genome size estimate for the samples in Mb. It is used to estimate the coverage and other assembly parameters andchecks (default: 0.012)"
        println "       --minCoverage_1_1           Minimum coverage for a sample to proceed. By default it's setto 0 to allow any coverage (default: 15)"
        println "       "
        println "       Component 'FASTQC_TRIMMOMATIC_1_2'"
        println "       ----------------------------------"
        println "       --adapters_1_2              Path to adapters files, if any. (default: 'None')"
        println "       --trimSlidingWindow_1_2     Perform sliding window trimming, cutting once the average quality within the window falls below a threshold. (default: '5:20')"
        println "       --trimLeading_1_2           Cut bases off the start of a read, if below a threshold quality. (default: 3)"
        println "       --trimTrailing_1_2          Cut bases of the end of a read, if below a threshold quality. (default: 3)"
        println "       --trimMinLength_1_2         Drop the read if it is below a specified length. (default: 55)"
        println "       "
        println "       Component 'FILTER_POLY_1_3'"
        println "       ---------------------------"
        println "       --adapter_1_3               Pattern to filter the reads. Please separate parametervalues with a space and separate new parameter sets with semicolon (;).Parameters are defined by two values: the pattern (any combination of theletters ATCGN), and the number of repeats or percentage of occurence. (default: 'A 50%; T 50%; N 50%')"
        println "       "
        println "       Component 'BOWTIE_1_4'"
        println "       ----------------------"
        println "       --reference_1_4             Specifies the reference genome to be provided to bowtie2-build. (default: '/ref/1_GenotypesDENV_14-05-18.fasta')"
        println "       --index_1_4                 Specifies the reference indexes to be provided to bowtie2. (default: null)"
        println "       "
        println "       Component 'REMOVE_HOST_1_6'"
        println "       ---------------------------"
        println "       --refIndex_1_6              Specifies the reference indexes to be provided to bowtie2. (default: '/index_hg19/hg19')"
        println "       "
        println "       Component 'VIRAL_ASSEMBLY_1_7'"
        println "       ------------------------------"
        println "       --minimumContigSize_1_7     Expected genome size in bases (default: 10000)"
        println "       --spadesMinCoverage_1_7     The minimum number of reads to consider an edge in the de Bruijn graph during the assembly (default: 2)"
        println "       --spadesMinKmerCoverage_1_7 Minimum contigs K-mer coverage. After assembly only keep contigs with reported k-mer coverage equal or above this value (default: 2)"
        println "       --spadesKmers_1_7           If 'auto' the SPAdes k-mer lengths will be determined from the maximum read length of each assembly. If 'default', SPAdes will use the default k-mer lengths.  (default: 'auto')"
        println "       --megahitKmers_1_7          If 'auto' the megahit k-mer lengths will be determined from the maximum read length of each assembly. If 'default', megahit will use the default k-mer lengths. (default: $params.megahitKmers) (default: 'auto')"
        println "       "
        println "       Component 'ASSEMBLY_MAPPING_1_8'"
        println "       --------------------------------"
        println "       --minAssemblyCoverage_1_8   In auto, the default minimum coverage for each assembled contig is 1/3 of the assembly mean coverage or 10x, if the mean coverage is below 10x (default: 'auto')"
        println "       --AMaxContigs_1_8           A warning is issued if the number of contigs is overthis threshold. (default: 1000)"
        println "       --genomeSize_1_8            Genome size estimate for the samples. It is used to check the ratio of contig number per genome MB (default: 0.01)"
        println "       "
        println "       Component 'SPLIT_ASSEMBLY_1_10'"
        println "       -------------------------------"
        println "       --size_1_10                 Minimum contig size (default: 10000)"
        println "       "
        println "       Component 'RAXML_3_13'"
        println "       ----------------------"
        println "       --substitutionModel_3_13    Substitution model. Option: GTRCAT, GTRCATI, ASC_GTRCAT, GTRGAMMA, ASC_GTRGAMMA etc  (default: 'GTRGAMMA')"
        println "       --seedNumber_3_13           Specify an integer number (random seed) and turn on rapid bootstrapping (default: 12345)"
        println "       --bootstrap_3_13            Specify the number of alternative runs on distinct starting trees (default: 500)"
        
    }

}

class CollectInitialMetadata {

    public static void print_metadata(nextflow.script.WorkflowMetadata workflow){

        def treeDag = new File(".treeDag.json").text
        def forkTree = new File(".forkTree.json").text

        def metadataJson = "{'nfMetadata':{'scriptId':'${workflow.scriptId}',\
'scriptName':'${workflow.scriptName}',\
'profile':'${workflow.profile}',\
'container':'${workflow.container}',\
'containerEngine':'${workflow.containerEngine}',\
'commandLine':'${workflow.commandLine}',\
'runName':'${workflow.runName}',\
'sessionId':'${workflow.sessionId}',\
'projectDir':'${workflow.projectDir}',\
'launchDir':'${workflow.launchDir}',\
'startTime':'${workflow.start}',\
'dag':${treeDag},\
'forks':${forkTree}}}"

        def json = metadataJson.replaceAll("'", '"')

        def jsonFile = new File(".metadata.json")
        jsonFile.write json
    }
}