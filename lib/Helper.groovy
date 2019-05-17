class Help {

    static def start_info(Map info, String time, String profile) {

        println ""
        println "============================================================"
        println "                     D E N - I M"
        println "============================================================"
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
        println "Built using flowcraft v1.4.0"
        println ""
        println ""
        println "Usage: "
        println "    nextflow run DEN-IM.nf"
        println ""
        println "       --fastq                     Path expression to paired-end fastq files. (default: $params.fastq) "
        println "       --genomeSize                Genome size estimate for the samples in Mb. It is used to estimate the coverage and other assembly parameters andchecks (integrity_coverage;check_coverage;assembly_mapping)"
        println "       --minCoverage               Minimum coverage for a sample to proceed. By default it's setto 0 to allow any coverage (integrity_coverage;check_coverage)"
        println "       --adapters                  Path to adapters files, if any. (fastqc_trimmomatic)"
        println "       --trimSlidingWindow         Perform sliding window trimming, cutting once the average quality within the window falls below a threshold. (fastqc_trimmomatic)"
        println "       --trimLeading               Cut bases off the start of a read, if below a threshold quality. (fastqc_trimmomatic)"
        println "       --trimTrailing              Cut bases of the end of a read, if below a threshold quality. (fastqc_trimmomatic)"
        println "       --trimMinLength             Drop the read if it is below a specified length. (fastqc_trimmomatic)"
        println "       --clearInput                Permanently removes temporary input files. This option is only useful to remove temporary files in large workflows and prevents nextflow's resume functionality. Use with caution. (fastqc_trimmomatic;filter_poly;bowtie;retrieve_mapped;viral_assembly;pilon)"
        println "       --pattern                   Pattern to filter the reads. Please separate parametervalues with a space and separate new parameter sets with semicolon (;). Parameters are defined by two values: the pattern (any combination of the letters ATCGN), and the number of repeats or percentage of occurence. (filter_poly)"
        println "       --reference                 Specifies the reference genome to be provided to bowtie2-build. (bowtie)"
        println "       --index                     Specifies the reference indexes to be provided to bowtie2. (bowtie)"
        println "       --minimumContigSize         Expected genome size in bases (viral_assembly)"
        println "       --spadesMinCoverage         The minimum number of reads to consider an edge in the de Bruijn graph during the assembly (viral_assembly)"
        println "       --spadesMinKmerCoverage     Minimum contigs K-mer coverage. After assembly only keep contigs with reported k-mer coverage equal or above this value (viral_assembly)"
        println "       --spadesKmers               If 'auto' the SPAdes k-mer lengths will be determined from the maximum read length of each assembly. If 'default', SPAdes will use the default k-mer lengths.  (viral_assembly)"
        println "       --megahitKmers              If 'auto' the megahit k-mer lengths will be determined from the maximum read length of each assembly. If 'default', megahit will use the default k-mer lengths. (default: $params.megahitKmers) (viral_assembly)"
        println "       --minAssemblyCoverage       In auto, the default minimum coverage for each assembled contig is 1/3 of the assembly mean coverage or 10x, if the mean coverage is below 10x (assembly_mapping)"
        println "       --AMaxContigs               A warning is issued if the number of contigs is overthis threshold. (assembly_mapping)"
        println "       --splitSize                 Minimum contig size (split_assembly)"
        println "       --typingReference           Typing database. (dengue_typing)"
        println "       --getGenome                 Retrieves the sequence of the closest reference. (dengue_typing)"
        println "       --substitutionModel         Substitution model. Option: GTRCAT, GTRCATI, ASC_GTRCAT, GTRGAMMA, ASC_GTRGAMMA etc  (raxml)"
        println "       --seedNumber                Specify an integer number (random seed) and turn on rapid bootstrapping (raxml)"
        println "       --bootstrap                 Specify the number of alternative runs on distinct starting trees (raxml)"
        println "       --simpleLabel               Simplify the labels in the newick tree (for interactive report only) (raxml)"
        
    }

}

class CollectInitialMetadata {

    public static void print_metadata(nextflow.script.WorkflowMetadata workflow){

        def treeDag = new File("${workflow.projectDir}/.treeDag.json").text
        def forkTree = new File("${workflow.projectDir}/.forkTree.json").text

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