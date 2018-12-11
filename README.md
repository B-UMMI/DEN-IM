# DEN-IM

A user-friendly, containerized and reproducible workflow for the analysis of dengue virus sequencing data, 
both from metagenomic or targeted sequencing methodologies.

As input DEN-IM accepts **raw paired-end Illumina sequencing reads**, both from metagenomic or targeted approaches,
and provides the user with an interactive HTML report. It's  implemented in Nextflow alongside Docker containers to 
facilitate installation and can be executed in local machines or high performance computing facilities.



## Installation

Before installing DEN-IM, a few dependencies must be installed in your system:

* **Nextflow**

Nextflow (version 0.26.x or higher) can be used on any POSIX compatible system (Linux, OS X, etc). It requires BASH and 
Java 8 (or higher) to be installed. More instructions are available [here](https://www.nextflow.io/docs/latest/getstarted.html).

* **Container Engine**

All components of DEN-IM are executed in docker containers, which means that you’ll need to have a container engine 
installed. The container engines available are the ones supported by Nextflow:
    - [Docker](https://www.nextflow.io/docs/latest/docker.html),
    - [Singularity](https://www.nextflow.io/docs/latest/singularity.html),
    - Shifter (undocumented)

If you already have any one of these installed, you are good to go as the provided docker containers are compatible 
with all engines available. If not, you’ll need to install one.

If you want a local installation of DEN-IM, you can clone this repository with 
`git clone https://github.com/assemblerflow/DEN-IM.git`, and all files will e cloned to your local machine.
Alternatively, you can run DEN-IM directly with nextflow with the command `nextflow run assemblerflow/DEN-IM`.



## Running DEN-IM

If you have a local installation of DEN-IM, you can get started using the workflow with:

`nextflow run DEN-IM.nf --fastq /path/to/input/files/*_{1,2}.fq.gz`

By default nextflow executes DEN-IM with singularity, but this can be changed by adding `-pofile docker` to the command.
Alternatively, you can run DEN-IM directly with nextflow with the command 

`nextflow run assemblerflow/DEN-IM --fastq /path/to/input/files/*_{1,2}.fq.gz`

Users can customize the workflow execution either by using command line options or by modifying a simple plain-text 
configuration file (`params.config`), where parameters are set as key-value pairs. The version of tools used can also 
be changed by providing new container tags in the appropriate configuration file (`containers.config`), as well as the 
resources for each process (`resources.config`).

More information on how to personalize the config files can be found in the [wiki](https://github.com/assemblerflow/DEN-IM/wiki/How-to-Run-DEN-IM).



## Output and Report

The output files are stored in the `results` folder in the directory where the workflow was executed. 
The nextflow log file for the execution of the pipeline can be found in the directory of execution. Log files for each
of the components in the workflow are stored inside the `results` folder.
DEN-IM creates an **interactive HTML report**, stored in the `pipeline_results` folder in the directory where the 
workflow was executed. To open the report simply click oh the **pipeline_report.html** file and the report will open on
 your default browser. 