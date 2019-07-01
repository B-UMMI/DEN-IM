# DEN-IM

![DEN-IM_logo](https://user-images.githubusercontent.com/15690332/60383483-47e5d900-9a69-11e9-9808-3a1d9023a32a.png)

DEN-IM is a  one-stop, user-friendly, containerised and reproducible workflow for the analysis of DENV sequencing data, 
both from shotgun and targeted metagenomics approaches.
Is is able to infer DENV coding sequence (CDS), identify serotype and genotype, and generate a phylogenetic tree.

As input DEN-IM accepts **raw paired-end Illumina sequencing reads**, and and informs the user with an interactive 
HTML report with information on the quality control, mapping, assembly typing and phylogenetic analysis, as well as 
all the output files of the whole pipeline.

It is implemented in Nextflow, a workflow management system that allows the effortless deployment and execution of 
complex distributed computational workflows. DEN-IM integrates Docker containerised images for all the tools necessary 
for its execution, ensuring reproducibility and the tracking of both software code and version, regardless of the 
operating system used.


## Installation

Before installing DEN-IM, a few dependencies must be installed in your system:

* **Nextflow**

Nextflow (version 0.31.0 or higher) can be used on any POSIX compatible system (Linux, OS X, etc). It requires BASH and 
Java 8 (or higher) to be installed. More instructions are available [here](https://www.nextflow.io/docs/latest/getstarted.html).

* **Container Engine**

All components of DEN-IM are executed in docker containers, which means that you’ll need to have a container engine 
installed. The container engines available are the ones supported by Nextflow:

* [Docker](https://www.nextflow.io/docs/latest/docker.html),
* [Singularity](https://www.nextflow.io/docs/latest/singularity.html),
* [Shifter](https://github.com/NERSC/shifter) (undocumented)

If you already have any one of these installed, you are good to go as the provided docker containers are compatible 
with all engines available. If not, you’ll need to install one.

* **Local Instalation** 

A local installation of DEN-IM is also required as Nextflow's pipeline sharing is not compatible with DEN-IM's 
architecture. You can clone this repository with `git clone https://github.com/B-UMMI/DEN-IM.git`, and all 
files will be in your local machine.


## Running DEN-IM

After you have a local installation of DEN-IM, you can get started using the workflow with:

`nextflow run DEN-IM.nf --fastq /path/to/input/files/*_{1,2}.fq.gz`

By default nextflow executes DEN-IM with singularity, but this can be changed by adding `-pofile docker` to the command.

Users can customize the workflow execution either by using command line options or by modifying a simple plain-text 
configuration file (`params.config`), where parameters are set as key-value pairs. The version of tools used can also 
be changed by providing new container tags in the appropriate configuration file (`containers.config`), as well as the 
resources for each process (`resources.config`).

More information on how to personalize the config files can be found in the [wiki](https://github.com/B-UMMI/DEN-IM/wiki/How-to-Run-DEN-IM).


## Output and Report

The output files are stored in the `results/` folder in the directory where the workflow was executed. 
The nextflow log file for the execution of the pipeline can be found in the directory of execution. Log files for each
of the components in the workflow are stored inside the `results/` folder.
DEN-IM creates an **interactive HTML report**, stored in the `pipeline_results/` folder in the directory where the 
workflow was executed. To open the report simply click oh the **pipeline_report.html** file and the report will open on 
your default browser. 


## Citation and Contacts

DEN-IM was developed by the [Molecular Microbiology and Infection Unit](http://darwin.phyloviz.net/wiki/doku.php), at 
Instituto de Medicina Molecular, in collaboration with [Personalised Microbiology Lab](https://rossenlab.com/) at UMCG. 
A pre-print is available at [biorxiv](https://www.biorxiv.org/content/10.1101/628073v1).

