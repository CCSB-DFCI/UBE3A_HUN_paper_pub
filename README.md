# UBE3A_HUN_paper_pub
This repository contains all code and files that accompany the following publication:
<citation>
<DOI>

This code was built and executed on a MacBook Pro under MacOS Sierra with Python 2.7.8 and Jupyter 4.3.1. All additional python packages that need to be installed in order to run all code can be found in requirements.txt.


The shell script run_all.sh can be run to execute all scripts in the right order, which will produce all files that were used in this publication to support the results. Almost all of these files are already present in folder "output_files" and folder "plots". Be aware, some scripts will run for several hours. If multiple processors are available, some of the especially time-consuming scripts can be run in parallel.


A few input and output files are not part of this repository. All input files necessary to build the QBCHL interactome are not provided in this repository, primarily for copyright and limits of file sizes. These input files can be downloaded from the following sources and copied into the folder input_files prior to running source_network_mapping.py or any script thereafter.

Hi-union: download all published, all unpublished, test space verified and test space validated datasets from interactome.baderlab.org

BioPlex: download the dataset with time stemp 06/12/2015 from http://bioplex.hms.harvard.edu/downloadInteractions.php

CoFrac: Supplementary table 2 from Wan et al (PubmedID: 26344197)

QUBIC: Supplementary table 2 from Hein et al (PubmedID: 26496610)

Lit-BM-13: Supplementary table 1 from Rolland et al (PubmedID: 25416956)

The second set of files that are not available in this repository (due to space limits) but can be created with the script random_graphs.py are the randomized QBCHL networks.
