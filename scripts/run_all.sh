# command line script to produce all data in this paper
#! /bin/bash

# generate prey files for every bait
python UBE3A_seed_gene_file_assembly.py
python HERC2_seed_gene_file_assembly.py
python seed_gene_file_assembly.py
# generate the source networks
python source_network_mapping.py HI-union,BioPlex,QUBIC,CoFrac,Lit-BM-13,DuB_AP,Howley_map
# generate the networks around the baits of interest
python UBE3A_network_analysis_pub.py
python HUN_network_analysis_pub.py
python CAMK2D_network_analysis_pub.py
# assess functional significance of networks and sets of preys
python significance_test_GO_network_pub.py 1000 QBCHL UBE3A
python significance_test_GO_network_pub.py 1000 QBCHL CAMK2D
python run_Funcassociate_pub.py HUN_network_QBCHL.node_attributes.txt 0
python run_Funcassociate_pub.py UBE3A_seed_file.txt 0
python run_Funcassociate_pub.py UBE3A_seed_file_with_proteasome.txt 0
python run_Funcassociate_pub.py CAMK2D_seed_file.txt 0
python run_Funcassociate_pub.py MAPK6_seed_file.txt 0
python run_Funcassociate_pub.py ECI2_seed_file.txt 0
python run_Funcassociate_pub.py ECH1_seed_file.txt 0
python run_Funcassociate_pub.py HERC2_seed_file.txt 0
python run_Funcassociate_pub.py NEURL4_seed_file.txt 0
python run_Funcassociate_pub.py HIF1AN_seed_file.txt 0
# generate degree-controlled randomized QBCHL networks
python random_graphs.py 1000
# compute significances on pulldown data and HUN network using
# randomized QBCHL networks
python significance_tests_pub.py seed_lcc 1000
python significance_tests_pub.py hun_lcc 1000
python significance_tests_pub.py camk2d_hun 1000
# produce plots
jupyter nbconvert --to notebook --execute UBE3A_HUN_paper_plots.ipynb
rm UBE3A_HUN_paper_plots.nbconvert.ipynb
