# script to build networks around CAMK2D

import config
import network_classes_igraph_pub


if __name__ == '__main__':

	bioplex_tup = ('BioPlex',config.output_path + 'BioPlex.txt')
	ccsb_tup = ('HI-union',config.output_path + 'HI-union.txt')
	litbm_tup = ('Lit-BM-13',config.output_path + 'Lit-BM-13.txt')
	qubic_tup = ('QUBIC',config.output_path + 'QUBIC.txt')
	cofrac_tup = ('CoFrac',config.output_path + 'CoFrac.txt')
	howley_tup = ('Howley_map',config.output_path + 'Howley_map.txt')
	harper_tup = ('DuB_AP',config.output_path + 'DuB_AP.txt')
	Y2H_UBE3A_tup = ('Y2H_UBE3A',config.output_path + 'Y2H_UBE3A.txt')
	recipIP_tup = ('recipIP',config.output_path + 'reciprocal_Howley_IPs.txt')
	AlHakim_tup = ('AlHakimIP',config.output_path + 'AlHakim_IPs.txt')

	entrez_file = config.input_path + 'entrez_gene.gene_info_human_20161113.tsv'

	# read in seed genes with entrez geneID in string format
	seed_gene_file = config.output_path + 'CAMK2D_seed_file.txt'
	column = 0
	seed_genes = network_classes_igraph_pub.read_seed_genes(seed_gene_file,column)
	print 'number of seed genes:', len(seed_genes)

	# get all the HUN complex preys
	HUN_seed_file = config.output_path + 'HUN_UBE3A_seed_file.txt'
	HUN_preys = set(network_classes_igraph_pub.read_seed_genes(HUN_seed_file,0))
	print 'number of HUN complex preys:', len(HUN_preys)

	# get a dict with prey and cell line info for UBE3A preys (because they are not
	# part of the HUN seed gene file)
	UBE3A_dict = {}
	UBE3A_file = config.output_path + 'UBE3A_seed_file.txt'
	file1 = open(UBE3A_file,'r')
	entries = file1.readlines()
	file1.close()
	for line in entries[1:]:
		tab_list = str.split(line[:-1],'\t')
		geneID = tab_list[0]
		cell_lines = tab_list[6]
		UBE3A_dict[geneID] = cell_lines

	### build a network that connects the CAMK2D preys to the preys from the HUN complex
	### using once only QBCHL and once QBCHL and the additional sources for edges
	for i in range(2):

		if i == 0:
			source_network_name = 'QBCHL'
			source_networks = [bioplex_tup,ccsb_tup,litbm_tup,qubic_tup,cofrac_tup]
		else:
			source_network_name = 'compl'
			source_networks = [bioplex_tup,ccsb_tup,litbm_tup,qubic_tup,cofrac_tup,Y2H_UBE3A_tup,recipIP_tup,AlHakim_tup]

		print source_network_name
		g = network_classes_igraph_pub.Graph(n=0,vertex_attrs={'name':[]})
		g.get_source_network(source_networks)
		g.get_entrez_gene_annotations_for_nodes(entrez_file)
		print 'number of nodes:', len(g.vs)
		print 'number of edges:', len(g.es)
		seed_genes_in_g = g.get_seed_genes_in_network(seed_genes)
		print 'number of seed genes in network:', len(seed_genes_in_g)
		gsub = g.get_network_with_neighbors_of_seed_genes(seed_genes_in_g,1)
		del gsub.vs['source']
		# delete all edges that are not between two CAMK2D preys or that do not involve 1 CAMK2D prey
		# and 1 HUN complex prey
		edges_to_remove = []
		for edge_tup in gsub.get_edgelist():
			gene_a = gsub.vs[edge_tup[0]]['name']
			gene_b = gsub.vs[edge_tup[1]]['name']
			if not (gene_a in seed_genes and gene_b in seed_genes) and \
			   not (gene_a in seed_genes and gene_b in HUN_preys) and \
			   not (gene_b in seed_genes and gene_a in HUN_preys):
				edges_to_remove.append(edge_tup)
		gsub.delete_edges(edges_to_remove)
		orphan_nodes = gsub.vs.select(_degree = 0)
		orphan_nodes_indices = [v.index for v in orphan_nodes]
		gsub.delete_vertices(orphan_nodes_indices)
		# add HCIP cutoff for CAMK2D preys and baits and cell lines info from the HUN
		# seel file to the nodes
		network_classes_igraph_pub.add_node_attributes(gsub,seed_gene_file,[7])
		HUN_file = config.output_path + 'HUN_complex_seed_file.txt'
		network_classes_igraph_pub.add_node_attributes(gsub,HUN_file,[4,6])
		# merge the node attribute info on the cell lines and baits
		for node in gsub.vs:
			geneID = node['name']
			if geneID in seed_genes_in_g:
				if node['baits'] is None:
					node['baits'] = 'CAMK2D'
				else:
					node['baits'] = node['baits'] + '|CAMK2D'
				if node['cell_lines'] is None:
					node['cell_lines'] = 's'
				elif node['cell_lines'] == 'h':
					node['cell_lines'] = 'b'
			if geneID in UBE3A_dict:
				if node['baits'] is None:
					node['baits'] = 'UBE3A'
				else:
					node['baits'] = node['baits'] + '|UBE3A'
				if node['cell_lines'] is None:
					node['cell_lines'] = UBE3A_dict[geneID]
				elif node['cell_lines'] == 'h' and UBE3A_dict[geneID] != 'h':
					node['cell_lines'] = 'b'
		print 'number of nodes:', len(gsub.vs)
		print 'number of edges:', len(gsub.es)
		# write out
		outfile_nodes = config.output_path + 'CAMK2D_network_' + source_network_name + '.node_attributes.txt'
		network_classes_igraph_pub.write_node_attributes(gsub,outfile_nodes)
		outfile_edges = config.output_path + 'CAMK2D_network_' + source_network_name + '.edge_attributes.txt'
		network_classes_igraph_pub.write_edge_attributes(gsub,outfile_edges)
