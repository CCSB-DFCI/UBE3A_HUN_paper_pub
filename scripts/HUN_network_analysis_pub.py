# script to build networks for HUN complex

import config
import network_classes_igraph_pub

# function that gets a graph and deletes all edges between nodes that don't share a given
# annotation
def delete_edges_between_preys_not_same_attribute(g,attr,separator):

	edges_to_remove = []
	for edge_tup in g.get_edgelist():
		n1_att = set(g.vs[edge_tup[0]][attr].split(separator))
		n2_att = set(g.vs[edge_tup[1]][attr].split(separator))
		if len(n1_att.intersection(n2_att)) == 0:
			edges_to_remove.append(edge_tup)

	g.delete_edges(edges_to_remove)

	orphan_nodes = g.vs.select(_degree = 0)
	orphan_nodes_indices = [v.index for v in orphan_nodes]
	g.delete_vertices(orphan_nodes_indices)

	return g


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
	seed_gene_file = config.output_path + 'HUN_complex_seed_file.txt'
	column = 0
	seed_genes = network_classes_igraph_pub.read_seed_genes(seed_gene_file,column)
	print 'number of seed genes:', len(seed_genes)


	### build network for HUN complex by linking the preys with each other if they
	### were pulled down with the same bait
	### use QBCHL as source
	source_networks = [bioplex_tup,ccsb_tup,litbm_tup,qubic_tup,cofrac_tup]
	g = network_classes_igraph_pub.Graph(n=0,vertex_attrs={'name':[]})
	g.get_source_network(source_networks)
	print 'number of nodes:', len(g.vs)
	print 'number of edges:', len(g.es)
	g.get_entrez_gene_annotations_for_nodes(entrez_file)
	seed_genes_in_g = g.get_seed_genes_in_network(seed_genes)
	print 'number of seed genes in network QBCHL:', len(seed_genes_in_g)
	gsub_s = g.get_network_connecting_seed_genes(seed_genes_in_g)
	print len(gsub_s.vs), len(gsub_s.es)
	network_classes_igraph_pub.add_node_attributes(gsub_s,seed_gene_file,[3,4,5,6])
	gsub_s = delete_edges_between_preys_not_same_attribute(gsub_s,'baits','|')
	print len(gsub_s.vs), len(gsub_s.es)
	# add node attribute that indicates which proteins were prey from UBE3A
	UBE3A_file = config.output_path + 'UBE3A_seed_file.txt'
	UBE3A_preys = set(network_classes_igraph_pub.read_seed_genes(UBE3A_file,0))
	gsub_s.vs['UBE3A_prey'] = [None for i in range(0,len(gsub_s.vs))]
	for node in gsub_s.vs:
		if node['name'] in UBE3A_preys:
			node['UBE3A_prey'] = '1'
	print 'number of nodes in network QBCHL:', len(gsub_s.vs)
	print 'number of edges in network QBCHL:', len(gsub_s.es)
	outfile_nodes = config.output_path + 'HUN_network_QBCHL.node_attributes.txt'
	network_classes_igraph_pub.write_node_attributes(gsub_s,outfile_nodes)
	outfile_edges = config.output_path + 'HUN_network_QBCHL.edge_attributes.txt'
	network_classes_igraph_pub.write_edge_attributes(gsub_s,outfile_edges)


	### build network for HUN complex by linking the preys with each other if they
	### were pulled down with the same bait
	### use QBCHL as source and additional smaller scale sources
	source_networks = [bioplex_tup,ccsb_tup,litbm_tup,qubic_tup,cofrac_tup,recipIP_tup,AlHakim_tup,Y2H_UBE3A_tup]
	g = network_classes_igraph_pub.Graph(n=0,vertex_attrs={'name':[]})
	g.get_source_network(source_networks)
	print 'number of nodes:', len(g.vs)
	print 'number of edges:', len(g.es)
	g.get_entrez_gene_annotations_for_nodes(entrez_file)
	seed_genes_in_g = g.get_seed_genes_in_network(seed_genes)
	print 'number of seed genes in network compl:', len(seed_genes_in_g)
	gsub_s = g.get_network_connecting_seed_genes(seed_genes_in_g)
	print len(gsub_s.vs), len(gsub_s.es)
	network_classes_igraph_pub.add_node_attributes(gsub_s,seed_gene_file,[3,4,5,6])
	gsub_s = delete_edges_between_preys_not_same_attribute(gsub_s,'baits','|')
	print len(gsub_s.vs), len(gsub_s.es)
	# add node attribute that indicates which proteins were prey from UBE3A
	UBE3A_file = config.output_path + 'UBE3A_seed_file.txt'
	UBE3A_preys = set(network_classes_igraph_pub.read_seed_genes(UBE3A_file,0))
	gsub_s.vs['UBE3A_prey'] = [None for i in range(0,len(gsub_s.vs))]
	for node in gsub_s.vs:
		if node['name'] in UBE3A_preys:
			node['UBE3A_prey'] = '1'
	print 'number of nodes in network compl:', len(gsub_s.vs)
	print 'number of edges in network compl:', len(gsub_s.es)
	outfile_nodes = config.output_path + 'HUN_network_compl.node_attributes.txt'
	network_classes_igraph_pub.write_node_attributes(gsub_s,outfile_nodes)
	outfile_edges = config.output_path + 'HUN_network_compl.edge_attributes.txt'
	network_classes_igraph_pub.write_edge_attributes(gsub_s,outfile_edges)
