# script to build the UBE3A networks

import config
import network_classes_igraph_pub

# function that deletes all edges and resulting orphan nodes from network that
# have less than 2 origins
def delete_edges_with_one_origin(g):

	edges_to_delete = []
	for i,edge in enumerate(g.get_edgelist()):
		if len(g.es[i]['origins']) < 2:
			edges_to_delete.append(edge)
	g.delete_edges(edges_to_delete)

	orphan_nodes = g.vs.select(_degree = 0)
	orphan_nodes_indices = [v.index for v in orphan_nodes]
	g.delete_vertices(orphan_nodes_indices)

	return g


# function that updates the node source info in a given network to
# highlight nodes as preys even if they were excluded from the seed gene file such as
# the proteasome subunits
def update_source_info_to_network_from_excluded_preys(g,attr_file,attributes,cols):

	preys = {}
	file1 = open(attr_file,'r')
	entries = file1.readlines()
	file1.close()
	for line in entries[1:]:
		tab_list = str.split(line[:-1],'\t')
		geneID = tab_list[0]
		preys[geneID] = tab_list

	for node in g.vs:
		if node['name'] in preys and node['source'] != 'seed':
			node['source'] = 'seed'
			for i,attr in enumerate(attributes):
				node[attr] = preys[node['name']][cols[i]]


# build UBE3A network
def get_UBE3A_network(source_networks,seed_genes,entrez_file,seed_gene_file,complete_prey_file,out_prefix):

	print out_prefix
	g = network_classes_igraph_pub.Graph(n=0,vertex_attrs={'name':[]})
	g.get_source_network(source_networks)
	print 'number of nodes in source network:', len(g.vs)
	print 'number of edges in source network:', len(g.es)
	g.get_entrez_gene_annotations_for_nodes(entrez_file)
	seed_genes_in_g = g.get_seed_genes_in_network(seed_genes)
	print 'number of seed genes in source network:', len(seed_genes_in_g)
	gsub_sn2s = g.get_network_with_neighbors_linked_to_2_seed_genes(seed_genes_in_g,1)
	network_classes_igraph_pub.add_node_attributes(gsub_sn2s,seed_gene_file,[3,4,5,6,7,8])
	update_source_info_to_network_from_excluded_preys(gsub_sn2s,complete_prey_file,['num_hek','num_sy5y','frac_occ','cell_lines','hcip_cutoff_sy5y'],[3,4,5,6,7])
	update_source_info_to_network_from_excluded_preys(gsub_sn2s,seed_gene_file,['Y2H_int'],[8])
	gsub_sn2s = delete_edges_with_one_origin(gsub_sn2s)
	seed_nodes = gsub_sn2s.vs.select(source_eq='seed')

	# update the source info to distinguish UBE3A seeds coming from Y2H or pulldown
	for seed_node in seed_nodes:
		# is this the bait?
		if seed_node['name'] == '7337':
			seed_node['source'] = 'bait'
		# is this seed node only been found as interactor of UBE3A from Y2H screen?
		elif seed_node['Y2H_int'] == 'y' and seed_node['cell_lines'] == 'yeast':
			seed_node['source'] = 'Y2H_seed'
		# is this seed node identified as interactor of UBE3A from Y2H screen and pulldown?
		elif seed_node['Y2H_int'] == 'y' and seed_node['cell_lines'].find('_yeast') > -1:
			seed_node['source'] = 'Y2H_APMS_seed'
		# is this seed node identified as interactor of UBE3A only in pulldown?
		elif (seed_node['Y2H_int'] == 'n' or seed_node['Y2H_int'] == None) and len(seed_node['cell_lines']) == 1:
			seed_node['source'] = 'APMS_seed'
		else:
			continue

	print 'number of nodes in UBE3A network:', len(gsub_sn2s.vs)
	print 'number of edges in UBE3A network:', len(gsub_sn2s.es)
	outfile_nodes = config.output_path + out_prefix + '.node_attributes.txt'
	network_classes_igraph_pub.write_node_attributes(gsub_sn2s,outfile_nodes)
	outfile_edges = config.output_path + out_prefix + '.edge_attributes.txt'
	network_classes_igraph_pub.write_edge_attributes(gsub_sn2s,outfile_edges)

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

	# read in UBE3A seed genes with proteasome subunits
	complete_prey_file = config.output_path + 'UBE3A_seed_file_with_proteasome.txt'
	prot_seed_genes = network_classes_igraph_pub.read_seed_genes(complete_prey_file,0)
	print 'Number of UBE3A seed genes with proteasome subunits:', len(prot_seed_genes)

	# read in UBE3A seed genes without proteasome subunits
	seed_gene_file = config.output_path + 'UBE3A_seed_file.txt'
	seed_genes = network_classes_igraph_pub.read_seed_genes(seed_gene_file,0)
	print 'Number of UBE3A seed genes without proteasome subunits:', len(seed_genes)

	# read in UBE3A seed genes without proteasome subunits and without Y2H seeds
	noY2H_seed_gene_file = config.output_path + 'UBE3A_seed_file_no_Y2H.txt'
	noY2H_seed_genes = network_classes_igraph_pub.read_seed_genes(noY2H_seed_gene_file,0)
	print 'Number of UBE3A seed genes without proteasome subunits and without Y2H interactors:', len(noY2H_seed_genes)

	### build network for preys and all neighbors linked to 2 seeds excluding proteasome related preys, QBCHL as source,
	### keep only edges with at least 2 evidences
	out_prefix = 'UBE3A_network_QBCHL'
	source_networks = [bioplex_tup,ccsb_tup,litbm_tup,qubic_tup,cofrac_tup]
	get_UBE3A_network(source_networks,seed_genes,entrez_file,seed_gene_file,complete_prey_file,out_prefix)

	### build network for preys and all neighbors linked to 2 seeds excluding proteasome related preys, QBCHL as source,
	### as well as additional more focussed and smaller scale sources
	### keep only edges with at least 2 evidences
	out_prefix = 'UBE3A_network_compl'
	source_networks = [bioplex_tup,ccsb_tup,litbm_tup,qubic_tup,cofrac_tup,howley_tup,Y2H_UBE3A_tup,harper_tup,recipIP_tup,AlHakim_tup]
	get_UBE3A_network(source_networks,seed_genes,entrez_file,seed_gene_file,complete_prey_file,out_prefix)
