# script to generate random graphs

import sys
import os
import igraph
import network_classes_igraph_pub
import config

def get_real_network(source):

	nodes = set()
	edge_list = []

	file1 = open(source,'r')
	entries = file1.readlines()
	file1.close()
	for line in entries[1:]:
		tab_list = str.split(line[:-1],'\t')
		gene_a = tab_list[0]
		gene_b = tab_list[1]
		if gene_a < gene_b:
			edge = (gene_a,gene_b)
		else:
			edge = (gene_b,gene_a)
		edge_list.append(edge)

	for edge in edge_list:
		nodes.add(edge[0])
		nodes.add(edge[1])

	g = igraph.Graph()
	# add all nodes to the network
	g.add_vertices(list(nodes))
	# add all edges to the network
	g.add_edges(edge_list)

	return g


def get_random_graph(real_network,outfile):

	degree_seq = real_network.degree()
	rand_graph = igraph.Graph.Degree_Sequence(degree_seq,method='vl')
	rand_graph.vs['name'] = real_network.vs['name']
	rand_edges = rand_graph.get_edgelist()
	rand_edges_mapped = []
	for edge in rand_edges:
		gene_a = int(real_network.vs[edge[0]]['name'])
		gene_b = int(real_network.vs[edge[1]]['name'])
		if gene_a < gene_b:
			rand_edges_mapped.append((gene_a,gene_b))
		else:
			rand_edges_mapped.append((gene_b,gene_a))

	target = open(outfile,'w')
	target.write('gene_a\tgene_b\n')
	for edge in rand_edges_mapped:
		target.write(str(edge[0]) + '\t' + str(edge[1]) + '\n')
	target.close()


def get_random_graphs(source,outpath,source_name,num_graphs):

	real_network = get_real_network(source)

	for i in range(0,num_graphs):
		if i%100 == 0:
			print i
		outfile = outpath + source_name + '_rand_network_' + str(i) + '.txt'
		get_random_graph(real_network,outfile)


if '__main__' == __name__:


	num_graphs = int(sys.argv[1])

	bioplex_tup = ('BioPlex',config.output_path + 'BioPlex.txt')
	ccsb_tup = ('HI-union',config.output_path + 'HI-union.txt')
	litbm_tup = ('Lit-BM-13',config.output_path + 'Lit-BM-13.txt')
	qubic_tup = ('QUBIC',config.output_path + 'QUBIC.txt')
	cofrac_tup = ('CoFrac',config.output_path + 'CoFrac.txt')

	source_networks = [bioplex_tup,ccsb_tup,litbm_tup,qubic_tup,cofrac_tup]
	g = network_classes_igraph_pub.Graph(n=0,vertex_attrs={'name':[]})
	g.get_source_network(source_networks)

	for i in range(0,num_graphs):
		if i%100 == 0:
			print i
		outfile = config.rand_path + config.rand_prefix + str(i) + '.txt'
		get_random_graph(g,outfile)
