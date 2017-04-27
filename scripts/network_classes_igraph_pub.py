
# script extending the Graph class from igraph python module

import igraph
import numpy
import random
import copy

class Graph(igraph.Graph):

	space = set()

	# function that adds edges and the respective nodes from the given
	# source networks
	def get_source_network(self,source_networks):

		# stores all geneIDs
		nodes = set()
		# stores for every edge a set of source names
		edge_dict = {}

		for tup in source_networks:
			source_name = tup[0]
			source = tup[1]
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
				if edge not in edge_dict:
					edge_dict[edge] = set()
				edge_dict[edge].add(source_name)
		for edge,origins in edge_dict.iteritems():
			nodes.add(edge[0])
			nodes.add(edge[1])
		# add all nodes to the network
		self.add_vertices(list(nodes))
		# add all edges to the network together with the origins attribute
		items = edge_dict.items()
		edges = [tup[0] for tup in items]
		ori_lists = [list(tup[1]) for tup in items]
		self.add_edges(edges)
		self.es['origins'] = ori_lists


	# function that tries to add to all nodes of the network annotations from
	# the entrez DB
	def get_entrez_gene_annotations_for_nodes(self,infile):

		annot_dict = {}
		file1 = open(infile,'r')
		entries = file1.readlines()
		file1.close()
		for row in entries[1:]:
			tab_list = str.split(row[:-1],'\t')
			geneID = str(tab_list[1])
			symbol = tab_list[2]
			synonyms = tab_list[4]
			descr = tab_list[8]
			annot_dict[geneID] = (symbol,synonyms,descr)

		for node in self.vs:
			geneID = node['name']
			if geneID in annot_dict:
				node['symbol'] = annot_dict[geneID][0]
				node['synonyms'] = annot_dict[geneID][1].split('|')
				node['descr'] = annot_dict[geneID][2]


	# get all seed genes that are part of the source network
	def get_seed_genes_in_network(self,seed_genes):

		in_network_sg = []
		for geneID in seed_genes:
			if geneID in self.vs['name']:
				in_network_sg.append(geneID)

		return in_network_sg


	# function that returns a graph object containing all the seed genes and edges between them
	# does not contain seed genes without an edge to any other seed gene
	# needs as input a list of seed genes that are all part of the source network
	def get_network_connecting_seed_genes(self,seed_genes):

		# get sub network with all edges between the seed genes that are part of the source
		# network
		gsub = self.subgraph(seed_genes)

		# delete from the subnetwork all nodes that are not linked
		orphan_nodes = gsub.vs.select(_degree = 0)
		orphan_nodes_indices = [v.index for v in orphan_nodes]
		gsub.delete_vertices(orphan_nodes_indices)

		return gsub


	# function that returns a network containing all seed genes and their neighbors with
	# links between them
	def get_network_with_neighbors_of_seed_genes(self,seed_genes,distance):

		# get the neighbors of the seed genes
		neighbor_lists = self.neighborhood(vertices=seed_genes,order=distance)

		# build a list of vertex names of seed genes and neighbors
		nodes = set()
		for i,neighbor_list in enumerate(neighbor_lists):
			for neighbor in neighbor_list:
				nodes.add(self.vs[neighbor]['name'])
		for seed_gene in seed_genes:
			nodes.add(seed_gene)

		# get the subgraph spanning these vertices
		gsub = self.subgraph(nodes)

		# add an attribute discriminating between neighbors and seed genes
		for seed_gene in seed_genes:
			gsub.vs.find(seed_gene)['source'] = 'seed'
		for node in nodes:
			if not gsub.vs.find(node)['source']:
				gsub.vs.find(node)['source'] = 'extended'

		return gsub


	# function that returns a network containing all seed genes and their neighbors if
	# these neighbors are linked to at least 2 seed genes
	def get_network_with_neighbors_linked_to_2_seed_genes(self,seed_genes,distance):

		gsub = self.get_network_with_neighbors_of_seed_genes(seed_genes,distance)

		# get all neighbors that are linked to less than 2 seed genes
		neighbors_to_remove = []
		for node in gsub.vs.select(source='extended'):
			NoN_nodes = gsub.neighbors(node)
			count_seed_genes = 0
			i = 0
			while count_seed_genes < 2 and i < len(NoN_nodes):
				if gsub.vs[NoN_nodes[i]]['source'] == 'seed':
					count_seed_genes = count_seed_genes + 1
				i = i + 1
			if count_seed_genes < 2:
				neighbors_to_remove.append(node.index)

		# remove these neighbors
		gsub.delete_vertices(neighbors_to_remove)

		# remove orphan nodes
		orphan_nodes = gsub.vs.select(_degree = 0)
		orphan_nodes_indices = [v.index for v in orphan_nodes]
		gsub.delete_vertices(orphan_nodes_indices)

		return gsub


	# function that updates the node attribute of given attr_name and list of genes names
	# with the given value
	def update_node_attribute(self,attr_name,genes,value):

		if attr_name not in self.vs.attributes():
			self.vs[attr_name] = [None for i in range(0,len(self.vs))]

		for gene in genes:
			node_list = self.vs.select(name=gene)
			for node in node_list:
				node[attr_name] = value


# function that writes out the nodes and their attributes from the given graph object in cytoscape readable format
def write_node_attributes(graph,outfile):

	target = open(outfile,'w')
	node_attributes = graph.vertex_attributes()
	node_attributes.sort()
	target.write('geneID')
	for attr in node_attributes:
		if attr != 'name':
			target.write('\t' + attr)
	target.write('\n')

	for node in graph.vs:
		target.write(node['name'])
		for attr in node_attributes:
			if attr != 'name':
				if isinstance(node[attr],list):
					outstring = ''
					for value in node[attr]:
						outstring = outstring + str(value) + '|'
					outstring = outstring[:-1]
				else:
					outstring = str(node[attr])
				target.write('\t' + outstring)
		target.write('\n')

	target.close()


# function that writes out the edges and their attributes from the given graph object in cytoscape readable format
def write_edge_attributes(graph,outfile):

	target = open(outfile,'w')
	edge_attributes = graph.edge_attributes()
	edge_attributes.sort()
	target.write('geneID_A\tgeneID_B')
	for attr in edge_attributes:
		target.write('\t' + attr)
	target.write('\n')

	for edge in graph.es:
		source_node = int(graph.vs[edge.source]['name'])
		target_node = int(graph.vs[edge.target]['name'])
		if source_node < target_node:
			target.write(str(source_node) + '\t' + str(target_node))
		else:
			target.write(str(target_node) + '\t' + str(source_node))
		for attr in edge_attributes:
			if isinstance(edge[attr],list):
				outstring = ''
				for value in edge[attr]:
					outstring = outstring + str(value) + '|'
				outstring = outstring[:-1]
			else:
				outstring = str(edge[attr])
			target.write('\t' + outstring)
		target.write('\n')

	target.close()


# function that gets a file from which it reads in the seed genes from the
# given column
def read_seed_genes(infile,column):

	genes = set()
	file1 = open(infile,'r')
	entries = file1.readlines()
	file1.close()

	for line in entries[1:]:
		tab_list = str.split(line[:-1],'\t')
		genes.add(tab_list[column])

	return list(genes)


# function that adds attributes to the nodes in the given graph from the file for
# the given columns
def add_node_attributes(graph,attribute_file,cols):

	file1 = open(attribute_file,'r')
	entries = file1.readlines()
	file1.close()
	header = entries[0][:-1].split('\t')

	for col in cols:
		attr_name = header[col]
		if attr_name not in graph.vs.attributes():
			graph.vs[attr_name] = [None for i in range(0,len(graph.vs))]

	for line in entries[1:]:
		tab_list = str.split(line[:-1],'\t')
		geneID = tab_list[0]
		values = tab_list
		for i in cols:
			node_list = graph.vs.select(name=geneID)
			for node in node_list:
				node[header[i]] = values[i]


# function that adds attributes to the edges in the given graph from the file for
# the given columns
def add_edge_attributes(graph,attribute_file,cols):

	file1 = open(attribute_file,'r')
	entries = file1.readlines()
	file1.close()
	header = entries[0][:-1].split('\t')

	for col in cols:
		attr_name = header[col]
		if attr_name not in graph.es.attributes():
			graph.es[attr_name] = [None for i in range(0,len(graph.es))]

	for line in entries[1:]:
		tab_list = str.split(line[:-1],'\t')
		gene_a = tab_list[0]
		gene_b = tab_list[1]
		node_list = graph.vs.select(name=gene_a)
		if len(node_list) == 0:
			continue
		index_a = node_list[0].index
		node_list = graph.vs.select(name=gene_b)
		if len(node_list) == 0:
			continue
		index_b = node_list[0].index
		if graph.are_connected(index_a,index_b):
			edge_id = graph.get_eid(index_a,index_b)
			values = tab_list
			for i in cols:
				graph.es[edge_id][header[i]] = values[i]
