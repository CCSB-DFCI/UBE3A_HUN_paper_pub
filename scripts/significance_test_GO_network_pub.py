# function that randomly chooses the same number of genes as UBE3A seed genes with
# the same degree distribution to then rebuilt the network the same way we did for the
# UBE3A network and to then test which GO terms are enriched in this network
# for every GO term that was seen to be enriched in the original UBE3A network, we need
# to count how often this term was also observed to be enriched in the random control
# and how often this term was at least as enriched as in the real network
# output a table with all GO terms that were enriched in the original network and how
# they do in the control networks
# a table that lists for every control network the number of seed genes that were in the
# final network, the total number of genes in the network, min, 25th percentil, median,
# 75th percentile, max number of genes for enriched GO terms, total number of enriched
# GO terms, LCC size, number of random seed genes that overlap with original seed genes

import random
import sys
import math
import pyjsonrpc
import network_classes_igraph_pub
import UBE3A_network_analysis_pub
import run_Funcassociate_pub
import config


class GOterm:

	def __init__(self,GOID):

		self.GOID = GOID
		self.descr = ''
		self.num_genes = None
		self.num_real_seed = None
		self.real_LOD = None
		self.real_pvalue = None
		self.real_adj_pvalue = None
		self.rand_LODs = []
		self.rand_pvalues = []
		self.nums_rand_seed = []


def get_subnetwork(mode,seeds,g):

	if mode == 'UBE3A':
		ube3a_network = g.get_network_with_neighbors_linked_to_2_seed_genes(seeds,1)
		ube3a_network = UBE3A_network_analysis_pub.delete_edges_with_one_origin(ube3a_network)
	elif mode == 'CAMK2D':
		HUN_preys = set()
		HUN_baits = ['UBE3A','HERC2','NEURL4','ECH1','ECI2','MAPK6']
		for HUN_bait in HUN_baits:
			HUN_preys = HUN_preys.union(set(network_classes_igraph_pub.read_seed_genes(config.output_path + HUN_bait + '_seed_file.txt',0)))
		ube3a_network = g.get_network_with_neighbors_of_seed_genes(seeds,1)
		# delete all edges that are not between two CAMK2D preys or that do not involve 1 CAMK2D prey
		# and 1 HUN complex prey
		edges_to_remove = []
		for edge_tup in ube3a_network.get_edgelist():
			gene_a = ube3a_network.vs[edge_tup[0]]['name']
			gene_b = ube3a_network.vs[edge_tup[1]]['name']
			if not (gene_a in seeds and gene_b in seeds) and \
			   not (gene_a in seeds and gene_b in HUN_preys) and \
			   not (gene_b in seeds and gene_a in HUN_preys):
				edges_to_remove.append(edge_tup)
		ube3a_network.delete_edges(edges_to_remove)
		orphan_nodes = ube3a_network.vs.select(_degree = 0)
		orphan_nodes_indices = [v.index for v in orphan_nodes]
		ube3a_network.delete_vertices(orphan_nodes_indices)
	else:
		print 'unknown input for mode (type of network)'
		sys.exit()

	return ube3a_network


if __name__ == '__main__':

	# define the source networks
	bioplex_tup = ('BioPlex',config.output_path + 'BioPlex.txt')
	ccsb_tup = ('HI-union',config.output_path + 'HI-union.txt')
	litbm_tup = ('Lit-BM-13',config.output_path + 'Lit-BM-13.txt')
	qubic_tup = ('QUBIC',config.output_path + 'QUBIC.txt')
	cofrac_tup = ('CoFrac',config.output_path + 'CoFrac.txt')

	# get a dict that stores for evey GO term all genes annotated with it as a set
	http_client = pyjsonrpc.HttpClient(url = "http://llama.mshri.on.ca/cgi/funcassociate/serv")
	go_dict = run_Funcassociate_pub.get_go_annot_dict(http_client)
	entrez_file = config.input_path + 'entrez_gene.gene_info_human_20161113.tsv'
	id_symbol_dict = run_Funcassociate_pub.get_geneID_symbol_dict(entrez_file)

	# read in the input arguments
	num_rand = int(sys.argv[1])
	network_name = sys.argv[2]
	mode = sys.argv[3]
	if network_name == 'QBCHL':
		source_networks = [bioplex_tup,ccsb_tup,litbm_tup,qubic_tup,cofrac_tup]
	else:
		print 'unknown set of source networks'
		sys.exit()

	# build the source network
	g = network_classes_igraph_pub.Graph(n=0,vertex_attrs={'name':[]})
	g.get_source_network(source_networks)

	# get the seed genes and build the network around the seed genes
	UBE3A_seed_file = config.output_path + mode + '_seed_file.txt'
	UBE3A_seeds = network_classes_igraph_pub.read_seed_genes(UBE3A_seed_file,0)
	UBE3A_seeds = g.get_seed_genes_in_network(UBE3A_seeds)
	ube3a_network = get_subnetwork(mode,UBE3A_seeds,g)
	clusters = ube3a_network.clusters(mode='WEAK')

	# build the input dictionary for running the functional enrichment test with the
	# Funcassociate web service
	query_nodes = ube3a_network.vs['name']
	genespace = [geneID for geneID in g.vs['name']]
	genespace_weights = {}
	for geneID in g.vs['name']:
		genespace_weights[geneID] = 1
	input_dict = {
		'query':query_nodes,
		'species':'Homo sapiens',
		'namespace':'entrezgene',
		'genespace':genespace,
		'genespace_weights':genespace_weights,
		'support':['EXP', 'IC', 'IDA', 'IEA', 'IEP', 'IGC', 'IGI', 'IMP', 'IPI', 'ISA', 'ISM', 'ISO', 'ISS', 'NAS', 'RCA', 'TAS'],
		'mode':'unordered',
		'which':'over',
		'reps':1000,
		'cutoff':0.05
	}
	# run GO term enrichment test for real network
	request_json = {"method": "functionate", "params": [input_dict], "id": "hello", "jsonrpc": "2.0"}
	response =  http_client.call(request_json)
	enriched_GO_terms = {}
	for enrichment in response['over']:
		GOID = str(enrichment[6])
		GO_obj = GOterm(GOID)
		GO_obj.descr = str(enrichment[7])
		GO_obj.num_genes = int(enrichment[2])
		GO_obj.num_real_seed = int(enrichment[0])
		GO_obj.real_LOD = float(enrichment[3])
		GO_obj.real_pvalue = float(enrichment[4])
		GO_obj.real_adj_pvalue = float(enrichment[5])
		enriched_GO_terms[GOID] = GO_obj

	num_genes_enriched_GO = [GO_obj.num_genes for GO_obj in enriched_GO_terms.values()]

	# build a matrix that stores for every degree a list of geneIDs that have that degree from
	# the source network
	max_degree = max(g.degree())
	source_net_degree_matrix = [[] for i in range(0,max_degree+1)]
	for node in g.vs:
		degree = node.degree()
		source_net_degree_matrix[degree].append(node['name'])

	# build a dict that stores the number of seed genes with a given degree
	UBE3A_degree_dict = {}
	for seed in UBE3A_seeds:
		nodes = g.vs.select(name=seed)
		degree = nodes[0].degree()
		if degree not in UBE3A_degree_dict:
			UBE3A_degree_dict[degree] = 0
		UBE3A_degree_dict[degree] += 1

	# generate a lookup list that contains tuples of set of genes and a count for indicating the
	# number of times genes need to be sampled from to build a random selection of seed genes
	# in a degree-controlled way
	sample_list = []
	for degree, count in UBE3A_degree_dict.iteritems():
		delta = math.floor(math.log(degree,2)) ** 2
		sample_genes = []
		i = 0
		while len(sample_genes) < num_rand and i <= delta:
			if i == 0:
				sample_genes = sample_genes + source_net_degree_matrix[degree]
			else:
				sample_genes = sample_genes + source_net_degree_matrix[degree+i]
				sample_genes = sample_genes + source_net_degree_matrix[degree-i]
			i += 1
		sample_list.append((count,sample_genes))

	# initiate output file that will contain information on the randomly built networks and their
	# GO term enrichments and in its first data line info on the real network
	target_graphs = open(config.output_path + 'summary_test_signif_' + mode + '_' + network_name + '_rand_graph_stats.txt','w')
	target_graphs.write('index\tnum_seeds_start\tnum_seeds_network\tnum_seeds_overlap_' + mode + '\t' + \
				 'num_nodes\tnum_nodes_LCC\t' + \
				 'num_enriched_GO\tnum_enriched_GO_overlap_' + mode + '\n')

	# initiate output file that stores the random seed genes and in its first data line the set of real seed genes
	target_rand_seeds = open(config.output_path + 'summary_test_signif_' + mode + '_' + network_name + '_rand_seeds.txt','w')
	target_rand_seeds.write('index\tgeneIDs\n')

	# write out specifics about the real network
	seeds_in_real_network = ube3a_network.vs.select(source='seed')
	real_lcc_size = len(clusters.giant().vs)
	target_graphs.write('real\t' + str(len(UBE3A_seeds)) + '\t' + str(len(seeds_in_real_network)) + '\t' + \
				 str(len(seeds_in_real_network)) + '\t' + str(len(ube3a_network.vs)) + '\t' + str(real_lcc_size) + '\t' + \
				 str(len(response['over'])) + '\t' + str(len(response['over'])) + '\n')
	target_rand_seeds.write('real\t' + '|'.join(UBE3A_seeds) + '\n')

	# determine for every random selection of seeds the network and its GO term enrichments
	for r in range(0,num_rand):

		# randomly choose the seed genes such that the number of random seeds is equivalent to the
		# number of the real seed genes
		rand_seeds = set()
		for tup in sample_list:
			rand_seeds = rand_seeds.union(set(random.sample(set(tup[1]).difference(rand_seeds),tup[0])))
		# write out the randomly picked seed genes
		target_rand_seeds.write(str(r) + '\t' + '|'.join(rand_seeds) + '\n')

		# build the network around the random seed genes
		rand_network = get_subnetwork(mode,rand_seeds,g)
		rand_clusters = rand_network.clusters(mode='WEAK')

		rand_query_nodes = rand_network.vs['name']
		if len(rand_query_nodes) > 0:
			# get the enriched GO terms for this randomly built network
			input_dict = {
				'query':rand_query_nodes,
				'species':'Homo sapiens',
				'namespace':'entrezgene',
				'genespace':genespace,
				'genespace_weights':genespace_weights,
				'support':['EXP', 'IC', 'IDA', 'IEA', 'IEP', 'IGC', 'IGI', 'IMP', 'IPI', 'ISA', 'ISM', 'ISO', 'ISS', 'NAS', 'RCA', 'TAS'],
				'mode':'unordered',
				'which':'over',
				'reps':1000,
				'cutoff':0.05
			}
			request_json = {"method": "functionate", "params": [input_dict], "id": network_name + str(r), "jsonrpc": "2.0"}
			response =  http_client.call(request_json)
		else:
			response['over'] = []

		# compare these enriched GO terms to the enriched GO terms from the real network
		count_GO_overlap = 0
		num_genes_enriched_GO = []
		for enrichment in response['over']:
			GOID = enrichment[6]
			num_genes_enriched_GO.append(int(enrichment[2]))
			if GOID in enriched_GO_terms:
				enriched_GO_terms[GOID].rand_LODs.append(float(enrichment[3]))
				enriched_GO_terms[GOID].rand_pvalues.append(float(enrichment[4]))
				enriched_GO_terms[GOID].nums_rand_seed.append(int(enrichment[0]))
				count_GO_overlap += 1

		# write out info on the randomly generated network and its enriched GO terms
		seeds_in_rand_network = rand_network.vs.select(source='seed')
		seeds_in_rand_network = set([v['name'] for v in seeds_in_rand_network])
		seeds_overlap = seeds_in_rand_network.intersection(UBE3A_seeds)
		if len(rand_query_nodes) == 0:
			rand_lcc_size = 0
		else:
			rand_lcc_size = len(rand_clusters.giant().vs)
		target_graphs.write(str(r) + '\t' + str(len(rand_seeds)) + '\t' + str(len(seeds_in_rand_network)) + '\t' + \
					 str(len(seeds_overlap)) + '\t' + str(len(rand_network.vs)) + '\t' + str(rand_lcc_size) + '\t' + \
					 str(len(response['over'])) + '\t' + str(count_GO_overlap) + '\n')

	target_graphs.close()

	# write out summary stats about the enriched GO terms from the real network and how often they were seen to
	# be enriched or at least or more enriched in the randomly built networks
	target_GO = open(config.output_path + 'summary_test_signif_' + mode + '_' + network_name + '_GO_stats.txt','w')
	target_GO.write('GO_ID\tGO_name\tnum_UBE3A_nodes\tnum_genes_annot\tLOD\tFuncAssociatePvalue\tFuncAssociatePvalueAdjust\t' + \
					'num_rand_networks_same_GO\tnum_rand_networks_as_enriched\tempirical_p_value\tsignificant\t' + \
					'gene_symbols_UBE3A_nodes\tgeneIDs_UBE3A_nodes\n')

	# order the GO terms from most to least significant empirical pvalue
	sorted_GOs = []
	for GOID, GO_obj in enriched_GO_terms.iteritems():
		empirical_p_value = len(filter(lambda x: x >= GO_obj.real_LOD, GO_obj.rand_LODs))/float(num_rand)
		sorted_GOs.append((empirical_p_value,GOID))
	sorted_GOs.sort()
	for tup in sorted_GOs:
		GO_obj = enriched_GO_terms[tup[1]]
		if tup[0] <= 0.05:
			signif = 'yes'
		else:
			signif = 'no'
		target_GO.write(tup[1] + '\t' + GO_obj.descr + '\t' + str(GO_obj.num_real_seed) + '\t' + str(GO_obj.num_genes) + \
					 '\t' + str(GO_obj.real_LOD) + '\t' + str(GO_obj.real_pvalue) + '\t' + str(GO_obj.real_adj_pvalue) \
					 + '\t' + str(len(GO_obj.rand_LODs)) + '\t' + \
					 str(len(filter(lambda x: x >= GO_obj.real_LOD, GO_obj.rand_LODs))) + '\t' + str(tup[0]) + '\t' + \
					 signif + '\t')
		# get all gene IDs of nodes annotated with this GO term
		seedIDs = list(go_dict[tup[1]].intersection(set(query_nodes)))
		seed_symbols = []
		for seedID in seedIDs:
			if seedID in id_symbol_dict:
				seed_symbols.append(id_symbol_dict[seedID])
			else:
				seed_symbols.append('None')
		target_GO.write(','.join(seed_symbols) + '\t' + ','.join(seedIDs) + '\n')


	target_GO.close()
	target_rand_seeds.close()
