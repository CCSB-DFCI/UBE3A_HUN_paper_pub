# script that computes significances of LCC of networks formed by preys from every baits,
# LCC of HUN network, and closeness of CAMK2D preys to HUN complex

import sys
import numpy
import network_classes_igraph_pub
import random_graphs
import HUN_network_analysis_pub
import config


# function that test for every given bait and random networks the signficance of the LCC
# formed by the preys of the bait
def test_signif_LCC_of_seed_genes(num_rand,seed_files,real_network,outfile_stats):

	# build dictionary that saves every set of seed genes that are in the real network
	seeds = {}
	for seed_file in seed_files:
		seed_genes = network_classes_igraph_pub.read_seed_genes(config.output_path + seed_file,0)
		seed_genes_in_g = real_network.get_seed_genes_in_network(seed_genes)
		seeds[seed_file] = seed_genes_in_g

	# for every random network, compute for every set of seed genes the lcc
	rand_lccs = numpy.zeros((len(seeds),num_rand))
	for r in range(0,num_rand):
		# for every random network
		rand_file = config.rand_path + config.rand_prefix + str(r) + '.txt'
		rg = random_graphs.get_real_network(rand_file)
		for b,seed_file in enumerate(seed_files):
			# determine the size of the LCC and save
			sg = rg.subgraph(seeds[seed_file])
			clusters = sg.clusters(mode='WEAK')
			lcc_size = len(clusters.giant().vs)
			rand_lccs[b,r] = lcc_size

	# save the statistics for every set of seed genes
	target = open(outfile_stats,'w')
	target.write('seed_file\tp_value\tmean\tstd\tz_score\tnum_nodes_real_network\tsize_real_LCC\tnum_equal_or_larger_rand_LCC\tnum_rand\n')

	# for every bait, compute the statistics
	for b,seed_file in enumerate(seed_files):
		# determine the real size of the LCC of its preys
		sg = real_network.subgraph(seeds[seed_file])
		clusters = sg.clusters(mode='WEAK')
		lcc_size = len(clusters.giant().vs)
		# determine the p-value, mean, std, z-score based on distribution of random LCCs
		b_rand_lccs = rand_lccs[b,:]
		num_larger_lcc = len(filter(lambda s: s >= lcc_size,b_rand_lccs))
		mean = b_rand_lccs.mean()
		std = b_rand_lccs.std()
		z_score = (lcc_size - mean)/std
		p_value = num_larger_lcc/float(num_rand)
		target.write(seed_file + '\t' + str(p_value) + '\t' + str(mean) + '\t' + str(std) + '\t' + \
					 str(z_score) + '\t' + str(len(sg.vs)) + '\t' + str(lcc_size) + '\t' + str(num_larger_lcc) + \
					 '\t' + str(num_rand) + '\n')

		# save the rand lccs
		target_lccs = open(config.output_path + seed_file[:-4] + '_rand_lccs.txt','w')
		target_lccs.write('index\tnum_nodes_lcc\n')
		target_lccs.write('real\t' + str(lcc_size) + '\n')
		for r,rand_lcc in enumerate(b_rand_lccs):
			target_lccs.write(str(r) + '\t' + str(rand_lcc) + '\n')
		target_lccs.close()

	target.close()



# function that tests the significance of the size of the LCC of the HUN network compared to random networks
def test_signif_LCC_HUN_network(node_file,real_network,num_rand,outfile_stats):

	seed_genes = network_classes_igraph_pub.read_seed_genes(node_file,0)
	seed_genes_in_g = real_network.get_seed_genes_in_network(seed_genes)

	sg = real_network.subgraph(seed_genes_in_g)
	network_classes_igraph_pub.add_node_attributes(sg,node_file,[2])
	sg = HUN_network_analysis_pub.delete_edges_between_preys_not_same_attribute(sg,'baits','|')
	clusters = sg.clusters(mode='WEAK')
	lcc_size = len(clusters.giant().vs)

	lccs = []
	for r in range(0,num_rand):
		# for every random network
		rand_file = config.rand_path + config.rand_prefix + str(r) + '.txt'
		rg = random_graphs.get_real_network(rand_file)
		srg = rg.subgraph(seed_genes_in_g)
		network_classes_igraph_pub.add_node_attributes(srg,node_file,[2])
		srg = HUN_network_analysis_pub.delete_edges_between_preys_not_same_attribute(srg,'baits','|')
		clusters = srg.clusters(mode='WEAK')
		rand_lcc_size = len(clusters.giant().vs)
		lccs.append(rand_lcc_size)

	target = open(node_file[:-4]  + '_rand_lccs.txt','w')
	target.write('index\tnum_nodes_lcc\n')
	target.write('real\t' + str(lcc_size) + '\n')
	for i,lcc in enumerate(lccs):
		target.write(str(i) + '\t' + str(lcc) + '\n')
	target.close()

	num_larger_lcc = len(filter(lambda s: s >= lcc_size,lccs))
	mean = numpy.mean(lccs)
	std = numpy.std(lccs)
	z_score = (lcc_size - mean)/std
	p_value = num_larger_lcc/float(num_rand)
	target = open(outfile_stats,'w')
	target.write('p_value\tmean\tstd\tz_score\tnum_nodes_real_network\tsize_real_LCC\tnum_larger_rand_LCC\tnum_rand\n')
	target.write(str(p_value) + '\t' + str(mean) + '\t' + str(std) + '\t' + str(z_score) + '\t' + str(len(sg.vs)) + \
				 '\t' + str(lcc_size) + '\t' + str(num_larger_lcc) + '\t' + str(num_rand) + '\n')
	target.close()


# test for closeness of CAMK2D preys to HUN complex
def test_signif_CAMK2D_HUNcomplex_link(bait_file,complex_file,core_complex_members,real_network,num_rand,out_prefix):

	bait_preys = network_classes_igraph_pub.read_seed_genes(bait_file,0)
	bait_preys_in_g = real_network.get_seed_genes_in_network(bait_preys)

	complex_preys = set(network_classes_igraph_pub.read_seed_genes(complex_file,0))
	complex_preys_in_g = real_network.get_seed_genes_in_network(complex_preys)

	neighbor_lists = real_network.neighborhood(bait_preys_in_g,1)
	HUN_neighbors = set()
	count_preys_connected_to_core_complex = 0
	count_preys_connected_to_complex_preys = 0
	for i,neighbor_list in enumerate(neighbor_lists):
		neighbors = set()
		for neighbor in neighbor_list:
			neighbors.add(real_network.vs[neighbor]['name'])
		if len(neighbors.intersection(core_complex_members)) > 0:
			count_preys_connected_to_core_complex += 1
		if len(neighbors.intersection(complex_preys_in_g)) > 0:
			count_preys_connected_to_complex_preys += 1
		HUN_neighbors = HUN_neighbors.union(neighbors.intersection(complex_preys_in_g))

	counts_preys_connected_to_core_complex = []
	counts_preys_connected_to_complex_preys = []
	counts_complex_interactors = []
	for r in range(0,num_rand):
		# for every random network
		rand_file = config.rand_path + config.rand_prefix + str(r) + '.txt'
		rg = random_graphs.get_real_network(rand_file)
		neighbor_lists = rg.neighborhood(bait_preys_in_g,1)
		rand_HUN_neighbors = set()
		rand_count_preys_connected_to_core_complex = 0
		rand_count_preys_connected_to_complex_preys = 0
		for i,neighbor_list in enumerate(neighbor_lists):
			neighbors = set()
			for neighbor in neighbor_list:
				neighbors.add(rg.vs[neighbor]['name'])
			if len(neighbors.intersection(core_complex_members)) > 0:
				rand_count_preys_connected_to_core_complex += 1
			if len(neighbors.intersection(complex_preys_in_g)) > 0:
				rand_count_preys_connected_to_complex_preys += 1
			rand_HUN_neighbors = rand_HUN_neighbors.union(neighbors.intersection(complex_preys_in_g))
		counts_preys_connected_to_core_complex.append(rand_count_preys_connected_to_core_complex)
		counts_preys_connected_to_complex_preys.append(rand_count_preys_connected_to_complex_preys)
		counts_complex_interactors.append(len(rand_HUN_neighbors))

	target = open(config.output_path + out_prefix + 'counts_preys_connected_to_core_complex_members_rand_distr.txt','w')
	target.write('index\tcount\n')
	target.write('real\t' + str(count_preys_connected_to_core_complex) + '\n')
	for i,v in enumerate(counts_preys_connected_to_core_complex):
		target.write(str(i) + '\t' + str(v) + '\n')
	target.close()
	target = open(config.output_path + out_prefix + 'counts_preys_connected_to_complex_preys_rand_distr.txt','w')
	target.write('index\tcount\n')
	target.write('real\t' + str(count_preys_connected_to_complex_preys) + '\n')
	for i,v in enumerate(counts_preys_connected_to_complex_preys):
		target.write(str(i) + '\t' + str(v) + '\n')
	target.close()
	target = open(config.output_path + out_prefix + 'counts_complex_preys_connected_to_preys_rand_distr.txt','w')
	target.write('index\tcount\n')
	target.write('real\t' + str(len(HUN_neighbors)) + '\n')
	for i,v in enumerate(counts_complex_interactors):
		target.write(str(i) + '\t' + str(v) + '\n')
	target.close()

	target = open(config.output_path + out_prefix + 'summary_signif_closeness.txt','w')
	target.write('test\treal_obs\tp_value\tmean\tstd\tz_score\tnum_rand\n')

	num_larger = len(filter(lambda s: s >= count_preys_connected_to_core_complex,counts_preys_connected_to_core_complex))
	mean = numpy.mean(counts_preys_connected_to_core_complex)
	std = numpy.std(counts_preys_connected_to_core_complex)
	z_score = (count_preys_connected_to_core_complex - mean)/std
	p_value = num_larger/float(num_rand)
	target.write('num_preys_connected_core_complex\t' + str(count_preys_connected_to_core_complex) + '\t' + \
				  str(p_value) + '\t' + str(mean) + '\t' + str(std) + '\t' + str(z_score) + '\t' + str(num_rand) + '\n')

	num_larger = len(filter(lambda s: s >= count_preys_connected_to_complex_preys,counts_preys_connected_to_complex_preys))
	mean = numpy.mean(counts_preys_connected_to_complex_preys)
	std = numpy.std(counts_preys_connected_to_complex_preys)
	z_score = (count_preys_connected_to_complex_preys - mean)/std
	p_value = num_larger/float(num_rand)
	target.write('num_preys_connected_complex_preys\t' + str(count_preys_connected_to_complex_preys) + '\t' + \
				  str(p_value) + '\t' + str(mean) + '\t' + str(std) + '\t' + str(z_score) + '\t' + str(num_rand) + '\n')

	num_larger = len(filter(lambda s: s >= len(HUN_neighbors),counts_complex_interactors))
	mean = numpy.mean(counts_complex_interactors)
	std = numpy.std(counts_complex_interactors)
	z_score = (len(HUN_neighbors) - mean)/std
	p_value = num_larger/float(num_rand)
	target.write('num_complex_interactors\t' + str(len(HUN_neighbors)) + '\t' + str(p_value) + '\t' + str(mean) + \
				 '\t' + str(std) + '\t' + str(z_score) + '\t' + str(num_rand) + '\n')

	target.close()



if __name__ == '__main__':

	bioplex_tup = ('BioPlex',config.output_path + 'BioPlex.txt')
	ccsb_tup = ('HI-union',config.output_path + 'HI-union.txt')
	litbm_tup = ('Lit-BM-13',config.output_path + 'Lit-BM-13.txt')
	qubic_tup = ('QUBIC',config.output_path + 'QUBIC.txt')
	cofrac_tup = ('CoFrac',config.output_path + 'CoFrac.txt')
	source_networks = [bioplex_tup,ccsb_tup,litbm_tup,qubic_tup,cofrac_tup]

	analysis = sys.argv[1]
	num_rand = int(sys.argv[2])

	real_network = network_classes_igraph_pub.Graph(n=0,vertex_attrs={'name':[]})
	real_network.get_source_network(source_networks)

	if analysis == 'seed_lcc':

		outfile_stats = config.output_path + 'summary_signif_LCC_seed_files.txt'
		seed_files = ['CAMK2D_seed_file.txt','ECH1_seed_file.txt','ECI2_seed_file.txt','HERC2_seed_file.txt',\
					  'HIF1AN_seed_file.txt','MAPK6_seed_file.txt','NEURL4_seed_file.txt','UBE3A_seed_file.txt',\
					  'UBE3A_seed_file_with_proteasome.txt']
		test_signif_LCC_of_seed_genes(num_rand,seed_files,real_network,outfile_stats)

	elif analysis == 'hun_lcc':

		node_file = config.output_path + 'HUN_network_QBCHL.node_attributes.txt'
		outfile_stats = config.output_path + 'summary_signif_LCC_HUN_network.txt'
		test_signif_LCC_HUN_network(node_file,real_network,num_rand,outfile_stats)

	elif analysis == 'camk2d_hun':

		bait_file = config.output_path + 'CAMK2D_seed_file.txt'
		complex_file = config.output_path + 'HUN_UBE3A_seed_file.txt'
		core_complex_members = set(['5597','1891','10455','84461','8924','7337'])
		out_prefix = 'CAMK2D_HUN_'
		test_signif_CAMK2D_HUNcomplex_link(bait_file,complex_file,core_complex_members,real_network,num_rand,out_prefix)

	else:

		print 'available analyses: seed_lcc, hun_lcc, camk2d_hun'
		sys.exit()
