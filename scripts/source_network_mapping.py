# script for mapping Ensembl and uniprot IDs to entrez gene IDs for the source networks
# HI-union, BioPlex, QUBIC, CoFrac, Lit-BM-13

import config
import itertools
import sys

# function that returns a dict to map from uniprot or ensembl ID to a list of entrez gene IDs
def get_map_dict(infile,ID_type):

	file1 = open(infile,'r')
	entries = file1.readlines()
	file1.close()

	map_dict = {}

	if ID_type == 'ensembl':
		for line in entries[1:]:
			tab_list = str.split(line[:-1],'\t')
			ensgID = tab_list[0]
			geneID = tab_list[1]
			if geneID != '':
				if ensgID not in map_dict:
					map_dict[ensgID] = []
				map_dict[ensgID].append(int(geneID))
	elif ID_type == 'uniprot':
		for line in entries:
			tab_list = str.split(line[:-1],'\t')
			uniID = tab_list[0]
			geneIDs = tab_list[2].split('; ')
			if len(geneIDs) > 0 and geneIDs[0] != '':
				geneIDs = map(int,geneIDs)
				if uniID not in map_dict:
					map_dict[uniID] = []
				map_dict[uniID] = map_dict[uniID] + geneIDs

	return map_dict

# function that receives a list with the original edges, maps them to the entrez Gene IDs
# with the given dict and sorts the interactors of each edge and returns a set of edges
def map_network(ori_edges,map_dict):

	new_edges = set()
	for edge in ori_edges:
		if edge[0] in map_dict and edge[1] in map_dict:
			geneIDsA = map_dict[edge[0]]
			geneIDsB = map_dict[edge[1]]
			pairs = list(itertools.product(geneIDsA,geneIDsB))
			sorted_pairs = [tuple(sorted(t)) for t in pairs]
			new_edges = new_edges.union(set(sorted_pairs))
	print len(ori_edges), len(new_edges)
	return new_edges

# function that gets a list of edges and returns the list of edges wtihout homodimers
def remove_homodimers(edges):

	new_edges = set()
	for edge in edges:
		if edge[0] != edge[1]:
			new_edges.add(edge)

	return new_edges

# function that writes out the given list of edges to the given outfile
def write_network(outfile,edges):

	target = open(outfile,'w')
	target.write('geneID_A\tgeneID_B\n')
	for edge in edges:
		target.write(str(edge[0]) + '\t' + str(edge[1]) + '\n')
	target.close()

if __name__ == '__main__':

	networks = sys.argv[1].split(',')

	uniprot_map_dict = get_map_dict(config.input_path + 'HUMAN_9606_idmapping_selected_20170303.tab','uniprot')
	ensembl_map_dict = get_map_dict(config.input_path + 'mart_export_20170303.txt','ensembl')

	for network in networks:
		print network

		if network == 'HI-union':
			# get the mapped HI-union dataset
			files = ['CCSB_test_space_validated_PSI_MI.2.7_20170306.psi','CCSB_test_space_verified_PSI_MI.2.7_20170306.psi',\
					 'CCSB_unpublished_PSI_MI_2.7_20170306.psi','all_published.psi']
			ori_edges = set()
			for file_name in files:
				file1 = open(config.input_path + file_name,'Ur')
				entries = file1.readlines()
				file1.close()
				for line in entries[1:]:
					tab_list = str.split(line[:-1],'\t')
					interactorA = tab_list[0].split(':')[1]
					interactorB = tab_list[1].split(':')[1]
					ori_edges.add((interactorA,interactorB))
			new_edges = map_network(ori_edges,uniprot_map_dict)
			new_edges = remove_homodimers(new_edges)
			print len(new_edges)

		elif network == 'BioPlex':
			# get the mapped BioPlex dataset
			files = ['BioPlex_interactionList_v4.tsv']
			new_edges = set()
			for file_name in files:
				file1 = open(config.input_path + file_name,'Ur')
				entries = file1.readlines()
				file1.close()
				for line in entries[1:]:
					tab_list = str.split(line[:-1],'\t')
					new_edges.add(tuple(sorted(map(int,tab_list[:2]))))
			print len(entries)-1, len(new_edges)
			new_edges = remove_homodimers(new_edges)
			print len(new_edges)

		elif network == 'QUBIC':
			# get the mapped QUBIC dataset
			file1 = open(config.input_path + 'QUBIC_TableS2.txt','Ur')
			entries = file1.readlines()
			file1.close()
			ori_edges = set()
			for line in entries[1:]:
				tab_list = str.split(line[:-1],'\t')
				baitIDs = tab_list[0].split(';')
				proc_baitIDs = []
				for baitID in baitIDs:
					pos = baitID.find('-')
					if pos > -1:
						proc_baitIDs.append(baitID[:pos])
					else:
						proc_baitIDs.append(baitID)
				preyIDs = tab_list[1].split(';')
				proc_preyIDs = []
				for preyID in preyIDs:
					pos = preyID.find('-')
					if pos > -1:
						proc_preyIDs.append(preyID[:pos])
					else:
						proc_preyIDs.append(preyID)
				ori_edges = ori_edges.union(set(itertools.product(baitIDs,preyIDs)))
			print len(entries)-1
			new_edges = map_network(ori_edges,uniprot_map_dict)
			new_edges = remove_homodimers(new_edges)
			print len(new_edges)

		elif network == 'CoFrac':
			# get the mapped CoFrac dataset
			file1 = open(config.input_path + 'Wan_et_al_TableS2.txt','Ur')
			entries = file1.readlines()
			file1.close()
			ori_edges = set()
			for line in entries[1:]:
				tab_list = str.split(line[:-1],'\t')
				interactorA = tab_list[0].strip()
				interactorB = tab_list[1].strip()
				ori_edges.add((interactorA,interactorB))
			print len(entries)-1
			new_edges = map_network(ori_edges,ensembl_map_dict)
			new_edges = remove_homodimers(new_edges)
			print len(new_edges)

		elif network == 'Lit-BM-13':
			# get the mapped Lit-BM-13 dataset
			file1 = open(config.input_path + 'Lit-BM-13_TableS1.txt','Ur')
			entries = file1.readlines()
			file1.close()
			new_edges = set()
			for line in entries[1:]:
				tab_list = str.split(line[:-1],'\t')
				new_edges.add(tuple(sorted(map(int,tab_list[:2]))))
			print len(entries)-1, len(new_edges)
			new_edges = remove_homodimers(new_edges)
			print len(new_edges)

		elif network == 'DuB_AP':
			# get the mapped deubiquitination PPI dataset
			files = ['IntAct_PMID_19615732_170424.txt','IntAct_PMID_20562859_170424.txt']
			ori_edges = set()
			for file_name in files:
				file1 = open(config.input_path + file_name,'Ur')
				entries = file1.readlines()
				file1.close()
				for line in entries[1:]:
					tab_list = str.split(line[:-1],'\t')
					interactorA = tab_list[0].split(':')[1]
					interactorB = tab_list[1].split(':')[1]
					ori_edges.add((interactorA,interactorB))
			new_edges = map_network(ori_edges,uniprot_map_dict)
			new_edges = remove_homodimers(new_edges)
			print len(new_edges)

		elif network == 'Howley_map':
			# get a network of all protein pairs from the pulldowns worked with in this project
			# from the Howley lab
			baits = [("CAMK2D",817),("MAPK6",5597),("ECH1",1891),("ECI2",10455),("HIF1AN",55662),\
			 		 ("NEURL4",84461),("HERC2",8924),("UBE3A",7337)]
			new_edges = set()
			for tup in baits:
				bait = tup[0]
				bait_geneID = tup[1]
				file1 = open(config.output_path + bait + '_seed_file.txt','r')
				entries = file1.readlines()
				file1.close()
				for line in entries[1:]:
					tab_list = str.split(line[:-1],'\t')
					geneID = int(tab_list[0])
					if geneID < bait_geneID:
						pair = (str(geneID),str(bait_geneID))
					else:
						pair = (str(bait_geneID),str(geneID))
					new_edges.add(pair)
			new_edges = remove_homodimers(new_edges)
			print len(new_edges)

		else:
			print 'unknown network name'
			sys.exit()

		outfile = config.output_path + network + '.txt'
		write_network(outfile,new_edges)
