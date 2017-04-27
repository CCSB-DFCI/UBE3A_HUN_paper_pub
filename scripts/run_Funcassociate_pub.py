# script that takes any input file and computes functional enrichments on the
# list of genes using the Funcassociate web service
# writes out the enriched GO terms with stats and a list of input genes for each
# GO term with gene IDs and symbols that are annotated with each GO term

import sys
import pyjsonrpc
import config

# get a dictionary to map between gene IDs and gene symbols
def get_geneID_symbol_dict(entrez_file):

	id_symbol_dict = {}
	file1 = open(entrez_file,'r')
	entries = file1.readlines()
	file1.close()
	for row in entries[1:]:
		tab_list = str.split(row[:-1],'\t')
		geneID = str(tab_list[1])
		symbol = tab_list[2]
		id_symbol_dict[geneID] = symbol

	return id_symbol_dict

# get a dict that stores for evey GO term all genes annotated with it as a set
def get_go_annot_dict(http_client):

	input_list = ['Homo sapiens','entrezgene',['EXP', 'IC', 'IDA', 'IEA', 'IEP', 'IGC', 'IGI', 'IMP', 'IPI', 'ISA', 'ISM', 'ISO', 'ISS', 'NAS', 'RCA', 'TAS']]
	request_json = {"method": "go_associations", "params": input_list, "id": 0, "jsonrpc": "2.0"}
	go_associations =  http_client.call(request_json)
	go_dict = {}
	for go_assoc in go_associations:
		go_dict[go_assoc[0]] = set(go_assoc[1:])

	return go_dict


if __name__ == '__main__':

	infile = sys.argv[1]
	col = int(sys.argv[2])

	# get the seed genes
	file1 = open(config.output_path + infile,'r')
	entries = file1.readlines()
	file1.close()
	genes = set()
	for line in entries[1:]:
		tab_list = str.split(line[:-1],'\t')
		genes.add(tab_list[col])

	http_client = pyjsonrpc.HttpClient(url = "http://llama.mshri.on.ca/cgi/funcassociate/serv")
	go_dict = get_go_annot_dict(http_client)

	entrez_file = config.input_path + 'entrez_gene.gene_info_human_20161113.tsv'
	id_symbol_dict = get_geneID_symbol_dict(entrez_file)

	# get the enriched GO terms
	input_dict = {
		'query':list(genes),
		'species':'Homo sapiens',
		'namespace':'entrezgene',
		'support':['EXP', 'IC', 'IDA', 'IEA', 'IEP', 'IGC', 'IGI', 'IMP', 'IPI', 'ISA', 'ISM', 'ISO', 'ISS', 'NAS', 'RCA', 'TAS'],
		'mode':'unordered',
		'which':'over',
		'reps':1000,
		'cutoff':0.05
	}
	request_json = {"method": "functionate", "params": [input_dict], "id": "hello", "jsonrpc": "2.0"}
	response =  http_client.call(request_json)

	outfile = config.output_path + infile[:-4] + '_Funcassociate_results.txt'
	target = open(outfile,'w')
	target.write('GO_ID\tGO_name\tnum_seed_genes\tnum_total_seed_genes\tnum_genes_genespace\tLOD\tp_value\tp_adjust\tseed_gene_symbols\tseed_gene_IDs\n')

	# order the enriched GO terms by their LOD from highest to lowest
	LODs = []
	for i,GO_result in enumerate(response['over']):
		LOD = GO_result[3]
		LODs.append((LOD,i))
	LODs.sort()
	LODs = LODs[::-1]

	for tup in LODs:
		GO_result = response['over'][tup[1]]
		target.write(GO_result[6] + '\t' + GO_result[7] + '\t' + str(GO_result[0]) + '\t' + str(GO_result[1]) + \
					'\t' + str(GO_result[2]) + '\t' + str(GO_result[3]) + '\t' + str(GO_result[4]) + '\t' + \
					str(GO_result[5]) + '\t')
		# get all gene IDs of seed genes annotated with this GO term
		seedIDs = list(go_dict[GO_result[6]].intersection(set(genes)))
		seed_symbols = []
		for seedID in seedIDs:
			if seedID in id_symbol_dict:
				seed_symbols.append(id_symbol_dict[seedID])
			else:
				seed_symbols.append('None')
		target.write(','.join(seed_symbols) + '\t' + ','.join(seedIDs) + '\n')
		print GO_result[0], len(seedIDs)

	target.close()
