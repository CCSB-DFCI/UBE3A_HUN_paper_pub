# script to assemble seed genes (preys) of UBE3A

import config

class Prey(object):
	def __init__(self, geneID):
		super(Prey,self).__init__()
		self.geneID = geneID
		self.symbol = None
		self.descr = None
		self.num_hek = 0
		self.num_sy5y_90 = 0
		self.num_sy5y_95 = 0
		self.num_sy5y = 0
		self.hcip_sy5y = 'NA'
		self.frac_occ = 0
		self.cell_lines = None
		self.num_cell_lines = 0


def get_symbol_geneID_dict(infile):

	symbol_gene_dict = {}
	file1 = open(infile,'r')
	entries = file1.readlines()
	file1.close()
	for row in entries[1:]:
		tab_list = str.split(row[:-1],'\t')
		symbol = tab_list[2]
		geneID = int(tab_list[1])
		if symbol not in symbol_gene_dict:
			symbol_gene_dict[symbol] = []
		symbol_gene_dict[symbol].append(geneID)

	return symbol_gene_dict


def get_geneID_descr_dict(infile):

	gene_descr_dict = {}
	file1 = open(infile,'r')
	entries = file1.readlines()
	file1.close()
	for row in entries[1:]:
		tab_list = str.split(row[:-1],'\t')
		descr = tab_list[8]
		geneID = int(tab_list[1])
		gene_descr_dict[geneID] = descr

	return gene_descr_dict

if __name__ == '__main__':

	entrez_file = config.input_path + 'entrez_gene.gene_info_human_20161113.tsv'
	symbol_gene_dict = get_symbol_geneID_dict(entrez_file)
	gene_descr_dict = get_geneID_descr_dict(entrez_file)

	ube3a_preys = {}

	# read in the HEK293T pulldown file
	file1 = open(config.input_path + "UBE3A_HEK293T_pulldown_summary.txt",'r')
	entries = file1.readlines()
	file1.close()
	for line in entries[1:]:
		tab_list = str.split(line[:-1],'\t')
		symbol = tab_list[0]
		num_hek = int(round(float(tab_list[19])))
		if symbol not in symbol_gene_dict:
			print symbol, "not in symbol_gene_dict"
		else:
			geneIDs = symbol_gene_dict[symbol]
			if len(geneIDs) > 1:
				print symbol, geneIDs
			else:
				geneID = geneIDs[0]
				if geneID not in ube3a_preys:
					ube3a_preys[geneID] = Prey(geneID)
				ube3a_preys[geneID].symbol = symbol
				ube3a_preys[geneID].num_hek = num_hek

	# read in th SY5Y pulldown files and get the right hcip cutoff for every prey and UBE3A isoform
	# the 90 file contains the 95 cutoff genes as well
	sy5y_preys = {}
	cutoffs = ['90','95']
	for cutoff in cutoffs:
		file1 = open(config.input_path + "UBE3A_SH-SY5Y_cutoff" + cutoff + ".txt",'Ur')
		entries = file1.readlines()
		file1.close()
		for line in entries[1:]:
			tab_list = str.split(line[:-1],'\t')
			geneID = int(tab_list[0])
			symbol = tab_list[1].strip()
			isoform = tab_list[8].strip()
			active_form = tab_list[9].strip()
			gene_tup = (geneID,isoform,active_form)
			if gene_tup not in sy5y_preys:
				sy5y_preys[gene_tup] = Prey(geneID)
				sy5y_preys[gene_tup].symbol = symbol
				sy5y_preys[gene_tup].hcip_sy5y = 90
			else:
				sy5y_preys[gene_tup].hcip_sy5y = 95

	# collapse the SY5Y data by counting how many occurrences of a prey there are at 90 and 95 cutoff
	# and only keep the 95 cutoff as value if at least half of all occurrences were found at this
	# cutoff
	for gene_tup,prey_obj in sy5y_preys.iteritems():
		geneID = gene_tup[0]
		if geneID not in ube3a_preys:
			ube3a_preys[geneID] = Prey(geneID)
			ube3a_preys[geneID].symbol = prey_obj.symbol
		if prey_obj.hcip_sy5y == 90:
			ube3a_preys[geneID].num_sy5y_90 += 1
		else:
			ube3a_preys[geneID].num_sy5y_95 += 1

	for geneID,prey_obj in ube3a_preys.iteritems():
		prey_obj.num_sy5y = prey_obj.num_sy5y_95 + prey_obj.num_sy5y_90
		if prey_obj.num_sy5y_90 > 0 or prey_obj.num_sy5y_95 > 0:
			if prey_obj.num_sy5y_95 >= prey_obj.num_sy5y_90:
				prey_obj.hcip_sy5y = 95
			else:
				prey_obj.hcip_sy5y = 90

	# get the description for every gene
	for geneID, prey_obj in ube3a_preys.iteritems():
		if geneID in gene_descr_dict:
			prey_obj.descr = gene_descr_dict[geneID]

	# delete all preys that occurred in only 1 out of all pulldowns from a given cell line unless
	# it is MCM-like or it was seen in a pulldown in the other cell line
	geneIDs_to_remove = set()
	for geneID, prey_obj in ube3a_preys.iteritems():
		if ("MCM" not in prey_obj.symbol and prey_obj.num_sy5y == 1 and prey_obj.num_hek == 0) or \
		   ("MCM" not in prey_obj.symbol and prey_obj.num_sy5y == 0 and prey_obj.num_hek == 1):
			geneIDs_to_remove.add(geneID)

	for geneID in geneIDs_to_remove:
		del ube3a_preys[geneID]

	# determine by which cell lines a prey has been found
	for geneID,prey_obj in ube3a_preys.iteritems():
		if prey_obj.num_hek > 0 and prey_obj.num_sy5y > 0:
			prey_obj.cell_lines = "b"
			prey_obj.frac_occ = 1.0
		else:
			if prey_obj.num_hek > 0:
				prey_obj.cell_lines = "h"
				prey_obj.frac_occ = prey_obj.num_hek/12.0
			else:
				prey_obj.cell_lines = "s"
				prey_obj.frac_occ = prey_obj.num_sy5y/6.0

	# mark UBE3A as bait
	ube3a_preys[7337].cell_lines = "bait"

	# write the seed file with the proteasome subunits
	target = open(config.output_path + "UBE3A_seed_file_with_proteasome.txt","w")
	target.write("geneID\tsymbol\tdescr\tnum_hek\tnum_sy5y\tfrac_occ\tcell_lines\thcip_cutoff_sy5y\n")
	for geneID,prey_obj in ube3a_preys.iteritems():
		target.write(str(geneID) + "\t" + prey_obj.symbol + "\t" + prey_obj.descr + "\t" + str(prey_obj.num_hek) + \
					 "\t" + str(prey_obj.num_sy5y) + "\t" + str(prey_obj.frac_occ) + "\t" + \
					 prey_obj.cell_lines + "\t" + str(prey_obj.hcip_sy5y) + "\n")
	target.close()

	# remove proteasome-related proteins
	geneIDs_to_remove = set()
	for geneID, prey_obj in ube3a_preys.iteritems():
		if prey_obj.symbol.find('PSM') > -1 and prey_obj.symbol != 'PSMD4':
			print geneID, prey_obj.descr
			geneIDs_to_remove.add(geneID)
	for geneID in geneIDs_to_remove:
		del ube3a_preys[geneID]

	# write the seed file without the proteasome subunits
	target = open(config.output_path + "UBE3A_seed_file.txt","w")
	target.write("geneID\tsymbol\tdescr\tnum_hek\tnum_sy5y\tfrac_occ\tcell_lines\thcip_cutoff_sy5y\n")
	for geneID,prey_obj in ube3a_preys.iteritems():
		target.write(str(geneID) + "\t" + prey_obj.symbol + "\t" + prey_obj.descr + "\t" + str(prey_obj.num_hek) + \
					 "\t" + str(prey_obj.num_sy5y) + "\t" + str(prey_obj.frac_occ) + "\t" + \
					 prey_obj.cell_lines + "\t" + str(prey_obj.hcip_sy5y) + "\n")
	target.close()
