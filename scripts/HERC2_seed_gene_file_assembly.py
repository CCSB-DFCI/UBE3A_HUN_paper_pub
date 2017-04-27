# script to assemble a seed gene file for HERC2

import config
import UBE3A_seed_gene_file_assembly

if __name__ == '__main__':

	entrez_file = config.input_path + 'entrez_gene.gene_info_human_20161113.tsv'
	gene_descr_dict = UBE3A_seed_gene_file_assembly.get_geneID_descr_dict(entrez_file)

	herc2_preys = {}

	# read in the HEK293T pulldown file
	file1 = open(config.input_path + "HERC2_HEK293T_pulldown_summary_1sheet.txt",'r')
	entries = file1.readlines()
	file1.close()
	for line in entries[1:]:
		tab_list = str.split(line[:-1],'\t')
		geneID = int(tab_list[2])
		symbol = tab_list[1]
		if geneID not in herc2_preys:
			herc2_preys[geneID] = UBE3A_seed_gene_file_assembly.Prey(geneID)
			herc2_preys[geneID].num_hek += 1
			herc2_preys[geneID].num_cell_lines += 1
			herc2_preys[geneID].symbol = symbol

	# read the SY5Y pulldown file
	file1 = open(config.input_path + "HERC2_SH-SY5Y_cutoff95-90.txt",'r')
	entries = file1.readlines()
	file1.close()
	for line in entries[1:]:
		tab_list = str.split(line[:-1],'\t')
		geneID = int(tab_list[0])
		symbol = tab_list[1]
		hcip_sy5y = int(tab_list[11][:-1])
		if geneID not in herc2_preys:
			herc2_preys[geneID] = UBE3A_seed_gene_file_assembly.Prey(geneID)
			herc2_preys[geneID].symbol = symbol
		if herc2_preys[geneID].num_sy5y == 0:
			herc2_preys[geneID].num_cell_lines += 1
			herc2_preys[geneID].num_sy5y += 1
		if herc2_preys[geneID].hcip_sy5y == 'NA':
			herc2_preys[geneID].hcip_sy5y = hcip_sy5y
		elif herc2_preys[geneID].hcip_sy5y < hcip_sy5y:
			herc2_preys[geneID].hcip_sy5y = hcip_sy5y

	# get the description for every gene
	for geneID, prey_obj in herc2_preys.iteritems():
		if geneID in gene_descr_dict:
			prey_obj.descr = gene_descr_dict[geneID]

	# remove proteasome-related proteins
	geneIDs_to_remove = set()
	for geneID, prey_obj in herc2_preys.iteritems():
		if "proteas" in prey_obj.descr:
			print geneID, prey_obj.descr
			geneIDs_to_remove.add(geneID)
	for geneID in geneIDs_to_remove:
		del herc2_preys[geneID]

	# determine by which cell lines a prey has been found
	for geneID,prey_obj in herc2_preys.iteritems():
		if prey_obj.num_cell_lines == 2:
			prey_obj.cell_lines = "b"
		else:
			if prey_obj.num_hek > 0:
				prey_obj.cell_lines = "h"
			else:
				prey_obj.cell_lines = "s"

	# write the seed file
	target = open(config.output_path + "HERC2_seed_file.txt","w")
	target.write("geneID\tsymbol\tdescr\tnum_hek\tnum_sy5y\tnum_cell_lines\tcell_lines\thcip_cutoff_sy5y\tbait\n")
	for geneID,prey_obj in herc2_preys.iteritems():
		target.write(str(geneID) + "\t" + prey_obj.symbol + "\t" + prey_obj.descr + "\t" + str(prey_obj.num_hek) + \
					 "\t" + str(prey_obj.num_sy5y) + "\t" + str(prey_obj.num_cell_lines) + "\t" + \
					 prey_obj.cell_lines + "\t" + str(prey_obj.hcip_sy5y) + "HERC2\n")
	target.close()
