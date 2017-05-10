# script to assemble seed gene file for a couple of baits

import config
import UBE3A_seed_gene_file_assembly

def get_seed_gene_file(gene_descr_dict,bait,geneID_bait):

	preys = {}

	# read the HEK293 pulldown file
	file1 = open(config.input_path + bait + "_T-REx293.txt",'Ur')
	entries = file1.readlines()
	file1.close()
	for line in entries[1:]:
		tab_list = str.split(line[:-1],'\t')
		geneID = int(tab_list[0])
		symbol = tab_list[1]
		if geneID not in preys:
			preys[geneID] = UBE3A_seed_gene_file_assembly.Prey(geneID)
			preys[geneID].symbol = symbol
			preys[geneID].num_cell_lines += 1
			preys[geneID].num_hek += 1

	# read the SY5Y pulldown file
	file1 = open(config.input_path + bait + "_SH-SY5Y_cutoff95-90.txt",'Ur')
	entries = file1.readlines()
	file1.close()
	for line in entries[1:]:
		tab_list = str.split(line[:-1],'\t')
		geneID = int(tab_list[0])
		symbol = tab_list[1]
		hcip_sy5y = int(tab_list[10][:-1])
		if geneID not in preys:
			preys[geneID] = UBE3A_seed_gene_file_assembly.Prey(geneID)
			preys[geneID].symbol = symbol
		if preys[geneID].num_sy5y == 0:
			preys[geneID].num_cell_lines += 1
			preys[geneID].num_sy5y += 1
		if preys[geneID].hcip_sy5y == 'NA':
			preys[geneID].hcip_sy5y = hcip_sy5y
		elif preys[geneID].hcip_sy5y < hcip_sy5y:
			preys[geneID].hcip_sy5y = hcip_sy5y

	# get the description for every gene
	for geneID, prey_obj in preys.iteritems():
		if geneID in gene_descr_dict:
			prey_obj.descr = gene_descr_dict[geneID]

	# remove proteasome-related proteins
	geneIDs_to_remove = set()
	for geneID, prey_obj in preys.iteritems():
		if "proteas" in prey_obj.descr:
			geneIDs_to_remove.add(geneID)
	for geneID in geneIDs_to_remove:
		del preys[geneID]

	# determine by which cell lines a prey has been found
	for geneID,prey_obj in preys.iteritems():
		if prey_obj.num_cell_lines == 2:
			prey_obj.cell_lines = "b"
		else:
			if prey_obj.num_hek > 0:
				prey_obj.cell_lines = "h"
			else:
				prey_obj.cell_lines = "s"

	preys[geneID_bait].cell_lines = "bait"

	# write the seed file
	target = open(config.output_path + bait + "_seed_file.txt","w")
	target.write("geneID\tsymbol\tdescr\tnum_hek\tnum_sy5y\tnum_cell_lines\tcell_lines\thcip_cutoff_sy5y\tbait\n")
	for geneID,prey_obj in preys.iteritems():
		target.write(str(geneID) + "\t" + prey_obj.symbol + "\t" + prey_obj.descr + "\t" + str(prey_obj.num_hek) + \
					 "\t" + str(prey_obj.num_sy5y) + "\t" + str(prey_obj.num_cell_lines) + "\t" + \
					 prey_obj.cell_lines + "\t" + str(prey_obj.hcip_sy5y) + "\t" + bait + "\n")
	target.close()


def build_HUN_seed_gene_file(baits,outfile):

	target = open(outfile,"w")
	target.write("geneID\tsymbol\tdescr\torigin\tbaits\tnum_baits\tcell_lines\tY2H_int\n")

	preys = {}

	for bait in baits:
		file1 = open(config.output_path + bait + '_seed_file.txt','r')
		entries = file1.readlines()
		file1.close()
		for line in entries[1:]:
			tab_list = str.split(line[:-1],'\t')
			geneID = tab_list[0]
			cell_lines = tab_list[6]
			if geneID not in preys:
				preys[geneID] = UBE3A_seed_gene_file_assembly.Prey(geneID)
				preys[geneID].symbol = tab_list[1]
				preys[geneID].cell_lines = []
				if preys[geneID].symbol in baits:
					preys[geneID].origin = "bait"
				else:
					preys[geneID].origin = "prey"
				preys[geneID].descr = tab_list[2]
				preys[geneID].baits = []
			preys[geneID].cell_lines.append(cell_lines)
			preys[geneID].baits.append(bait)
			# get information on Y2H UBE3A interactors
			if bait == 'UBE3A':
				Y2H_int = tab_list[8]
				preys[geneID].Y2H_int = Y2H_int

	for geneID, prey_obj in preys.iteritems():
		# decide on the final cell line annotation
		yeast = None
		bait = 1
		b = None
		s = None
		h = None
		cl_list = prey_obj.cell_lines
		for cl in prey_obj.cell_lines:
			if cl == 'bait':
				bait = 1
			else:
				if cl.find('yeast') > -1:
					yeast = 1
				if cl[0] == 's':
					s = 1
				elif cl[0] == 'b':
					b = 1
				elif cl[0] == 'h':
					h = 1
		if (yeast and b) or (yeast and s and h):
	  		prey_obj.cell_lines = 'b_yeast'
		elif yeast and h:
	  		prey_obj.cell_lines = 'h_yeast'
		elif yeast and s:
	  		prey_obj.cell_lines = 's_yeast'
		elif yeast:
			prey_obj.cell_lines = 'yeast'
		elif (h and s) or b:
			prey_obj.cell_lines = 'b'
		elif h:
			prey_obj.cell_lines = 'h'
		elif b:
			prey_obj.cell_lines = 'b'
		elif s:
			prey_obj.cell_lines = 's'
		else:
			prey_obj.cell_lines = 'bait'

		target.write(geneID + '\t' + prey_obj.symbol + '\t' + prey_obj.descr + '\t' + prey_obj.origin + '\t' + \
					 '|'.join(prey_obj.baits) + '\t' + str(len(prey_obj.baits)) + '\t' + prey_obj.cell_lines + '\n')

	target.close()


if __name__ == '__main__':

	entrez_file = config.input_path + 'entrez_gene.gene_info_human_20161113.tsv'
	gene_descr_dict = UBE3A_seed_gene_file_assembly.get_geneID_descr_dict(entrez_file)

	baits = [("CAMK2D",817),("MAPK6",5597),("ECH1",1891),("ECI2",10455),("HIF1AN",55662),("NEURL4",84461)]
	for bait_tup in baits:
		bait = bait_tup[0]
		geneID_bait = bait_tup[1]
		get_seed_gene_file(gene_descr_dict,bait,geneID_bait)

	baits = ['HERC2','NEURL4','MAPK6','ECH1','ECI2']
	outfile = config.output_path + "HUN_complex_seed_file.txt"
	build_HUN_seed_gene_file(baits,outfile)

	baits = ['HERC2','NEURL4','MAPK6','ECH1','ECI2','UBE3A']
	outfile = config.output_path + "HUN_UBE3A_seed_file.txt"
	build_HUN_seed_gene_file(baits,outfile)
