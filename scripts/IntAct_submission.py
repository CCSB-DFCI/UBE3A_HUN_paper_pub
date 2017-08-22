# script that generates the file for submitting the interaction data to IntAct for the UBE3A paper

import sys
import database_utils

# add the Y2H data to the submission file
# read in the Y2H file
# map the geneIDs to the corresponding new orf ID
# get the ensembl gene ID and Uniprot AC for the new ORF ID
# get a mapping of the E6AP proteoform
# write out
def format_Y2H_data(target,Y2H_file,cursor,exp_number_screen,exp_number_pt,int_count,fig_ref):

	gene_orf_dict = {}
	query = 'select distinct entrez_gene_id,orf_id from hi_ref.master_ref'
	cursor.execute(query)
	for row in cursor:
		gene_id = int(row[0])
		orf_id = int(row[1])
		if gene_id not in gene_orf_dict:
			gene_orf_dict[gene_id] = []
		gene_orf_dict[gene_id].append(orf_id)

	orf_outside_id_dict = {}
	query = 'select distinct orf_id,ensembl_gene_id,ensembl_protein_id,uniprot_AC_iso from bioinfo_katja.orf_class_map_ensg'
	cursor.execute(query)
	for row in cursor:
		orf_id = int(row[0])
		ensg_id = row[1]
		ensp_id = row[2]
		uni_id = row[3]
		if orf_id not in orf_outside_id_dict:
			orf_outside_id_dict[orf_id] = []
		orf_outside_id_dict[orf_id].append((ensg_id,ensp_id,uni_id))

	ube3a_dict = {'WT1':'Q05086-2','WT2':'Q05086-1','WT3':'Q05086-3','INA1':'Q05086-2','INA2':'Q05086-1','INA3':'Q05086-3'}

	file1 = open(Y2H_file,'r')
	entries = file1.readlines()
	file1.close()
	for line in entries[1:]:
		tab_list = str.split(line[:-1],'\t')
		ube3a_id = ube3a_dict[tab_list[3]]
		gene_id = int(tab_list[0])
		if gene_id in gene_orf_dict:
			orf_ids = gene_orf_dict[gene_id]
			if len(orf_ids) > 1:
				print gene_id,orf_ids
				sys.exit()
			else:
				orf_id = orf_ids[0]
				if orf_id in orf_outside_id_dict:
					outside_ids = orf_outside_id_dict[orf_id]
					if len(outside_ids) > 1:
						print gene_id,orf_id,outside_ids
						sys.exit()
					else:
						## screening step
						# bait
						target.write(str(exp_number_screen) + '\ttwo hybrid prey pooling approach\tpartial DNA sequence identification\tSaccharomyces cerevisiae / 4932\tCCSB ORFeome collection as prey library\t' + \
									str(int_count) + '\t' + fig_ref + '\t\t' + \
									ube3a_id + '\tUniProtKB\t\tbait\tunspecified role\tprotein\t\t\t\t\tHomo sapiens / 9606\tGal4 DNA-binding domain N-terminal fusion,pQZ212 vector,Y8930c yeast strain')
						if tab_list[3] == 'INA1':
							target.write('C820A\t')
						elif tab_list[3] == 'INA2':
							target.write('C843A\t')
						elif tab_list[3] == 'INA3':
							target.write('C840A\t')
						else:
							target.write('\t')
						target.write('\t\n')
						# prey
						target.write(str(exp_number_screen) + '\ttwo hybrid prey pooling approach\tpartial DNA sequence identification\tSaccharomyces cerevisiae / 4932\tCCSB ORFeome collection as prey library\t' + \
									str(int_count) + '\t' + fig_ref + '\t\t' + \
									outside_ids[0][2] + '\tUniProtKB\t\tprey\tunspecified role\tprotein\t' + str(orf_id) + '\tCCSB ORFeome collection\t' + outside_ids[0][1] + '\tEnsembl\tHomo sapiens / 9606\tGal4 activation domain N-terminal fusion,pDEST-AD vector,Y8800c yeast strain\t\t\n')
						int_count += 1

						## pairwise test step
						# bait
						target.write(str(exp_number_pt) + '\ttwo hybrid array\tpartial DNA sequence identification\tSaccharomyces cerevisiae / 4932\t\t' + \
									str(int_count) + '\t' + fig_ref + '\t\t' + \
									ube3a_id + '\tUniProtKB\t\tbait\tunspecified role\tprotein\t\t\t\t\tHomo sapiens / 9606\tGal4 DNA-binding domain N-terminal fusion,pQZ212 vector,Y8930c yeast strain')
						if tab_list[3] == 'INA1':
							target.write('C820A\t')
						elif tab_list[3] == 'INA2':
							target.write('C843A\t')
						elif tab_list[3] == 'INA3':
							target.write('C840A\t')
						else:
							target.write('\t')
						target.write('\t\n')
						# prey
						target.write(str(exp_number_pt) + '\ttwo hybrid array\tpartial DNA sequence identification\tSaccharomyces cerevisiae / 4932\t\t' + \
									str(int_count) + '\t' + fig_ref + '\t\t' + \
									outside_ids[0][2] + '\tUniProtKB\t\tprey\tunspecified role\tprotein\t' + str(orf_id) + '\tCCSB ORFeome collection\t' + outside_ids[0][1] + '\tEnsembl\tHomo sapiens / 9606\tGal4 activation domain N-terminal fusion,pDEST-AD vector,Y8800c yeast strain\t\t\n')
						int_count += 1

	return target, int_count


# reads in the seed files, reads in the corresponding file from CompPASS, extracts the gene IDs and maps them to
# the Uniprot ACs from the CompPASS file
# write out the exp info
def format_AP_MS_data(baits,comppass_path,seed_file_path,exp_number_hek,exp_number_sy5y,target,int_count,herc2_dict,ube3a_mut_dict,ube3a_uni_dict,fig_ref):

	celllines = {'SY5Y':'SH-SY5Y','293':'T-REx 293'}

	for bait in baits:

		bait_symbol = bait[0].split('-')[0]
		file1 = open(seed_file_path + bait_symbol + '_seed_file.txt','r')
		entries = file1.readlines()
		file1.close()
		seed_gene_dict = {}
		for line in entries[1:]:
			tab_list = str.split(line[:-1],'\t')
			gene_id = int(tab_list[0])
			symbol = tab_list[1]
			seed_gene_dict[gene_id] = symbol

		for cellline_short, cellline_long in celllines.items():

			try:
				if cellline_short == 'SY5Y':
					file1 = open(comppass_path + bait[0] + '_' + cellline_short + '_90.txt','Ur')
				else:
					file1 = open(comppass_path + bait[0] + '_' + cellline_short + '_95.txt','Ur')
			except IOError:
				print bait[0], cellline_short
				continue

			print bait, cellline_short
			entries = file1.readlines()
			file1.close()
			comp_dict = {}
			for line in entries[1:]:
				tab_list = str.split(line[:-1],'\t')
				if tab_list[1].isdigit():
					gene_id = int(tab_list[1])
				else:
					continue
				uni_ac = tab_list[3].split('|')[1]
				comp_dict[gene_id] = uni_ac

			# write out the bait
			if cellline_short == 'SY5Y':
				exp_number = exp_number_sy5y
			else:
				exp_number = exp_number_hek
			target.write(str(exp_number) + '\taffinity purification\tidentification by mass spectrometry\tHomo sapiens / 9606\t' + cellline_short + ' cell line\t' + \
						str(int_count) + '\t' + fig_ref + '\t\t')
			if bait_symbol == 'UBE3A':
				target.write(ube3a_uni_dict[bait[0].split('-')[1]])
			else:
				target.write(comp_dict[bait[1]])
			target.write('\tUniProtKB\t' + bait_symbol + '\tbait\tunspecified role\tprotein\t\t\t\t\tHomo sapiens / 9606\tN-terminal HA tag')
			if bait_symbol == 'HERC2':
				target.write(', fragment\t' + herc2_dict[bait[0]] + '\t\n')
			elif len(bait[0].split('-')) == 3 and bait[0].split('-')[2] == 'CA':
				target.write(', inactivating mutation\t' + ube3a_mut_dict[bait[0]] + '\t\n')
			else:
				target.write('\t\t\n')

			for seed_gene_id,symbol in seed_gene_dict.items():

				if seed_gene_id in comp_dict and symbol != bait_symbol:
					target.write(str(exp_number) + '\taffinity purification\tidentification by mass spectrometry\tHomo sapiens / 9606\t' + cellline_short + ' cell line\t' + \
								str(int_count) + '\t' + fig_ref + '\t\t' + \
								comp_dict[seed_gene_id] + '\tUniProtKB\t' + symbol + '\tprey\tunspecified role\tprotein\t\t\t\t\tHomo sapiens / 9606\t\t\t\n')

				else:
					print seed_gene_id,symbol


			int_count += 1

	return target,int_count


def format_follow_up_exps(target,int_count,exp_number,fig_refs):

	exp_info = str(exp_number) + '\tanti tag coimmunoprecipitation\tanti tag western blot\tHomo sapiens / 9606\tT-REx 293 cell line\t'

	target.write(exp_info + str(int_count) + '\t' + fig_refs[2] + '\t\t' + \
				'Q9UJU6-2\tUniProtKB\tDBNL\tbait\tunspecified role\tprotein\t\t\t\t\tHomo sapiens / 9606\tN-terminal V5 tag, vector pHAGE-P CMVt N-V5 GAW\t\t\n')
	target.write(exp_info + str(int_count) + '\t' + fig_refs[2] + '\t\t' + \
				'Q13557\tUniProtKB\tCAMK2D\tprey\tunspecified role\tprotein\t\t\t\t\tHomo sapiens / 9606\tN-terminal HA tag, vector pHAGE-P CMVt N-HA GAW\t\t\n')
	int_count += 1

	target.write(exp_info + str(int_count) + '\t' + fig_refs[2] + '\t\t' + \
				'Q9UJU6-2\tUniProtKB\tDBNL\tbait\tunspecified role\tprotein\t\t\t\t\tHomo sapiens / 9606\tN-terminal V5 tag, vector pHAGE-P CMVt N-V5 GAW\t\t\n')
	target.write(exp_info + str(int_count) + '\t' + fig_refs[2] + '\t\t' + \
				'Q9UQM7\tUniProtKB\tCAMK2A\tprey\tunspecified role\tprotein\t\t\t\t\tHomo sapiens / 9606\tN-terminal HA tag, vector pHAGE-P CMVt N-HA GAW\t\t\n')
	int_count += 1

	target.write(exp_info + str(int_count) + '\t' + fig_refs[2] + '\t\t' + \
				'Q9UJU6-1\tUniProtKB\tDBNL\tbait\tunspecified role\tprotein\t\t\t\t\tHomo sapiens / 9606\tN-terminal V5 tag, vector pHAGE-P CMVt N-V5 GAW\t\t\n')
	target.write(exp_info + str(int_count) + '\t' + fig_refs[2] + '\t\t' + \
				'Q13557\tUniProtKB\tCAMK2D\tprey\tunspecified role\tprotein\t\t\t\t\tHomo sapiens / 9606\tN-terminal HA tag, vector pHAGE-P CMVt N-HA GAW\t\t\n')
	int_count += 1

	target.write(exp_info + str(int_count) + '\t' + fig_refs[2] + '\t\t' + \
				'Q9UJU6-1\tUniProtKB\tDBNL\tbait\tunspecified role\tprotein\t\t\t\t\tHomo sapiens / 9606\tN-terminal V5 tag, vector pHAGE-P CMVt N-V5 GAW\t\t\n')
	target.write(exp_info + str(int_count) + '\t' + fig_refs[2] + '\t\t' + \
				'Q9UQM7\tUniProtKB\tCAMK2A\tprey\tunspecified role\tprotein\t\t\t\t\tHomo sapiens / 9606\tN-terminal HA tag, vector pHAGE-P CMVt N-HA GAW\t\t\n')
	int_count += 1

	target.write(exp_info + str(int_count) + '\t' + fig_refs[0] + '\t\t' + \
				'Q14566\tUniProtKB\tMCM6\tbait\tunspecified role\tprotein\t\t\t\t\tHomo sapiens / 9606\tN-terminal V5 tag, vector pHAGE-P CMVt N-V5 GAW\t\t\n')
	target.write(exp_info + str(int_count) + '\t' + fig_refs[0] + '\t\t' + \
				'Q05086-3\tUniProtKB\tUBE3A\tprey\tunspecified role\tprotein\t\t\t\t\tHomo sapiens / 9606\tN-terminal HA tag, vector pHAGE-P CMVt N-HA GAW, inactivating mutation\tC840A\t\n')
	int_count += 1

	target.write(exp_info + str(int_count) + '\t' + fig_refs[0] + '\t\t' + \
				'Q9Y2Z0\tUniProtKB\tSUGT1\tbait\tunspecified role\tprotein\t\t\t\t\tHomo sapiens / 9606\tN-terminal V5 tag, vector pHAGE-P CMVt N-V5 GAW\t\t\n')
	target.write(exp_info + str(int_count) + '\t' + fig_refs[0] + '\t\t' + \
				'Q05086-3\tUniProtKB\tUBE3A\tprey\tunspecified role\tprotein\t\t\t\t\tHomo sapiens / 9606\tN-terminal HA tag, vector pHAGE-P CMVt N-HA GAW, inactivating mutation\tC840A\t\n')
	int_count += 1

	target.write(exp_info + str(int_count) + '\t' + fig_refs[0] + '\t\t' + \
				'Q99613\tUniProtKB\tEIF3C\tbait\tunspecified role\tprotein\t\t\t\t\tHomo sapiens / 9606\tN-terminal V5 tag, vector pHAGE-P CMVt N-V5 GAW\t\t\n')
	target.write(exp_info + str(int_count) + '\t' + fig_refs[0] + '\t\t' + \
				'Q05086-3\tUniProtKB\tUBE3A\tprey\tunspecified role\tprotein\t\t\t\t\tHomo sapiens / 9606\tN-terminal HA tag, vector pHAGE-P CMVt N-HA GAW, inactivating mutation\tC840A\t\n')
	int_count += 1

	target.write(exp_info + str(int_count) + '\t' + fig_refs[0] + '\t\t' + \
				'P55036\tUniProtKB\tPSMD4\tbait\tunspecified role\tprotein\t\t\t\t\tHomo sapiens / 9606\tN-terminal V5 tag, vector pHAGE-P CMVt N-V5 GAW\t\t\n')
	target.write(exp_info + str(int_count) + '\t' + fig_refs[0] + '\t\t' + \
				'Q05086-3\tUniProtKB\tUBE3A\tprey\tunspecified role\tprotein\t\t\t\t\tHomo sapiens / 9606\tN-terminal HA tag, vector pHAGE-P CMVt N-HA GAW, inactivating mutation\tC840A\t\n')
	int_count += 1

	target.write(exp_info + str(int_count) + '\t' + fig_refs[1] + '\t\t' + \
				'Q05086-2\tUniProtKB\tUBE3A\tbait\tunspecified role\tprotein\t\t\t\t\tHomo sapiens / 9606\tN-terminal HA tag, vector pHAGE-P CMVt N-HA GAW, inactivating mutation\tC820A\t\n')
	target.write(exp_info + str(int_count) + '\t' + fig_refs[1] + '\t\t' + \
				'Q13625\tUniProtKB\tTP53BP2\tprey\tunspecified role\tprotein\t\t\t\t\tHomo sapiens / 9606\tN-terminal HA tag, vector pHAGE-P CMVt N-HA GAW\t\t\n')
	int_count += 1

	return target,int_count


if __name__ == '__main__':

	connect = database_utils.get_connection()
	cursor = connect.cursor()

	path = '/Users/luck/CCSB/Howley/2016/'
	Y2H_file = path + 'source_networks/Y3H_data.txt'

	outfile = path + 'manuscript/IntAct_submission_20170822.txt'

	target = open(outfile,'w')
	target.write('Experiment number\tInteraction detection method\tParticipant identification method\tHost organism\t' + \
				'Experiment annotations, comments\tInteraction number\tFigure legend\tInteraction annotations, comments\t' + \
				'Participant ID\tDatabase\tParticipant name\tExperimental role\tBiological role\tParticipant type\t' + \
				'Participant xreference\txref Database\tParticipant xreference\txref Database\tParticipant organism\tParticipant feature(s)\t' + \
				'Participant feature(s) range(s)\tParticipant expressed in\n')

	# write out the Y2H data
	exp_number_screen = 1
	exp_number_pt = 2
	int_count = 1
	fig_ref = 'Table 3'
	target, int_count = format_Y2H_data(target,Y2H_file,cursor,exp_number_screen,exp_number_pt,int_count,fig_ref)

	# write out the AP-MS data
	fig_ref = 'GitHub repository'
	baits = [("CAMK2D",817),("MAPK6",5597),("ECH1",1891),("ECI2",10455),("HIF1AN",55662),("NEURL4",84461),\
			 ("HERC2-Frag1",8924),("HERC2-Frag2",8924),("HERC2-Frag3",8924),("HERC2-Frag4",8924),("HERC2-Frag5",8924),("HERC2-Frag6",8924),\
			 ("UBE3A-ISO1-WT",7337),("UBE3A-ISO1-CA",7337),("UBE3A-ISO2-WT",7337),("UBE3A-ISO2-CA",7337),("UBE3A-ISO3-WT",7337),("UBE3A-ISO3-CA",7337)]
	comppass_path = path + 'input_files/CompPASS_files/'
	seed_file_path = path + 'output_files/'
	exp_number_hek = 4
	exp_number_sy5y = 3
	herc2_dict = {'HERC2-Frag1':'1-1000','HERC2-Frag2':'950-1750','HERC2-Frag3':'1700-2700','HERC2-Frag4':'2600-3600','HERC2-Frag5':'3550-4500','HERC2-Frag6':'4421-4834'}
	ube3a_mut_dict = {'UBE3A-ISO1-CA':'C820A','UBE3A-ISO2-CA':'C843A','UBE3A-ISO3-CA':'C840A'}
	ube3a_uni_dict = {'ISO1':'Q05086-2','ISO2':'Q05086-1','ISO3':'Q05086-3'}
	target, int_count = format_AP_MS_data(baits,comppass_path,seed_file_path,exp_number_hek,exp_number_sy5y,target,int_count,herc2_dict,ube3a_mut_dict,ube3a_uni_dict,fig_ref)

	# write out the follow-up experiments
	exp_number = 5
	fig_refs = ['Figure 2','Figure 5','Figure 7']
	target,int_count = format_follow_up_exps(target,int_count,exp_number,fig_refs)
	target.close()
