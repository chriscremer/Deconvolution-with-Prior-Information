
import os
from os.path import expanduser
home = expanduser("~")

import csv

import numpy as np

import pickle

import return_gene_data as rgd

chr_lens = {
			'1': 248956422,
			'2': 242193529,
			'3': 198295559,
			'4': 190214555,
			'5': 181538259,
			'6': 170805979,
			'7': 159345973,
			'8': 145138636,
			'9': 138394717,
			'10': 133797422,
			'11': 135086622,
			'12': 133275309,
			'13': 114364328,
			'14': 107043718,
			'15': 101991189,
			'16': 90338345,
			'17': 83257441,
			'18': 80373285,
			'19': 58617616,
			'20': 64444167,
			'21': 46709983,
			'22': 50818468,
			'X': 156040895,
			'Y': 57227415
			}

centromere_begin = {
			'1': 121500000,
			'2': 90500000,
			'3': 87900000,
			'4': 48200000,
			'5': 46100000,
			'6': 58700000,
			'7': 58000000,
			'8': 43100000,
			'9': 47300000,
			'10': 38000000,
			'11': 51600000,
			'12': 33300000,
			'13': 16300000,
			'14': 16100000,
			'15': 15800000,
			'16': 34600000,
			'17': 22200000,
			'18': 15400000,
			'19': 24400000,
			'20': 25600000,
			'21': 10900000,
			'22': 12200000,
			'X': 58100000,
			'Y': 11600000
			}

centromere_end = {
			'1': 128900000,
			'2': 96800000,
			'3': 93900000,
			'4': 52700000,
			'5': 50700000,
			'6': 63300000,
			'7': 61700000,
			'8': 48100000,
			'9': 50700000,
			'10': 42300000,
			'11': 55700000,
			'12': 38200000,
			'13': 19500000,
			'14': 19100000,
			'15': 20700000,
			'16': 38600000,
			'17': 25800000,
			'18': 19000000,
			'19': 28600000,
			'20': 29400000,
			'21': 14300000,
			'22': 17900000,
			'X': 63000000,
			'Y': 13400000
			}

#GET IDS OF PBCAs
sample_ids = []
sample_type = []
with open(home + '/../jawinter/scratch/pwgs/share/summary_tables/summary_table.consprelim.sample-5000.psub.txt', 'rb') as csvfile:
	reader = csv.reader(csvfile, delimiter='\t')
	for row in reader:
		# print row
		# if row[0] == 'BRCA':
		# if row[0] == 'BRCA':
		# if row[0] == 'PACA':

			# print row
			sample_ids.append(row[1])
			sample_type.append(row[0])


# def analyze_chroms(stuff):
# 	'''
# 	Returns a list of the percent each chromosome has changed.  
# 	'''

# 	#First get the clonal frequecy, ie how much of the sample is tumour, aka purity
# 	clonal_freq = 0
# 	for i in range(len(stuff)):

# 		# print stuff[i][6]
# 		if float(stuff[i][6]) > clonal_freq:
# 			clonal_freq = float(stuff[i][6])

# 	#Only care about cnvs that are in atleast 80% of the tumour cells. ie .8 * clonal_freq

# 	#weight the cnvs by the length of the cnv compared to the length of the chromosome
# 	#so it could cnv of 3 on half, so 3*.5 then the rest will be .5*2 because its normal
# 	# now i have to keep track of which chromosome im on, since I need the total not changed chromosome length


# 	percent_change = []


# 	#Put all copy numbers in a dict
# 	copy_numb_dict = {}

# 	#for each chromosome
# 	for chrom in range(1,23):

# 		#keep track of how much of the chrom isnt changed
# 		normal_chrom_frac = 1

# 		#look at every entry about the sample, even though its in order
# 		for i in range(len(stuff)):
# 			#if its the current chromosome
# 			if stuff[i][0] == str(chrom):
# 				#sometimes its NA, ignore it
# 				if stuff[i][3] != 'NA':
# 					#sometimes its a large number, ignore it.
# 					if float(stuff[i][3]) < 11.:
# 						#only if its a significant amount of the clonal population
# 						if float(stuff[i][6]) > clonal_freq*.8:

# 							#fraction of chromosome
# 							len_cnv = int(stuff[i][2]) - int(stuff[i][1])
# 							frac_cnv = float(len_cnv) / chr_lens[str(chrom)]
# 							# print frac_cnv, float(stuff[i][3])
# 							#update how much is left normal
# 							normal_chrom_frac = normal_chrom_frac - frac_cnv

# 							#if its the first one for this chrom into the dict
# 							if stuff[i][0] not in copy_numb_dict:
# 								copy_numb_dict[stuff[i][0]] = [float(stuff[i][3]) * frac_cnv]
# 							else:
# 								copy_numb_dict[stuff[i][0]].append(float(stuff[i][3]) * frac_cnv)

# 		#now that weve seen every entry, add the frac that hasnt changed
# 		#if its the first one for this chrom into the dict
# 		if str(chrom) not in copy_numb_dict:
# 			copy_numb_dict[str(chrom)] = [2. * normal_chrom_frac]
# 		else:
# 			copy_numb_dict[str(chrom)].append(2. * normal_chrom_frac)


# 	#Take avergae copy number
# 	for i in range(1,23):

# 		if str(i) in copy_numb_dict:
# 			percent_change.append(sum(copy_numb_dict[str(i)]))
# 		else:
# 			percent_change.append(2.)


# 	return percent_change




# def run_analyze_chroms_on_all():


# 	all_cnv_files =  os.listdir(home + '/../jawinter/scratch/pwgs/cnvs/vanloo_wedge')

# 	n_found =0
# 	n_sig_change = 0
# 	for s_id in sample_ids:
# 		# found = 0

# 		for f in all_cnv_files:

# 			if s_id in f:

# 				print s_id


# 				stuff = []
# 				with open(home + '/../jawinter/scratch/pwgs/cnvs/vanloo_wedge/' + f, 'rb') as csvfile:
# 					reader = csv.reader(csvfile, delimiter='\t')
# 					for row in reader:
# 						if row[0] == 'chromosome':
# 							continue
# 						stuff.append(row)

# 				chrom_changes = analyze_chroms(stuff)
# 				print chrom_changes
# 				dif_from_normal = np.mean(abs(np.array(chrom_changes) - np.array([2.]*len(chrom_changes))))
# 				print dif_from_normal
# 				if dif_from_normal > .05:
# 					n_sig_change +=1
# 				print

# 				n_found +=1

# 				break


# 	# print len(sample_ids)
# 	# print len(all_cnv_files)
# 	print n_found
# 	print n_sig_change








#Now I need to look for ones that have a subclonal change that is seen in multiple samples
# First look for how many samples have a subclonal change

def analyze_chroms2(stuff):

	#First get the clonal frequecy, ie how much of the sample is tumour, aka purity
	clonal_freq = 0
	for i in range(len(stuff)):
		if float(stuff[i][6]) > clonal_freq:
			clonal_freq = float(stuff[i][6])


	#this will keep track of if a chrom has a subclonal cnv
	sig_subclones = [False]*23

	subpop_freq = -1

	#for each chromosome
	for chrom in range(1,23):
		#look at every entry about the sample, even though its in order
		for i in range(len(stuff)):
			#if its the current chromosome
			if stuff[i][0] == str(chrom):
				#sometimes its NA, ignore it
				if stuff[i][3] != 'NA':
					#sometimes its a large number, ignore it.
					if float(stuff[i][3]) < 11.:

						#make sure there is a cnv, because it does report normal areas too
						if float(stuff[i][3]) < 2: # or float(stuff[i][3]) > 2:
						# if float(stuff[i][3]) == 0.:

							# if it ends before the centromere
							#this is to make sure there realtively in same location
							# if int(stuff[i][2]) < centromere_begin[str(chrom)]:

							#if begins in second half
							# if int(stuff[i][1]) > centromere_end[str(chrom)]:

								#only if its a significant amount of the clonal population
								if float(stuff[i][6]) > clonal_freq*.2 and float(stuff[i][6]) < clonal_freq*.9:

									#fraction of chromosome
									len_cnv = int(stuff[i][2]) - int(stuff[i][1])
									frac_cnv = float(len_cnv) / chr_lens[str(chrom)]

									if frac_cnv > .3:

										if chrom == 15:# or chrom == 7:
											# 	print stuff[i], clonal_freq
											print stuff[i], clonal_freq
											subpop_freq = float(stuff[i][6])

										sig_subclones[chrom] = True
										break
	proportions = []
	#Normal
	proportions.append(1.-clonal_freq)
	#IF subpop
	if sig_subclones[15]:
		proportions.append(subpop_freq)
		proportions.append(clonal_freq - subpop_freq)
	else:
		proportions.append(clonal_freq)



	return sig_subclones, proportions




def run_analyze_chroms2_on_all():

	all_cnv_files =  os.listdir(home + '/../jawinter/scratch/pwgs/cnvs/vanloo_wedge')

	print 'Length of sample_ids', len(sample_ids)


	total = np.zeros(23)

	all_proportions = []
	cnv_ids_with_gene_exp = []
	rna_ids_with_gene_exp = []

	n_found =0
	n_found_with_gene_exp = 0
	for s_id in range(len(sample_ids)):
		if s_id %100 == 0:
			print s_id, n_found_with_gene_exp
		for f in all_cnv_files:
			if sample_ids[s_id] in f:

				# print
				# print s_id
				n_found +=1

				gene_name = rgd.find_gene_exp_name(sample_ids[s_id])
				if gene_name == '':
					# print 'Not here'
					break
				else:
					# print gene_name
					n_found_with_gene_exp +=1
					# print 'with gene exp:', n_found_with_gene_exp


				stuff = []
				with open(home + '/../jawinter/scratch/pwgs/cnvs/vanloo_wedge/' + f, 'rb') as csvfile:
					reader = csv.reader(csvfile, delimiter='\t')
					for row in reader:
						if row[0] == 'chromosome':
							continue
						stuff.append(row)

				chrom_changes, proportions = analyze_chroms2(stuff)

				all_proportions.append(proportions)
				cnv_ids_with_gene_exp.append(sample_ids[s_id])
				rna_ids_with_gene_exp.append(gene_name)

				if chrom_changes[15]:
					print sample_type[s_id] 
				# print s_id
				# print chrom_changes

				chrom_changes_ints = np.array([int(x) for x in chrom_changes])

				total = total + chrom_changes_ints
				break


	for i in range(len(total)):
		print str(i) + ': ' + str(total[i])

	print n_found
	print n_found_with_gene_exp
	print len(all_proportions)
	# print all_proportions[:10]

	return all_proportions, rna_ids_with_gene_exp, cnv_ids_with_gene_exp







if __name__ == '__main__':


	# all_proportions, rna_ids_with_gene_exp, cnv_ids_with_gene_exp = run_analyze_chroms2_on_all()

	# with open('proportions_and_ids.pkl', 'wb') as f:
	# 	pickle.dump( [all_proportions, rna_ids_with_gene_exp, cnv_ids_with_gene_exp], f)

	with open('proportions_and_ids.pkl', 'rb') as f:
		stuff = pickle.load(f)
	all_proportions, rna_ids_with_gene_exp, cnv_ids_with_gene_exp = stuff

	print 'Getting data'
	gene_exp_data, samp_names = rgd.get_gene_expression_matrix(cnv_ids_with_gene_exp)

	print gene_exp_data.shape
	print len(samp_names)

	with open('gene_exps_and_names.pkl', 'wb') as f:
		pickle.dump( [gene_exp_data, samp_names], f)

	print 'saved. All done'












