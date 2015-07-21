

import os

import numpy as np


def read_data_folder():

	data = []

	directory = '/data1/morrislab/ccremer/TCGA_data/rnaseqv2_clinical_onlytumours_25nov2014/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/'

	count =0

	for file123 in os.listdir(directory):

		if 'genes.normalized_results' in file123:
			samp = []
			#print file123
			with open(directory + file123) as f2:
				firstLine = 1
				for line2 in f2:
					if firstLine == 1:
						firstLine = 0
						continue
					firstSplit = line2.split()
					samp.append(firstSplit[1])

			data.append(samp)

			count += 1
			#if count > 40:
			#	break

	data = np.array(data)
	print 'real data shape ' + str(data.shape)

	return data