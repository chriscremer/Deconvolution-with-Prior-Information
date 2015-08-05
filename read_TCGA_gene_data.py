

import os

import numpy as np
import csv

def read_data_folder():

	data = []
	file_names = []

	directory = '/data1/morrislab/ccremer/ISOpure_PCAWG/batch_108_prostate/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/'

	manifest_file = '../../ISOpure_PCAWG/batch_108_prostate/file_manifest.txt'


	for file123 in os.listdir(directory):

		if 'genes.normalized_results' in file123:

			#see if its a normal or tumour
			with open(manifest_file, 'rU') as f:
				reader = csv.reader(f, delimiter='\t')
				for row in reader:
					#if the file names match
					if file123 == row[6]:
						#sample name
						sample_name = row[4]
						#if tumour 01, normal = 11
						if sample_name.endswith('01'):
							tumour = 1
						else: 
							tumour = 0
			f.close()

			#allow tumour or both
			if tumour == 1 or tumour == 0:

				samp = []
				with open(directory + file123) as f2:
					firstLine = 1
					for line2 in f2:
						if firstLine == 1:
							firstLine = 0
							continue
						firstSplit = line2.split()
						samp.append(float(firstSplit[1]))

				data.append(samp)
				file_names.append(file123)
				f2.close()

	data = np.array(data)

	return data, file_names




def read_cellularities():

	cellularity_file = '/data1/morrislab/ccremer/ISOpure_PCAWG/batch_108_prostate/Clinical/Biotab/nationwidechildrens.org_biospecimen_tumor_sample_prad.txt'

	biotab_cellularities = []

	with open(cellularity_file, 'rU') as f:
		reader = csv.reader(f, delimiter='\t')
		b = 0 #used to skip the first row
		for row in reader:

			if b <= 1:
				b += 1
				continue

			biotab_cellularities.append([row[1], row[4]])

			b += 1
	f.close()

	return biotab_cellularities




def match_cellularities_to_files(file_names, cellularities):

	manifest_file = '../../ISOpure_PCAWG/batch_108_prostate/file_manifest.txt'

	ordered_cellularities = []

	for file_name in file_names:
		with open(manifest_file, 'rU') as f:
			reader = csv.reader(f, delimiter='\t')
			
			for row in reader:
				#if the file names match
				if file_name == row[6]:
					#append sample name
					sample_name = row[4]
					#if tumour 01, normal = 11
					if sample_name.endswith('01'):

						for samp in cellularities:
							if samp[0] in sample_name:
								ordered_cellularities.append(samp[1])
								break
								
					if sample_name.endswith('11'):
						ordered_cellularities.append('100')
						break

					else:
						break

	return ordered_cellularities


def convert_cells_to_freqs(ordered_cellularities):


	freqs = []
	for cell in ordered_cellularities:

		cell_float = float(cell) / 100
		freqs.append(np.array([cell_float, 1.0-cell_float]))

	return np.array(freqs)