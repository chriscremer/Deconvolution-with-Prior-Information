

import os

import numpy as np
import csv

def read_data_folder(data_directory):

	data = []
	file_names = []

	directory = data_directory + 'RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/'

	manifest_file = data_directory + 'file_manifest.txt'


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




def read_cellularities(data_directory):

	cellularity_file = data_directory + 'Clinical/Biotab/nationwidechildrens.org_biospecimen_tumor_sample_brca.txt'

	biotab_cellularities = []

	with open(cellularity_file, 'rU') as f:
		reader = csv.reader(f, delimiter='\t')
		b = 0 #used to skip the first row
		for row in reader:

			if b <= 1:
				b += 1
				continue

			if 'prad' in cellularity_file:
				biotab_cellularities.append([row[1], row[4]])
			elif 'brca' in cellularity_file:
				biotab_cellularities.append([row[1], row[5]])

			b += 1
	f.close()

	return biotab_cellularities




def match_cellularities_to_files(file_names, cellularities, data_directory):

	manifest_file = data_directory + 'file_manifest.txt'

	ordered_cellularities = []
	sample_order = []

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
								sample_order.append(sample_name)
								
						break

					elif sample_name.endswith('11'):
						ordered_cellularities.append('100')
						sample_order.append(sample_name)
						break

					else:
						print 'HEYHEYHEY'
						break

			#print file_name
			#print ordered_cellularities

	return ordered_cellularities, sample_order


def convert_cells_to_freqs(ordered_cellularities):


	freqs = []
	for cell in ordered_cellularities:

		cell_float = float(cell) / 100
		freqs.append(np.array([cell_float, 1.0-cell_float]))

	return np.array(freqs)



def get_subtype(sample_order, data_directory):

	ordered_subtypes = []

	patient_file = data_directory + 'Clinical/Biotab/nationwidechildrens.org_clinical_patient_brca.txt' 

	for samp_name in sample_order:

		patient = samp_name[:12]

		with open(patient_file, 'rU') as f:
			reader = csv.reader(f, delimiter='\t')
			count = 0
			for row in reader:
				#print count
				if count < 3:
					count += 1
					continue


				if row[1] != patient:
					#print row[1] + ' looking for ' + patient
					continue


				#print row[1] + ' ==  ' + patient

				er_status = row[43]
				her2_status = row[55]

				#print 'ER status ' + str(er_status)
				#print 'HER2 status ' + str(her2_status)

				subtype = 'none'
				if samp_name.endswith('11'):
					subtype = 'normal'
				elif er_status == 'Positive':
					if her2_status == 'Positive':
						subtype = 'Lum B'
					elif her2_status == 'Negative':
						subtype = 'Lum A'
				elif er_status == 'Negative':
					if her2_status == 'Positive':
						subtype = 'Her2'
					elif her2_status == 'Negative':
						subtype = 'Basal'

				#print 'Subtype ' + str(subtype) + '\n'

				ordered_subtypes.append(subtype)

				count +=1


		#print patient
		#print len(ordered_subtypes)

	return ordered_subtypes
