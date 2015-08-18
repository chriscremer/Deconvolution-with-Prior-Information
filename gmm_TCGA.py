
import numpy as np

import read_TCGA_gene_data as rtgd

from sklearn import mixture


def main():

	data_directory = '/data1/morrislab/ccremer/TCGA_data/breast_batch61/'

	X, file_names = rtgd.read_data_folder(data_directory)
	print 'real data shape ' + str(X.shape)
	print 'file names len ' + str(len(file_names))

	cellularities = rtgd.read_cellularities(data_directory)
	ordered_cellularities, sample_order = rtgd.match_cellularities_to_files(file_names, cellularities, data_directory)
	ordered_subtypes = rtgd.get_subtype(sample_order, data_directory)
	possible_subtypes = list(set(ordered_subtypes))

	min_components = 2
	max_components = 3

	for numb_components in range(min_components, max_components+1):

		gmm = mixture.GMM(n_components=numb_components)
		gmm.fit(X)


		bic = gmm.bic(X)
		print str(numb_components) + ' components. BIC ' + str(bic)

		print gmm.weights_
		#print gmm.means_
		#print gmm.covars_
		print gmm.converged_

		predictions = gmm.predict(X)

		for i in range(len(possible_subtypes)):

			counts = [0]*numb_components

			for j in range(len(ordered_subtypes)):

				if ordered_subtypes[j] == possible_subtypes[i]:

					counts[predictions[j]] = counts[predictions[j]] + 1

			print possible_subtypes[i]
			print counts
			print





		# for i in range(len(predictions)):
		# 	print str(predictions[i]) + '  ' + ordered_subtypes[i]

		# for i in range(len(ordered_subtypes)):
		# 	if ordered_subtypes[i] == 'Lum A':
		# 		print str(predictions[i]) + '  ' + ordered_subtypes[i]

		# for i in range(len(ordered_subtypes)):
		# 	if ordered_subtypes[i] == 'Lum B':
		# 		print str(predictions[i]) + '  ' + ordered_subtypes[i]

		# for i in range(len(ordered_subtypes)):
		# 	if ordered_subtypes[i] == 'Basal':
		# 		print str(predictions[i]) + '  ' + ordered_subtypes[i]

		# for i in range(len(ordered_subtypes)):
		# 	if ordered_subtypes[i] == 'Her2':
		# 		print str(predictions[i]) + '  ' + ordered_subtypes[i]

		# for i in range(len(ordered_subtypes)):
		# 	if ordered_subtypes[i] == 'normal':
		# 		print str(predictions[i]) + '  ' + ordered_subtypes[i]





if __name__ == "__main__":

	main()
