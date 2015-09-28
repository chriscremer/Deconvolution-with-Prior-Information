



import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np

from os.path import expanduser
home = expanduser("~")

import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
sys.path.insert(0,'..')
sys.path.insert(0,'../simulated_data')


import make_real_simulated_data as mrsd

from sklearn.decomposition import NMF

from deconvol_model import Deconvol


if __name__ == "__main__":



	models = ['nmf', 'nnls']
	k = 10

	W_L1_error = [[] for x in models]
	Z_L1_error = [[] for x in models]
	noise_amount = [0.0, .2, .4, .6, .8, 1.]
	for noise1 in noise_amount:
		print '\n\n\nNoise ' + str(noise1)
		#Make data
		print 'Making data...'
		subpops, fractions, X = mrsd.make_and_return(n_subpops=k, n_samps=100, probability_of_zero=.7, noise=noise1)
		print 'X ' + str(X[0][:10])
		print 'subpops ' + str(subpops[0][:10])
		print 'fractions ' + str(fractions[0])
		print 'Preprocessing ...'
		#remove bottom 25% of lowest mean expressed genes
		genes_mean_exps = [np.mean(x) for x in X.T]
		#largest to smallest
		sorted_indexes = np.argsort(genes_mean_exps)[::-1]
		#keep top 75%
		sorted_indexes = sorted_indexes[:len(sorted_indexes)*.8]
		back_in_order = sorted(sorted_indexes)
		X = X.T[back_in_order].T
		subpops = subpops.T[back_in_order].T
		print X.shape
		print 'X ' + str(X[0][:10])
		print 'subpops ' + str(subpops[0][:10])

		#scale by mean
		mean_list = [np.mean(x) for x in X.T]
		for gene in range(len(X.T)):
			X.T[gene] = X.T[gene]/mean_list[gene]
		#scale subpops too
		for gene in range(len(subpops.T)):
			subpops.T[gene] = subpops.T[gene]/mean_list[gene]
		print 'X ' + str(X[0][:10])
		print 'subpops ' + str(subpops[0][:10])

		#remove bottom 25% of lowest var expressed genes
		genes_var_exps = [np.std(x) for x in X.T]
		#largest to smallest
		sorted_indexes = np.argsort(genes_var_exps)[::-1]
		#keep top 75%
		sorted_indexes = sorted_indexes[:len(sorted_indexes)*.1]
		back_in_order = sorted(sorted_indexes)
		X = X.T[back_in_order].T
		subpops = subpops.T[back_in_order].T
		print X.shape
		print 'X ' + str(X[0][:10])
		print 'subpops ' + str(subpops[0][:10])



		# X_minus_low_genes = []
		# for gene in range(len(X.T)):
		# 	if np.mean(X.T[gene]) > .5:
		# 		X_minus_low_genes.append(X.T[gene])
		# X = np.array(X_minus_low_genes).T
		# print X.shape

		# for gene in range(len(X.T)):
		# 	X.T[gene] = (X.T[gene])/np.mean(X.T[gene])

		# X_minus_low_var_genes = []
		# for gene in range(len(X.T)):
		# 	if np.std(X.T[gene]) > .5:
		# 		X_minus_low_var_genes.append(X.T[gene])
		# X = np.array(X_minus_low_var_genes).T
		# print X.shape
		

		for model in range(len(models)):

			#Select model
			if model == 0:

				print 'NMF'
				decomposer = NMF(n_components=k, sparseness='components', max_iter=1000, tol=.001)
				decomposer.fit(X)
				W = decomposer.transform(X)
				Z = decomposer.components_ 
				print 'Z shape ' + str(Z.shape)

			if model == 1:

				print 'NNLS'
				decomposer = Deconvol(n_components=k, rand_inits=1)
				decomposer.fit(X)
				W = decomposer.transform(X)
				Z = decomposer.components_ 
				print 'Z shape ' + str(Z.shape)

			#scale W for each sample so that sum = 1
			for i in range(len(W)):
				W[i] = W[i]/sum(W[i])

			print 'Matching...'
			#Hungarian Algorithm
			norm_matrix = []
			for learned_profile in range(len(Z)):
				this_samp = []
				for real_profile in range(len(subpops)):
					# this_samp.append(np.linalg.norm(Z[learned_profile] - subpops[real_profile]))
					# print 'Learned ' + str(learned_profile) + ' ' + str(Z[learned_profile][:6])
					# print 'Real ' + str(real_profile) + ' ' + str(subpops[real_profile][:6])
					this_samp.append(sum(abs(Z[learned_profile] - subpops[real_profile])))
				norm_matrix.append(this_samp)
			# for list1 in norm_matrix:
			# 	print str(['%.2f' % elem for elem in list1])
			from munkres import Munkres, print_matrix
			m = Munkres()
			indexes = m.compute(norm_matrix)
			indexes2 = [x[1] for x in indexes]
			print indexes2
			# print_matrix(norm_matrix, msg='Lowest cost through this matrix:')
			total = 0
			for row, column in indexes:
				value = norm_matrix[row][column]
				# value2 = np.linalg.norm(TZ[row] - real_profiles[column])
				total += value
				print '(%d, %d) -> %f' % (row, column, value)
			print 'total cost: %f' % total


			# for learned_profile in range(len(Z)):
			# 	for real_profile in range(len(subpops)):
			# 		print '(%d, %d) -> %f' % (learned_profile, real_profile, norm_matrix[learned_profile][real_profile])

			# print 'Fraction: ' + str(fractions[0]) + ' W ' + str(W[0])
			# print 'Fraction: ' + str(fractions[1]) + ' W ' + str(W[1])

			# W = W.T[indexes2].T
			rearranged_fractions = fractions.T[indexes2].T

			# print 'Fraction: ' + str(rearranged_fractions[0]) + ' W ' + str(W[0])
			# print 'Fraction: ' + str(rearranged_fractions[1]) + ' W ' + str(W[1])

			W_error = sum(sum(abs(rearranged_fractions - W)))

			W_L1_error[model].append(W_error)

			print 'Error = ' + str(W_error)



	#plot the errors
	print noise_amount
	for model in range(len(models)):

		print W_L1_error[model]

		plt.plot(noise_amount, W_L1_error[model], label=models[model])


	plt.xlabel('Noise')
	plt.ylabel('Fraction Error')
	plt.legend()

	plt.savefig('error_with_noise2.png')
	print 'Saved plot'










