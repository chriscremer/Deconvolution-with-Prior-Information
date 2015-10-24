



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
from deconvol_model import DIFI_strict
from deconvol_model import DIFI_w_deviates
from deconvol_model import Deconvol_normalized
from deconvol_model import DIFI_strict_top5
from deconvol_model import DIFI_deviate_top5
from deconvol_model import DIFI_nmf_deviate_top5
from deconvol_model import DIFI_nmf_deviate_matchW
from deconvol_model import DIFI_nmf_top5_nmf
from deconvol_model import DIFI_match_dev_optZ


from perturb_fractions import distribute_evenly




if __name__ == "__main__":

	# models = ['DIFI_match_dev_optZ']
	models = ['Deconvol_normalized', 'DIFI_match_dev_optZ', 'DIFI_match_dev_optZ_off_by_25', 'DIFI_match_dev_optZ_off_by_100']
	# models = ['Deconvol_normalized', 'DIFI_nmf_deviate_matchW', 'Difi_ndm_off_by_25', 'Difi_ndm_off_by_50', 'Difi_ndm_off_by_100']
	# models = ['Deconvol_normalized', 'DIFI_nmf_deviate_matchW', 'DIFI_nmf_top5_nmf']
	# models = ['Deconvol_normalized', 'DIFI_nmf_deviate_top5', 'DIFI_off_by_25', 'DIFI_off_by_50', 'DIFI_off_by_100']
	# models = ['Deconvol_normalized', 'DIFI_nmf_deviate_top5', 'DIFI_off_by_50', 'DIFI_off_by_75', 'DIFI_off_by_100']
	# models = ['Deconvol_normalized', 'DIFI_strict_top5', 'DIFI_deviate_top5', 'negative_control', 'negative_control2',]
	# models = ['DIFI_strict', 'DIFI_w_deviates', 'Deconvol_normalized', 'DIFI_strict_top5', 'negative_control']
	# models = ['nmf', 'nnls', 'DIFI_strict', 'DIFI_w_deviates']
	# models = ['']
	k = 10
	n_rand_inits =1
	n_samps = 30
	# p_of_zero=.7
	tol=1e-1
	average_over_x_iters = 1
	noise1 = .5

	variables = [.1, .3, .5, .7, .9]
	# variables = [0., .5, 1.]
	# variables = [0.0, 1.]

	


	W_L1_error = [[] for x in models]
	Z_L1_error = [[] for x in models]
	results_W = []
	results_Z = []

	for iter_to_avg in range(average_over_x_iters):

		for var in variables:
			p_of_zero = var
			print '\n\n\nIter ' + str(iter_to_avg)
			print 'var ' + str(var)
			#Make data
			print 'Making data...'
			subpops, fractions, X = mrsd.make_and_return(n_subpops=k, n_samps=n_samps, probability_of_zero=p_of_zero, noise=noise1)
			# print 'fractions shape ' + str(fractions.shape)
			# print 'X ' + str(X[0][:10])
			# print 'subpops ' + str(subpops[0][:10])
			# print 'fractions ' + str(fractions[0])
			# print 'Preprocessing ...'
			#remove bottom 25% of lowest mean expressed genes
			genes_mean_exps = [np.mean(x) for x in X.T]
			#largest to smallest
			sorted_indexes = np.argsort(genes_mean_exps)[::-1]
			#keep top 75%
			sorted_indexes = sorted_indexes[:len(sorted_indexes)*.5]
			back_in_order = sorted(sorted_indexes)
			X = X.T[back_in_order].T
			subpops = subpops.T[back_in_order].T
			print X.shape
			# print 'X ' + str(X[0][:10])
			# print 'subpops ' + str(subpops[0][:10])

			#scale by mean
			mean_list = [np.mean(x) for x in X.T]
			for gene in range(len(X.T)):
				X.T[gene] = X.T[gene]/mean_list[gene]
			#scale subpops too
			for gene in range(len(subpops.T)):
				subpops.T[gene] = subpops.T[gene]/mean_list[gene]
			# print 'X ' + str(X[0][:10])
			# print 'subpops ' + str(subpops[0][:10])

			#remove bottom 25% of lowest var expressed genes
			genes_var_exps = [np.std(x) for x in X.T]
			#largest to smallest
			sorted_indexes = np.argsort(genes_var_exps)[::-1]
			#keep top 75%
			sorted_indexes = sorted_indexes[:len(sorted_indexes)*.03]
			back_in_order = sorted(sorted_indexes)
			X = X.T[back_in_order].T
			subpops = subpops.T[back_in_order].T
			print X.shape
			# print 'fractions[0] ' + str(fractions[0])

			#remove 0s from fractions
			cleaned_fractions = []
			for samp in range(len(fractions)):
				frac_list = []
				for ji in range(len(fractions[samp])):
					if fractions[samp][ji] > 0:
						frac_list.append(fractions[samp][ji])
				cleaned_fractions.append(frac_list)

			# print cleaned_fractions[0]
			print

			



			# print 'X ' + str(X[0][:10])
			# print 'subpops ' + str(subpops[0][:10])



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
				if models[model] == 'nmf':

					print 'NMF'
					# decomposer = NMF(n_components=k, sparseness='components', max_iter=1000, tol=.001)
					decomposer = NMF(n_components=k)
					decomposer.fit(X)
					W = decomposer.transform(X)
					Z = decomposer.components_ 
					# print 'Z shape ' + str(Z.shape)

				if models[model] == 'nnls':

					print 'NNLS'
					decomposer = Deconvol(n_components=k)
					decomposer.fit(X)
					W = decomposer.transform(X)
					Z = decomposer.components_ 
					# print 'Z shape ' + str(Z.shape)

				if models[model] == 'DIFI_strict':

					print 'DIFI_strict'
					decomposer = DIFI_strict(n_components=k, rand_inits=n_rand_inits, tol=tol)
					decomposer.fit(X, fractions)
					W = decomposer.transform(X)
					Z = decomposer.components_ 
					# print 'Z shape ' + str(Z.shape)


				if models[model] == 'DIFI_w_deviates':
					print 'DIFI_w_deviate'
					decomposer = DIFI_w_deviates(n_components=k)
					decomposer.fit(X, fractions)
					W = decomposer.transform(X)
					Z = decomposer.components_ 

				if models[model] == 'Deconvol_normalized':
					print 'Deconvol_normalized'
					decomposer = Deconvol_normalized(n_components=k, rand_inits=n_rand_inits, tol=tol)
					decomposer.fit(X)
					W = decomposer.transform(X)
					Z = decomposer.components_ 

				if models[model] == 'DIFI_strict_top5':
					print 'DIFI_strict_top5'
					decomposer = DIFI_strict_top5(n_components=k, rand_inits=n_rand_inits, tol=1e-1)
					decomposer.fit(X, fractions)
					W = decomposer.transform(X, fractions)
					Z = decomposer.components_ 

				if models[model] == 'negative_control':
					print 'negative_control'
					W = np.zeros((len(X), k))
					Z = np.zeros((k,len(X[0])))

				if models[model] == 'negative_control2':
					print 'negative_control2'
					W = np.array([[1./k]*k]*len(X))
					Z = np.array([[np.mean(X.T[gene]) for gene in range(len(X[0]))]]*k)

				if models[model] == 'DIFI_deviate_top5':
					print 'DIFI_deviate_top5'
					decomposer = DIFI_deviate_top5(n_components=k, rand_inits=n_rand_inits, tol=1e-1)
					decomposer.fit(X, fractions)
					W = decomposer.transform(X, fractions)
					Z = decomposer.components_ 

				if models[model] == 'DIFI_nmf_deviate_top5':
					print 'DIFI_nmf_deviate_top5'
					decomposer = DIFI_nmf_deviate_top5(n_components=k, rand_inits=n_rand_inits, tol=tol)
					decomposer.fit(X, fractions)
					W = decomposer.transform(X, fractions)
					Z = decomposer.components_ 

				if models[model] == 'DIFI_off_by_10':
					print 'DIFI_off_by_10'
					decomposer = DIFI_nmf_deviate_top5(n_components=k, rand_inits=n_rand_inits, tol=tol)
					off_by_10 = distribute_evenly(fractions, .1)
					decomposer.fit(X, off_by_10)
					W = decomposer.transform(X, off_by_10)
					Z = decomposer.components_ 


				if models[model] == 'DIFI_off_by_25':
					print 'DIFI_off_by_25'
					decomposer = DIFI_nmf_deviate_top5(n_components=k, rand_inits=n_rand_inits, tol=tol)
					off_by_10 = distribute_evenly(fractions, .25)
					decomposer.fit(X, off_by_10)
					W = decomposer.transform(X, off_by_10)
					Z = decomposer.components_ 


				if models[model] == 'DIFI_off_by_50':
					print 'DIFI_off_by_50'
					decomposer = DIFI_nmf_deviate_top5(n_components=k, rand_inits=n_rand_inits, tol=tol)
					off_by_10 = distribute_evenly(fractions, .5)
					decomposer.fit(X, off_by_10)
					W = decomposer.transform(X, off_by_10)
					Z = decomposer.components_ 

				if models[model] == 'DIFI_off_by_75':
					print 'DIFI_off_by_75'
					decomposer = DIFI_nmf_deviate_top5(n_components=k, rand_inits=n_rand_inits, tol=tol)
					off_by_10 = distribute_evenly(fractions, .75)
					decomposer.fit(X, off_by_10)
					W = decomposer.transform(X, off_by_10)
					Z = decomposer.components_ 

				if models[model] == 'DIFI_off_by_100':
					print 'DIFI_off_by_100'
					decomposer = DIFI_nmf_deviate_top5(n_components=k, rand_inits=n_rand_inits, tol=tol)
					off_by_10 = distribute_evenly(fractions, 1.)
					decomposer.fit(X, off_by_10)
					W = decomposer.transform(X, off_by_10)
					Z = decomposer.components_ 

				if models[model] == 'DIFI_nmf_deviate_matchW':
					print 'DIFI_nmf_deviate_matchW'
					decomposer = DIFI_nmf_deviate_matchW(n_components=k, rand_inits=n_rand_inits, tol=tol)
					# off_by_10 = distribute_evenly(fractions, 1.)
					decomposer.fit(X, cleaned_fractions)
					W = decomposer.W
					Z = decomposer.Z
					# W = decomposer.transform(X, off_by_10)
					# Z = decomposer.components_ 

				if models[model] == 'Difi_ndm_off_by_25':
					print 'Difi_ndm_off_by_25'
					decomposer = DIFI_nmf_deviate_matchW(n_components=k, rand_inits=n_rand_inits, tol=tol)
					decomposer.fit(X, distribute_evenly(cleaned_fractions, .25)) 
					W = decomposer.W
					Z = decomposer.Z

				if models[model] == 'Difi_ndm_off_by_50':
					print 'Difi_ndm_off_by_50'
					decomposer = DIFI_nmf_deviate_matchW(n_components=k, rand_inits=n_rand_inits, tol=tol)
					decomposer.fit(X, distribute_evenly(cleaned_fractions, .5)) 
					W = decomposer.W
					Z = decomposer.Z

				if models[model] == 'Difi_ndm_off_by_100':
					print 'Difi_ndm_off_by_100'
					decomposer = DIFI_nmf_deviate_matchW(n_components=k, rand_inits=n_rand_inits, tol=tol)
					decomposer.fit(X, distribute_evenly(cleaned_fractions, 1.)) 
					W = decomposer.W
					Z = decomposer.Z

				if models[model] == 'DIFI_nmf_top5_nmf':
					print 'DIFI_nmf_top5_nmf'
					decomposer = DIFI_nmf_top5_nmf(n_components=k, rand_inits=n_rand_inits, tol=tol)
					# off_by_10 = distribute_evenly(fractions, 1.)
					decomposer.fit(X, cleaned_fractions)
					W = decomposer.W
					Z = decomposer.Z


				if models[model] == 'DIFI_match_dev_optZ':
					print 'DIFI_match_dev_optZ'
					decomposer = DIFI_match_dev_optZ(n_components=k, rand_inits=n_rand_inits, tol=tol)
					# off_by_10 = distribute_evenly(fractions, 1.)
					decomposer.fit(X, cleaned_fractions)
					W = decomposer.W
					Z = decomposer.Z

				if models[model] == 'DIFI_match_dev_optZ_off_by_25':
					print 'DIFI_match_dev_optZ_off_by_25'
					decomposer = DIFI_match_dev_optZ(n_components=k, rand_inits=n_rand_inits, tol=tol)
					# off_by_10 = distribute_evenly(fractions, 1.)
					decomposer.fit(X, distribute_evenly(cleaned_fractions, .25))
					W = decomposer.W
					Z = decomposer.Z

				if models[model] == 'DIFI_match_dev_optZ_off_by_100':
					print 'DIFI_match_dev_optZ_off_by_100'
					decomposer = DIFI_match_dev_optZ(n_components=k, rand_inits=n_rand_inits, tol=tol)
					# off_by_10 = distribute_evenly(fractions, 1.)
					decomposer.fit(X, distribute_evenly(cleaned_fractions, 1.))
					W = decomposer.W
					Z = decomposer.Z

				#scale W for each sample so that sum = 1
				# for i in range(len(W)):
				# 	W[i] = W[i]/sum(W[i])

				# print 'Matching...'
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
				# print indexes2
				# print_matrix(norm_matrix, msg='Lowest cost through this matrix:')
				# total = 0
				# for row, column in indexes:
				# 	value = norm_matrix[row][column]
				# 	total += value
					# print '(%d, %d) -> %f' % (row, column, value)
				# print 'total cost: %f' % total


				# for learned_profile in range(len(Z)):
				# 	for real_profile in range(len(subpops)):
				# 		print '(%d, %d) -> %f' % (learned_profile, real_profile, norm_matrix[learned_profile][real_profile])

				# print 'Fraction: ' + str(fractions[0]) + ' W ' + str(W[0])
				# print 'Fraction: ' + str(fractions[1]) + ' W ' + str(W[1])

				# W = W.T[indexes2].T
				rearranged_fractions = fractions.T[indexes2].T
				rearranged_profiles = subpops[indexes2]

				# print 'Fraction: ' + str(rearranged_fractions[0]) + ' W ' + str(W[0])
				# print 'Fraction: ' + str(rearranged_fractions[1]) + ' W ' + str(W[1])

				# print 'real frac ' + str(rearranged_fractions[0]) 
				# print 'W[0] ' + str(W[0])
				# print 'real frac ' + str(rearranged_fractions[1]) 
				# print 'W[1] ' + str(W[1])

				W_error = sum(sum(abs(rearranged_fractions - W))) / var
				Z_error = sum(sum(abs(rearranged_profiles - Z))) / sum(sum(abs(rearranged_profiles)))

				# print 'sum of W ' + str(sum(sum(abs(rearranged_fractions))))
				# print 'sum of Z ' + str(sum(sum(abs(rearranged_profiles))))

				W_L1_error[model].append(W_error)
				Z_L1_error[model].append(Z_error)

				print 'W Error = ' + str(W_error)
				print 'Z Error = ' + str(Z_error)

				print

		results_W.append(W_L1_error)
		results_Z.append(Z_L1_error)

	#average the results
	#results is a list of lists (one for each iter)
	#each list has a list for each model
	#each list has an entry for the errors at some noise
	W_L1_error_avg = [[] for x in models]
	Z_L1_error_avg = [[] for x in models]

	#so for each model, get average error
	for mod in range(len(models)):
		for noi in range(len(variables)):

			noises_for_each_iter_W = []
			noises_for_each_iter_Z = []
			for ite in range(average_over_x_iters):

				noises_for_each_iter_W.append(results_W[ite][mod][noi])
				noises_for_each_iter_Z.append(results_Z[ite][mod][noi])

			W_L1_error_avg[mod].append(np.mean(noises_for_each_iter_W))
			Z_L1_error_avg[mod].append(np.mean(noises_for_each_iter_Z))


	#plot the errors
	plt.figure(1)

	ax = plt.subplot(211)
	# plt.subplot(211)
	for model in range(len(models)):
		# print W_L1_error[model]
		plt.plot(variables, W_L1_error_avg[model], label=models[model])
	plt.xlabel('Sparsity')
	plt.ylabel('|Wp - Wr| Per Sample')


	# Shrink current axis by 20%
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
	ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size':8})
	plt.xlim( .1, .9 )  


	# plt.legend(prop={'size':8}, loc=1)

	ax = plt.subplot(212)
	for model in range(len(models)):
		# print Z_L1_error[model]
		plt.plot(variables, Z_L1_error_avg[model], label=models[model])
	plt.ylabel('|Zp - Zr| / |Zr|')
	plt.xlabel('k=' + str(k) + ' | noise= ' + str(noise1) + ' | n_samps= ' + str(n_samps) + ' | tol= ' + str(tol) + ' | iters= ' + str(average_over_x_iters) + ' | lambda= ' + str(len(X[0])))



	# Shrink current axis by 20%
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
	ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size':8})
	plt.xlim( .1, .9 )  
	# plt.legend(prop={'size':8}, loc=1)


	plt.savefig('Performances_with_sparsity_oct16_1.png')
	print 'Saved plot'




# ax = plt.subplot(211)
# for i in xrange(5):
#     ax.plot(x, i * x, label='$y = %ix$'%i)

# # Shrink current axis by 20%
# box = ax.get_position()
# ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# # Put a legend to the right of the current axis
# ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))









