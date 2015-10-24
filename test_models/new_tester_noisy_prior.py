
import numpy as np

from os.path import expanduser
home = expanduser('~')

import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
sys.path.insert(0,'..')
sys.path.insert(0,'../simulated_data')
sys.path.insert(0, home+'/plotting')

import make_real_simulated_data as mrsd
import nice_plot
from perturb_fractions import add_noise

from sklearn.decomposition import NMF
from deconvol_model import ALternate_NNLS
from difi_model import Difi

def preprocess(X, subpops):

	#Simple way
	#Reduce size to speed up testing
	# X = X.T[:1000].T
	# subpops = subpops.T[:1000].T
	# print 'X ' + str(X.shape)
	# #Scale by mean
	# mean_list = [np.mean(x) for x in X.T]
	# for gene in range(len(X.T)):
	# 	X.T[gene] = X.T[gene]/mean_list[gene]
	# #Scale subpops too
	# for gene in range(len(subpops.T)):
	# 	subpops.T[gene] = subpops.T[gene]/mean_list[gene]


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

	return X,subpops


if __name__ == "__main__":


	models = ['NMF', 'Difi']

	#Data Parameters
	n_samps = 50
	n_comps = 5
	p_of_zero=.7
	noise1 = .7

	#Model Parameters
	n_rand_inits = 1
	tol=1e-2
	lambda1 = None

	#Testing Parameters
	average_over_x_iters = 10

	#Variable to test
	# x_values = [0.01,.2,.4,.6,.8,.99]
	# x_values = [1000,1,10000]
	x_values = [0., .2, .4, .6, .8, 1.]

	#Recording results
	this_iter_results_W = []
	this_iter_results_Z = []

	for iter_to_avg in range(average_over_x_iters):

		W_L2_error = [[] for x in models]
		Z_L2_error = [[] for x in models]

		for var in x_values:

			print '\n\n\nIter ' + str(iter_to_avg)
			print 'Var ' + str(var) 
			print 'Making data...'
			subpops, fractions, X = mrsd.make_and_return(n_subpops=n_comps, n_samps=n_samps, probability_of_zero=p_of_zero, noise=noise1)
			
			#Remove 0s from proportions
			cleaned_proportions = []
			for samp in range(len(fractions)):
				frac_list = []
				for ji in range(len(fractions[samp])):
					if fractions[samp][ji] > 0:
						frac_list.append(fractions[samp][ji])
				cleaned_proportions.append(frac_list)


			#Preprocess
			X, subpops = preprocess(X, subpops)

			#Run method
			for model in range(len(models)):

				if models[model] == 'NMF2':
					print models[model]
					decomposer = NMF(n_components=n_comps, tol=tol)
					decomposer.fit(X)
					W = decomposer.transform(X)
					Z = decomposer.components_ 

				if models[model] == 'NMF':
					print models[model]
					decomposer = ALternate_NNLS(n_components=n_comps, tol=tol)
					decomposer.fit(X)
					W = decomposer.W
					Z = decomposer.Z

				if models[model] == 'Difi':
					print models[model]
					decomposer = Difi(n_components=n_comps, tol=tol, lambda1=None)
					decomposer.fit(X, add_noise(cleaned_proportions, var))
					W = decomposer.W
					Z = decomposer.Z

				if models[model] == 'Difi_uniform':
					print models[model]
					decomposer = Difi(n_components=n_comps, tol=tol, lambda1=None)
					decomposer.fit(X, distribute_evenly(cleaned_proportions, 1.))
					W = decomposer.W
					Z = decomposer.Z


				#Hungarian Algorithm to match predicted to real Zi
				norm_matrix = []
				for learned_profile in range(len(Z)):
					this_samp = []
					for real_profile in range(len(subpops)):
						# this_samp.append(sum(abs(Z[learned_profile] - subpops[real_profile])))
						this_samp.append(np.linalg.norm(Z[learned_profile] - subpops[real_profile]))
					norm_matrix.append(this_samp)
				from munkres import Munkres, print_matrix
				m = Munkres()
				indexes = m.compute(norm_matrix)
				indexes2 = [x[1] for x in indexes]

				rearranged_fractions = fractions.T[indexes2].T
				rearranged_profiles = subpops[indexes2]


				#Calculate error
				# W_error = sum(sum(abs(rearranged_fractions - W)))
				# Z_error = sum(sum(abs(rearranged_profiles - Z))) / sum(sum(abs(rearranged_profiles)))
				W_error = np.linalg.norm(rearranged_fractions - W)
				Z_error = np.linalg.norm(rearranged_profiles - Z) / np.linalg.norm(rearranged_profiles)
				print 'W Error = ' + str(W_error)
				print 'Z Error = ' + str(Z_error)
				print

				#Store results
				W_L2_error[model].append(W_error)
				Z_L2_error[model].append(Z_error)

		#Store this iterations results
		this_iter_results_W.append(W_L2_error)
		this_iter_results_Z.append(Z_L2_error)


	#Get the average and standard deviation of the iterations for each model/variable 
	W_error_avg = [[] for x in models]
	Z_error_avg = [[] for x in models]
	W_error_std = [[] for x in models]
	Z_error_std = [[] for x in models]
	for mod in range(len(models)):
		for var in range(len(x_values)):
			var_for_each_iter_W = []
			var_for_each_iter_Z = []

			for ite in range(average_over_x_iters):
				var_for_each_iter_W.append(this_iter_results_W[ite][mod][var])
				var_for_each_iter_Z.append(this_iter_results_Z[ite][mod][var])

			W_error_avg[mod].append(np.mean(var_for_each_iter_W))
			Z_error_avg[mod].append(np.mean(var_for_each_iter_Z))
			W_error_std[mod].append(np.std(var_for_each_iter_W))
			Z_error_std[mod].append(np.std(var_for_each_iter_Z))


	# for i in range(len(this_iter_results_W)):
	# 	print this_iter_results_W[i]
	# 	print


	# print
	# print W_error_std
	# print Z_error_std
	# print var_for_each_iter_W
	# print var_for_each_iter_Z
	# nice_plot.plot_lines_2graphs_errorbars(x_values, W_error_avg, W_error_std, models, 'Regularization Strength', 'W Reconstruction Error', 'na', 1, 0)
	# nice_plot.plot_lines_2graphs_errorbars(x_values, Z_error_avg, Z_error_std, models, 'Regularization Strength', 'Z Reconstruction Error', 'noise_plot.png', 0, 1)

	nice_plot.plot_lines_2graphs_errorbars3(x_values, W_error_avg, W_error_std, models, 'Prior Noise', 'W Reconstruction Error', 'na', 1, 0, xlim=[0,1.01])
	nice_plot.plot_lines_2graphs_errorbars3(x_values, Z_error_avg, Z_error_std, models, 'Prior Noise', 'Z Reconstruction Error', 'noisy_prior_plot.png', 0, 1, xlim=[0,1.01])

