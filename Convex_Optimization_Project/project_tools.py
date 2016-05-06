


import numpy as np

import scipy

import itertools



def preprocess_data(X, subpops, dimensions_to_keep=100, normalize=True):

	#remove bottom 25% of lowest mean expressed genes
	genes_mean_exps = [np.mean(x) for x in X.T]
	#largest to smallest
	sorted_indexes = np.argsort(genes_mean_exps)[::-1]
	#keep top 75%
	sorted_indexes = sorted_indexes[:len(sorted_indexes)*.2]
	back_in_order = sorted(sorted_indexes)
	X = X.T[back_in_order].T
	subpops = subpops.T[back_in_order].T
	# print X.shape
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
	# print X.shape




	#Reduce size to speed up testing
	X = X.T[:dimensions_to_keep].T
	subpops = subpops.T[:dimensions_to_keep].T
	# print 'X ' + str(X.shape)

	# # #Scale by mean
	# mean_list = [np.mean(x) for x in X.T]
	# for gene in range(len(X.T)):
	# 	if mean_list[gene] == 0:
	# 		continue
	# 	X.T[gene] = X.T[gene]/mean_list[gene]
	# # #Scale subpops too
	# for gene in range(len(subpops.T)):
	# 	if mean_list[gene] == 0:
	# 		continue
	# 	subpops.T[gene] = subpops.T[gene]/mean_list[gene]

	return X, subpops


def get_avg_std(results_W, results_Z, models, var_list, n_iters):

	W_error_avg = [[] for x in models]
	Z_error_avg = [[] for x in models]
	W_error_std = [[] for x in models]
	Z_error_std = [[] for x in models]
	for m in range(len(models)):
		for var in range(len(var_list)):
			var_for_each_iter_W = []
			var_for_each_iter_Z = []

			for ite in range(n_iters):
				var_for_each_iter_W.append(results_W[ite][m][var])
				var_for_each_iter_Z.append(results_Z[ite][m][var])

			W_error_avg[m].append(np.mean(var_for_each_iter_W))
			Z_error_avg[m].append(np.mean(var_for_each_iter_Z))
			W_error_std[m].append(np.std(var_for_each_iter_W))
			Z_error_std[m].append(np.std(var_for_each_iter_Z))

	return W_error_avg, Z_error_avg, W_error_std, Z_error_std



def match_Z_error(W, Z, fractions, subpops):

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

	# W_error = np.linalg.norm(rearranged_fractions - W)
	# Z_error = np.linalg.norm(rearranged_profiles - Z) / np.linalg.norm(rearranged_profiles)


	#Spearman Cor
	# W_error = sum([abs(scipy.stats.spearmanr(rearranged_fractions[x], W[x])[0]) for x in range(len(W))]) / len(W)
	# Z_error = sum([abs(scipy.stats.spearmanr(rearranged_profiles[x], Z[x])[0]) for x in range(len(Z))]) / len(Z)

	#Pearson Cor
	# W_error = sum([abs(scipy.stats.pearsonr(rearranged_fractions[x], W[x])[0]) for x in range(len(W))]) / len(W)
	# Z_error = sum([abs(scipy.stats.pearsonr(rearranged_profiles[x], Z[x])[0]) for x in range(len(Z))]) / len(Z)

	#L1 error
	# W_error = []
	# for i in range(len(W)):
	# 	W_error.append(sum(abs(rearranged_fractions[i] - W[i])))
	# W_error = sum(W_error) / len(W)
	# Z_error = []
	# for i in range(len(Z)):
	# 	Z_error.append(sum(abs(rearranged_profiles[i] - Z[i])))
	# Z_error = sum(Z_error) / len(Z)

	#L2 error 
	M = rearranged_fractions - W
	sum1 = 0
	for i in range(len(M)):
		for j in range(len(M.T)):
			sum1 += (M[i][j])**2
	W_error = sum1 / len(W)
	M = rearranged_profiles - Z
	sum1 = 0
	for i in range(len(M)):
		for j in range(len(M.T)):
			sum1 += (M[i][j])**2
	Z_error = sum1 / len(Z)

	return W_error, Z_error




def assign_W_to_top_5_components(X, freqs, Z):

	new_W = []

	n_comps = len(Z)
	if n_comps >= 5:
		closest_x = 5
	else:
		closest_x = n_comps

	for samp in range(len(X)):

		#get 5 closest profiles to the sample
		dist_to_profiles = []
		for p in range(len(Z)):
			# dist_to_profiles.append(sum(abs(X[samp] - Z[p])))
			dist_to_profiles.append(np.linalg.norm(X[samp] - Z[p]))
		index_of_closest_profiles = np.argsort(dist_to_profiles)[:closest_x]
		Z_reduced = Z[index_of_closest_profiles]

		#get non zero freqs
		non_zero_freqs = []
		for freq in range(len(freqs[samp])):
			if freqs[samp][freq] > 0:
				non_zero_freqs.append(freqs[samp][freq])

		#if more than 5 freqs, remove lower freqs and re-normalize
		if len(non_zero_freqs) > closest_x:
			sorted_freqs = sorted(non_zero_freqs)[::-1]
			top_5 = np.array(sorted_freqs[:closest_x])
			top_5 = top_5 / sum(top_5)
			non_zero_freqs = list(top_5)

		#if less than 5, add some zeros
		if len(non_zero_freqs) < closest_x:
			while len(non_zero_freqs) < closest_x:
				non_zero_freqs.append(0.0)

		#find permuations and remove the duplicates
		perms = list(set(itertools.permutations(non_zero_freqs)))
		best_perm = -1
		best_perm_error = -1
		for perm in range(len(perms)):

			perm1 = np.array(perms[perm])
			X_i_hat = np.dot(perm1, Z_reduced)
			# perm_error = sum(abs(X_i_hat - X[samp]))
			perm_error = np.linalg.norm(X_i_hat - X[samp])

			if perm_error < best_perm_error or best_perm_error == -1:
				best_perm_error = perm_error
				best_perm = perm1

		#convert indexes back to the W vector indexes
		W_i = np.zeros(len(Z))
		for index in range(len(index_of_closest_profiles)):
			W_i[index_of_closest_profiles[index]] = best_perm[index]
		new_W.append(W_i)

	W = np.array(new_W)

	# W = use_all_components(X, W, Z)

	return W
