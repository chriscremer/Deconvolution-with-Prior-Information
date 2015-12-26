


import numpy as np

from scipy.optimize import nnls

from nnls import nnls_with_prior

import itertools

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt



def init_Z(X, n_components):

	#Initialize Z
	# Z = X[random.sample(range(len(X)), self.n_components)]
	Z = []
	used = []
	for i in range(n_components):
		if len(Z) == 0:
			#so that they start at the same place
			rand = 0
			# rand = random.randint(0,len(X)-1)
			Z.append(X[rand])
			used.append(rand)
		else:
			#select sample that is furthest from mean of already selected samples
			Z_array = np.array(Z)
			mean1 = [np.mean(x) for x in Z_array.T]
			furthest_samp = -1
			furthest_dist = -1
			for samp in range(len(X)):
				if samp in used:
					continue
				else:
					distance = sum(abs(X[samp] - mean1))
				if distance > furthest_dist or furthest_dist == -1:
					furthest_dist = distance
					furthest_samp = samp
			Z.append(X[furthest_samp])
			used.append(furthest_samp)
	Z = np.array(Z)
	return Z
	

def use_all_components(X, W, Z):

	#If no assignments to a profile, give it the one with the worst norm
	count = 0 #Numb of zero profiles found
	calculated = 0 #To know if I've aleady calculated the difference matrix
	for i in range(len(W.T)):
		if sum(W.T[i]) == 0:
			if calculated == 0:
				calculated = 1
				dif = X - np.dot(W, Z)
				norms = np.linalg.norm(dif, axis=1)
				sort_indexes = np.argsort(norms)[::-1] #descending order
			for j in range(len(W[sort_indexes[count]])):
				if j == i:
					W[sort_indexes[count]][j] = 1.0
				else:
					W[sort_indexes[count]][j] = 0.0
			count += 1

	return W


def optimize_z_using_nnls(X, W):

	#Optimize Z
	Z = []
	for d in range(len(X.T)):
		Z_d = nnls(W, X.T[d])[0]
		Z.append(Z_d)
	Z = np.array(Z).T

	return Z


def optimize_z_and_w_using_nnls_given_X_and_Z(X, Z, max_iters, tol):

	to_plot = []
	converged = 0
	for iter1 in range(max_iters):

		#Optimize W
		W = []
		for i in range(len(X)):
			W_i = nnls(Z.T, X[i])[0]
			W.append(W_i)
		W = np.array(W)

		W = use_all_components(X, W, Z)

		#scale W for each sample so that sum = 1
		for i in range(len(W)):
			W[i] = W[i]/sum(W[i])

		norm = np.linalg.norm(X - np.dot(W, Z))
		# print 'norm ' + str(norm)

		#Optimize Z
		Z = []
		for d in range(len(X.T)):
			Z_d = nnls(W, X.T[d])[0]
			Z.append(Z_d)
		Z = np.array(Z).T

		new_norm = np.linalg.norm(X - np.dot(W, Z))
		# print 'new_norm ' + str(new_norm)

		to_plot.append(new_norm)

		if (norm - new_norm) < tol:
			# print '# iters until optimized= ' + str(iter1)
			converged = 1
			break

	if converged == 0:
		print 'Did not converge1'

	return W, Z, to_plot


def assign_W_to_top_5_components(X, freqs, Z):

	new_W = []
	for samp in range(len(X)):

		#get 5 closest profiles to the sample
		dist_to_profiles = []
		for p in range(len(Z)):
			# dist_to_profiles.append(sum(abs(X[samp] - Z[p])))
			dist_to_profiles.append(np.linalg.norm(X[samp] - Z[p]))
		index_of_closest_profiles = np.argsort(dist_to_profiles)[:5]
		Z_reduced = Z[index_of_closest_profiles]

		#get non zero freqs
		non_zero_freqs = []
		for freq in range(len(freqs[samp])):
			if freqs[samp][freq] > 0:
				non_zero_freqs.append(freqs[samp][freq])

		#if more than 5 freqs, remove lower freqs and re-normalize
		if len(non_zero_freqs) > 5:
			sorted_freqs = sorted(non_zero_freqs)[::-1]
			top_5 = np.array(sorted_freqs[:5])
			top_5 = top_5 / sum(top_5)
			non_zero_freqs = list(top_5)

		#if less than 5, add some zeros
		if len(non_zero_freqs) < 5:
			while len(non_zero_freqs) < 5:
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

	W = use_all_components(X, W, Z)

	return W


def fit_prior_freqs(W, freqs):


	#TODO
	#match the Ws

	new_W = []
	for samp in range(len(W)):
		norm_matrix = []
		for freq in range(len(freqs[samp])):
			this_samp = []
			for w in range(len(W[samp])):
				this_samp.append(abs(freqs[samp][freq] - W[samp][w]))
			norm_matrix.append(this_samp)
		# for list1 in norm_matrix:
		# 	print str(['%.2f' % elem for elem in list1])
		from munkres import Munkres, print_matrix
		m = Munkres()
		indexes = m.compute(norm_matrix)
		indexes2 = [x[1] for x in indexes]

		new_W_i = [0]*len(W[samp])
		for i in range(len(indexes2)):
			new_W_i[indexes2[i]] = freqs[samp][i]
		new_W.append(new_W_i)

	W = np.array(new_W)

	return W


def nnls_with_prior_for_each_samp(X,W,Z,lambda1):

	W_new = []
	for i in range(len(X)):
		# print i
		# print W[i]
		W_i = nnls_with_prior(Z.T,X[i],W[i],lambda1)
		W_new.append(W_i)
		# print W_i

	W = np.array(W_new)

	W = use_all_components(X, W, Z)

	#scale W for each sample so that sum = 1
	for i in range(len(W)):
		W[i] = W[i]/sum(W[i])

	return W




class Difi2():
	"""
	Deconvolution Incorporating Frequency Information (Difi)
	A simple linear model with latent variables.
	"""


	def __init__(self, n_components=3, tol=1e-3, max_iter=200, rand_inits=1, lambda1=None):
		self.n_components = n_components
		self.tol = tol
		self.max_iter = max_iter
		self.rand_inits = rand_inits
		self.lambda1 = lambda1
		self.components_ = None
		self.norm = None
		self.W = None
		self.Z = None


	def fit(self, X, freqs):
		"""
		Fit the DIFI model to X

		Parameters
		----------
		X : array-like, shape (n_samples, n_features)
			Training data.

		Returns
		-------
		self
		"""

		to_plot = []

		best_norm = -1
		best_components = -1
		for inits in range(self.rand_inits):

			#Initialize Z
			Z = init_Z(X, self.n_components)

			#Optimize model with alternating NNLS
			W, Z, to_plot2 = optimize_z_and_w_using_nnls_given_X_and_Z(X, Z, self.max_iter, self.tol)

			to_plot.extend(to_plot2)

			to_beat = -1
			count = 0
			for i in range(20):
				print count
				#Fit prior to W
				W = fit_prior_freqs(W, freqs)
				to_plot.append(np.linalg.norm(X - np.dot(W, Z)))

				norm = np.linalg.norm(X - np.dot(W, Z))
				print norm
				# if to_beat - norm < self.tol and to_beat != -1:
				# 	break
				# to_beat = norm
				count += 1

				#Optimize Z
				Z = optimize_z_using_nnls(X, W)
				to_plot.append(np.linalg.norm(X - np.dot(W, Z)))

				#Optimize model with alternating NNLS
				W, Z, to_plot2 = optimize_z_and_w_using_nnls_given_X_and_Z(X, Z, self.max_iter, self.tol)
				to_plot.extend(to_plot2)


				

				




		#plot
		plt.plot(range(len(to_plot)), to_plot)
		plt.savefig("error.png")


		self.W = W
		self.Z = Z

		return self


	def transform(self, X, freqs):
		"""Solve for W given Z and X

		Parameters
		----------
		X : array-like, shape (n_samples, n_features)

		Returns
		-------
		X_new : array-like, shape (n_samples, n_components)
			The latent variables of X.
		"""

		print 'This needs to be implemented'

		return None
