
"""
Deconvolution
A model for deconvolution of mixed samples.
"""

# Author: Chris Cremer <ccremer@cs.toronto.edu>

import numpy as np
import random
from scipy.optimize import nnls
import itertools
from numpy.linalg import pinv

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


def regression_with_prior(X, W, Z):

	#get the error of each samp using current W
	#dont change unless it becomes better
	W_errors = []
	for i in range(len(W)):
		W_errors.append(sum(abs(X[i] - np.dot(W[i], Z))))

	#allow W to deviate
	W_new = []
	lambda1 = 1.
	ZZT = np.dot(Z, Z.T)
	denominator = pinv(lambda1*np.identity(len(Z)) + ZZT)
	for samp in range(len(X)):
		W_i = np.dot(np.dot(X[samp], Z.T) + lambda1*W[samp], denominator)
		W_new.append(W_i)
	W_new = np.array(W_new)

	#remove negatives
	for i in range(len(W_new)):
		for j in range(len(W_new[i])):
			if W_new[i][j] < 0.0:
				W_new[i][j] = 0.0


	W_errors2 = []
	for i in range(len(W_new)):
		W_errors2.append(sum(abs(X[i] - np.dot(W_new[i], Z))))


	W_best = []
	for i in range(len(X)):
		if W_errors[i] < W_errors2[i]:
			W_best.append(W[i])
		else:
			W_best.append(W_new[i])

	W_best = np.array(W_best)

	return W_best


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

	return W


def optimize_z_using_nnls(X, W):

	#Optimize Z
	Z = []
	for d in range(len(X.T)):
		Z_d = nnls(W, X.T[d])[0]
		Z.append(Z_d)
	Z = np.array(Z).T

	return Z

class Deconvol():
	"""
	Deconvol
	A simple linear model with latent variables.
	"""


	def __init__(self, n_components=3, tol=1e-3, max_iter=200, rand_inits=1):
		self.n_components = n_components
		self.tol = tol
		self.max_iter = max_iter
		self.rand_inits = rand_inits
		self.components_ = None
		self.norm = None


	def fit(self, X):
		"""
		Fit the Deconvol model to X

		Parameters
		----------
		X : array-like, shape (n_samples, n_features)
			Training data.

		Returns
		-------
		self
		"""

		best_norm = -1
		lowest_W = -1
		best_components = -1
		for inits in range(self.rand_inits):
			# print 'Rand Init ' + str(inits)

			#Initialize Z
			Z = init_Z(X, self.n_components)

			#Optimize model
			converged = 0
			for iter1 in range(self.max_iter):

				#Optimize W
				W = []
				for i in range(len(X)):
					W_i = nnls(Z.T, X[i])[0]
					W.append(W_i)
				W = np.array(W)

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


				dif = X - np.dot(W, Z)
				norm = np.linalg.norm(dif)

				#Optimize Z
				Z = []
				for d in range(len(X.T)):
					Z_d = nnls(W, X.T[d])[0]
					Z.append(Z_d)
				Z = np.array(Z).T

				dif = X - np.dot(W, Z)
				new_norm = np.linalg.norm(dif)


				if (norm - new_norm) < self.tol:
					print '# iters until optimized= ' + str(iter1)
					converged = 1
					break

			if converged == 0:
				print 'Did not converge. Went to ' + str(self.max_iter) + ' iterations'

			if new_norm < best_norm or best_norm == -1:
				# W_sum = sum(sum(abs(W)))
				# if W_sum < lowest_W or lowest_W == -1:
				best_norm = new_norm
				best_components = Z
				# lowest_W = W_sum

		# self.W = W
		self.components_ = best_components
		self.norm = best_norm

		return self


	def transform(self, X):
		"""Solve for W given Z and X

		Parameters
		----------
		X : array-like, shape (n_samples, n_features)

		Returns
		-------
		X_new : array-like, shape (n_samples, n_components)
			The latent variables of X.
		"""

		W = []
		for i in range(len(X)):
			W_i = nnls(self.components_.T, X[i])[0]
			W.append(W_i)
		W = np.array(W)

		return W



class DIFI_strict():
	"""
	Deconvolution Including Frequency Information (DIFI)
	A simple linear model with latent variables.
	"""


	def __init__(self, n_components=3, tol=1e-3, max_iter=200, rand_inits=1):
		self.n_components = n_components
		self.tol = tol
		self.max_iter = max_iter
		self.rand_inits = rand_inits
		self.components_ = None
		self.norm = None


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

		best_norm = -1
		lowest_W = -1
		best_components = -1
		for inits in range(self.rand_inits):
			# print 'Rand Init ' + str(inits)

			#Initialize Z
			Z = init_Z(X, self.n_components)

			#Optimize model
			converged = 0
			for iter1 in range(self.max_iter):

				#Optimize W

				#Check all permutations: #profiles P #freqs
				#Instead
				#Check permutations of closest x profiles: (#freqs + x) P #freqs


				#For each sample
					#Get distance of each profile to sample
					#Try all permutations of the freqs on closest (#freqs+x) profiles 

				W = []
				for samp in range(len(X)):
					dist_to_profiles = []
					for p in range(len(Z)):
						dist_to_profiles.append(sum(abs(X[samp] - Z[p])))
					index_of_closest_profiles = np.argsort(dist_to_profiles)

					#find permuations and remove the duplicates
					perms = list(set(itertools.permutations(freqs[samp])))
					best_perm = -1
					best_norm = -1
					for perm in range(len(perms)):


						# print perms[perm]
						perm1 = np.array(perms[perm])

						# print perm1.shape
						# print Z.shape

						X_i_hat = np.dot(perm1, Z)

						# print X_i_hat.shape

						norm = sum(abs(X_i_hat - X[samp]))

						# print norm

						if norm < best_norm or best_norm == -1:
							best_norm = norm
							best_perm = perm1

					W.append(best_perm)
				W = np.array(W)

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

				#scale W for each sample so that sum = 1
				for i in range(len(W)):
					W[i] = W[i]/sum(W[i])


				dif = X - np.dot(W, Z)
				norm = np.linalg.norm(dif)

				#Optimize Z
				Z = []
				for d in range(len(X.T)):
					Z_d = nnls(W, X.T[d])[0]
					Z.append(Z_d)
				Z = np.array(Z).T

				dif = X - np.dot(W, Z)
				new_norm = np.linalg.norm(dif)


				if (norm - new_norm) < self.tol:
					print '# iters until optimized= ' + str(iter1)
					converged = 1
					break

			if converged == 0:
				print 'Did not converge. Went to ' + str(self.max_iter) + ' iterations'

			if new_norm < best_norm or best_norm == -1:
				# W_sum = sum(sum(abs(W)))
				# if W_sum < lowest_W or lowest_W == -1:
				best_norm = new_norm
				best_components = Z
				# lowest_W = W_sum

		# self.W = W
		self.components_ = best_components
		self.norm = best_norm

		return self


	def transform(self, X):
		"""Solve for W given Z and X

		Parameters
		----------
		X : array-like, shape (n_samples, n_features)

		Returns
		-------
		X_new : array-like, shape (n_samples, n_components)
			The latent variables of X.
		"""

		W = []
		for i in range(len(X)):
			W_i = nnls(self.components_.T, X[i])[0]
			W.append(W_i)
		W = np.array(W)

		return W



class DIFI_w_deviates():
	"""
	Deconvolution Including Frequency Information (DIFI)
	A simple linear model with latent variables.
	"""


	def __init__(self, n_components=3, tol=1e-3, max_iter=200, rand_inits=1):
		self.n_components = n_components
		self.tol = tol
		self.max_iter = max_iter
		self.rand_inits = rand_inits
		self.components_ = None
		self.norm = None


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

		best_norm = -1
		lowest_W = -1
		best_components = -1
		for inits in range(self.rand_inits):
			# print 'Rand Init ' + str(inits)

			#Initialize Z
			Z = init_Z(X, self.n_components)

			#Optimize model
			converged = 0
			norm = -1
			new_norm = -1
			for iter1 in range(self.max_iter):

				#Optimize W

				#Check all permutations: #profiles P #freqs
				#Instead
				#Check permutations of closest x profiles: (#freqs + x) P #freqs


				#For each sample
					#Get distance of each profile to sample
					#Try all permutations of the freqs on closest (#freqs+x) profiles 

				W = []
				for samp in range(len(X)):
					dist_to_profiles = []
					for p in range(len(Z)):
						dist_to_profiles.append(sum(abs(X[samp] - Z[p])))
					index_of_closest_profiles = np.argsort(dist_to_profiles)

					#find permuations and remove the duplicates
					perms = list(set(itertools.permutations(freqs[samp])))
					best_perm = -1
					best_lin_error = -1
					for perm in range(len(perms)):


						# print perms[perm]
						perm1 = np.array(perms[perm])

						# print perm1.shape
						# print Z.shape

						X_i_hat = np.dot(perm1, Z)

						# print X_i_hat.shape

						linear_error = sum(abs(X_i_hat - X[samp]))

						# print norm

						if linear_error < best_lin_error or best_lin_error == -1:
							best_lin_error = linear_error
							best_perm = perm1

					W.append(best_perm)
				W = np.array(W)

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


				#allow W to deviate
				W = []
				lambda1 = 1.
				ZZT = np.dot(Z, Z.T)
				denominator = pinv(lambda1*np.identity(len(Z)) + ZZT)
				for samp in range(len(X)):
					W_i = np.dot(np.dot(X[samp], Z.T), denominator)
					W.append(W_i)
				W = np.array(W)

				#remove negatives
				for i in range(len(W)):
					for j in range(len(W[i])):
						if W[i][j] < 0.0:
							W[i][j] == 0.0


				#scale W for each sample so that sum = 1
				for i in range(len(W)):
					W[i] = W[i]/sum(W[i])

				dif = X - np.dot(W, Z)
				norm = np.linalg.norm(dif)

				# print '\nNorm ' + str(norm)

				if norm > new_norm and new_norm != -1:
					print '# iters until optimized= ' + str(iter1)
					converged = 1
					break

				#Optimize Z
				Z = []
				for d in range(len(X.T)):
					Z_d = nnls(W, X.T[d])[0]
					Z.append(Z_d)
				Z = np.array(Z).T

				dif = X - np.dot(W, Z)
				new_norm = np.linalg.norm(dif)

				# print 'New Norm ' + str(new_norm)


				if (norm - new_norm) < self.tol:
					print '# iters until optimized= ' + str(iter1)
					converged = 1
					break

			if converged == 0:
				print 'Did not converge. Went to ' + str(self.max_iter) + ' iterations'

			if new_norm < best_norm or best_norm == -1:
				# W_sum = sum(sum(abs(W)))
				# if W_sum < lowest_W or lowest_W == -1:
				best_norm = new_norm
				best_components = Z
				# lowest_W = W_sum

		# self.W = W
		self.components_ = best_components
		self.norm = best_norm

		return self


	def transform(self, X):
		"""Solve for W given Z and X

		Parameters
		----------
		X : array-like, shape (n_samples, n_features)

		Returns
		-------
		X_new : array-like, shape (n_samples, n_components)
			The latent variables of X.
		"""

		W = []
		for i in range(len(X)):
			W_i = nnls(self.components_.T, X[i])[0]
			W.append(W_i)
		W = np.array(W)

		return W



class Deconvol_normalized():
	"""
	Deconvol
	A simple linear model with latent variables.
	"""


	def __init__(self, n_components=3, tol=1e-3, max_iter=200, rand_inits=1):
		self.n_components = n_components
		self.tol = tol
		self.max_iter = max_iter
		self.rand_inits = rand_inits
		self.components_ = None
		self.norm = None


	def fit(self, X):
		"""
		Fit the Deconvol model to X

		Parameters
		----------
		X : array-like, shape (n_samples, n_features)
			Training data.

		Returns
		-------
		self
		"""

		best_norm = -1
		lowest_W = -1
		best_components = -1
		for inits in range(self.rand_inits):
			# print 'Rand Init ' + str(inits)

			#Initialize Z
			Z = init_Z(X, self.n_components)

			#Optimize model
			converged = 0
			for iter1 in range(self.max_iter):

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

				dif = X - np.dot(W, Z)
				norm = np.linalg.norm(dif)

				#Optimize Z
				Z = []
				for d in range(len(X.T)):
					Z_d = nnls(W, X.T[d])[0]
					Z.append(Z_d)
				Z = np.array(Z).T

				dif = X - np.dot(W, Z)
				new_norm = np.linalg.norm(dif)


				if (norm - new_norm) < self.tol:
					print '# iters until optimized= ' + str(iter1)
					converged = 1
					break

			if converged == 0:
				print 'Did not converge. Went to ' + str(self.max_iter) + ' iterations'

			if new_norm < best_norm or best_norm == -1:
				# W_sum = sum(sum(abs(W)))
				# if W_sum < lowest_W or lowest_W == -1:
				best_norm = new_norm
				best_components = Z
				# lowest_W = W_sum

		# self.W = W
		self.components_ = best_components
		self.norm = best_norm

		return self


	def transform(self, X):
		"""Solve for W given Z and X

		Parameters
		----------
		X : array-like, shape (n_samples, n_features)

		Returns
		-------
		X_new : array-like, shape (n_samples, n_components)
			The latent variables of X.
		"""

		W = []
		for i in range(len(X)):
			W_i = nnls(self.components_.T, X[i])[0]
			W.append(W_i)
		W = np.array(W)

		return W



class DIFI_strict_top5():
	"""
	Deconvolution Incorporating Frequency Information (Difi)
	A simple linear model with latent variables.
	"""


	def __init__(self, n_components=3, tol=1e-3, max_iter=200, rand_inits=1):
		self.n_components = n_components
		self.tol = tol
		self.max_iter = max_iter
		self.rand_inits = rand_inits
		self.components_ = None
		self.norm = None


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

		best_norm = -1
		lowest_W = -1
		best_components = -1
		for inits in range(self.rand_inits):
			# print 'Rand Init ' + str(inits)

			#Initialize Z
			Z = init_Z(X, self.n_components)

			#Optimize model
			converged = 0
			norm = -1
			new_norm = -1
			for iter1 in range(self.max_iter):

				#Optimize W

				#Check all permutations: #profiles P #freqs
				#Instead
				#Check permutations of closest x profiles: (#freqs + x) P #freqs


				#For each sample
					#Get distance of each profile to sample
					#Try all permutations of the freqs on closest (#freqs+x) profiles 

				W = assign_W_to_top_5_components(X, freqs, Z)

				W = use_all_components(X, W, Z)

				#scale W for each sample so that sum = 1
				for i in range(len(W)):
					W[i] = W[i]/sum(W[i])

				dif = X - np.dot(W, Z)
				norm = np.linalg.norm(dif)
				# norm = sum(sum(abs(dif)))

				print 'norm ' + str(norm)


				# if norm > new_norm and new_norm != -1:
				# 	print '# iters until optimized= ' + str(iter1)
				# 	converged = 1
				# 	break


				Z = optimize_z_using_nnls(X, W)

				dif = X - np.dot(W, Z)
				new_norm = np.linalg.norm(dif)
				# new_norm = sum(sum(abs(dif)))

				print 'new_norm ' + str(new_norm)

				if (norm - new_norm) < self.tol:
					print '# iters until optimized= ' + str(iter1)
					converged = 1
					break

			if converged == 0:
				print 'Did not converge. Went to ' + str(self.max_iter) + ' iterations'

			if new_norm < best_norm or best_norm == -1:
				# W_sum = sum(sum(abs(W)))
				# if W_sum < lowest_W or lowest_W == -1:
				best_norm = new_norm
				best_components = Z
				# lowest_W = W_sum

		# self.W = W
		self.components_ = best_components
		self.norm = best_norm

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

		Z = self.components_
		W = assign_W_to_top_5_components(X, freqs, Z)
		W = use_all_components(X, W, Z)
		#scale W for each sample so that sum = 1
		for i in range(len(W)):
			W[i] = W[i]/sum(W[i])

		return W



class DIFI_deviate_top5():
	"""
	Deconvolution Incorporating Frequency Information (Difi)
	A simple linear model with latent variables.
	"""


	def __init__(self, n_components=3, tol=1e-3, max_iter=200, rand_inits=1):
		self.n_components = n_components
		self.tol = tol
		self.max_iter = max_iter
		self.rand_inits = rand_inits
		self.components_ = None
		self.norm = None


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

		best_norm = -1
		best_components = -1
		for inits in range(self.rand_inits):
			# print 'Rand Init ' + str(inits)

			#Initialize Z
			Z = init_Z(X, self.n_components)

			#Optimize model
			converged = 0
			norm = -1
			new_norm = -1
			for iter1 in range(self.max_iter):

				#Optimize W

				#Check all permutations: #profiles P #freqs
				#Instead
				#Check permutations of closest x profiles: (#freqs + x) P #freqs


				#For each sample
					#Get distance of each profile to sample
					#Try all permutations of the freqs on closest (#freqs+x) profiles 


				W = assign_W_to_top_5_components(X, freqs, Z)


				W = use_all_components(X, W, Z)


				W = regression_with_prior(X, W, Z)

				#scale W for each sample so that sum = 1
				for i in range(len(W)):
					W[i] = W[i]/sum(W[i])

				dif = X - np.dot(W, Z)
				norm = np.linalg.norm(dif)

				# print 'norm ' + str(norm)

				if norm > new_norm and new_norm != -1:
					print '# stoped after ' + str(iter1)
					converged = 1
					break

				#Optimize Z
				Z = []
				for d in range(len(X.T)):
					Z_d = nnls(W, X.T[d])[0]
					Z.append(Z_d)
				Z = np.array(Z).T

				dif = X - np.dot(W, Z)
				new_norm = np.linalg.norm(dif)

				# print 'new_norm ' + str(new_norm)


				if (norm - new_norm) < self.tol:
					print '# iters until optimized= ' + str(iter1)
					converged = 1
					break

			if converged == 0:
				print 'Did not converge. Went to ' + str(self.max_iter) + ' iterations'

			if new_norm < best_norm or best_norm == -1:
				# W_sum = sum(sum(abs(W)))
				# if W_sum < lowest_W or lowest_W == -1:
				best_norm = new_norm
				best_components = Z
				# lowest_W = W_sum

		# self.W = W
		self.components_ = best_components
		self.norm = best_norm

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

		Z = self.components_
		W = assign_W_to_top_5_components(X, freqs, Z)
		W = use_all_components(X, W, Z)
		W = regression_with_prior(X, W, Z)
		#scale W for each sample so that sum = 1
		for i in range(len(W)):
			W[i] = W[i]/sum(W[i])

		return W





class DIFI_nmf_deviate_top5():
	"""
	Deconvolution Incorporating Frequency Information (Difi)
	A simple linear model with latent variables.
	"""


	def __init__(self, n_components=3, tol=1e-3, max_iter=200, rand_inits=1):
		self.n_components = n_components
		self.tol = tol
		self.max_iter = max_iter
		self.rand_inits = rand_inits
		self.components_ = None
		self.norm = None


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

		best_norm = -1
		best_components = -1
		for inits in range(self.rand_inits):
			# print 'Rand Init ' + str(inits)

			#Initialize Z
			Z = init_Z(X, self.n_components)

			print 'Step 1'
			#Optimize model with alternating NNLS
			#This is to get a good Z
			converged = 0
			for iter1 in range(self.max_iter):

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

				if (norm - new_norm) < self.tol:
					# print '# iters until optimized= ' + str(iter1)
					converged = 1
					break

			print 'Step 2'
			#Optimize model using the known frequencies
			#Can start Z where the NMF ended
			converged = 0
			norm = -1
			new_norm = -1
			for iter1 in range(self.max_iter):

				#Optimize W
				W = assign_W_to_top_5_components(X, freqs, Z)
				W = use_all_components(X, W, Z)
				W = regression_with_prior(X, W, Z)

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


				if (norm - new_norm) < self.tol:
					print '# iters until optimized= ' + str(iter1)
					converged = 1
					break

			if converged == 0:
				print 'Did not converge. Went to ' + str(self.max_iter) + ' iterations'

			if new_norm < best_norm or best_norm == -1:
				best_norm = new_norm
				best_components = Z


		# self.W = W
		self.components_ = best_components
		self.norm = best_norm

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

		Z = self.components_
		W = assign_W_to_top_5_components(X, freqs, Z)
		W = use_all_components(X, W, Z)
		W = regression_with_prior(X, W, Z)
		#scale W for each sample so that sum = 1
		for i in range(len(W)):
			W[i] = W[i]/sum(W[i])

		return W


