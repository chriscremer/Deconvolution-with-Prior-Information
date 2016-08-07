

import numpy as np
import make_data as md

from scipy.optimize import nnls
import itertools

from os.path import expanduser
home = expanduser("~")


from nnls import nnls_with_prior

class regNMF():
	"""
	Deconvolution Incorporating Frequency Information (DIFI)
	"""


	def __init__(self, n_components, tol, max_iters):

		self.n_components = n_components
		self.tol = tol
		self.max_iters = max_iters


	def deconvolve(self, X, proportions):


		N = X.shape[0] #number of samps
		M = X.shape[1] #number of genes
		
		#Preprocess
		X_preprocessed = preprocess(X)

		#Initialize Z
		Z = init_Z(X_preprocessed, self.n_components)

		#Optimize model with alternating NNLS, where rows of W sum to 1
		W, Z = optimize_z_and_w_using_nnls_given_X_and_Z(X_preprocessed, Z, self.max_iters, self.tol)


		#Undo preprocessing - Transform Z back to full gene expressions
		Z = np.dot(W.T,X)

		return W,Z






#Functions 

def preprocess(X):

	#Remove low expression genes
	genes_mean_exps = [np.mean(x) for x in X.T]
	#largest to smallest
	sorted_indexes = np.argsort(genes_mean_exps)[::-1]
	#keep top 50%
	sorted_indexes = sorted_indexes[:len(sorted_indexes)/2]
	back_in_order = sorted(sorted_indexes)
	X = X.T[back_in_order].T


	#Scale by mean
	mean_list = [np.mean(x) for x in X.T]
	for gene in range(len(X.T)):
		X.T[gene] = X.T[gene]/mean_list[gene]

	#Remove bottom 25% of lowest variance genes
	#Instead keep top 1000 genes
	genes_var_exps = [np.std(x) for x in X.T]
	#largest to smallest
	sorted_indexes = np.argsort(genes_var_exps)[::-1]

	sorted_indexes = sorted_indexes[:1000]
	back_in_order = sorted(sorted_indexes)
	X = X.T[back_in_order].T

	return X




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



def optimize_z_and_w_using_nnls_given_X_and_Z(X, Z, max_iters, tol):

	print 'Optimizing...'
	converged = 0
	for iter1 in range(max_iters):

		#Optimize W
		W = []
		for i in range(len(X)):
			# W_i = nnls(Z.T, X[i])[0]
			W_i = nnls_with_prior(Z.T, X[i], np.array([0]*len(Z)), 1.)
			W.append(W_i)
		W = np.array(W)

		W = use_all_components(X, W, Z)

		#scale W for each sample so that sum = 1
		for i in range(len(W)):
			W[i] = W[i]/sum(W[i])

		norm = np.linalg.norm(X - np.dot(W, Z))
		# print 'norm ' + str(norm)

		#Optimize Z
		Z = optimize_z_using_nnls(X, W)

		new_norm = np.linalg.norm(X - np.dot(W, Z))
		# print 'new_norm ' + str(new_norm)

		if (norm - new_norm) < tol:
			print '# iters until optimized= ' + str(iter1) + ' last norm: ' + str(new_norm) 
			converged = 1
			break

	if converged == 0:
		print 'Did not converge1'

	return W, Z




def optimize_z_using_nnls(X, W):

	#Optimize Z
	Z = []
	for d in range(len(X.T)):
		Z_d = nnls(W, X.T[d])[0]
		Z.append(Z_d)
	Z = np.array(Z).T

	return Z




#For testing
if __name__ == "__main__":


	n_subpops = 10
	n_samps = 100
	probability_of_zero = .20
	noise = .05
	n_components = n_subpops
	tol = .001
	max_iters = 400

	X, subpops, fractions = md.make_and_return(n_subpops, n_samps, probability_of_zero, noise)

	print X.shape
	print subpops.shape
	print fractions.shape

	model = NMF(n_components=n_components, tol=tol, max_iters=max_iters)
	model.deconvolve(X, fractions)

	print 'All Done'






