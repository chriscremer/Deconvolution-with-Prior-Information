

import numpy as np
import make_data as md

from scipy.optimize import nnls
import itertools


from nnls import nnls_with_prior

class nsubpopsDIFI():
	"""
	Deconvolution Incorporating Frequency Information (DIFI)
	"""


	def __init__(self, n_components, tol, max_iters):

		self.n_components = n_components
		self.tol = tol
		self.max_iters = max_iters

	def __str__(self):
		return 'nsubpopsDIFI'

	def deconvolve(self, X, proportions):


		N = X.shape[0] #number of samps
		M = X.shape[1] #number of genes
		
		#Preprocess
		X_preprocessed = preprocess(X)

		#Initialize Z
		Z = init_Z(X_preprocessed, self.n_components)

		#Optimize model with alternating NNLS, where rows of W sum to 1
		W, Z = optimize_z_and_w_using_nnls_given_X_and_Z(X_preprocessed, Z, self.max_iters, self.tol)

		#Optimize model using the proportions
		W, Z = optimize_using_proportions(X_preprocessed, W, Z, proportions, self.max_iters, self.tol)

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


def optimize_using_proportions(X, W, Z, proportions, max_iters, tol):
	'''
	Assume there are up to 5 proportions for each sample

	'''

	print 'Optimizing...'
	converged = 0
	for iter1 in range(max_iters):

		#Fit proportions to nnls W, tries permuatations
		W = fit_n_proportions_to_W(X,W,Z,proportions)

		W = use_all_components(X, W, Z)

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
		print 'Did not converge2'

	return W, Z


def optimize_z_using_nnls(X, W):

	#Optimize Z
	Z = []
	for d in range(len(X.T)):
		Z_d = nnls(W, X.T[d])[0]
		Z.append(Z_d)
	Z = np.array(Z).T

	return Z

def fit_n_proportions_to_W(X,W,Z,proportions):

	#allow it to be less than n_subpops, but not more
	#so were assuming that same DNA = same RNA, but dif DNA != dif RNA
	new_W =[]
	for i in range(len(proportions)):
		# print W[i]
		# print proportions[i]

		n_subpops = len(proportions[i])
		#set W to best of using 1 to n_subpops
		best_norm = -1
		best_w = []
		best_n = -1
		# print W[i]
		for n in range(1,n_subpops+1):

			# print
			# print n

			#indexes of the top n_subpops weights
			top_indexes = np.argsort(W[i])[::-1][:n]
			# print top_indexes
			w = np.zeros(len(W[i]))
			for index_ in top_indexes:
				# w[index_] = W[i][index_]
				# w[index_] = -500.
				w[index_] = 10.
			# w = w + 500.
			# w = np.array(w)

			# print w

			W_i = nnls_with_prior(Z.T, X[i], w, 5000.)

			
			# print W_i
			# print w

			#sum to 1
			# print w
			# print sum(w)
			W_i = W_i / sum(W_i)

			# print W_i

			# print w
			#try it
			norm = np.linalg.norm(X[i] - np.dot(W_i, Z))
			# print norm

			if norm < best_norm or best_norm == -1:
				best_norm = norm
				best_w = W_i
				best_n = n

		

		new_W.append(best_w)
		# print 'actual n', n_subpops, 'chosen n', best_n 
		# print best_w

	# fsdfa


	return np.array(new_W)









# def fit_proportions_to_W(W,proportions):

# 	#so if there are 3 proportions, try matching it to the top 3  in W
# 	#then try all combos of 2, then match them to top 2 in W
# 	#then put them all together, and match to top in W
# 	#so there are 3, it checks 5 different assignments

# 	#if 1, = 1
# 	#if 2, = 2
# 	#if 3, = 5
# 	#if 4, 4C4 (1) , 4C3 (6), 2 (3), 4C1 (1) = 11
# 	#if 5, 5C5 (1) , 5C4 (10),  5C3 (15), 2 (10), 451 (1) = 37

# 	new_W =[]
# 	for i in range(len(proportions)):

# 		n_subpops = len(proportions[i])
# 		props = np.array(proportions[i])

# 		to_try = []
# 		#enumerate the possible combinatioms
# 		if n_subpops == 1:
# 			to_try.append(props) #Combos of 1
# 		if n_subpops == 2:
# 			to_try.append(props) #Combos of 2
# 			to_try.append([sum(props)]) #Combos of 1
# 		if n_subpops == 3:
# 			to_try.append(props) #Combos of 3
# 			#Combos of 2
# 			for j in range(n_subpops):
# 				group1 = props[j]
# 				group2_indexes = range(n_subpops)
# 				group2_indexes.pop(j)
# 				group2 = sum(props[group2_indexes])
# 				to_try.append(np.array([group1, group2])) 
# 			to_try.append(np.array([sum(props)])) #Combos of 1
# 		if n_subpops == 4:
# 			to_try.append(props) #Combos of 4
# 			#Combos of 2
# 			for j in range(n_subpops):
# 				group1 = props[j]
# 				group2_indexes = range(n_subpops)
# 				group2_indexes.pop(j)
# 				group2 = sum(props[group2_indexes])
# 				to_try.append(np.array([group1, group2]))


# 			pairs_of_2 = list(itertools.combinations(range(4), 2))
# 			# to_try.append(np.array([sum(proportions[i]])) #Combos of 1

# 		# if n_subpops == 5:
# 		# 	to_try.append(proportions[i]) #Combos of 5
# 		# 	to_try.append([sum(proportions[i]])) #Combos of 1


# 		#try them all, select the best

# 	return to_try


# def all_combos_of_2(props):
# 	#since len of props is max of 5, we only need to look at 2 indexes
# 	props = np.array(props)

# 	to_try = []

# 	n_subpops = len(props)

# 	#Leave one out


# 	#Combine only 2
# 	for i in range(n_subpops):
# 		for j in range(i+1,n_subpops):
# 				group1=props[i] + props[j]
# 				group2_indexes = range(n_subpops)
# 				group2_indexes.pop(i)
# 				group2_indexes.pop(j)
# 				group2 = sum(props[group2_indexes])
# 				to_try.append([group1, group2])





# 	# for i in range(n_subpops):

# 	# 	for j in range(i,n_subpops):

# 	# 		print i ,j

# 	# 		if i==j:
# 	# 			group1 = props[i]
# 	# 			print range(n_subpops)
# 	# 			group2_indexes = range(n_subpops)
# 	# 			del group2_indexes[i]
# 	# 			print group2_indexes
# 	# 			group2 = sum(props[group2_indexes])
# 	# 		else:
# 	# 			print props
# 	# 			group1=sum(props[i,j])
# 	# 			group2_indexes = range(n_subpops)
# 	# 			del group2_indexes[i]
# 	# 			del group2_indexes[j]
# 	# 			group2 = sum(props[group2_indexes])

# 	# 		to_try.append([group1, group2])
# 	# 		print to_try

# 	return to_try








	










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

	model = DIFI(n_components=n_components, tol=tol, max_iters=max_iters)
	model.deconvolve(X, fractions)

	print 'All Done'


	#MAkE SURE TO CHECKC THAT THE PREPROCESSING WORKS CPRRECTLY



	# print list(itertools.permutations([.3,.2,.1,.2,.2], 5))

	# fdsfa

	# props = [[.3,.2,.5]]
	# # props = [.3,.7]

	# print props
	# print fit_proportions_to_W([],props)


	# dfasdfa




