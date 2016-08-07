
import numpy as np
import random
import os
from os.path import expanduser
home = expanduser("~")
import math as m

'''
Simulated Data Maker


Call run_and_return with params:
	n_subpops, 
	n_samps, 
	probability_of_zero, 
	noise


It will take k breast cancer samples and mix them based on their proportions

Sample noise: multiplicative noise. ie add exp(X) to each gene expression where 
				X is sampled from a Gaussian.


'''



def make_subpopulations(n_subpops):
	'''
	Get real samples.
	Select the first x breast samples
	'''

	directory = home + '/TCGA_data/BRCA_rest/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/'

	count = 0
	subpops = []
	for f in os.listdir(directory):
		if 'genes.normalized_results' in f:
			# print(f)
			
			gene_expressions = []
			with open(directory + f) as f2:
				b = 0
				for line in f2:
					if b == 0:
						b+=1
						continue
					firstSplit = line.split()
					secondSplit = firstSplit[0].split('|')

					if secondSplit[0] == '?':
						continue

					gene_expressions.append(float(firstSplit[1]))
			subpops.append(gene_expressions)
			count += 1
			if count == n_subpops:
				break

	subpops = np.array(subpops)
	return subpops




def subpop_fractions(n_samps, n_subpops ,probability_of_zero):
	'''
	Define the frequency of each subpop in each sample
	'''

	freqs = []
	for i in range(n_samps):
		this_samp = []
		for j in range(n_subpops):
			rand = random.random()
			if rand < probability_of_zero:
				this_samp.append(0.0)
			else:
				this_samp.append(random.random())
		freqs.append(this_samp)
	freqs = np.array(freqs)

	#each vector should sum to 1
	for i in range(len(freqs)):
		sum_of_row = sum(freqs[i])
		if sum_of_row == 0.0:
			#avoid sample with all zeros
			a = random.randint(0, n_subpops-1)
			freqs[i][a] = 1.0
			continue
		#divide row by sum
		freqs[i] = freqs[i] / sum_of_row

	freqs = np.array(freqs)
	return freqs



def make_samps(subpops, fractions, noise):


	X = np.dot(fractions, subpops)






	'''
	# additive noise

	gene_variance = [np.std(x) for x in X.T]
	for i in range(len(X)):
		
		# print amount
		# print np.random.normal(0,gene_variance)
		for j in range(len(X[i])):
			if gene_variance[j] <= 0:
				continue
			X[i][j] = X[i][j] + (np.random.normal(0,gene_variance[j])*noise)
			if X[i][j] < 0:
				X[i][j] = 0.
	'''
	
	


	
	#This is multiplicative noise
	#cibersort does this
	global_std = np.std(X)
	log_glocal_std = np.log(global_std)

	for i in range(len(X)):
		for j in range(len(X.T)):

			# print X[i][j]
			X[i][j] = X[i][j] + m.exp(np.random.normal(0, (noise*log_glocal_std)))
			# print X[i][j]
			# print
	


	return X



def make_and_return(n_subpops, n_samps, probability_of_zero, noise):

	subpops = make_subpopulations(n_subpops)
	fractions = subpop_fractions(n_samps, n_subpops, probability_of_zero)
	X = make_samps(subpops, fractions, noise)

	return X, subpops, fractions


def split_fractions(fractions):
	'''
	Fractions will be a list of lists
	The lists are the fractions of each sample without the zeros
	'''

	#Allow fractions to come from same subpops
	#ie. dif DNA, but same RNA

	new_fractions = []
	for i in range(len(fractions)):

		this_samp = []

		#theres a 30% chance any fraction will be split in half

		for j in range(len(fractions[i])):

			rand = random.random()
			if rand < .3:

				half = fractions[i][j] / 2.
				this_samp.append(half)
				this_samp.append(half)

			else:
				this_samp.append(fractions[i][j])

		new_fractions.append(this_samp)

	return new_fractions








if __name__ == "__main__":

	n_subpops = 10
	n_samps = 10
	probability_of_zero = .20
	noise = .5

	X, subpops, fractions = make_and_return(n_subpops, n_samps, probability_of_zero, noise)

	print X.shape
	print subpops.shape
	print fractions.shape








