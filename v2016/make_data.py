
import numpy as np
import random
import os
import math as m

'''
Simulated Data Maker


Call run_and_return with params: ...


It will take k breast cancer samples and mix them based on their proportions

Sample noise: multiplicative noise. ie add exp(X) to each gene expression where 
				X is sampled from a Gaussian.


'''



def make_subpopulations(n_subpops):
	'''
	Get real samples.
	Select the first x breast samples
	'''

	directory = '../../../TCGA_data/BRCA_rest/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/'

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

	#noise is added based on the variance of the gene
	#noise 0 is nothing changed
	#noise 1 is change by a full standard deviation
	#noise between 0 and 1 is intermediate

	'''
	#type 1 noise
	gene_variance = [np.std(x) for x in X.T]
	# print X[0][:10]
	
	noise_matrix = []
	for samp in X:
		# samp_noise = [random.random()*noise*std*2. for std in gene_variance]
		samp_noise = [random.random()*noise*std for std in gene_variance]

		noise_matrix.append(samp_noise)
	noise_matrix = np.array(noise_matrix)
	X = X + noise_matrix
	'''

	'''
	additive noise
	'''
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
	
	

	#multiplicative noise
	# print np.log(0)
	# print np.std(np.log(0))
	# fafsd

	# std_log_gene = [np.std(np.log(x)) for x in X.T[:100]]
	#since i cant take log of zero

	# std_gene = [np.std(x) for x in X.T[:100]]
	# log_std_gene = []
	# for i in range(len(std_gene)):
	# 	if std_gene[i] == 0:
	# 		log_std_gene.append(0.)
	# 	else:
	# 		log_std_gene.append(np.log(std_gene[i]))


	'''
	#This is cibersort noise
	global_std = np.std(X)
	log_glocal_std = np.log(global_std)

	for i in range(len(X)):
		for j in range(len(X.T)):

			# print X[i][j]
			X[i][j] = X[i][j] + m.exp(np.random.normal(0, (noise*log_glocal_std)))
			# print X[i][j]
			# print
	'''


	return X



def make_and_return(n_subpops, n_samps, probability_of_zero, noise):

	# n_subpops = 10
	# n_samps = 10
	# probability_of_zero = .20
	# noise = .5

	subpops = make_subpopulations(n_subpops)
	fractions = subpop_fractions(n_samps, n_subpops, probability_of_zero)
	X = make_samps(subpops, fractions, noise)

	print subpops.shape
	print fractions.shape
	print X.shape

	return subpops, fractions, X



if __name__ == "__main__":

	n_subpops = 10
	n_samps = 10
	probability_of_zero = .20
	noise = .5

	subpops = make_subpopulations(n_subpops)
	fractions = subpop_fractions(n_samps, n_subpops, probability_of_zero)
	X = make_samps(subpops, fractions, noise)

	# print subpops.shape
	# print fractions.shape
	# print X.shape








