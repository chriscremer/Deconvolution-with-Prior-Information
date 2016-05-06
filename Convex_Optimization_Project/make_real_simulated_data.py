



import numpy as np
import random
import math as m


import os,sys,inspect
home = os.path.expanduser('~')
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
# sys.path.insert(0,'..')
# sys.path.insert(0,'../simulated_data')
# sys.path.insert(0, home+'/plotting')


def make_subpopulations(n_subpops):
	'''
	Get real samples.
	Select the first x breast samples
	'''

	directory = home +'/TCGA_data/BRCA_rest/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/'

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



def make_samps(subpops, fractions, noise, n_outliers, n_subpops):

	X = np.dot(fractions, subpops)

	# print 'clean X'
	# print X[0][:4]


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
	gene_std = [np.std(x) for x in X.T]
	# print gene_std[:4]
	# gene_mean = [np.mean(x) for x in X.T]
	# gene_std = [1.]*len(X.T)
	for i in range(len(X)):
		
		# print amount
		# print np.random.normal(0,gene_variance)
		for j in range(len(X[i])):
			if gene_std[j] <= 0:
				continue


			# rand = random.random()
			# if rand > .5:
			# 	#Note, this normal function takes in std, whereas regualar N(m,s) notation, the s is variance
			# 	#so s in this case is sqrt(std)*noise
			# 	X[i][j] = X[i][j] + (np.random.normal(0,gene_std[j]*noise))
			# 	if X[i][j] < 0:
			# 		X[i][j] = 0.
			# else:
			# 	X[i][j] = X[i][j] + (np.random.normal(0,gene_std[j]*noise*2))
			# 	if X[i][j] < 0:
			# 		X[i][j] = 0.

			
			X[i][j] = X[i][j] + (np.random.normal(0,gene_std[j]*noise))

			if X[i][j] < 0:
				X[i][j] = 0.


	# print 'noise X'
	# print X[0][:4]

	# asdf

	#ADD OTHER BREAST SAMPLES, THESE ARE OUTLIERS
	if n_outliers > 0:
		directory = home +'/TCGA_data/BRCA_rest/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/'

		count = 0
		outliers = []
		#skip the number of ones that are already used for components
		for f in os.listdir(directory):
			if 'genes.normalized_results' in f:

				if count <= n_subpops:
					count +=1
					continue
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
				outliers.append(gene_expressions)
				count += 1
				if count == n_subpops+n_outliers+1:
					break

		outliers = np.array(outliers)

		for i in range(len(outliers)):

			X[i] = outliers[i]

		# X = np.concatenate((X,outliers), axis=0)


	

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



def make_and_return(n_subpops, n_samps, probability_of_zero, noise, n_outliers):

	# n_subpops = 10
	# n_samps = 10
	# probability_of_zero = .20
	# noise = .5

	subpops = make_subpopulations(n_subpops)
	fractions = subpop_fractions(n_samps, n_subpops, probability_of_zero)
	X = make_samps(subpops, fractions, noise, n_outliers, n_subpops)

	# print subpops.shape
	# print fractions.shape
	# print X.shape

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


