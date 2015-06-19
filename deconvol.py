


#add  directory to path to get my packages there
import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
#parentdir = os.path.dirname(currentdir)
#print sys.path
#sys.path.insert(0,parentdir) 
sys.path.insert(0,currentdir+'/simulated_data')

#for many things
import numpy as np

from numpy.linalg import inv

import make_convoluted_data

from sklearn.decomposition import PCA
from sklearn.decomposition import NMF





def read_data_folder():

	data = []

	directory = '/data1/morrislab/ccremer/TCGA_data/rnaseqv2_clinical_onlytumours_25nov2014/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/'

	count =0

	for file123 in os.listdir(directory):

		if 'genes.normalized_results' in file123:
			samp = []
			#print file123
			with open(directory + file123) as f2:
				firstLine = 1
				for line2 in f2:
					if firstLine == 1:
						firstLine = 0
						continue
					firstSplit = line2.split()
					samp.append(firstSplit[1])

			data.append(samp)

			count += 1
			if count > 40:
				break

	data = np.array(data)
	print 'real data shape ' + str(data.shape)

	return data


def get_possible_Ws(freqs):

	import itertools

	possible_Ws = []

	for i in range(len(freqs)):

		#find permuations and remove the duplicates
		perms = list(set(itertools.permutations(freqs[i])))

		#now find all combinations of sums of freqs (combine 2)
		combs = list(itertools.combinations(range(len(freqs[i])), 2))
		for j in range(len(combs)):
			#see if this combo has any zeros
			#if it does skip it
			skip = 0
			for index in combs[j]:
				if freqs[i][index] == 0.0:
					skip =1 
					break
			if skip == 1:
				continue
			#combine the freqencies and replace one with a zero
			new_freq = list(freqs[i])
			new_freq[combs[j][0]] = new_freq[combs[j][0]] + new_freq[combs[j][1]]
			new_freq[combs[j][1]] = 0.0
			new_perms = list(set(itertools.permutations(new_freq)))
			perms.extend(new_perms)
			
		#now find all combinations of sums of freqs (combine 3)
		combs = list(itertools.combinations(range(len(freqs[i])), 3))
		for j in range(len(combs)):
			#see if this combo has any zeros
			#if it does skip it
			skip = 0
			for index in combs[j]:
				if freqs[i][index] == 0.0:
					skip =1 
					break
			if skip == 1:
				continue
			#combine the freqencies and replace one with a zero
			new_freq = list(freqs[i])
			new_freq[combs[j][0]] = new_freq[combs[j][0]] + new_freq[combs[j][1]] + new_freq[combs[j][2]]
			new_freq[combs[j][1]] = 0.0
			new_freq[combs[j][2]] = 0.0
			new_perms = list(set(itertools.permutations(new_freq)))
			perms.extend(new_perms)

		#print len(combs)
		#print combs

		perms = list(set(perms))

		possible_Ws.append(perms)


def select_w(X, freqs, TZ):


	import itertools

	W = []

	for i in range(len(freqs)):

		#find permuations and remove the duplicates
		perms = list(set(itertools.permutations(freqs[i])))

		#now find all combinations of sums of freqs (combine 2)
		combs = list(itertools.combinations(range(len(freqs[i])), 2))
		for j in range(len(combs)):
			#see if this combo has any zeros
			#if it does skip it
			skip = 0
			for index in combs[j]:
				if freqs[i][index] == 0.0:
					skip =1 
					break
			if skip == 1:
				continue
			#combine the freqencies and replace one with a zero
			new_freq = list(freqs[i])
			new_freq[combs[j][0]] = new_freq[combs[j][0]] + new_freq[combs[j][1]]
			new_freq[combs[j][1]] = 0.0
			new_perms = list(set(itertools.permutations(new_freq)))
			perms.extend(new_perms)
		#now find all combinations of sums of freqs (combine 3)
		combs = list(itertools.combinations(range(len(freqs[i])), 3))
		for j in range(len(combs)):
			#see if this combo has any zeros
			#if it does skip it
			skip = 0
			for index in combs[j]:
				if freqs[i][index] == 0.0:
					skip =1 
					break
			if skip == 1:
				continue
			#combine the freqencies and replace one with a zero
			new_freq = list(freqs[i])
			new_freq[combs[j][0]] = new_freq[combs[j][0]] + new_freq[combs[j][1]] + new_freq[combs[j][2]]
			new_freq[combs[j][1]] = 0.0
			new_freq[combs[j][2]] = 0.0
			new_perms = list(set(itertools.permutations(new_freq)))
			perms.extend(new_perms)

		#print len(combs)
		#print combs

		perms = list(set(perms))

		#print len(perms)

		best_perm = []
		best_norm = -1

		for perm in perms:

			perm = np.array(perm)
			X_hat = np.dot(TZ.T, perm)
			norm = np.linalg.norm(X[i] - X_hat)

			if norm < best_norm or best_norm == -1:
				best_norm = norm
				best_perm = np.array(perm)

		W.append(best_perm)

	return np.reshape(np.array(W), (len(W),len(W[0]))) 




if __name__ == "__main__":


	real_data = read_data_folder()

	#params: #subpops, #feats, #samps
	numb_samps = len(real_data)
	numb_feats = len(real_data[0])
	numb_subpops = 3
	samps, freqs, subpops = make_convoluted_data.run_and_return(numb_subpops, numb_feats, numb_samps)

	X = samps

	print 'samps shape ' + str(samps.shape)
	print 'freqs shape ' + str(freqs.shape)
	print 'subpops shape ' + str(subpops.shape)

	numb_model_subpops= 4
	pca = PCA(n_components=numb_model_subpops)
	pca.fit(samps.T)
	#print(pca.explained_variance_ratio_)



	#########################################################
	#make all frequencies have same number of entries
	new_freqs = []
	for j in range(len(freqs)):
		freq = list(freqs[j])
		#if less than number of model, add zeros
		while len(freq) < numb_model_subpops:
			freq.append(0.0)
		if len(freq) > numb_model_subpops:
			sorted_freq = sorted(freq, reverse=True)
			while len(freq) > numb_model_subpops:
				freq.pop(len(freq) - 1)
			#scale so sums to 1
			freq_sum = sum(freq)
			for i in range(len(freq)):
				freq[i] = freq[i] / float(freq_sum)
		new_freqs.append(freq)

	freqs = np.array(new_freqs)
	print 'new_freqs shape ' + str(freqs.shape)
	#########################################################


	Z = pca.transform(samps.T).T
	print 'Z shape ' + str(Z.shape)

	T = np.identity(len(Z))
	print 'T shape ' + str(T.shape)

	TZ = np.dot(T, Z)
	print 'TZ shape ' + str(TZ.shape)

	X_hat = np.dot(freqs, TZ)
	print 'X_hat shape ' + str(X_hat.shape)
	norm = np.linalg.norm(samps - X_hat)
	print 'Initial norm ' + str(norm)


	possible_Ws = get_possible_Ws(freqs)

	for i in range(20):

		print 'Iter ' + str(i)

		print 'selecting W'
		W = select_w(X, freqs, TZ)
		X_hat = np.dot(W, TZ)
		norm = np.linalg.norm(X - X_hat)
		print '         Norm ' + str(norm)

		print 'optimizig TZ'
		TZ = np.dot(inv(np.dot(W.T,W)), np.dot(W.T, X))
		new_X_hat = np.dot(W, TZ)
		new_norm = np.linalg.norm(X - new_X_hat)
		print '         Norm ' + str(new_norm)

		if norm == new_norm:
			break


	#print W

	print 'DONE'




