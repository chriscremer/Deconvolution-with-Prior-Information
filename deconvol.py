


#add  directory to path to get my packages there
import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
#parentdir = os.path.dirname(currentdir)
#print sys.path
#sys.path.insert(0,parentdir) 
sys.path.insert(0,currentdir+'/simulated_data')

import numpy as np

import itertools

from numpy.linalg import inv

import make_convoluted_data

from sklearn.decomposition import PCA
from sklearn.decomposition import NMF
from sklearn.decomposition import fastica

import multiprocessing as mp



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
			#if count > 40:
			#	break

	data = np.array(data)
	print 'real data shape ' + str(data.shape)

	return data


def get_possible_Ws(freqs):

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

	return possible_Ws


def select_w(X, possible_Ws, TZ):

	W = []

	for i in range(len(possible_Ws)):

		best_perm = []
		best_norm = -1

		for perm in possible_Ws[i]:

			perm1 = np.array(perm)
			X_hat = np.dot(TZ.T, perm1)
			norm = np.linalg.norm(X[i] - X_hat)

			if norm < best_norm or best_norm == -1:
				best_norm = norm
				best_perm = np.array(perm)

		W.append(best_perm)

	return np.reshape(np.array(W), (len(W),len(W[0]))) 



def doWork(samp_index):

	#print 'Doing sample ' + str(samp_index)

	best_perm = []
	best_norm = -1
	for perm in possible_Ws[samp_index]:

		perm1 = np.array(perm)
		X_hat = np.dot(TZ.T, perm1)
		norm = np.linalg.norm(X[samp_index] - X_hat)

		if norm < best_norm or best_norm == -1:
			best_norm = norm
			best_perm = np.array(perm)

	return best_perm

	

def select_w_parallel():

	try:

		samp_indexes = [i for i in range(len(possible_Ws))]
		numb_cpus = mp.cpu_count() - 4
		pool = mp.Pool(processes=numb_cpus)

		W = pool.map(doWork, samp_indexes)


		#print 'W ' + str(len(W)) + ' ' + str(len(W[0]))

	except KeyboardInterrupt:
		print 'keyboard interruption'

	#print np.array(W).shape

	#return np.reshape(np.array(W), (len(W)))

	return np.array(W)


if __name__ == "__main__":

	#real_data = read_data_folder()

	#params: #subpops, #feats, #samps
	#numb_samps = len(real_data)
	#numb_feats = len(real_data[0])
	numb_samps = 50
	numb_feats = 10000
	numb_subpops = 5
	samps, freqs, subpops = make_convoluted_data.run_and_return(numb_subpops, numb_feats, numb_samps)

	global X
	X = samps

	print 'samps shape ' + str(samps.shape)
	print 'freqs shape ' + str(freqs.shape)
	print 'subpops shape ' + str(subpops.shape)

	numb_model_subpops= 5
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
		#if  more than number of model,then take largest
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

	print 'Finding all weight permutations..'
	global possible_Ws
	possible_Ws = get_possible_Ws(freqs)
	print 'Completed.'


	print 'Optimizing model..'

	for i in range(20):

		print 'Iter ' + str(i)

		print 'selecting W'
		# W = select_w(X, possible_Ws, TZ)
		W = select_w_parallel()
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

	print 
	#figure out how close the model components are to the actual profiles
	#and match the components to their profiles so printing makes sense
	# so for each actual profile, find the row of TZ that is most similar to it

	'''
	#this list keeps track of which profile correspond to which compoents
	#also this list is so that I dont assign a profiles to the same components
	component_order = []
	for profile_index in range(len(subpops)):

		#what happens when # subpops != #components??
		best_norm = -1
		best_component = -1

		for component_index in range(len(TZ)):

			if component_index in component_order:
				continue

			norm = np.linalg.norm(subpops[profile_index] - TZ[component_index])
			if norm < best_norm or best_norm == -1:
				best_norm = norm
				best_component = component_index


		component_order.append(best_component)

		#print how close each component is to each matching profile
		print 'Profile ' + str(profile_index) + ' norm ' + str(best_norm)
	'''

	#that is one way of doing it but I dont think its optimal
	#a better way is too find all permutations and pick the one that minizes the sum of norms
	possible_component_order = list(set(itertools.permutations(range(len(TZ)))))
	best_norm_sum = -1
	best_order = -1
	for i in range(len(possible_component_order)):

		norm_sum = 0
		for profile_index in range(len(subpops)):

			norm_sum += np.linalg.norm(subpops[profile_index] - TZ[possible_component_order[i][profile_index]])

		if norm_sum < best_norm_sum or best_norm_sum == -1:
			best_norm_sum = norm_sum
			best_order = i

	component_order = possible_component_order[best_order]

	print 
	#print how close each component is to each matching profile
	for profile_index in range(len(subpops)):

		print 'Profile ' + str(profile_index) + ' norm ' + str(np.linalg.norm(subpops[profile_index] - TZ[component_order[profile_index]]))




	#print actual frequencies of samples and predicted assignment of frequencies
	print 'Actual - Predicted'
	for i in range(len(freqs)):

		#print str(freqs[i]) + '  ' + str(W[i])
		#print str(freqs[i]) + '  ' + str([W[i][x] for x in component_order])

		#print("%.2f" % a)
		#myFormattedList = [ '%.2f' % elem for elem in myList ]

		print str(['%.2f' % elem for elem in freqs[i]]) + '  ' + str(['%.2f' % elem for elem in [W[i][x] for x in component_order]])

	#print freqs

	#print W

	print 'DONE'




