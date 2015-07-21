


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
from numpy.linalg import pinv


import make_convoluted_data

from sklearn.decomposition import PCA
from sklearn.decomposition import NMF
from sklearn.decomposition import FastICA

import multiprocessing as mp

#to plot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


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

def plot_bar_chart(types, means, stds, xaxis_labels):

	# the x locations for the groups
	ind = np.arange(len(xaxis_labels))
	# the width of the bars
	width = 0.20

	fig, ax = plt.subplots()
	random = ax.bar(ind-width, means[0], width, color='purple', yerr=stds[0])
	pca = ax.bar(ind, means[1], width, color='red', yerr=stds[1])
	#ica = ax.bar(ind+width, means[2], width, color='blue', yerr=stds[2])



	# add some text for labels, title and axes ticks
	ax.set_ylabel('Error')
	ax.set_xlabel('Number of Components')
	#ax.set_title('Scores by group and gender')
	ax.set_xticks(ind+width)
	ax.set_xticklabels(xaxis_labels)
	#ax.legend( (random[0], pca[0], ica[0]), ('Random', 'PCA', 'ICA') )
	ax.legend( (random[0], pca[0]), ('Random', 'PCA') )


	plt.savefig('bar_plot_different_initializations.pdf')
	print 'Saved plot'

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

		#samp_indexes = [i for i in range(len(possible_Ws))]
		samp_indexes = range(len(possible_Ws))
		numb_cpus = mp.cpu_count()
		#print 'numb of cpus' + str(numb_cpus)
		pool = mp.Pool()

		W = pool.map(doWork, samp_indexes)


		#print 'W ' + str(len(W)) + ' ' + str(len(W[0]))

	except KeyboardInterrupt:
		print 'keyboard interruption'

	#print np.array(W).shape

	#return np.reshape(np.array(W), (len(W)))

	return np.array(W)


if __name__ == "__main__":

	rand_means = []
	pca_means = []
	ica_means = []
	rand_stds = []
	pca_stds = []
	ica_stds = []

	for numb_components in range(2, 7):

		print 'Numb of Components =' + str(numb_components)

		means_for_these_components = []
		stds_for_these_components = []

		#########################################################
		#Define data and subpops

		#real_data = read_data_folder()

		#numb_samps = len(real_data)
		#numb_feats = len(real_data[0])
		numb_samps = 50
		numb_feats = 10000

		#this is the number of profiles that exist in the simulated data
		numb_subpops = numb_components
		#this is the number of components that are kept in pca
		numb_model_subpops = numb_components
		#########################################################



		#########################################################
		#Evaluating the different Z initializations
		# 1) Random 2) PCA 3) ICA
		#types = ['Random', 'PCA', 'ICA']
		types = ['Random', 'PCA']

		profile_norm_scores = []
		profile_norm_stds = []
		freqs_norm_scores = []

		#for each initialization type
		for init_type in range(len(types)):

			print 'Testing ' + str(types[init_type])

			profile_norm_store = []
			freqs_norm_store = []

			#iterate x times and get average
			for iteration in range(10):

				print 'Iter ' + str(iteration)

				#########################################################
				#Make data
				#params: #subpops, #feats, #samps
				samps, freqs, subpops = make_convoluted_data.run_and_return(numb_subpops, numb_feats, numb_samps)

				global X
				X = samps

				#print 'samps shape ' + str(samps.shape)
				#print 'freqs shape ' + str(freqs.shape)
				#print 'subpops shape ' + str(subpops.shape)
				#########################################################

				#########################################################
				#make all frequencies have same number of entries
				new_freqs = []
				#for each sample
				for j in range(len(freqs)):
					freq = list(freqs[j])
					#if less than number of model, add zeros
					while len(freq) < numb_model_subpops:
						freq.append(0.0)
					#if  more than number of model,then take largest
					if len(freq) > numb_model_subpops:
						freq = sorted(freq, reverse=True)
						while len(freq) > numb_model_subpops:
							freq.pop(len(freq) - 1)
					#scale so sums to 1
					freq_sum = sum(freq)
					if freq_sum != 1.0:
						for i in range(len(freq)):
							freq[i] = freq[i] / float(freq_sum)
					np.random.shuffle(freq)
					new_freqs.append(freq)

				new_freqs = np.array(new_freqs)
				#print 'new_freqs shape ' + str(new_freqs.shape)
				#########################################################

				#########################################################
				#make all frequencies have same number of entries, WITHOUT SHUFFLING, used for comparing at the end
				start_freqs = []
				#for each sample
				for j in range(len(freqs)):
					freq = list(freqs[j])
					#if less than number of model, add zeros
					while len(freq) < numb_model_subpops:
						freq.append(0.0)
					#if  more than number of model,then take largest
					if len(freq) > numb_model_subpops:
						freq = sorted(freq, reverse=True)
						while len(freq) > numb_model_subpops:
							freq.pop(len(freq) - 1)
					#scale so sums to 1
					freq_sum = sum(freq)
					if freq_sum != 1.0:
						for i in range(len(freq)):
							freq[i] = freq[i] / float(freq_sum)
					#np.random.shuffle(freq)
					start_freqs.append(freq)

				start_freqs = np.array(start_freqs)
				#print 'start_freqs shape ' + str(start_freqs.shape)
				#########################################################


				#########################################################
				#Initializing model

				if init_type == 0:
					Z = np.random.rand(numb_model_subpops, len(samps[0]))
				elif init_type == 1:
					pca = PCA(n_components=numb_model_subpops)
					pca.fit(samps.T)
					Z = pca.transform(samps.T).T
				else:
					pca = FastICA(n_components=numb_model_subpops)
					pca.fit(samps.T)
					Z = pca.transform(samps.T).T		

				#print 'Z shape ' + str(Z.shape)

				T = np.identity(len(Z))
				#print 'T shape ' + str(T.shape)

				TZ = np.dot(T, Z)
				#print 'TZ shape ' + str(TZ.shape)

				X_hat = np.dot(new_freqs, TZ)
				#print 'X_hat shape ' + str(X_hat.shape)
				norm = np.linalg.norm(samps - X_hat)
				#print 'Initial norm ' + str(norm)

				#print 'Finding all weight permutations..'
				global possible_Ws
				possible_Ws = get_possible_Ws(new_freqs)
				#print 'Completed.'
				#########################################################


				#########################################################
				print 'Optimizing model..'

				for i in range(50):

					#print 'Iter ' + str(i)

					#print 'selecting W'
					W = select_w(X, possible_Ws, TZ)
					#W = select_w_parallel()
					X_hat = np.dot(W, TZ)
					norm = np.linalg.norm(X - X_hat)
					#print '         	Norm ' + str(norm)


					#print 'optimizig TZ'
					TZ = np.dot(pinv(np.dot(W.T,W)), np.dot(W.T, X))
					new_X_hat = np.dot(W, TZ)
					new_norm = np.linalg.norm(X - new_X_hat)
					#print '         	Norm ' + str(new_norm)

					if norm == new_norm:
						break
				#########################################################

				#print

				######################################################### 
				# match the components to their profiles so printing makes sense
				# so for each actual profile, find the row of TZ that is most similar to it
				possible_component_order = list(set(itertools.permutations(range(len(TZ)))))
				best_norm_sum = -1
				best_order = -1
				#for each possible ordering of the components
				for i in range(len(possible_component_order)):
					#keep track of the sum of the norms
					norm_sum = 0
					#for each profile 
					for profile_index in range(len(subpops)):
						#add to the norm sum
						#the norm of the difference betweem the profile and the corresponding component given this order
						norm_sum += np.linalg.norm(subpops[profile_index] - TZ[possible_component_order[i][profile_index]])

					if norm_sum < best_norm_sum or best_norm_sum == -1:
						best_norm_sum = norm_sum
						best_order = i

				component_order = possible_component_order[best_order]
				#########################################################

				#print

				#########################################################
				#print how close (norm) each component is to each matching profile
				sum1=0
				for profile_index in range(len(subpops)):
					norm1 = np.linalg.norm(subpops[profile_index] - TZ[component_order[profile_index]])
					#print 'Profile ' + str(profile_index) + ' norm ' + str(norm1)
					sum1 += norm1
				print 'Sum profile norm ' + str(sum1)
				profile_norm_store.append(sum1)
				#########################################################

				#print 

				#########################################################
				#this prints the actual frequencies of samples and predicted assignment of frequencies
				# print 'Actual - Predicted'
				# for i in range(len(freqs)):
				# 	#only print first 10
				# 	if i > 9:
				# 		break
				# 	print str(['%.2f' % elem for elem in start_freqs[i]]) + '  ' + str(['%.2f' % elem for elem in [W[i][x] for x in component_order]])
				#########################################################


				#########################################################
				#print average norm of assigned frequencies vs actual
				sum1 = 0
				for i in range(len(freqs)):
					dif_array = start_freqs[i] - [W[i][x] for x in component_order]
					sum1 += np.linalg.norm(dif_array)
				print 'Sum freq assignemtn norm ' + str(sum1)
				freqs_norm_store.append(sum1)
				#########################################################

			profile_norm_scores.append(np.mean(profile_norm_store))
			profile_norm_stds.append(np.std(profile_norm_store))
			freqs_norm_scores.append(np.mean(freqs_norm_store))

			if init_type == 0:
				rand_means.append(np.mean(profile_norm_store))
				rand_stds.append(np.std(profile_norm_store))
			elif init_type == 1:
				pca_means.append(np.mean(profile_norm_store))
				pca_stds.append(np.std(profile_norm_store))
			else:
				ica_means.append(np.mean(profile_norm_store))
				ica_stds.append(np.std(profile_norm_store))


		#for i in range(len(types)):
		#	print types[i] + ' Profiles_norm ' + str(profile_norm_scores[i]) + ' Freqs_norm ' + str(freqs_norm_scores[i])

		#means_for_these_components.append()

	means = []
	stds = []
	means.append(rand_means)
	means.append(pca_means)
	means.append(ica_means)
	stds.append(rand_stds)
	stds.append(pca_stds)
	stds.append(ica_stds)
	plot_bar_chart(types, means, stds, ['2', '3', '4', '5', '6'])

	print '\nDONE'




