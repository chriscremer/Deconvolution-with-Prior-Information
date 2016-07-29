


#add  directory to path to get my packages there
import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
#parentdir = os.path.dirname(currentdir)
#print sys.path
#sys.path.insert(0,parentdir) 
sys.path.insert(0,currentdir+'/simulated_data')
import make_convoluted_data

import numpy as np
import scipy
import csv

import itertools

import plot_bar_chart as pbc
import read_TCGA_gene_data as rtgd
import weight_tools as wt
import model_tools as mt

import global_variables as gv

import time

from scipy.stats import multivariate_normal as mn

def main():


	k = 2
	min_components = k
	max_components = k
	numb_of_iterations = 1
	numb_of_iters_to_remove_local_minima = 1
	init_types = ['Random_Samples']

	start = time.time()

	data_file = '../affymetrix_data2/data/data.txt'
	X = np.genfromtxt(data_file, delimiter='\t', skip_header=1, usecols=range(1,34))
	X = X.T

	#iterate x times and get average
	for iteration in range(numb_of_iterations):

		print 'Iter ' + str(iteration)

		#Make data
		# numb_subpops = k
		# numb_feats = 2
		# numb_samps = 10
		# percent_hidden = 0.8
		# X, freqs, real_profiles, real_freqs = make_convoluted_data.run_and_return(numb_subpops, numb_feats, numb_samps, percent_hidden)


		#this for testing the all permutations way
		# new_freqs = wt.same_numb_of_entries(k, freqs)
		# possible_Ws = wt.get_possible_Ws(new_freqs, allow_combos=False)
		# gv.set_Ws_global(possible_Ws)

		gv.set_X_global(X)
		# gv.set_freqs_global(freqs)

		print 'X shape ' + str(X.shape)
		# print 'W len ' + str(len(freqs))
		# print 'Proportions shape ' + str(freqs.shape)
		# print 'Example freq ' + str(freqs[0])
		# print 'Real TZ shape ' + str(real_profiles.shape)

		#for each initialization type
		# for init_type in range(len(init_types)):

		# print '-------Testing ' + init_types[init_type] + '--------'

		#try the same initilization type multiple times to prevent really bad local minimums
		best_norm = -1
		for try1 in range(numb_of_iters_to_remove_local_minima):

			#Initializing model
			TZ = mt.init_model('Random_Samples', k, X)
			# print TZ.shape

			#just to show it gets perfect
			# TZ = real_profiles
			# print TZ.shape

			# initial_TZ = np.copy(TZ)
			# print 'W init'
			# print W
			# gv.set_current_W(W)
			gv.set_current_TZ(TZ)
			print 'Optimizing model..'
			# W, TZ = mt.optimize_model_with_loglike(TZ)
			W, TZ = mt.optimize_model()
			X_hat = np.dot(W, TZ)
			norm = np.linalg.norm(X - X_hat)
			if norm < best_norm or best_norm == -1:
				best_norm = norm
				best_W = W
				best_TZ = TZ
				# best_init = initial_TZ

		W = best_W
		TZ = best_TZ

		print W


		# X_hat = np.dot(W, TZ)


		# print '\nEvaluate Performance'
		# indexes = mt.match_profiles(TZ, real_profiles)

		# for learned_profile in range(len(TZ)):
		# 	hidden_norms = []
		# 	for hidden_profile in range(len(real_profiles)):
		# 		hidden_norms.append(np.linalg.norm(TZ[learned_profile] - real_profiles[hidden_profile]))
		# 	val, idx = min((val, idx) for (idx, val) in enumerate(hidden_norms))
		# 	print 'Predicted k ' + str(learned_profile) + ' avg norm with hidden profiles ' + str(np.mean(hidden_norms)) + ' min norm ' + str(val) + ' with real profile' + str(idx)


		# print 'LEARNED'
		# for list1 in TZ:
		# 	print str(['%.2f' % elem for elem in list1])
		# print 'REAL'
		# for list1 in real_profiles[indexes]:
		# 	print str(['%.2f' % elem for elem in list1])
		# print 'W'
		# for list1 in W:
		# 	print str(['%.2f' % elem for elem in list1])
		# print 'Real freqs'
		# r_f = real_freqs.T[indexes].T
		# for list1 in r_f:
		# 	print str(['%.2f' % elem for elem in list1])



	# pbc.plot_visualize_learning_iters('../plots/visualize_learning/nnls.png', X, real_profiles[indexes], W, real_freqs.T[indexes].T, best_init, TZ)



	with open(data_file, 'rU') as f:
		reader = csv.reader(f, delimiter='\t')
		for row in reader:
			header = row
			break
	header.pop(0)
	for samp in range(len(header)):
		if 'mix1' in header[samp]:
			print str(W[samp]) + ' = mix1 (0, 1)'
			# cellularities.append([0.0, 1.0])
		elif 'mix2' in header[samp]:
			print str(W[samp]) + ' = mix2 (0.05, 0.95)'
			# cellularities.append([0.05, 0.95])
		elif 'mix3' in header[samp]:
			print str(W[samp]) + ' = mix3 (0.10, 0.90)'
			# cellularities.append([0.10, 0.90])
		elif 'mix4' in header[samp]:
			print str(W[samp]) + ' = mix4 (0.25, 0.75)'
			# cellularities.append([0.25, 0.75])
		elif 'mix5' in header[samp]:
			print str(W[samp]) + ' = mix5 (0.5, 0.5)'
			# cellularities.append([0.5, 0.5])
		elif 'mix6' in header[samp]:
			print str(W[samp]) + ' = mix6 (0.75, 0.25)'
			# cellularities.append([0.75, 0.25])
		elif 'mix7' in header[samp]:
			print str(W[samp]) + ' = mix7 (0.9, 0.1)'
			# cellularities.append([0.9, 0.1])
		elif 'mix8' in header[samp]:
			print str(W[samp]) + ' = mix8 (0.95, 0.05)'
			# cellularities.append([0.95, 0.05])
		elif 'mix9' in header[samp]:
			print str(W[samp]) + ' = mix9 (1.0, 0.0)'


	print '\nDONE'
	end = time.time()
	print end - start






if __name__ == "__main__":

	main()



