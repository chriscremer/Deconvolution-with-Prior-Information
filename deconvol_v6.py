


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


	k = 20
	min_components = k
	max_components = k
	numb_of_iterations = 1
	numb_of_iters_to_remove_local_minima = 2
	init_types = ['Random_Samples']

	start = time.time()

	#iterate x times and get average
	for iteration in range(numb_of_iterations):

		print 'Iter ' + str(iteration)

		#Make data
		numb_subpops = k
		numb_feats = 1000
		numb_samps = 100
		percent_hidden = 0.5
		X, freqs, real_profiles = make_convoluted_data.run_and_return(numb_subpops, numb_feats, numb_samps, percent_hidden)

		gv.set_X_global(X)
		gv.set_freqs_global(freqs)

		print X.shape
		print freqs.shape
		# print 'Example freq ' + str(freqs[0])
		print real_profiles.shape


		#for each initialization type
		for init_type in range(len(init_types)):

			# print '-------Testing ' + init_types[init_type] + '--------'

			#try the same initilization type multiple times to prevent really bad local minimums
			best_norm = -1
			for try1 in range(numb_of_iters_to_remove_local_minima):

				#Initializing model
				TZ = mt.init_model(init_types[init_type], numb_subpops, X)
				gv.set_current_TZ(TZ)
				print 'Optimizing model..'
				W, TZ = mt.optimize_model_with_loglike(TZ)

				X_hat = np.dot(W, TZ)
				norm = np.linalg.norm(X - X_hat)
				if norm < best_norm or best_norm == -1:
					best_norm = norm
					best_W = W
					best_TZ = TZ

			W = best_W
			#print W
			TZ = best_TZ
			X_hat = np.dot(W, TZ)


		for learned_profile in range(len(TZ)):
			hidden_norms = []
			for hidden_profile in range(len(real_profiles)):
				hidden_norms.append(np.linalg.norm(TZ[learned_profile] - real_profiles[hidden_profile]))
			val, idx = min((val, idx) for (idx, val) in enumerate(hidden_norms))
			print 'Predicted k ' + str(learned_profile) + ' avg norm with hidden profiles ' + str(np.mean(hidden_norms)) + ' min norm ' + str(val) + ' with real profile' + str(idx)


		#LIKELIHODD
		# X=WZ
		# log likelihood = sum( ln (P(x|x_hat, cov))
		# cov = []
		# for 



	print '\nDONE'
	end = time.time()
	print end - start






if __name__ == "__main__":

	main()



