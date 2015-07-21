


#add  directory to path to get my packages there
import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
#parentdir = os.path.dirname(currentdir)
#print sys.path
#sys.path.insert(0,parentdir) 
sys.path.insert(0,currentdir+'/simulated_data')
import make_convoluted_data

import numpy as np

import itertools

import plot_bar_chart as pbc
import read_TCGA_gene_data as rtgd
import weight_tools as wt
import model_tools as mt

import global_variables as gv

def main():

	#PARAMETERS
	plot_file_name = '../plots/bar_plot_different_initializations_test_first_x_samps_as_init.pdf'
	numb_samps = 50
	numb_feats = 10000



	rand_means = []
	pca_means = []
	ica_means = []
	rand_stds = []
	pca_stds = []
	ica_stds = []

	for numb_components in range(2, 6):

		print 'Numb of Components =' + str(numb_components)

		means_for_these_components = []
		stds_for_these_components = []

		#this is the number of profiles that exist in the simulated data
		numb_subpops = numb_components
		#this is the number of components that are kept in pca
		numb_model_subpops = numb_components
		
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
			for iteration in range(6):

				print 'Iter ' + str(iteration)

				#########################################################
				#Make data
				#params: #subpops, #feats, #samps
				samps, freqs, subpops = make_convoluted_data.run_and_return(numb_subpops, numb_feats, numb_samps)

				#global X
				X = samps
				gv.set_X_global(X)

				#make all frequencies have same number of entries
				new_freqs = wt.same_numb_of_entries(numb_model_subpops, freqs)
				#make all frequencies have same number of entries, WITHOUT SHUFFLING, used for comparing at the end
				start_freqs = wt.same_numb_of_entries_no_shuffle(numb_model_subpops, freqs)

				#########################################################
				#global possible_Ws
				possible_Ws = wt.get_possible_Ws(new_freqs)
				gv.set_Ws_global(possible_Ws)


				#try the same initilization type multiple times to prevent really bad local minimums
				best_norm = -1
				for try1 in range(5):

					#Initializing model
					TZ = mt.init_model(init_type, numb_model_subpops, X)
					gv.set_current_TZ(TZ)
					#print 'Optimizing model..'
					W, TZ = mt.optimize_model(possible_Ws, TZ)

					X_hat = np.dot(W, TZ)
					#print 'X_hat shape ' + str(X_hat.shape)
					norm = np.linalg.norm(X - X_hat)
					#print 'Initial norm ' + str(norm)
					if norm < best_norm or best_norm == -1:
						best_norm = norm
						best_W = W
						best_TZ = TZ

				W = best_W
				TZ = best_TZ

				#print 'Finding all weight permutations..'
				######################################################### 
				# match the components to their profiles so printing makes sense
				# so for each actual profile, find the row of TZ that is most similar to it
				component_order = wt.match_components_to_profiles(W, TZ, subpops)

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
	pbc.plot_bar_chart(types, means, stds, ['2', '3', '4', '5'], plot_file_name)

	print '\nDONE'






if __name__ == "__main__":

	main()



