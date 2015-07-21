


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



def main():

	#PARAMETERS
	plot_file_name = 'xxxx'



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

				#make all frequencies have same number of entries
				new_freqs = wt.same_numb_of_entries(numb_model_subpops, freqs)
				#make all frequencies have same number of entries, WITHOUT SHUFFLING, used for comparing at the end
				start_freqs = wt.same_numb_of_entries_no_shuffle(numb_model_subpops, freqs)

				#########################################################
				#Initializing model
				TZ = mt.init_model(init_type, numb_model_subpops, samps)

				X_hat = np.dot(new_freqs, TZ)
				#print 'X_hat shape ' + str(X_hat.shape)
				norm = np.linalg.norm(samps - X_hat)
				#print 'Initial norm ' + str(norm)

				#print 'Finding all weight permutations..'
				global possible_Ws
				possible_Ws = wt.get_possible_Ws(new_freqs)
				#print 'Completed.'
				#########################################################


				#########################################################
				print 'Optimizing model..'
				W, TZ = mt.optimize_model(possible_Ws, TZ)


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
	pbc.plot_bar_chart(types, means, stds, ['2', '3', '4', '5', '6'])

	print '\nDONE'






if __name__ == "__main__":

	main()



