
import numpy as np
import random


#for now, make samples have 1000 features




def make_subpopulations(numb_of_subpops, numb_of_feats):
	'''
	Define possible subpopulations/subtypes profiles
	'''

	subpops = np.random.rand(numb_of_subpops, numb_of_feats)

	return subpops


def subpop_frequencies(numb_samps, numb_of_subpops):
	'''
	Define the frequency of each subpop in each sample
	'''

	#freqs = np.random.rand(numb_samps, numb_of_subpops)
	#I want some samples to have 0 frequency for some profiles
	freqs = []
	for i in range(numb_samps):
		this_samp = []
		for j in range(numb_of_subpops):
			rand = random.random()
			if rand > 0.5:
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
			a = random.randint(0, numb_of_subpops-1)
			freqs[i][a] = 1.0
			continue
		#divide row by sum
		freqs[i] = freqs[i] / sum_of_row

	
	return freqs


def make_samps(numb_samps, subpops, freqs):

	samps = np.zeros((numb_samps, len(subpops[0])))


	for i in range(len(freqs)):

		for j in range(len(freqs[0])):

			# print 'samps1 shape ' + str(samps[i].shape)
			# print 'freqs1 shape ' + str(freqs[i][j].shape)
			# print 'subpops1 shape ' + str(subpops[j].shape)
			samps[i] = samps[i] + freqs[i][j]*subpops[j]

	
	#add some noise, divide to make it smaller
	noise = np.random.rand(numb_samps, len(subpops[0])) / 100.0
	samps = samps + noise
	#need to make sure that no cell is negative
	#cant be negative cuz the random is from 0 to 1 not -1 to 1
	#so could remove this
	for i in range(len(samps)):
		for j in range(len(samps[0])):
			if samps[i][j] < 0.0:
				samps[i][j] = 0.0



	return samps




def make_csv(data):
	'''
	Writes the data to a csv file
	'''

	f = open('simulated_data_' + str(numb_samps) + '_' + str(numb_of_feats) + '.csv', 'w')

	for samp in range(len(data)):
		for feature in range(len(data[samp])):
			f.write(str(data[samp][feature]) + ', ')
		f.write('\n')
	f.close()


	print 'File has been made.'



def run_and_return(numb_of_subpops, numb_of_feats, numb_samps):
	'''
	Run the methods to create the data and return it
	'''

	subpops = make_subpopulations(numb_of_subpops, numb_of_feats)
	#print 'subpops shape ' + str(subpops.shape)
	freqs = subpop_frequencies(numb_samps, numb_of_subpops)
	#print 'freqs shape ' + str(freqs.shape)
	samps = make_samps(numb_samps, subpops, freqs)
	#print 'samps shape ' + str(samps.shape)

	return samps, freqs, subpops









if __name__ == "__main__":

	numb_of_subpops = 10
	numb_of_feats = 100
	numb_samps = 50

	# subpops = make_subpopulations(numb_of_subpops, numb_of_feats)
	# print 'subpops shape ' + str(subpops.shape)
	# freqs = subpop_frequencies(numb_samps, numb_of_subpops)
	# print 'freqs shape ' + str(freqs.shape)
	# samps = make_samps(numb_samps, subpops, freqs)
	# print 'samps shape ' + str(samps.shape)
	
	#make_csv(samps)
	#run_and_return(numb_of_subpops, numb_of_feats, numb_samps)

