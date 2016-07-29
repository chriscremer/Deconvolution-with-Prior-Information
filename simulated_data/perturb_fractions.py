



import numpy as np
import math



def distribute_evenly(fractions, amount):
	'''
	Multiply fractions by amount
	Distribute that amount evenly among fractions
	So 100 percent will result in an even prior
	'''

	new_freqs = []
	for samp in fractions:
		new_samp = []
		for freq in samp:
			new_samp.append(freq - (freq * amount) + (sum(samp)*amount/len(samp)))
		new_freqs.append(new_samp)

	return new_freqs


def add_noise(fractions, amount):
	'''
	Add some noise
	'''
	if amount == 0.:
		return fractions

	new_freqs = []
	for samp in fractions:
		new_samp = []
		for freq in samp:
			# new_freq = freq + np.random.normal(0,amount)
			new_freq = freq + (np.random.normal(0,amount**2))
			#no negatives
			if new_freq < 0:
				new_freq = 0.001
			new_samp.append(new_freq)
		#sum to 1
		# print new_samp
		# print new_samp[0]
		# print sum(new_samp[i])
		for i in range(len(new_samp)):
			# print i
			new_samp[i] = new_samp[i]/sum(new_samp)

		new_freqs.append(new_samp)

	return new_freqs




# #test
# if __name__ == "__main__":
# 	freqs = [[.7, .2, .1],
# 				[1.],
# 				[.5, .5]]

# 	print freqs

# 	# print distribute_evenly(freqs, .1)

# 	aaa =  add_noise(freqs, 0.5)

# 	for i in aaa:
# 		print i