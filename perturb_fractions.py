



import numpy as np




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



#test

# freqs = [[.7, .2, .1],
# 			[1.],
# 			[.5, .5]]

# print freqs

# print distribute_evenly(freqs, .1)