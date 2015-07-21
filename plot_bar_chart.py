


#to plot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np


def plot_bar_chart(types, means, stds, xaxis_labels, save_as_name):

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


	plt.savefig(save_as_name)
	print 'Saved plot'