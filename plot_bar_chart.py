


#to plot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np


def plot_bar_chart(types, means, stds, xaxis_labels, save_as_name):

	# the x locations for the groups
	ind = np.arange(len(xaxis_labels))
	# the width of the bars
	width = 0.1

	fig, ax = plt.subplots()

	#random = ax.bar(ind-width, means[0], width, color='purple', yerr=stds[0])
	#pca = ax.bar(ind, means[1], width, color='red', yerr=stds[1])
	#ica = ax.bar(ind+width, means[2], width, color='blue', yerr=stds[2])

	colors = ['red', 'blue', 'green', 'purple']

	for i in range(len(types)):
		ax.bar(ind+(float(width*i)), means[i], width, yerr=stds[i], label=types[i], color=colors[i])

	# add some text for labels, title and axes ticks
	ax.set_ylabel('Average Norm Error')
	ax.set_xlabel('Number of Components')
	#ax.set_title('Scores by group and gender')
	ax.set_xticks(ind+width)
	ax.set_xticklabels(xaxis_labels)
	#ax.legend( (random[0], pca[0], ica[0]), ('Random', 'PCA', 'ICA') )
	#ax.legend( (random[0], pca[0]), types, fontsize=5 )
	ax.legend(fontsize=7)


	plt.savefig(save_as_name)
	print 'Saved plot'


def plot_bar_chart_no_stds(types, means, xaxis_labels, save_as_name):

	# the x locations for the groups
	ind = np.arange(len(xaxis_labels))
	# the width of the bars
	width = 0.1

	fig, ax = plt.subplots()

	#random = ax.bar(ind-width, means[0], width, color='purple', yerr=stds[0])
	#pca = ax.bar(ind, means[1], width, color='red', yerr=stds[1])
	#ica = ax.bar(ind+width, means[2], width, color='blue', yerr=stds[2])

	colors = ['red', 'blue', 'green', 'purple']

	for i in range(len(types)):
		ax.bar(ind+(float(width*i)), means[i], width, label=types[i], color=colors[i])

	# add some text for labels, title and axes ticks
	ax.set_ylabel('Average Norm Error')
	ax.set_xlabel('Number of Components')
	#ax.set_title('Scores by group and gender')
	ax.set_xticks(ind+width)
	ax.set_xticklabels(xaxis_labels)
	#ax.legend( (random[0], pca[0], ica[0]), ('Random', 'PCA', 'ICA') )
	#ax.legend( (random[0], pca[0]), types, fontsize=5 )
	ax.legend(fontsize=7)


	plt.savefig(save_as_name)
	print 'Saved plot'


def plot_line(y_values, x_values, y_label, x_label, save_as_name):

	plt.plot(x_values, y_values, 'bo')
	plt.xlabel(x_label)
	plt.ylabel(y_label)


	plt.savefig(save_as_name)
	print 'Saved plot'



def plot_line_with_text(y_values, x_values, y_label, x_label, save_as_name, text):

	fig = plt.figure()
	ax1 = fig.add_axes((.1,.4,.8,.5))
	ax1.plot(x_values, y_values, 'b.-')
	fig.text(.1,.15,text, size='x-small')
	plt.xlabel(x_label)
	plt.ylabel(y_label)


	plt.savefig(save_as_name)
	print 'Saved plot'

def plot_two_lines_with_text(y_values, y2_values, x_values, y_label, x_label, save_as_name, text):

	fig = plt.figure()
	ax1 = fig.add_axes((.1,.4,.8,.5))
	ax1.plot(x_values, y_values, 'b.-')
	ax1.plot(x_values, y2_values, 'r.-')
	fig.text(.1,.15,text, size='x-small')

	plt.xlabel(x_label)
	plt.ylabel(y_label)


	plt.savefig(save_as_name)
	print 'Saved plot'




def plot_three_lines_with_text(y_values, y2_values, y3_values, x_values, y_label, x_label, save_as_name, text):

	fig = plt.figure()
	ax1 = fig.add_axes((.1,.4,.8,.5))
	ax1.plot(x_values, y_values, 'b.-')
	ax1.plot(x_values, y2_values, 'r.-')
	ax1.plot(x_values, y3_values, 'g.-')
	fig.text(.1,.15,text, size='x-small')

	plt.xlabel(x_label)
	plt.ylabel(y_label)


	plt.savefig(save_as_name)
	print 'Saved plot'


	


