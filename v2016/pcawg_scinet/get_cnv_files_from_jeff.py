
#Plan: read the type file, then copy all the cnv files of one type here. 

from shutil import copyfile


from os.path import expanduser
home = expanduser("~")

import csv

import os

# CNV: /scratch/q/qmorris/jawinter/pwgs/cnvs/vanloo_wedge
#NO vanloo_wedge_segs has more than just wedge 

# Types: /scratch/q/qmorris/jawinter/pwgs/share/summary_tables/summary_table.consprelim.sample-5000.psub.txt


# print os.listdir(home + '/../jawinter/scratch/pwgs/share/summary_tables')


#GET IDS OF PBCAs
sample_ids = []
with open(home + '/../jawinter/scratch/pwgs/share/summary_tables/summary_table.consprelim.sample-5000.psub.txt', 'rb') as csvfile:
	reader = csv.reader(csvfile, delimiter='\t')
	for row in reader:
		# print row
		if row[0] == 'PBCA':

			# print row
			sample_ids.append(row[1])



#COPY THEIR CNV FILES TO THIS DIRECTORY

all_cnv_files =  os.listdir(home + '/../jawinter/scratch/pwgs/cnvs/vanloo_wedge_segs')

n_found =0
for s_id in sample_ids:

	found = 0

	for f in all_cnv_files:

		if s_id in f:
			# print 'found'
			found = 1
			n_found +=1

			#COPY THE FILE
			copyfile(home + '/../jawinter/scratch/pwgs/cnvs/vanloo_wedge_segs/' + f, home + '/Deconvolution_stuff/cnv_files/'+f)

			break

	if found == 0:
		print 'not found'
	else: 
		print 'yessssss'




print len(sample_ids)
print len(all_cnv_files)

print n_found










