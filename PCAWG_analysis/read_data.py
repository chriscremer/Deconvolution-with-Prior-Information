
import csv
import os

from os.path import expanduser
home = expanduser("~")


gene_expression_file = home+'/PCAWG_data_17aug2016/joint_fpkm_uq.tsv'

#get gene expression sample names
with open(gene_expression_file, 'rU') as f:
	reader = csv.reader(f, delimiter='\t')
	count = 0
	for row in reader:
		# print count
		print row#[:10]
		# if count == 0:
		# 	sample_names = row[1:]
		# if count>0:
		# 	break
		count += 1
		if count > 0:
			break
# print 'Number of Samples: ' + str(count)
print len(row)

# print row[:10]
if '0292a46f-a282-4b7f-a7d6-ac55cc7324fb' in row:
	print 'yes'

else:
	print len(row)
	print 'no'



