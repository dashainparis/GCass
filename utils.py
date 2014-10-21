import logging
from scipy.stats import norm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import numpy as np
import bisect

logger = logging.getLogger('main_logger')

####################################################################
# Calculate different statistics
####################################################################
def get_chromosome_stats(fr, rf, config,historam_flag):
	# Build confidence interval. Take fragment lengths for most popular type (it may be fr or rf)
	if len(fr) > len(rf):
		fragment_lengths = [frag.length for frag in fr]
		flag_direction='fr'
	else:
		fragment_lengths = [frag.length for frag in rf]
		flag_direction='rf'
	
	fragment_lengths = sorted(fragment_lengths)
	normal_lengths = fragment_lengths[0:bisect.bisect_left(fragment_lengths,10000)]
	h = int(round(config['alpha'] * len(fragment_lengths)))	# Number of segments that should be excluded from each side
	normal_lengths = normal_lengths[int(round(h))+1:-int(round(h)) -1]
	#normal_lengths =[i for i in fragment_lengths if fragment_lengths[0]<=i<= fragment_lengths[-1]]
	#normal_lengths =[i for i in fragment_lengths if 900<=i<= 4000]	
	smallest_normal = normal_lengths[0]
	biggest_normal = normal_lengths[-1]
	
	# Mean, standart deviation and Radius
	mean = sum(normal_lengths) / len(normal_lengths)
	q_sum = 0
	for i in range(0, len(normal_lengths) - 1):
		q_sum += (normal_lengths[i] - mean)**2
	sigma = (q_sum / len(normal_lengths))**0.5
	normal_lengths.sort()
	if len(normal_lengths) % 2 == 0:
		n = len(normal_lengths)
		mediana = (normal_lengths[n/2-1]+ normal_lengths[n/2] )/2
	else:
		mediana =normal_lengths[len(normal_lengths)/2] 
	R = (mediana) + 3*sigma/(2**0.5)
	epsilon = mediana / 5.0

	#Histogram calculation
	logger.info('I started calculate the histogram')
	logger.info('normal_lengths = ['+str(normal_lengths[0])+str('; ')+str(normal_lengths[-1])+str(']'))
	a=0
	b=0
	p=0
	if not historam_flag:
		#bins = set(normal_lengths)
		historam_file = open(config['working_dir']+config['length_histogram_file'],'w')
		n, bins, patches = plt.hist(normal_lengths, len(set(normal_lengths)), normed=1, facecolor='green', alpha=0.75)
		print len(n)
		print len(bins)
		for i in range(len(bins)-1):
			historam_file.write(str(round(bins[i])) + ' '+str(round(n[i],9))+'\n')
		historam_file.close()
		plt.close()
		#(mu, sigma) = norm.fit(normal_lengths)
		#logger.info('mu = '+str(mu))
		#logger.info('sigma = '+str(sigma))
		#a=mu-(12**0.5)*0.5*sigma
		#b=mu+(12**0.5)*0.5*sigma
		#p = 1/(b-a)
		#f=open(config['working_dir']+'length.txt','w')
		#f.write(str(fragment_lengths))
		#f.close()
	return (flag_direction, smallest_normal, biggest_normal, mean, mediana, R, epsilon,a,b,p)
