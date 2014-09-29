import yaml
import optparse
import glob
import itertools
import logging
import random
import matplotlib.pyplot as plt
import fragment
import utils
import lambda_damien
# Parse command line arguments to get config file name
parser = optparse.OptionParser()
parser.add_option('-c', '--config', dest='config_file_name', help='Name of configuration file to use')

(options, args) = parser.parse_args()

# Load configuration
config_file = open(options.config_file_name)
config = yaml.load(config_file)
config_file.close()

# Initialize logging
# Levels of logging:
# DEBUG - all messages
# INFO - basic information about execution phases etc.
# WARNING - warnings
# ERROR - errors
# CRYTICAL - crytical errors (not actually used)
logger = logging.getLogger('main_logger')
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
#handler = logging.FileHandler('main_clust_Colo.log')

formatter = logging.Formatter('%(filename)-20sline:%(lineno)-5d%(levelname)-8s [%(asctime)s]   %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)

# Dict in which we will save differnet stats to serialize it to file
# Stats from this file will be than used by probabilities.py
stats_to_serialize = dict()
stats_to_serialize['per_chr_stats'] = dict()


#Flag that shows if the histogram was already calculated ir not.
histogram_flag = 0 
# Get list of sam files in sam_files_dir
data_files = glob.glob(config['working_dir'] + config['sam_files_dir'] + '*.sam')
logger.info('Processing bam files from ' + config['working_dir'] + config['sam_files_dir'])
logger.info('Files: ' + str(data_files))
# Global list for translocations. Accumulated in memory while processing chromosomes and processed separately in the end.
all_translocations = []
normal_fragments_all = []
# Process data_files one by one. 
# Translocations are written to global list, so using joblib parallel is unsafe
# For adding parallel translocations should be instead dumped to sepaprate files
for data_file in data_files:

	fragments = fragment.load_sam_file(data_file)
	
	# Found tarnslocations and append them to global list
	chr_translocations = [frag for frag in fragments if frag.first_read_chr != frag.second_read_chr]
	all_translocations.extend(chr_translocations)
	fragments = list(set(fragments).difference(chr_translocations))

	# Split other records by direction
	ff = [frag for frag in fragments if frag.direction == 'ff' and frag.first_read_chr == frag.second_read_chr]
	rr = [frag for frag in fragments if frag.direction == 'rr' and frag.first_read_chr == frag.second_read_chr]
	fr = [frag for frag in fragments if frag.direction == 'fr' and frag.first_read_chr == frag.second_read_chr]
	rf = [frag for frag in fragments if frag.direction == 'rf' and frag.first_read_chr == frag.second_read_chr]

	# Chromosome name (same for all fragments except for translocation, pick from one)
	chromosome = rf[0].first_read_chr
	# Calculate different chromosome statistics

	(flag_direction, smallest_normal, biggest_normal, mean, median, R, epsilon,a,b,p) = utils.get_chromosome_stats(fr, rf, config,histogram_flag)
	logger.info('The appoximation normal distribution of the fragment lenght by the uniform distribution.')
	if not histogram_flag:
		logger.info('!!!!!!!!!!!!!!!!!!!!')
		logger.info('a: ' + str(a))
		logger.info('b: ' + str(b))
		logger.info('p: ' + str(p))



	logger.info('-----------------------------------------------------------------------------------')
	logger.info('Processing ' + data_file + '...')
	logger.info('-----------------------------------------------------------------------------------')

	logger.info('Fragments totally: ' + str(len(fragments)))
	logger.info('Number of translocations: ' + str(len(chr_translocations)))
	logger.info('Number of pair ended (fr) fragments: ' + str(len(fr)))
	logger.info('Number of mate ended (rf) fragments ' + str(len(rf)))
	logger.info('Number of double forward fragments: ' + str(len(ff)))
	logger.info('Number of double reverse fragments: ' + str(len(rr)))
	logger.info('-----------------------------------------------------------------------------------')

	# Print chromosome statistics
	logger.info('Direction flag: ' + flag_direction)
	logger.info('Confidence interval: [' + str(smallest_normal) + ';' + str(biggest_normal) + ']')
	logger.info('Mean: ' + str(mean))
	logger.info('Mediana: ' + str(median))
	logger.info('R: ' + str(R))
	logger.info('Epsilon: ' + str(epsilon))

	logger.info('-----------------------------------------------------------------------------------')

	# Save fragment lenghts to file for histogram building
	# logger.info('Saving fragment lengths...')
	# file_length = open(config['working_dir'] + config['results_dir'] + chromosome + '_fragment_lengths' + '.txt', 'w')
	# for i in fragment_lengths:
	# 	file_length.write(str(i) + ' ')
	# file_length.close()	

	# Sort lists with different directions by fragment begin
	ff.sort(key = lambda frag: frag.middle)
	rr.sort(key = lambda frag: frag.middle)
	fr.sort(key = lambda frag: frag.middle)
	rf.sort(key = lambda frag: frag.middle)
	
	# Get abnormal arrays, write abnormal and normal to files
	abn=0
	norm=0
	ff_abn=[]
	fr_abn=[]
	rf_abn=[]
	rr_abn=[]

	logger.info('Saving abnormal fragments...')
	# Saving abnormal fragments to file is not actually needed
	# fr_wr = open(config['working_dir'] + config['results_dir'] + chromosome + '_abnormal_frag.txt','w')
	normal_fragments = []
	for frag in fragments:
		if not frag.is_abnormal(smallest_normal, biggest_normal) or frag.direction == flag_direction:
			if frag.mapp_qul_flag or frag.unique_flag:
				normal_fragments.append([frag.begin, frag.first_read_chr, frag.unique_flag,frag.mapp_qul_flag,frag.name])
				norm += 1
	# fr_wr.close()

	normal_fragments.sort()
	normal_fragments_all.extend(normal_fragments)
	#Write abnornal fragments to the file
	fr_norm = open(config['working_dir'] + config['normal_fragments_dir'] + 'normal_fragments_' + chromosome + '.txt', 'w')
	fr_norm.write( '| chr | begin | unique_flag | mapp_qul_flag | name | \n')
	for frag in normal_fragments:
		fr_norm.write(str(frag[1]) + ' ' + str(frag[0]) +' '+str(frag[2]) + ' '+str(frag[3])+' '+frag[4]+'\n')
	fr_norm.close()
	logger.info('Number of normal fragments is ' + str(norm))
	# Write begin-end pairs and middles to separate files for visualization
	# write_begin_end_mid_files (fr_abn, rf_abn, rr_abn, ff_abn, chromosome, config)
	M = biggest_normal 	
	# Add stats for this chromosome for serialization
	stats_to_serialize['per_chr_stats'][chromosome] = dict()
	stats_to_serialize['per_chr_stats'][chromosome]['num_all_abn'] = abn
	stats_to_serialize['per_chr_stats'][chromosome]['R'] = R
	stats_to_serialize['per_chr_stats'][chromosome]['smallest_normal'] = smallest_normal
	stats_to_serialize['per_chr_stats'][chromosome]['biggest_normal'] = biggest_normal
	stats_to_serialize['per_chr_stats'][chromosome]['median'] = median
	if not histogram_flag:
		stats_to_serialize['a'] = float(a)
		stats_to_serialize['b'] = float(b)
		stats_to_serialize['p'] = float(p)
	histogram_flag+=1
lambda_gc = lambda_damien.LambdaCalculation(chromosome,config,normal_fragments_all,median)
# Serialize stats
serialized_stats_file = open(config['working_dir'] + config['serialized_stats_file'], 'w')
serialized_stats_file.write(yaml.dump(stats_to_serialize))
serialized_stats_file.close()
