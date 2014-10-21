import glob
import logging
from cluster import *
from numpy import arange,array,ones,linalg
from pylab import plot,show
from pylab import *
from operator import is_not
from functools import partial
class InputData:
	###########################
	# Several dictionaries of data related to one chromosome
	# With chromosome name as a key
	###########################
	# value - list of normal fragments begins
	chr_normal_fragments_dict = dict()
	chr_normal_fragments_dict_unique = dict()
	# value - list of centromer (begin, end) pairs
	chr_centrom_dict = dict()
	# value - string with mappability info from gem file
	chr_gem_dict = dict()
	# value - chromosome string from fa file
	chr_line_dict = dict()
	# value - list of sub_links, sorted by begin
	chr_sub_links_dict = dict()
	# CN information about each chromosome
	chr_cnv_dict = dict()
	###########################
	# Chromosomes for which data is loaded
	chromosomes = []
	# List of lambda data pairs (gc share, probability)
	lambda_gen = []
	# Lambda gen values for possible percents of gc (0-100)
	#lambda_gen_index = []
	lambda_norm_index = []
	lambda_abnorm_index =[]
	lambda_abnorm = []
	lambda_norm =[]
	# List of pairs (fragment length, probability of this length)
	length_probabilities = []
	# List of pairs (fragment length, probability of length <= than this)
	cummulative_length_probabilities = []
	# All sub links for all chromosomes stored here
	# TODO remove this, beacuse it duplicates chr_sub_links_dict
	sub_links = []
	# Sorted list of element numbers observed in clusters
	numb_elem = []

	def __init__(self, config, stats, chromosomes):
		logger.info('Loading data...')
		self.chromosomes = chromosomes
		# Load centrom, normal fragments and sublinks for all chromosomes
		# Do not load fa and gem line for all chromosomes - it requires too much memory,
		# so we load it for separate chromosomes on demand
		#self.__LoadFreec(config)
		self.__LoadNormalFragments(config)
		self.__LoadCentrom(config)
		self.__LoadLinksAndConstructSubLinks(config, stats)
		# Load common data
		self.__LoadLambda(config, stats)
		self.__LoadLengthProbabilities(config, stats)
		logger.info('Finished loading data')

	# Load fa and gem data required to process chromosome
	def LoadChrom(self, config, chrom):
		self.__LoadGem(config, chrom)
		self.__LoadChromLine(config, chrom)
		self.__LoadFreec(config,chrom)

	def UnloadChrom(self, chrom):
		self.chr_gem_dict.pop(chrom)
		self.chr_line_dict.pop(chrom)

	# def LogLoadedDataInfo(self):
	# 	logger.info('Loaded data for chromosomes ' + str(self.chromosomes))
	# 	logger.info('Number of elements in lambda_gen: ' + str(len(self.lambda_gen)))
	# 	logger.info('First: ' + str(self.lambda_gen[0]), ', last: ' + str(self.lambda_gen[-1]))
	# 	logger.info('Number of elements in length_probabilities: ' + str(len(self.length_probabilities)))
	# 	logger.info('First: ' + str(self.length_probabilities[0]), ', last: ' + str(self.length_probabilities[-1]))
	# 	logger.info('Number of elements in cummulative_length_probabilities: ' + str(len(self.cummulative_length_probabilities)))
	# 	logger.info('First: ' + str(self.cummulative_length_probabilities[0]), ', last: ' + str(self.cummulative_length_probabilities[-1]))
	# 	logger.info('Number of elements in sub_links: ' + str(len(self.sub_links)))
	# 	logger.info('First: ' + str(self.sub_links[0]), ', last: ' + str(self.sub_links[-1]))
	# 	logger.info('Number of elements in numb_elem: ' + str(len(self.numb_elem)))
	# 	logger.info('First: ' + str(self.numb_elem[0]), ', last: ' + str(self.numb_elem[-1]))

	#########################################################################
	# Functions below should not be used directly from outside the class
	# So names are mungled
	#########################################################################
	def __LoadGem(self, config, chrom):
		#f = open(config['working_dir'] + config['gem_files_dir'] + chrom + '_hg38')
		f = open(config['working_dir'] + config['gem_files_dir'] + chrom + '_gem.txt')
		gem_line = ''
		for line in f:
			if '~' not in line:
				if line[-1] == '\n':
					gem_line += line[:-1]
				else:
					gem_line += line
		f.close()
		self.chr_gem_dict[chrom] = gem_line

	def __LoadChromLine(self, config, chrom):
		f = open(config['working_dir'] + config['fa_files_dir'] + chrom + '.fa')
		chrom_line = ''
		for line in f:
			chrom_line += line.rstrip()
		f.close()
		self.chr_line_dict[chrom] = chrom_line.lower()