import glob
import bisect
import optparse
import yaml
import logging
from numpy import arange,array,ones,linalg
from pylab import plot,show
from pylab import *
import random
import matplotlib.pyplot as plt

logger = logging.getLogger('main_logger')
parser = optparse.OptionParser()
parser.add_option('-c', '--config', dest='config_file_name', help='Name of configuration file to use')

(options, args) = parser.parse_args()

def LoadChromLine(config, chrom):
	f = open(config['working_dir'] + config['fa_files_dir'] + chrom + '.fa')
	chrom_line = ''
	for line in f:
		chrom_line += line.rstrip()
	f.close()
	return(chrom_line.lower())
def LoadGem(config, chrom):
	f = open(config['working_dir'] + config['gem_files_dir'] + chrom + '_gem.txt')
	gem_line = ''
	for line in f:
		if '~' not in line:
			if line[-1] == '\n':
				gem_line += line[:-1]
			else:
				gem_line += line
	f.close()
	return(gem_line)

def LoadNormalFragments(chrom, config):
	chr_nf_list = []
	f = open(config['working_dir'] + config['normal_fragments_dir'] + 'normal_fragments_' + chrom + '.txt')
	for line in f:
		(c, nf) = line.split()
		chr_nf_list.append(int(nf))
	f.close()
	return(chr_nf_list)



# Load serialized stats, calculated by main_clustering.py

#chromosomes = config['chromosomes']


def LambdaCalculation(chrom,config,normal_fragments,median):
	# Load configuration
	config_file = open(options.config_file_name)
	config = yaml.load(config_file)
	config_file.close()
	Fgc_Ngc_mapp = 	dict()
	step = 21
	
	gem_val_mapp = []
	gc_rate = []
	gc_gem = []
	#Fgc_Ngc[chrom] = {'Fgc':0,'Ngc':0}
	chrom_line = LoadChromLine(config,chrom)
	gem_line = LoadGem(config,chrom)
	normal_fragments.sort()
	print 'Now for this chromosome ', chrom
	print 'Uploading process is finished'
	for i in xrange(0,len(gem_line)-median,step):
		if gem_line[i]=='!':
			gem_val_mapp.append(1)
		else:
			gem_val_mapp.append(-1)
	print 'Finished wth gem'
	for i in xrange(0,len(chrom_line)-median,step):
		gc_rate.append(round((chrom_line[i:i+median].count('g')+chrom_line[i:i+median].count('c'))/float(median),2))
	gc_gem_mapp = [gc_rate[i]*gem_val_mapp[i] for i in xrange(len(gem_val_mapp)-1)]
	print 'All windows are counted, now Fgc'
	for i in xrange(len(gc_gem_mapp)):
		if gc_gem_mapp[i]<0:
			gc_gem_mapp[i] = 'n'
	for i in set(gc_gem_mapp):
		if i not in Fgc_Ngc_mapp.keys():
			Fgc_Ngc_mapp[i] = {'Fgc':0,'Ngc':0}
		Fgc_Ngc_mapp[i]['Ngc'] += gc_gem_mapp.count(i)
	norm_ind=0
	ind_left=0
	print 'Biggest normal frag = ',normal_fragments[-1]
	print 'len chrom_line =',len(chrom_line)
	print 'lenght gc_gem_mapp = ',len(gc_gem_mapp)
	print 'lenght gem len mapp = ',len(gem_val_mapp)
	for i in normal_fragments:
		if i[0]>len(gem_line)-median:
			break
		ind_left+=1	
		if ind_left%50000==0:
			print 'Already processed ',ind_left,' fragments, left ',len(normal_fragments)-ind_left
		if i[0]%step == 0:
			if gc_gem_mapp[i[0]/step]!='n':
				Fgc_Ngc_mapp[gc_gem_mapp[i[0]/step]]['Fgc'] +=1
	lam_mapp = []		

	if 0.7 not in Fgc_Ngc_mapp.keys():
		Fgc_Ngc_mapp[0.66] = {'Fgc':0,'Ngc':0}
	if 0.3 not in Fgc_Ngc_mapp.keys():
		Fgc_Ngc_mapp[0.3] = {'Fgc':0,'Ngc':0}
	for key in Fgc_Ngc_mapp.keys():
		if key>0.66:
			Fgc_Ngc_mapp[0.66]['Fgc']+=Fgc_Ngc_mapp[key]['Fgc']
			Fgc_Ngc_mapp[0.66]['Ngc']+=Fgc_Ngc_mapp[key]['Ngc']
			print key, 'delete'
			del Fgc_Ngc_mapp[key]
		if key<0.3:
			Fgc_Ngc_mapp[0.3]['Fgc']+=Fgc_Ngc_mapp[key]['Fgc']
			Fgc_Ngc_mapp[0.3]['Ngc']+=Fgc_Ngc_mapp[key]['Ngc']
			print key, 'delete'
			del Fgc_Ngc_mapp[key]
	print 'Mappable lambda'
	print 'rate 	Fgc 	Ngc'
	for i in Fgc_Ngc_mapp.keys():
		print i,  Fgc_Ngc_mapp[i]['Fgc'], Fgc_Ngc_mapp[i]['Ngc']
		lam_mapp.append([i, float(Fgc_Ngc_mapp[i]['Fgc'])/Fgc_Ngc_mapp[i]['Ngc']])

	print 'lambda Mappable'
	print lam_mapp
	f = open(config['working_dir'] + config['lambda_file'],'w')

	return(lam_mapp)










