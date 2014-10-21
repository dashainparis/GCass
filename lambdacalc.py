import glob
import yaml
import optparse
import logging
import itertools
import pylab 
from operator import is_not
from functools import partial

logger = logging.getLogger('main_logger')
parser = optparse.OptionParser()
parser.add_option('-c', '--config', dest='config_file_name', help='Name of configuration file to use')



def SortLambda(l):
	return(l[0])
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


def LambdaCalculation(config,chr_normal_fragments_dict,median):
	# Load configuration
	#config_file = open(options.config_file_name)
	#config = yaml.load(config_file)
	#config_file.close()
	Fgc_Ngc_mapp = 	dict()
	step = 1
	contigs_lambda = []
	for contig in chr_normal_fragments_dict.keys():
		chrom_line = LoadChromLine(config,contig)
		if chrom_line>6000:
			contigs_lambda.append(contig)
	logger.info('We will use this contigs for lambda calculation')
	logger.info(contigs_lambda)
	for contig in contigs_lambda:
		gem_val_mapp = []
		gc_rate = []
		gc_gem = []
		#Fgc_Ngc[chrom] = {'Fgc':0,'Ngc':0}
		chrom_line = LoadChromLine(config,contig)
		gem_line = LoadGem(config,contig)
		normal_fragments = chr_normal_fragments_dict[contig]
		logger.info('Now for this chromosome ' + contig)
		logger.info('Uploading process is finished')
		for i in xrange(0,len(gem_line)-median,step):
			if gem_line[i]=='!':
				gem_val_mapp.append(1)
			else:
				gem_val_mapp.append(-1)
		logger.info('Finished wth gem')
		for i in xrange(0,len(chrom_line)-median,step):
			gc_rate.append(round((chrom_line[i:i+median].count('g')+chrom_line[i:i+median].count('c'))/float(median),2))
		gc_gem_mapp = [gc_rate[i]*gem_val_mapp[i] for i in xrange(len(gem_val_mapp)-1)]
		logger.info('All windows are counted, now Fgc')
		for i in xrange(len(gc_gem_mapp)):
			if gc_gem_mapp[i]<0:
				gc_gem_mapp[i] = 'n'
		for i in set(gc_gem_mapp):
			if i not in Fgc_Ngc_mapp.keys():
				Fgc_Ngc_mapp[i] = {'Fgc':0,'Ngc':0}
			Fgc_Ngc_mapp[i]['Ngc'] += gc_gem_mapp.count(i)
		norm_ind=0
		ind_left=0
		logger.info('Biggest normal frag = '+ str(normal_fragments[-1]))
		logger.info('len chrom_line ='+str(len(chrom_line)))
		logger.info('lenght gc_gem_mapp = '+str(len(gc_gem_mapp)))
		logger.info('lenght gem len mapp = '+str(len(gem_val_mapp)))
		for i in normal_fragments:
			if i>len(gc_gem_mapp)-median:
				break
			ind_left+=1	
			if ind_left%50000==0:
				logger.info('Already processed '+str(ind_left)+' fragments, left '+str(len(normal_fragments)-ind_left))
			if i%step == 0:
				#if i<=len(gc_gem_mapp):
				if gc_gem_mapp[i/step]!='n':
					Fgc_Ngc_mapp[gc_gem_mapp[i/step]]['Fgc'] +=1
		lam_mapp = []		

	if 0.7 not in Fgc_Ngc_mapp.keys():
		Fgc_Ngc_mapp[0.66] = {'Fgc':0,'Ngc':0}
	if 0.4 not in Fgc_Ngc_mapp.keys():
		Fgc_Ngc_mapp[0.4] = {'Fgc':0,'Ngc':0}
	for key in Fgc_Ngc_mapp.keys():
		if key>0.66:
			Fgc_Ngc_mapp[0.66]['Fgc']+=Fgc_Ngc_mapp[key]['Fgc']
			Fgc_Ngc_mapp[0.66]['Ngc']+=Fgc_Ngc_mapp[key]['Ngc']
			logger.info(str(key)+ 'delete')
			del Fgc_Ngc_mapp[key]
		if key<0.4:
			Fgc_Ngc_mapp[0.4]['Fgc']+=Fgc_Ngc_mapp[key]['Fgc']
			Fgc_Ngc_mapp[0.4]['Ngc']+=Fgc_Ngc_mapp[key]['Ngc']
			logger.info(str(key)+ 'delete')
			del Fgc_Ngc_mapp[key]
		
	logger.info('Mappable lambda')
	logger.info('rate 	Fgc 	Ngc')
	for i in Fgc_Ngc_mapp.keys():
		logger.info(str(i)+' ' + str(Fgc_Ngc_mapp[i]['Fgc'])+' '+str(Fgc_Ngc_mapp[i]['Ngc']))
		lam_mapp.append([i, float(Fgc_Ngc_mapp[i]['Fgc'])/Fgc_Ngc_mapp[i]['Ngc']])
	lam_mapp.sort(key = SortLambda)
	x = []
	y = []
	for k in lam_mapp:
		x.append(float(k[0]))
		y.append(float(k[1]))
	logger.info('lambda Mappable')
	logger.info(lam_mapp)
	f = open(config['working_dir'] + config['lambda_file'],'w')
	for i in lam_mapp:
		f.write(str(i)+'\n')
	f.close()
	
	logger.info(x)
	logger.info(y)
	pylab.plot(x,y)
	pylab.xlabel('gc-rate')
	pylab.ylabel('Expected number of fragment at the position')
	pylab.title('Estimation of lambda correction')
	pylab.grid(True)
	pylab.savefig(config['working_dir'] +'lambda.png')
	pylab.show()

	return(lam_mapp)










