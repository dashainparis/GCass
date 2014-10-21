import glob
#files=['sampe_1.sam','sampe_2_.sam','sampe_3_.sam','sampe_4_.sam','sampe_6_.sam','sampe_7_.sam','sampe_8_.sam','081224_EAS188_0073_FC30TTCAAXX_1_new.sam','081224_EAS188_0073_FC30TTCAAXX_3_new.sam']
chrom_names=[]
#for i in range(1,23,1):
#	chrom_names.append('chr'+str(i))

#chrom_names.append('chrY')
#chrom_names.append('chrX')
#chrom_names.remove('chr20')
#chrom_names.remove('chr7')
#chrom_names =['chr20']
#print chrom_names
chrom_names = ['NODE_101_length_14188_cov_4.233930', 'NODE_179_length_9479_cov_3.446671', 'NODE_250_length_7033_cov_4.151287', 'NODE_104_length_5184_cov_3.857446', 'NODE_17_length_1080_cov_3.637963', 'NODE_252_length_2736_cov_3.790936', 'NODE_106_length_3199_cov_4.113473', 'NODE_182_length_6780_cov_3.785988', 'NODE_255_length_4437_cov_3.717827', 'NODE_107_length_3905_cov_4.000000', 'NODE_187_length_15556_cov_3.934237', 'NODE_256_length_13513_cov_3.902316', 'NODE_108_length_3400_cov_3.532647', 'NODE_188_length_4981_cov_3.468380', 'NODE_258_length_5201_cov_3.787349', 'NODE_110_length_5089_cov_4.254863', 'NODE_18_length_2280_cov_4.506579', 'NODE_25_length_6296_cov_4.088310', 'NODE_113_length_7208_cov_3.724334', 'NODE_190_length_16726_cov_3.944338', 'NODE_261_length_3732_cov_3.874866', 'NODE_114_length_6840_cov_4.206433', 'NODE_192_length_3773_cov_3.499603', 'NODE_263_length_2881_cov_3.549809', 'NODE_115_length_2449_cov_3.543895', 'NODE_193_length_1949_cov_3.398666', 'NODE_265_length_6579_cov_3.935857', 'NODE_116_length_5490_cov_3.986521', 'NODE_194_length_21946_cov_3.779869', 'NODE_268_length_6040_cov_3.714570', 'NODE_117_length_7944_cov_4.199647', 'NODE_195_length_8677_cov_3.738734', 'NODE_272_length_1616_cov_3.925124', 'NODE_120_length_1808_cov_3.623894', 'NODE_197_length_11657_cov_3.711161', 'NODE_276_length_2039_cov_4.203531', 'NODE_123_length_3625_cov_4.095724', 'NODE_198_length_12557_cov_4.068089', 'NODE_279_length_1545_cov_3.806473', 'NODE_124_length_15482_cov_4.010529', 'NODE_200_length_8502_cov_3.978711', 'NODE_283_length_4163_cov_3.801826', 'NODE_125_length_1085_cov_3.274654', 'NODE_201_length_22760_cov_3.767311', 'NODE_286_length_1562_cov_3.814341', 'NODE_127_length_2906_cov_3.872677', 'NODE_202_length_2151_cov_4.605300', 'NODE_292_length_10666_cov_3.957716', 'NODE_128_length_2511_cov_3.477897', 'NODE_204_length_8292_cov_3.959841', 'NODE_293_length_30040_cov_3.986119', 'NODE_130_length_7360_cov_3.967391', 'NODE_206_length_19598_cov_3.938718', 'NODE_31_length_1236_cov_3.881877', 'NODE_133_length_3213_cov_3.886399', 'NODE_210_length_13118_cov_4.030111', 'NODE_33_length_1944_cov_4.741769', 'NODE_135_length_7663_cov_3.897429', 'NODE_211_length_13153_cov_3.510226', 'NODE_35_length_1841_cov_4.451385', 'NODE_137_length_3537_cov_3.956743', 'NODE_212_length_1439_cov_3.748436', 'NODE_36_length_1245_cov_5.691566', 'NODE_138_length_6411_cov_3.954453', 'NODE_213_length_11016_cov_3.904503', 'NODE_38_length_2075_cov_3.667952', 'NODE_13_length_1039_cov_3.559191', 'NODE_215_length_13407_cov_3.911017', 'NODE_48_length_2200_cov_3.869091', 'NODE_140_length_2211_cov_4.131162', 'NODE_217_length_8074_cov_3.765791', 'NODE_49_length_1394_cov_4.038737', 'NODE_144_length_2576_cov_4.541537', 'NODE_219_length_15920_cov_3.584234', 'NODE_50_length_1812_cov_3.320088', 'NODE_145_length_5865_cov_4.159591', 'NODE_21_length_1378_cov_3.547170', 'NODE_51_length_1301_cov_4.241353', 'NODE_146_length_5343_cov_3.642710', 'NODE_220_length_2541_cov_3.456513', 'NODE_52_length_2777_cov_3.598488', 'NODE_147_length_2112_cov_4.052557', 'NODE_223_length_2622_cov_3.490465', 'NODE_53_length_3956_cov_3.967897', 'NODE_148_length_4498_cov_3.676523', 'NODE_225_length_11656_cov_4.006692', 'NODE_56_length_6765_cov_3.903326', 'NODE_150_length_21950_cov_3.830433', 'NODE_226_length_2294_cov_4.039669', 'NODE_57_length_1229_cov_5.803905', 'NODE_155_length_3398_cov_3.577693', 'NODE_228_length_16792_cov_3.992020', 'NODE_66_length_1372_cov_3.696793', 'NODE_157_length_3189_cov_3.599247', 'NODE_229_length_6655_cov_3.675282', 'NODE_67_length_1754_cov_3.173318', 'NODE_159_length_10054_cov_4.023076', 'NODE_231_length_19195_cov_3.853816', 'NODE_72_length_1116_cov_2.767025', 'NODE_160_length_9937_cov_3.930562', 'NODE_232_length_7373_cov_3.636647', 'NODE_73_length_3813_cov_3.444007', 'NODE_161_length_10937_cov_3.858279', 'NODE_233_length_6170_cov_3.877958', 'NODE_7_length_3652_cov_4.288883', 'NODE_163_length_18198_cov_3.827564', 'NODE_234_length_24373_cov_3.784352', 'NODE_80_length_5278_cov_4.184729', 'NODE_165_length_11268_cov_3.766596', 'NODE_235_length_6281_cov_4.256010', 'NODE_87_length_9569_cov_3.779078', 'NODE_166_length_4273_cov_3.088931', 'NODE_237_length_2879_cov_3.267801', 'NODE_8_length_1243_cov_5.734513', 'NODE_167_length_7819_cov_4.150531', 'NODE_238_length_25248_cov_3.834799', 'NODE_91_length_10939_cov_3.958132', 'NODE_168_length_1375_cov_3.506909', 'NODE_242_length_15521_cov_4.142452', 'NODE_93_length_1373_cov_3.665696', 'NODE_170_length_14269_cov_4.055295', 'NODE_245_length_14977_cov_4.291046', 'NODE_94_length_4752_cov_4.045455', 'NODE_176_length_4231_cov_4.103049', 'NODE_246_length_10015_cov_3.979830', 'NODE_96_length_1036_cov_4.551158', 'NODE_178_length_10768_cov_3.818722', 'NODE_248_length_32660_cov_4.043356']
data_files = glob.glob('/Users/dasha/Dropbox/GCAss_DB/sample1/' + '*.sam')
#data_files = glob.glob('/home/dasha/Boeva/data/sam_files/files_sam_chrom/' + '*7.sam')
for chrom in chrom_names:
	print 'Now for this chromosome = ', chrom
	new_file = open('/Users/dasha/Dropbox/GCAss_DB/sample1/by_contigs/'+chrom+'.sam','w')
	for file_curr in data_files:
		print 'Now for this file = ', file_curr
		flag=0
		j=0
		n=1
		for line in open(file_curr,'r'):
			j+=1
			if line[0]=='@':
				#print 'I write!!!'
				new_file.write(line)
			elif not flag:
				line1=line
				flag=1
			elif flag:
				
				line2=line
				temp1=line1.split()
				temp2=line2.split()
				
				if temp2[2]== chrom or temp1[2]== chrom:
					new_file.write(line1)
					new_file.write(line2)
				flag=0
			if j == 100000*n:
				print 'prossed ', 100000*n,' lines'
				n+=1	
	new_file.close()			
print chrom_names
