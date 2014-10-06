f = open('/Users/dasha/PhD/data/chrom_files_hg38/chr10.fa','r')
median = 5000
contig = ''
for i in f:
	contig+=i[:-1]
f.close()
contig = contig.lower()
print contig[1000000:1000100]
windows = []
wind_init = contig[:median-1]
gc_count_init = wind_init.count('c')+wind_init.count('g')
gc_pers_init = float(gc_count_init)/len(wind_init)
ind_init = 0
summ = 0
contig = contig[median:]
k=1000000
length= len(contig)
f = open('regions.txt','w')
print 'gc_pers_init = ',gc_pers_init
while len(contig)>length*0.95:
	gc_pers_init = float(gc_count_init)/median
	summ += len(wind_init)
	if  summ/k>0:
		print 'Alredy did ',summ
		k+=1000000
	if (wind_init.count('n'))/float(len(wind_init))>=0.5:
		print 'skip'
		wind_init = contig[:median]
		contig = contig[len(wind_init):]
		gc_count_init = wind_init.count('c')+wind_init.count('g')
		gc_pers_init = float(gc_count_init)/len(wind_init)
	else:
		print 'I came here'
		wind_curr = contig[0:median]
		print len(wind_init)
		print len(wind_curr)
		if wind_curr==wind_init:
			print 'Whaaaatttt???'
		gc_count_curr = wind_curr.count('g')+wind_curr.count('c')
		gc_pers_curr = float(gc_count_curr)/(len(wind_init))
		#gc_pers_curr = float(gc_count_curr+gc_count_init)/(len(wind_init)+median)
		print 'gc_pers_curr = ',gc_pers_curr
		print 'gc_pers_init = ',gc_pers_init
		if abs(gc_pers_curr - gc_pers_init)>0.1:
			print 'I am here!'
			windows.append([wind_init,gc_pers_init])
			if len(wind_init)>5000:
				f.write(str(len(wind_init))+' '+str(gc_pers_init)+'\n')
				f.write(wind_init+(wind_curr)+'\n')
			contig = contig[len(wind_curr):]
			wind_init = (wind_curr)
			gc_count_init = gc_count_curr
		else:
			print 'I was here!'
			contig = contig[len(wind_curr):]
			wind_init+=(wind_curr)
			gc_count_init= gc_count_curr
			#print 'len(wind_init) = ',len(wind_init)
			
			
for i in windows:
	print 'wind len = ',len(i[0]), 'gc  = ', i[1]

f.close()
