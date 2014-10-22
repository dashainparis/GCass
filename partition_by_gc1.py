#f = open('/Users/dasha/PhD/data/chrom_files_hg38/chr10.fa','r')
f= open('/Users/dasha/Dropbox/GCAss_DB/sample1/ref_seq.fa','r')
median = 483
contig = ''
j=0
for i in f:
	j+=1
	#if j>200000:
	#	break
	contig+=i[:-1]
f.close()
contig = contig.lower()
wind_gc= []
wind_gc_pos =[]
for i in xrange(0,len(contig),median):
	#wind_gc.append(round((contig[i:i+median].count('g')+contig[i:i+median].count('c'))/float(median),2))
	wind_gc.append((contig[i:i+median].count('g')+contig[i:i+median].count('c'))/float(median))
#print wind_gc

curr_perc=0.0
for i in xrange(len(wind_gc)):
	if curr_perc ==0:
		curr_perc = wind_gc[i]
		wind_gc_pos.append([[i*median,(i+1)*median],wind_gc[i]])
	else:
		if abs(curr_perc-wind_gc[i])<0.01:
			wind_gc_pos[-1][0][1] = (i+1)*median
		else:
			curr_perc = wind_gc[i]
			wind_gc_pos.append([[i*median,(i+1)*median],wind_gc[i]])
#print wind_gc_pos

count = 0
for i in wind_gc_pos:
	if i[0][1]-i[0][0]>= median*2:
		count+=1
		#print i , (i[0][1]-i[0][0])/median
print 'len contig = ',len(contig)
print 'median = ',median
print 'Number of different regions = ',len(wind_gc_pos)
print 'Number of all nonoverlapping windows with lenght median = ',len(wind_gc)
print 'Number of all windows with lenght greater than 2*median = ',count
