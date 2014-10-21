f = open('/Users/dasha/Dropbox/GCAss/sample1/contigs23.fasta')
flag=0
for line in f:
	
	if '>' in line:
		print line
		if flag:
			new.close()
		new = open('/Users/dasha/Dropbox/GCAss/sample1/contigs/'+line[1:-1]+'.fa','w')
		flag = 1
	else:
		
		new.write(line)
new.close()