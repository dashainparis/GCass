#input_f='/home/dasha/Boeva/bwa/bwa-0.6.1/hg18.fa'
input_f='/home/dasha/PhD/Boeva/hg18.fa'
flag=0
name=''
arr=[]
g=open('result','a')
num_N=0
threshold=5
for line in open(input_f):
	if 'chr' in line:
		name=''
		for i in line:
			if i!='>' and i!='\n':
				name+=i
		arr=[]
		arr.append(name)
		pos_in_chr=0
		num_N=0
	else:
		for loc_pos in line:
			if loc_pos=='N' or 'n':
				pos_in_chr+=1
				if num_N==0:
					pos_in_gen_beg=pos_in_chr 
				num_N+=1
			elif loc_pos!='N' and num_N==0:
				pos_in_chr+=1
			elif loc_pos=='\n':	
				hhh=1	
			elif loc_pos!='N' and num_N!=0:
				#print num_N
				if num_N>=threshold:
					pos_in_gen_end=pos_in_chr
					arr.append(pos_in_gen_beg)
					arr.append(pos_in_gen_end)
					if arr:
						print arr
					num_N=0
					del arr[1:]
					pos_in_chr+=1
				else: 
					num_N=0
					pos_in_gen_end=0
					pos_in_gen_beg=0

