import sys
import os
import logging
import pysam

logger = logging.getLogger('main_logger')

####################################################################
# Fragment class
####################################################################
class Fragment(object):
	__slots__ = ('first_read_chr', 'second_read_chr', 'begin', 'end', 'length', 'middle', 'direction', 'name','unique_flag','mapp_qul_flag')

	# Constructor from all fields values
	def all_arg(self, first_read_chr, second_read_chr, begin, end, length, middle, direction, name):
		self.first_read_chr = first_read_chr # First read chromosome
		self.second_read_chr = second_read_chr # Second read chromosome
		self.begin = begin # Fragment begin
		self.end = end # Fragment end
		self.length = length # Length
		self.middle = middle # Middle
		# Direction. Supported types: 
		# rf (reverse - forward)
		# fr (forward - reverse)
		# ff (forward - forward)
		# rr (reverse - reverse)
		# tr (translocation)
		self.direction = direction
		self.name = name
		self.unique_flag = unique_flag
		self.mapp_qul_flag = mapp_qul_flag
	def part_arg(self, first_read_chr, second_read_chr,begin,end,direction,name,read_len,unique_flag, mapp_qul_flag):
		self.first_read_chr = first_read_chr # First read chromosome
		self.second_read_chr = second_read_chr # Second read chromosome
		self.begin = begin # Fragment begin
		self.end = end # Fragment end
		self.length = self.end - self.begin
		self.middle = (begin + end + read_len) / 2
		self.direction = direction
		self.name = name
		self.unique_flag = unique_flag
		self.mapp_qul_flag = mapp_qul_flag
		
	def __str__(self):
		return 'Fragment ' + self.name + ': [' + str(self.begin) + ';' + str(self.end) + ']' + \
			   ', direction ' + self.direction + ', chromosomes (' + \
			   self.first_read_chr + ',' + self.second_read_chr +' unique_flag ='+str(self.unique_flag)+ ','\
			   + 'mapp_qul_flag = '+str(self.mapp_qul_flag)+ ')'

	# Constructor from 2 lines of .sam file corresponding to 2 reads
	def from_reads(self, read1, read2, sam_in):
		self.first_read_chr = sam_in.getrname(read1.tid)
		self.second_read_chr = sam_in.getrname(read2.tid)
		self.middle = (read1.pos + read2.pos + read1.qlen) / 2
		self.name = read1.qname

		if read1.is_reverse and read2.is_reverse:
			self.direction = 'rr'
		elif (not read1.is_reverse and not read2.is_reverse):
			self.direction = 'ff'
		else:
			if read1.tlen > 0:
				(left_read, right_read) = (read1, read2)
			else:
				(left_read, right_read) = (read2, read1)

			if not left_read.is_reverse and right_read.is_reverse:
				self.direction = 'fr'
			else:
				self.direction = 'rf'

		if read1.pos > read2.pos:
			self.begin = read2.pos
			self.end = read1.pos + read1.qlen
		else:
			self.begin = read1.pos
			self.end = read2.pos + read2.qlen
		self.length = self.end - self.begin

	def is_abnormal(self, smallest_normal, biggest_normal):
		return self.length < smallest_normal or self.length > biggest_normal

	#@classmethod
	def from_string(self, str):
		(name, first_read_chr, second_read_chr, begin, end, length, middle, direction) = str.rstrip('\r\n').split(';')
		self.name = name
		self.first_read_chr = first_read_chr
		self.second_read_chr = second_read_chr
		self.begin = int(begin)
		self.end = int(end)
		self.length = int(length)
		self.middle = int(middle)
		self.direction = direction
		#return cls(first_read_chr, second_read_chr, int(begin), int(end), int(length), int(middle), direction)

####################################################################
# Load .sam file into array of fragments (pysam-based)
####################################################################
def load_sam_file(sam_file, max_lines_to_read=sys.maxint):
	num_unmapped_pairs = 0
	num_translocations = 0
	num_bad_quality = 0
	reads_processed = 0
	repeats = 0
	k=0
	logger.info('Reading file, creating a Read class from every line...')
	(name, extension) = os.path.splitext(sam_file)
	if extension == '.sam':
		sam_in = pysam.Samfile(sam_file, 'r')
	elif extension == '.bam':
		sam_in = pysam.Samfile(sam_file, 'rb')
	else:
		logger.info('Unsupported extension in load_sam_file')
		return
	
	#mate_reads = set()	
	short_keys =set()
	mate_reads = dict()
	# observed_extended_keys = set()
	fragments = []
	short_rf = 0
	short_fr = 0
	for read in sam_in.fetch():

		if reads_processed > max_lines_to_read:
			break
		reads_processed += 1
		if reads_processed % 1000000 == 0:
			logger.info('Reads processed: ' + str(reads_processed))

		# We use tuple (read 1 tid, read 1 pos, read 2 tid, read 2 pos) as a key to combine 2 reads
		# If a read with this key not observed yet, save it to a dictionaty
		# Otherwise combine 2 reads and remove key from dictionary
		unique_flag=0
		mapp_qul_flag =0
		if read.is_unmapped or read.mate_is_unmapped:
			num_unmapped_pairs += 1
			continue
		if read.opt('AM'):
			if int(read.opt('AM'))>= 20:
				mapp_qul_flag = 1
		if read.opt('XT')=='U':
			unique_flag = 1
		if read.is_read1:
			short_key = (read.tid, read.pos, read.rnext, read.pnext)
			key = (read.qname,read.tid, read.pos, read.rnext, read.pnext)
		else:
			short_key = (read.rnext, read.pnext, read.tid, read.pos)
			key = (read.qname,read.rnext, read.pnext, read.tid, read.pos)

		# Protection against wrong duplicates processing
		# More correct way, but causes 2x slow down

		# extended_key = (key, read.is_read1)
		# if extended_key in observed_extended_keys:
		# 	continue
		# else:
		# 	observed_extended_keys.add(extended_key)

		if key not in mate_reads:
			#mate_reads.add(key)
			mate_reads[key] = [unique_flag, mapp_qul_flag]
		else:
		#	mate = mate_reads[key]
			# Protection against wrong duplicates processing (simplified)
			#if read.is_read1 == mate.is_read1:
			#	continue
			if short_key in short_keys:
				continue
			else:
				short_keys.add(short_key)
			# Pairs that should be dropped
			
			unique_flag = mate_reads[key][0]*unique_flag
			mapp_qul_flag = mate_reads[key][1]*mapp_qul_flag
			if read.tid != read.rnext:
				num_translocations += 1
			# Create fragment from pair
			#frag = Fragment()
			direction = ''
			if read.pnext<read.pos:
				begin = read.pnext
				end = read.pos+read.qlen
				if read.mate_is_reverse:
					direction+='r'
				else:
					direction+='f'
				if read.is_reverse:
					direction+='r'
				else:
					direction+='f'
			else:
				begin = read.pos
				end = read.pnext+read.qlen
				if read.is_reverse:
					direction+='r'
				else:
					direction+='f'
				if read.mate_is_reverse:
					direction+='r'
				else:
					direction+='f'
			frag = Fragment()
			if unique_flag or mapp_qul_flag:
				frag.part_arg(sam_in.getrname(read.tid), sam_in.getrname(read.rnext), begin,end,direction,read.qname,read.qlen,unique_flag, mapp_qul_flag)
			#if read.is_read1:

			#	frag = Fragment()
			#	frag.from_reads(read, mate, sam_in)
			#else:
			#	frag.from_reads(mate, read, sam_in)
				
				fragments.append(frag)
			else:
				continue
				#if read.is_unmapped or read.mate_is_unmapped:
				#	print 'yes it is Unmapped!'
			#del mate_reads[key]
	sam_in.close()
	logger.info('Finished reading ' + sam_file)
	logger.info('In mate_reads dictionary: ' + str(len(mate_reads)))
	logger.info('Unmapped pairs: ' + str(num_unmapped_pairs))
	logger.info('Translocations: ' + str(num_translocations))
	logger.info('Repeats : ' + str(repeats))
	logger.info('Bad quality fragments : ' + str(num_bad_quality))
	return fragments
