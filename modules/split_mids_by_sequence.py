'''Splits a FASTQ file by MID sequence at the beginning of each read. The tag must match exactly
and must start at the first position. 454 Rapid Librray tags are hard-coded in, but any others must be
provided as a file in the format "Number\tSequence". EG:
1	ACGTCGATCA
2	ACGATACGAC
'''

# Copyright 2010, 2011 Simon Watson
#
# This file is part of QUASR.
#
# QUASR is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# QUASR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with QUASR.  If not, see <http://www.gnu.org/licenses/>.

import os
from fastq import *

def custom_tags_to_dict(infh):
	custom_dict = {}
	for line in custom_fh:
		line = line.rstrip("\r\n")
		split_line = line.split('\t')
		if len(split_line) != 2:
			raise RuntimeError('Incorrect formatting in line:\n[ERROR]: "%s"\n[ERROR]: Must be formatted as "MID num \\t MID sequence"' % line)
		try:
			mid_num = int(split_line[0])
		except ValueError:
			print('[ERROR]: Unable to convert "%s" into a MID number')
			raise
		mid_seq = split_line[1]
		if mid_seq in custom_dict.values():
			raise RuntimeError('MID sequence "%s" is present more than once')
		custom_dict[mid_num] = mid_seq
	return custom_dict

def main(infile, outprefix, mid_list, customfile=None):
	outhandles = {}
	out_nums = {}
	with open(infile, 'r') as infh:
		
		if customfile is not None:
			try:
				customfh = open(customfile, 'r')
			except IOError as err:
				print('[ERROR]: %s' % err)
				sys.exit(2)
			else:
				tags = custom_tags_to_dict(customfh)
				customfh.close()
		else:
			# Hard-coded 454 Rapid library MID tags
			tags = {
			1: 'ACACGACGACT', 2: 'ACACGTAGTAT', 3: 'ACACTACTCGT', 4: 'ACGACACGTAT',
			5: 'ACGAGTAGACT', 6: 'ACGCGTCTAGT' ,7: 'ACGTACACACT', 8: 'ACGTACTGTGT',
			9: 'ACGTAGATCGT', 10:'ACTACGTCTCT', 11:'ACTATACGAGT', 12:'ACTCGCGTCGT',
			13:'AGACTCGACGT', 14:'AGTACGAGAGT', 15:'AGTACTACTAT', 16:'AGTAGACGTCT',
			17:'AGTCGTACACT', 18:'AGTGTAGTAGT', 19:'ATAGTATACGT', 20:'CAGTACGTACT',
			21:'CGACGACGCGT', 22:'CGACGAGTACT', 23:'CGATACTACGT', 24:'CGTACGTCGAT',
			25:'CTACTCGTAGT', 26:'GTACAGTACGT', 27:'GTCGTACGTAT', 28:'GTGTACGACGT',
			29:'ACACAGTGAGT', 30:'ACACTCATACT', 31:'ACAGACAGCGT', 32:'ACAGACTATAT',
			33:'ACAGAGACTCT', 34:'ACAGCTCGTGT', 35:'ACAGTGTCGAT', 36:'ACGAGCGCGCT',
			37:'ACGATGAGTGT', 38:'ACGCGAGAGAT', 39:'ACGCTCTCTCT', 40:'ACGTCGCTGAT',
			41:'ACGTCTAGCAT', 42:'ACTAGTGATAT', 43:'ACTCACACTGT', 44:'ACTCACTAGCT',
			45:'ACTCTATATAT', 46:'ACTGATCTCGT', 47:'ACTGCTGTACT', 48:'ACTGTAGCGCT'
			}
		# "tags" now contains the MID and sequence as a dictionary
		# Open the output filehandles. 0 is for sequences which can't be assigned a MID
		try:
			for mid in mid_list:
				outfile = '%s.%d.fq' % (outprefix, mid)
				outhandles[mid] = open(outfile, 'w')
				out_nums[mid] = 0
		except IOError as err:
			raise
				
		# Loop through the FASTQ file and check if the sequence starts with the tag
		# If it doesn't match, it is assigned to 0.
		for header, sequence, quality in fastq_iterator(infh):
			read = FastqRecord(header, sequence, quality)
			for mid in mid_list:
				if sequence.startswith(tags[mid]):
					read.write_to_file(outhandles[mid], start=len(tags[mid]))
					out_nums[mid] += 1
					continue
			
	for k, v in outhandles.items():
		print('[INFO]: Sequences with MID %d: %d' % (k, out_nums[k]))
		v.close()
		if out_nums[k] == 0:
			os.unlink('%s.%d.fq' % (outprefix, k))