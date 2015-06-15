#! /software/bin/python3

'''Script written for removing N's from a FASTQ file. Will get the location of N's in a sequence. If multiple, discards read. If 1, trims from closest end.
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

import sys
sys.path.append('/nfs/users/nfs_s/sw10/QUASR6/modules/')
try:
	import fastq
except ImportError:
	print("[ERROR]: Path to QUASR modules not set in script. See 'docs/INSTALL' for more info")
	sys.exit(1)

if sys.version_info < (3,0):
	print("[ERROR]: QUASR requires Python3 to run. Please read 'docs/INSTALL' for more info")
	sys.exit(1)

if len(sys.argv) != 3:
	print('[USAGE]: %s infile.fastq outprefix')
	sys.exit()

infile = sys.argv[1]
outprefix = sys.argv[2]

outfile = outprefix + '.Ntrimmed.fastq'
passed_seqs = 0
failed_seqs = 0
total_seqs = 0
unparseable_seqs = 0
with open(infile, 'r') as infh, open(outfile, 'w') as outfh:
	print('[INFO]: Parsing %s' % infile)
	for header, sequence, quality in fastq.fastq_iterator(infh):
		total_seqs += 1
		try:
			record = fastq.FastqRecord(header, sequence, quality)
		except IOError as err:
			print('[ERROR]: Ignoring "%s": %s' % (header, err))
			unparseable_seqs += 1
			continue
		
		if 'N' in record.sequence:
			if record.sequence.count('N') > 1:
				failed_seqs += 1
				continue
			else:
				index = record.sequence.find('N')
				# If N is towards end then store beginning of read
				if len(record.sequence[:index]) > len(record.sequence[index+1:]):
					record.write_to_file(outfh, end=index)
					passed_seqs += 1
				# Or if N is towards beginning or in middle of sequence
				elif len(record.sequence[:index]) <= len(record.sequence[index+1:]):
					record.write_to_file(outfh, start=index+1)
					passed_seqs += 1
		else:
			record.write_to_file(outfh)
			passed_seqs += 1
	
	assert total_seqs == passed_seqs + failed_seqs + unparseable_seqs, '[ERROR]: Unparsed sequence in file'			
	print('[TOTAL]: Input sequences: %d' % total_seqs)
	print('[PASS]: Sequences written to "%s": %d' % (outfile, passed_seqs))
	print('[FAIL]: Sequences with >1 N: %d' % failed_seqs)
	print('[FAIL]: Unparseable sequences: %d' % unparseable_seqs)