'''Splits a FASTQ file by information in the header. By default, tries to match the
regex for EG #12/1 or #4/1 or #5/2
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

import re, os
from fastq import *

def main(infile, outprefix, mid_list):
	outhandles = {}
	out_nums = {}
	with open(infile, 'r') as infh:
		for mid in mid_list:
			outfile = '%s.%d.fq' % (outprefix, mid)
			try:
				outhandles[mid] = open(outfile, 'w')
				out_nums[mid] = 0
			except IOError as err:
				raise
				
		query = re.compile(r'#\d+/\d{1}$')
		for header, sequence, quality in fastq_iterator(infh):
			read = FastqRecord(header, sequence, quality)
			match = query.findall(read.get_header())
			if not match:
				print('[INFO]: MID value not found in "%s"' % header)
				continue
			m = match[0].split('/') # EG m = ['#3', '1']
			num = int(m[0][1:])
			if num in mid_list:
				read.write_to_file(outhandles[num])
				out_nums[num] += 1
				
	for k, v in outhandles.items():
		print('[INFO]: Sequences with MID %d: %d' % (k, out_nums[k]))
		v.close()
		if out_nums[k] == 0:
			os.unlink('%s.%d.fq' % (outprefix, k))