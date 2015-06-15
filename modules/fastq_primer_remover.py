'''This script will remove primer sequences from the beginning and end of a sequence'''

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

import re
from fastq import *

def parse_primerfile_to_regex(primerfh):
	ambiguity_codes = {
	'M':'(A|C)', 'R':'(A|G)', 'W':'(A|T)',
	'S':'(C|G)', 'Y':'(C|T)', 'K':'(G|T)',
	'V':'(A|C|G)', 'H':'(A|C|T)', 'D':'(A|G|T)',
	'B':'(C|G|T)', 'N':'(A|C|G|T)'
	}
	primers_regex = []
	for line in primerfh:
		line = line.rstrip()
		split_line = line.split('\t')
		length = len(split_line)
		if length == 1:
			if split_line[0] != '':
				reg = split_line[0]
				for key, value in ambiguity_codes.items():
					reg = reg.replace(key, value)
			primers_regex.append((reg, split_line[0]))
		elif length == 2:
			if split_line[0] != '':
				reg = split_line[0]
				for key, value in ambiguity_codes.items():
					reg = reg.replace(key, value)
				primers_regex.append((reg, split_line[0]))
			if split_line[1] != '':
				reg = split_line[1]
				for key, value in ambiguity_codes.items():
					reg = reg.replace(key, value)
				primers_regex.append((reg, split_line[1]))
		else:
			print('[WARN]: Unable to parse primer sequence from "%s"' % line)
	return primers_regex
	
def main(infile, outprefix, primerfile):
	outfile = outprefix + '.trim.fq'
	with open(infile, 'r') as infh, open(outfile, 'w') as outfh, open(primerfile, 'r') as primerfh:
		primers = parse_primerfile_to_regex(primerfh)
		primers_count = {}
		for header, sequence, quality in fastq_iterator(infh):
			record = FastqRecord(header, sequence, quality)
			length = record.get_sequence_length()
			for p in primers:
				primer = p[0]
				original = p[1]
				match = re.search(primer, sequence)
				if match:
					dist_from_end = length - match.end()
					if dist_from_end+1 < match.start():
						if dist_from_end <= 4:
							record.remove_bases(match.start(), length-1)
							primers_count[original] = primers_count.get(original, 0) + 1
							break
					else:
						if match.start() <= 4:
							record.remove_bases(0, match.end())
							primers_count[original] = primers_count.get(original, 0) + 1
							break
			if record.get_sequence_length() > 0:
				record.write_to_file(outfh)
			
	for key, value in primers_count.items():
		print('[INFO]: "%s" removed from %d sequences in %s' % (key, value, infile))
