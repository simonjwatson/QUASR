'''Represents a pileup file. Takes in a line and splits out the columns
into separate instance variables to allow for easier manipulation
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
from fasta import *

class PileupFile:
	
	AsciiTable = { 	'!': 33, '"': 34, '#': 35, '$': 36,
					'%': 37, '&': 38, "'": 39, '(': 40,
					')': 41, '*': 42, '+': 43, ',': 44,
					'-': 45, '.': 46, '/': 47, '0': 48,
					'1': 49, '2': 50, '3': 51, '4': 52,
					'5': 53, '6': 54, '7': 55, '8': 56,
					'9': 57, ':': 58, ';': 59, '<': 60,
					'=': 61, '>': 62, '?': 63, '@': 64,
					'A': 65, 'B': 66, 'C': 67, 'D': 68,
					'E': 69, 'F': 70, 'G': 71, 'H': 72,
					'I': 73, 'J': 74, 'K': 75, 'L': 76,
					'M': 77, 'N': 78, 'O': 79, 'P': 80,
					'Q': 81, 'R': 82, 'S': 83, 'T': 84,
					'U': 85, 'V': 86, 'W': 87, 'X': 88,
					'Y': 89, 'Z': 90, '[': 91, '\\':92,
					']': 93, '^': 94, '_': 95, '`': 96,
					'a': 97, 'b': 98, 'c': 99, 'd':100,
					'e':101, 'f':102, 'g':103, 'h':104,
					'i':105, 'j':106, 'k':107, 'l':108,
					'm':109, 'n':110, 'o':111, 'p':112,
					'q':113, 'r':114, 's':115, 't':116,
					'u':117, 'v':118, 'w':119, 'x':120,
					'y':121, 'z':122, '{':123, '|':124,
				 	'}':125, '~':126 }
	
	@classmethod
	def _convert_ascii_to_phred(cls, qualities, ascii_offset=33):
		return [PileupFile.AsciiTable[q]-ascii_offset for q in qualities]
	
	def __init__(self, pileup_fh, reference_fh=None):
		self._reference_name = []
		self._reference_base = []
		self._consensus_base = []
		self._consensus_quality = []
		self._snp_quality = []
		self._mapping_quality = []
		self._read_depth = []
		self._read_bases = []
		self._read_qualities = []
		self._reference_sizes = {} # this contains reference segment sizes only for comparison
		self._reference_seqs = {}
		
		if reference_fh is not None:
			for header, sequence in fasta_iterator(reference_fh):
				self._reference_sizes[header] = len(sequence)
				self._reference_seqs[header] = sequence
		
		count = 0
		segment = None
		for line in pileup_fh:
			if len(line) == 0:
				continue
			split_line = line.split()
			if split_line[2] == '*':
				continue
			if split_line[0] != segment:
				if reference_fh is not None and segment is not None:
					if count != self._reference_sizes[segment]:
						assert self._reference_sizes[segment] > count, 'Internal counter (%d) greater than reference position counter (%d) for %s' % (count, ref_pos, segment)
						while count < self._reference_sizes[segment]:
							count += 1
							self._reference_name.append(segment)
							self._reference_base.append(self._reference_seqs[segment][count-1])
							self._consensus_base.append('-')
							self._consensus_quality.append(0)
							self._snp_quality.append(0)
							self._mapping_quality.append(0)
							self._read_depth.append(0)
							self._read_bases.append('')
							self._read_qualities.append('')
				count = 1
			else:
				count += 1
			segment = split_line[0]
			ref_pos = int(split_line[1])
			if ref_pos != count:
				assert ref_pos > count, 'Internal counter (%d) greater than reference position counter (%d) for %s' % (count, ref_pos, segment)
				while count < ref_pos:
					self._reference_name.append(segment)
					self._reference_base.append(self._reference_seqs[segment][count-1])
					self._consensus_base.append('-')
					self._consensus_quality.append(0)
					self._snp_quality.append(0)
					self._mapping_quality.append(0)
					self._read_depth.append(0)
					self._read_bases.append('')
					self._read_qualities.append('')
					count += 1

			self._reference_name.append(split_line[0])
			self._reference_base.append(split_line[2].upper())
			columns = len(split_line)
			if columns == 10 or columns == 11: # -c flag or -c/-s flags
				self._consensus_base.append(split_line[3])
				self._consensus_quality.append(split_line[4]) # Phred-scaled
				self._snp_quality.append(split_line[5]) # Phred-scaled probability of consensus being identical to the reference
				self._mapping_quality.append(split_line[6]) # RMS mapping quality
				self._read_depth.append(int(split_line[7]))
				self._read_bases.append(split_line[8])
				self._read_qualities.append(split_line[9])
				# Ignore last column which is mapping quality per base
			elif columns == 6 or columns == 7: # -s flag or neither
				self._read_depth.append(int(split_line[3]))
				self._read_bases.append(split_line[4])
				self._read_qualities.append(split_line[5])
			elif columns == 4 and int(split_line[3]) == 0:
				self._read_depth.append(0)
				self._read_bases.append('')
				self._read_qualities.append('')
			else:
				raise IOError("Unidentifiable columns in pileup file: line %d" % count)
			
		if reference_fh is not None and count != self._reference_sizes[segment]:
			assert self._reference_sizes[segment] > count, 'Internal counter (%d) greater than reference position counter (%d) for %s' % (count, ref_pos, segment)
			while count < self._reference_sizes[segment]:
				self._reference_name.append(split_line[0])
				self._reference_base.append(self._reference_seqs[split_line[0]][count-1])
				self._consensus_base.append('-')
				self._consensus_quality.append(0)
				self._snp_quality.append(0)
				self._mapping_quality.append(0)
				self._read_depth.append(0)
				self._read_bases.append('')
				self._read_qualities.append('')
				count += 1
	
	def parse_read_bases(self, phred_cutoff=0, ascii_offset=33):
		import re
		transformed_bases = []
		for i in range(len(self._read_bases)):
			reads = self._read_bases[i]
			if reads == '':
				transformed_bases.append(self._reference_base[i])
				continue
			if '*' in reads:
				reads = re.sub('\*', '-', reads)
			if '^' in reads:
				reads = re.sub('\^.', '', reads) # remove read start token
			if '$' in reads:
				reads = reads.replace('$', '') # remove read end token
			for symbol in '.,':
				reads = reads.replace(symbol, self._reference_base[i])
			for match in re.finditer('\d+', reads): # remove the indels eg +1A, -2gg, +2at etc.
				regex = re.compile('(\+|\-)' + match.group() + '(A|C|G|T|N|-|Y|R|M|W|S|K|H|D|B|V|X){' + match.group() + '}', re.IGNORECASE)
				reads = re.sub(regex, '', reads)
			phred = [p for p in PileupFile._convert_ascii_to_phred(self._read_qualities[i], ascii_offset)]
			str_length = len(reads)
			assert str_length == len(phred), 'Unequal base and quality lengths at position %d:\n%s (%d)\n%s (%d)' % (i+1, reads, str_length, self._read_qualities[i], len(phred))
			passed = [reads[j] for j in range(str_length) if phred[j] >= phred_cutoff]
			transformed_bases.append(''.join(passed))
		
		return transformed_bases	
	
	def calc_segment_sizes(self):
		seg_counter = []
		prev_seg = None
		count = 0
		for r in range(0, len(self._reference_name)):
			count += 1
			if self._reference_name[r] != prev_seg:
				if prev_seg is None:
					prev_seg = self._reference_name[r]
					continue
				seg_counter.append((prev_seg, count))
				count = 0
				prev_seg = self._reference_name[r]
		seg_counter.append((self._reference_name[-1], count))
		return seg_counter
				
	def return_reference_names(self):
		return self._reference_name
	
	def return_read_depths(self, phred_cutoff=0, ascii_offset=33):
		passed_phred_lengths = []
		for pos in self._read_qualities:
			passed_phred = [p for p in PileupFile._convert_ascii_to_phred(pos, ascii_offset) if p >= phred_cutoff]
			passed_phred_lengths.append(len(passed_phred))
		return passed_phred_lengths		
