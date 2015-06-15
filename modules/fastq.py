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

from fasta import *

class FastqRecord(FastaRecord):
	
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
	
	def __init__(self, header, sequence, quality):
		# check that sequence and quality strings are the same length
		if len(sequence) != len(quality):
			raise IOError('Sequence and quality strings are unequal in length')
		# check that quality is in valid ASCII format
		invalid_ascii = [qual for qual in quality if qual not in FastqRecord.AsciiTable.keys()]
		if invalid_ascii:
			raise IOError('Contains invalid qualities: ' + str(invalid_ascii))
		else:
			FastaRecord.__init__(self, header, sequence)
			self.quality = quality
	
	def _convert_ascii_to_phred(self, qualities, ascii_offset=33):
		return [FastqRecord.AsciiTable[q]-ascii_offset for q in qualities]
		
	def _calculate_median(self, phred):
		phred = sorted(phred)
		count = len(phred)
		half = count//2
		if count % 2 == 0:
			return( (phred[half] + phred[half+1]) / 2)
		else:
			return phred[half]
			
	def _calculate_mean(self, phred):
		return ( float(sum(phred)) / len(phred) )
			
	def calculate_median_quality(self, start=None, end=None, ascii_offset=33):
		phred = self._convert_ascii_to_phred(self.quality[start:end], ascii_offset)
		return self._calculate_median(phred)
	
	def calculate_mean_quality(self, start=None, end=None, ascii_offset=33):
		phred = self._convert_ascii_to_phred(self.quality[start:end], ascii_offset)
		return self._calculate_mean(phred)
		
	def remove_nth_base(self, n):
		super().remove_nth_base(n)
		self.quality = self.quality[:n-1] + self.quality[n:]
	
	def remove_bases(self, start, end):
		super().remove_bases(start, end)
		self.quality = self.quality[:start] + self.quality[end:]
	
	def return_phred_scores(self, start=None, end=None, ascii_offset=33):
		return self._convert_ascii_to_phred(self.quality[start:end], ascii_offset)
			
	def write_to_file(self, filehandle, start=None, end=None):
		filehandle.write('@%s\n%s\n+\n%s\n' % (self.header, self.sequence[start:end], self.quality[start:end]))
			
# --- END OF CLASS --- #

def convert_phred_to_ascii(phred_scores, ascii_offset=33):
	invalid_phred = [q for q in phred_scores if not 0 <= q <= 93]
	if invalid_phred:
		raise IOError('Contains an out-of-range Phred score: ' + str(invalid_phred))
	else:
		reverse_ascii_table = { v: k for (k, v) in FastqRecord.AsciiTable.items() }
		output = [reverse_ascii_table.get(phred+ascii_offset, '!') for phred in phred_scores]
		return ''.join([ascii for ascii in output])
		
def fastq_iterator(filehandle):
	#Generator taken from BioPython. Yields the header, sequence and quality for a record'''
	while True:
		line = filehandle.readline()
		if line[0] == '@': break
	while True:
		header = line[1:].rstrip()
		sequence = filehandle.readline().rstrip()
		# This next loop is if there are any extra sequence lines
		while True:
			line = filehandle.readline()
			if line[0] == '+': break
			else: sequence += line.rstrip()	
		sequence_length = len(sequence) # For quality length comparison
		quality = filehandle.readline().rstrip()
		while True:
			line = filehandle.readline()
			if not line: break
			if line[0] == '@':
				if len(quality) >= sequence_length: break
			quality += line.rstrip()
		if len(quality) != sequence_length:
			print('[ERROR]: Length mismatch in quality and base sequence for \"%s\". Record ignored.' % header)
			if not line: return
			else: continue
		yield (header, sequence, quality)
		if not line: return
