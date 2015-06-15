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

class FastaRecord:	
	
	AcceptedBases = 'ACGTUNMRWSYKVHDBN-'
	AmbiguityCodes = { 	'AC': 'M', 'AG': 'R', 'AT': 'W',
						'CG': 'S', 'CT': 'Y', 'GT': 'K',
						'ACG': 'V', 'ACT': 'H', 'AGT': 'D',
						'CGT': 'B', 'ACGT': 'N' }
		
	def __init__(self, header, sequence):
		#check that the sequence only contains valid characters
		sequence = sequence.upper()
		invalid_chars = [char for char in sequence if char not in FastaRecord.AcceptedBases]
		if invalid_chars:
			raise IOError('Contains invalid bases: ' + str(invalid_chars))
		else:
			if header.startswith('>') or header.startswith('@'):
				self.header = header[1:]
			else:
				self.header = header
			self.sequence = sequence
			
	def get_header(self):
		return self.header
	
	def get_sequence(self):
		return self.sequence
	
	def get_sequence_length(self):
		return len(self.sequence)
		
	def calculate_gc_percentage(self):
		return ( ( (self.sequence.count('G') + self.sequence.count('C')) / self.get_sequence_length() ) * 100 )
		
	def remove_nth_base(self, n):
		self.sequence = self.sequence[:n-1] + self.sequence[n:]
	
	def remove_bases(self, start, end):
		self.sequence = self.sequence[:start] + self.sequence[end:]
			
	def write_to_file(self, filehandle):
		filehandle.write('>%s\n%s\n' % (header, sequence))
		
# --- END OF CLASS --- #

def fasta_iterator(fh):
	while True:
		line = fh.readline()
		if line[0] == '>': break
	while True:
		header = line[1:-1].rstrip()
		sequence = fh.readline().rstrip()
		while True:
			line = fh.readline()
			if not line: break
			if line.startswith('>'): break
			sequence += line.rstrip()
		yield(header, sequence)
		if not line: return