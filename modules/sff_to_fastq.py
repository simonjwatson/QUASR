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
#
#copyright Jose Blanca and Bastien Chevreux
#COMAV institute, Universidad Politecnica de Valencia (UPV)
#Valencia, Spain

# additions to handle paired end reads by Bastien Chevreux
# bugfixes for linker specific lengths: Lionel Guy

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

import struct
import fastq as fastq_module

def read_bin_fragment(struct_def, fileh, offset=0, data=None, byte_padding=None):
	if data is None:
		data = {}

	#we read each item
	bytes_read = 0
	for item in struct_def:
		#we go to the place and read
		fileh.seek(offset + bytes_read)
		n_bytes = struct.calcsize(item[1])
		buffer = fileh.read(n_bytes)
		read = struct.unpack('>' + item[1], buffer)
		if len(read) == 1:
			read = read[0]
		data[item[0]] = read
		bytes_read += n_bytes

	#if there is byte_padding the bytes_to_read should be a multiple of the
	#byte_padding
	if byte_padding is not None:
		pad = byte_padding
		bytes_read = ((bytes_read + pad - 1) // pad) * pad

	return (bytes_read, data)

def check_magic(magic):
	# It checks that the magic number of the file matches the sff magic.
	if not magic == 779314790:
		raise RuntimeError('This file does not seems to be an sff file.')

def check_version(version):
	# It checks that the version is supported, otherwise it raises an error.
	supported = (b'\x00', b'\x00', b'\x00', b'\x01')
	if not version == supported:
		raise RuntimeError('SFF version not supported. Please contact the author of the software.')

def read_header(fileh):
	# It reads the header from the sff file and returns a dict with the information
	#first we read the first part of the header
	head_struct = [
		('magic_number', 'I'),
		('version', 'cccc'),
		('index_offset', 'Q'),
		('index_length', 'I'),
		('number_of_reads', 'I'),
		('header_length', 'H'),
		('key_length', 'H'),
		('number_of_flows_per_read', 'H'),
		('flowgram_format_code', 'B'), ]
	data = {}
	first_bytes, data = read_bin_fragment(struct_def=head_struct, fileh=fileh, offset=0, data=data)
	check_magic(data['magic_number'])
	check_version(data['version'])
	#now that we know the number_of_flows_per_read and the key_length
	#we can read the second part of the header
	struct2 = [
		('flow_chars', str(data['number_of_flows_per_read']) + 'c'),
		('key_sequence', str(data['key_length']) + 'c')
	]
	read_bin_fragment(struct_def=struct2, fileh=fileh, offset=first_bytes,data=data)
	return data

def sequences(fileh, header):
	'''It returns a generator with the data for each read.'''
	#now we can read all the sequences
	fposition = header['header_length']	   #position in the file
	reads_read = 0
	while True:
		if fposition == header['index_offset']:
			#we have to skip the index section
			fposition += index_length
			continue
		else:
			bytes_read, seq_data = read_sequence(header=header, fileh=fileh, fposition=fposition)
			yield seq_data
			fposition += bytes_read
			reads_read += 1
			if reads_read >= header['number_of_reads']: break

def read_sequence(header, fileh, fposition):
	'''It reads one read from the sff file located at the fposition and
	returns a dict with the information.'''
	header_length = header['header_length']
	index_offset = header['index_offset']
	index_length = header['index_length']

	#the sequence struct
	read_header_1 = [
		('read_header_length', 'H'),
		('name_length', 'H'),
		('number_of_bases', 'I'),
		('clip_qual_left', 'H'),
		('clip_qual_right', 'H'),
		('clip_adapter_left', 'H'),
		('clip_adapter_right', 'H'),
	]
	def read_header_2(name_length):
		'''It returns the struct definition for the second part of the header'''
		return [('name', str(name_length) +'c')]
	def read_data(number_of_bases):
		'''It returns the struct definition for the read data section.'''
		if header['flowgram_format_code'] == 1:
			flow_type = 'H'
		else:
			raise Error('file version not supported')
		number_of_bases = str(number_of_bases)
		return [
			('flowgram_values', str(header['number_of_flows_per_read']) +
																	 flow_type),
			('flow_index_per_base', number_of_bases + 'B'),
			('bases', number_of_bases + 'c'),
			('quality_scores', number_of_bases + 'B'),
		]

	data = {}
	#we read the first part of the header
	bytes_read, data = read_bin_fragment(struct_def=read_header_1,
									fileh=fileh, offset=fposition, data=data)

	read_bin_fragment(struct_def=read_header_2(data['name_length']),
						  fileh=fileh, offset=fposition + bytes_read, data=data)
	#we join the letters of the name
	data['name'] = [c.decode("utf-8") for c in data['name']]
	data['name'] = ''.join(data['name'])
	offset = data['read_header_length']
	#we read the sequence and the quality
	read_data_st = read_data(data['number_of_bases'])
	bytes_read, data = read_bin_fragment(struct_def=read_data_st, fileh=fileh, offset=fposition + offset, data=data, byte_padding=8)
	#we join the bases
	data['bases'] = [c.decode("utf-8") for c in data['bases']]
	data['bases'] = ''.join(data['bases'])

	return data['read_header_length'] + bytes_read, data

def main(sff_file, outfile):
	with open(sff_file, 'rb') as sff_fh, open(outfile, 'w') as outfh:
		print('[INFO]: Processing SFF file "%s"' % sff_file)
		header_data = read_header(fileh=sff_fh)
		key_seq = [h.decode("utf-8") for h in header_data['key_sequence']]
		key_seq = ''.join(key_seq)
		key_len = len(key_seq)
		total_seqs = 0
		num_with_key = 0
		
		for seq_data in sequences(fileh=sff_fh, header=header_data):
			header = seq_data['name']
			sequence = seq_data['bases']
			try:
				quality = fastq_module.convert_phred_to_ascii(seq_data['quality_scores'], 33)
			except IOError as err:
				print('[ERROR]: Ignoring sequence "%s": %e' % (header, err))
				
			if sequence.startswith(key_seq):
				sequence = sequence[key_len:]
				quality = quality[key_len:]
				num_with_key += 1
			fastq = fastq_module.FastqRecord(header, sequence, quality)
			fastq.write_to_file(outfh)
			total_seqs += 1
			
	print('[INFO]: %d total sequences in SFF file' % total_seqs )
	print('[INFO]: %d had key sequence "%s" removed' % (num_with_key, key_seq) )
	print('[INFO]: Sequences written to "%s"' % outfile)