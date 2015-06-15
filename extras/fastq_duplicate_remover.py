#! /software/bin/python3

'''Takes in a FASTQ file and stores each record in a dictionary with the sequence 
as the key to output only unique sequences
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

import sys, os.path
sys.path.append('/nfs/users/nfs_s/sw10/QUASR6/modules/')
try:
	from fastq import *
except ImportError:
	print("[ERROR]: Path to QUASR modules not set in script. See 'docs/INSTALL' for more info")
	sys.exit(1)

if sys.version_info < (3,0):
	print("[ERROR]: QUASR requires Python3 to run. Please read 'docs/INSTALL' for more info")
	sys.exit(1)

def main(infile, outprefix, rev_file=None, offset=33):
	outfile = '%s.unique.f.fq' % outprefix
	try:
		infh = open(infile, 'r')
		if rev_file is not None:
			outfile_r = '%s.unique.r.fq' % outprefix
			revfh = open(rev_file, 'r')
			outfh_r = open(outfile_r, 'w')
		outfh = open(outfile, 'w')
	except IOError as e:
		raise
	else:
		print('[INFO]: Input and output files successfully opened')
		unique_seqs = {}
		duplicates = 0
		skipped = 0
	
	print('[INFO]: Parsing FASTQ records from "%s"' % infile)
	forward_reads = {header: FastqRecord(header, sequence, quality) for header, sequence, quality in fastq_iterator(infh)}
	infh.close()
	total = len(forward_reads)
	if rev_file is not None:
		print('[INFO]: Parsing reverse reads from "%s"' % rev_file)
		reverse_reads = {header: FastqRecord(header, sequence, quality) for header, sequence, quality in fastq_iterator(revfh)}
		revfh.close()
		assert len(reverse_reads) == total, "%d reads in forward file and %d in reverse file" % (total, len(reverse_reads))
	
	for header, record in forward_reads.items():
		f_mean_qual = record.calculate_mean_quality(ascii_offset=offset)
		f_sequence = record.get_sequence()
		if rev_file is not None:
			# PE duplicate check
			if not header in reverse_reads.keys():
				reverse_header = '%s2' % header[:-1]
				if not reverse_header in reverse_reads.keys():
					print('[ERROR]: Ignoring read "%s"; no matching reverse read' % header)
					skipped += 1
					continue
			else:
				reverse_header = header
			r_record = reverse_reads[reverse_header]
			comb_mean_qual = (f_mean_qual + r_record.calculate_mean_quality(ascii_offset=offset))/2
			comb_sequence = f_sequence + r_record.get_sequence()
			if comb_sequence not in unique_seqs.keys():
				unique_seqs[comb_sequence] = (record, r_record, comb_mean_qual)
				continue
			else:
				duplicates += 1
				if comb_mean_qual <= unique_seqs[comb_sequence][2]:
					continue
				else:
					unique_seqs[comb_sequence] = (record, r_record, comb_mean_qual)
			
		else:
			# SE duplicate check
			if f_sequence not in unique_seqs.keys():
				unique_seqs[f_sequence] = (record, f_mean_qual)
				continue
			else:
				duplicates += 1
				if f_mean_qual <= unique_seqs[f_sequence][1]:
					continue
				else:
					unique_seqs[f_sequence] = (record, f_mean_qual)
	
	if rev_file is not None:
		for f_fastq, r_fastq, ignore in unique_seqs.values():
			f_fastq.write_to_file(outfh)
			r_fastq.write_to_file(outfh_r)
		outfh_r.close()
	else:
		for fastq, ignore in unique_seqs.values():
			fastq.write_to_file(outfh)
	outfh.close()
	
	uniques = len(unique_seqs)
	assert	(duplicates + uniques + skipped == total), 'Unparsed reads in file!'
	print('[STATS]: %d inital sequences' % total)
	print('[STATS]: %d duplicate sequenes' % duplicates)
	print('[STATS]: %d unique sequences' % len(unique_seqs))

if __name__ == "__main__":
	import getopt
	
	prog = sys.argv[0]
	examples = '''
[EG]: %s -f input.fastq -o outDir/outPrefix
[EG]: %s -f input.fastq -r reverse.fastq -o outDir/outPrefix
[EG]: %s -f input.fastq -r reverse.fastq -o outDir/outPrefix -i''' % (prog, prog, prog)

	usage = '''[USAGE]: %s <options>
-h   [None]\tDisplay this usage message with examples (--help)
-f * [File]\tInput FASTQ file (--forward)
-r   [File]\tFASTQ file containing reverse reads if PE sequenced (--reverse)
-o * [String]\tOutput directory and file prefix (--outprefix)
-i   [None]\tIllumina ASCII offset (+64) used to encode quality (--illumina)''' % prog

	try:
		opts, args = getopt.getopt(sys.argv[1:], "hf:r:o:i", ["help", "forward=", "reverse=", "outprefix=", "illumina"])
	except getopt.GetoptError as err:
	        # print help information and exit:
		print(str(err)) # will print something like "option -a not recognized"
		print(usage)
		sys.exit(2)
	
	infile = None
	reverse = None
	outprefix = None
	offset = 33
	
	for o, a in opts:
		if o in ("-h", "--help"):
			print(usage)
			print(examples)
			sys.exit()
		elif o in ("-f", "--forward"):
			infile = a
		elif o in ("-r", "--reverse"):
			reverse = a
		elif o in ("-o", "--outprefix"):
			outprefix = a
		elif o in ("-i", "--illumina"):
			offset = 64
	
	if len(sys.argv) == 1:
		print(usage)
		sys.exit()
	
	if infile is None:
		print('[ERROR]: Input FASTQ file must be specified with the "-f" flag')
		sys.exit(2)
	elif outprefix is None:
		print('[ERROR]: Output directory and file prefix must be specified with the "-o" flag')
		sys.exit(2)
	
	try:
		main(infile, outprefix, reverse, offset)
	except IOError as e:
		print('[ERROR]: %s' % e)
	except AssertionError as e:
		print('[ERROR]: %s' % e)
	else:
		print('[INFO]: Unique sequences written to "%s.unique.f.fq"' % outprefix)
		if reverse is not None:
			print('[INFO]: Unique reverse reads written to "%s.unique.r.fq"' % outprefix)
