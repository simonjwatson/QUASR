#! /software/bin/python3

'''This script will perform the following processing steps on the data:
1) Convert SFF to FASTQ
2) Split FASTQ by MID
3) Remove primer or BAC sequences
4) Create QA graphs
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

import sys, getopt
sys.path.append('/nfs/users/nfs_s/sw10/QUASR6/modules/')
try:
	import sff_to_fastq, qa, split_mids_by_header, split_mids_by_sequence, fastq_primer_remover
except ImportError:
	print("[ERROR]: Path to QUASR modules not set in script. See 'docs/INSTALL' for more info")
	sys.exit(1)

if sys.version_info < (3,0):
	print("[ERROR]: QUASR requires Python3 to run. Please read 'docs/INSTALL' for more info")
	sys.exit(1)

prog = sys.argv[0]
examples = '''
[EG]: %s -f reads1.fastq -o example_output
[EG]: %s -f reads1.fastq -r reads2.fastq -o example_output
[EG]: %s -f reads1.fastq -o ../example_output''' % (prog, prog, prog)

usage = '''[USAGE]: %s <options>
-h   [None]\tDisplay help message with examples (--help)

--GENERAL--
-f * [File]\tInput FASTQ or SFF file (--infile)
-o * [String]\tOutput directory and filename prefix (--outprefix)
-i   [None]\tIllumina ASCII offset (+64) used to encode quality (--illumina)

--MID-SPLITTING--
-m + [CSV]\tComma-separated MID numbers to be extracted (--mids)
-d + [None]\tExtract MIDs by parsing header (--header)
OR
-s + [None]\tExtract MIDs by parsing sequence (--sequence)
-c   [File]\tFile containing tab-separated custom MID & tag sequence (--customfile)

--PRIMER-REMOVAL--
-r   [None]\tRemove primer sequences from readset (--remove)
-l   [CSV]\tComma-separated subset of MIDs to perform primer-removal [all -m] (--primerlist)
-t + [File]\tFile containing 5' primers in 1st column and 3' in 2nd (--trimfile)

--QUALITY-ASSURANCE--
-g   [None]\tGenerate QA graphs of output readset (--graphs)
-a   [CSV]\tComma-separated subset of MIDs to create QA graphs [all -m] (--graphlist)
-p   [File]\tPath to R binary for stats generation and graphing (--path)
-e   [Integer]\tNumber of bases to graph 3' mean quality dropoff [15] (--end)

[NOTE]: Options with * are mandatory. Those with + are mandatory for that optional section.''' % prog

try:
	opts, args = getopt.getopt(sys.argv[1:], "hf:o:im:c:dsrl:t:ga:p:e:", ["help", "infile=", "outprefix=", "illumina", "mids=", "customfile=", "header", "sequence", "remove", "primerlist=", "trimfile=", "graphs", "graphlist=", "path=", "end="])
except getopt.GetoptError as err:
	print(str(err))
	print(usage)
	sys.exit(2)

infile = None
outprefix = None

split_by_header = False
split_by_sequence = False
mid_list = None
customfile = None

remove_primers = False
primer_list = None
trimfile = None

perform_qa = False
graph_list = None
r_path = None
ascii_offset = 33
end_length = 15

for o,a in opts:
	if o in ("-h", "--help"):
		print(usage)
		print(examples)
		sys.exit()
	elif o in ("-f", "--infile"):
		infile = a
	elif o in ("-o", "--outprefix"):
		outprefix = a
	elif o in ("-i", "--illumina"):
		ascii_offset = 64
	elif o in ("-m", "--mids"):
		mid_list = a
	elif o in ("-c", "--customfile"):
		customfile = a
	elif o in ("-d", "--header"):
		split_by_header = True
	elif o in ("-s", "--sequence"):
		split_by_sequence = True
	elif o in ("-r", "--remove"):
		remove_primers = True
	elif o in ("-l", "--primerlist"):
		primer_list = a
	elif o in ("-t", "--trimfile"):
		trimfile = a
	elif o in ("-g", "--graphs"):
		perform_qa = True
	elif o in ("-a", "--graphlist"):
		graph_list = a
	elif o in ("-e", "--end"):
		end_length = int(a)
	elif o in ("-p", "--path"):
		r_path = a
	else:
		assert False, "[ERROR]: Unhandled option"
	
if len(sys.argv) == 1:
		print(usage)
		sys.exit()

# Check input arguments are sound
if infile is None:
	print('[ERROR]: Input FASTQ or SFF file must be specified')
	sys.exit(2)
elif outprefix is None:
	print('[ERROR]: Output directory and filename prefix must be specified')
	sys.exit(2)
if split_by_header is True and split_by_sequence is True:
	print('[ERROR]: Can split MIDs either by header information or by sequence, not both')
	sys.exit(2)
elif mid_list is not None and (split_by_header is False and split_by_sequence is False):
	print('[ERROR]: Where to find MID information must be specified (header or sequence)')
	sys.exit(2)
elif mid_list is None and (split_by_header is True or split_by_sequence is True):
	print('[ERROR]: One or more MIDs must be specified for extraction')
	sys.exit(2)
elif remove_primers is True and trimfile is None:
	print('[ERROR]: To remove primers, a file containing primer sequences must be given')
	sys.exit(2)
elif customfile is not None and split_by_sequence is False:
	print('[ERROR]: Cannot use custom MID file while parsing read header for MID information')
	sys.exit(2)
elif primer_list is not None and mid_list is None:
	print('[ERROR]: Primer list cannot be parsed without a MID list')
	sys.exit(2)
elif graph_list is not None and mid_list is None:
	print('[ERROR]: Graph list cannot be parsed without a MID list')
	sys.exit(2)


def convert_csv_to_list(csv):
	if csv == 'all':
		csv = '1,2,3,4,5,6,7,8,9,10,11,12'
	try:
		values = [int(mid) for mid in csv.split(',')]
	except ValueError as err:
		raise
	else:
		return values


# 1) Convert SFF to FASTQ
if infile.endswith('.sff') or infile.endswith('.SFF'):
	outfile = '%s.fq' % outprefix
	print('[INFO]: Converting "%s" to "%s"' % (infile, outfile))
	sff_to_fastq.main(infile, outfile)
	infile = outfile

# 2) Split by MID
if split_by_header is True or split_by_sequence is True:
	try:
		mid_list = convert_csv_to_list(mid_list)
	except ValueError as err:
		print('[ERROR]: Unable to parse MID list: %s' % err)
		sys.exit(2)
	else:
		print('[INFO]: Extracting MIDs ' + str(mid_list) + ' from "%s"' % infile)
	
	try:
		if split_by_header is True:
			pass
			split_mids_by_header.main(infile, outprefix, mid_list)
		elif split_by_sequence is True:
			split_mids_by_sequence.main(infile, outprefix, mid_list, customfile)
	except IOError as err:
		print('[ERROR]: %s' % err)
		sys.exit(2)
		
# 3) Remove primer or BAC sequences
if remove_primers is True:
	if mid_list is not None:
		if primer_list is not None:
			try:
				primer_list = convert_csv_to_list(graph_list)
			except ValueError as err:
				print('[ERROR]: Unable to parse primer list: %s' % err)
				sys.exit(2)
			primer_list = [g for g in primer_list if g in mid_list]
		else:
			primer_list = mid_list
		print('[INFO]: Removing primer sequences from MIDs ' + str(primer_list))
		for a in primer_list:
			infile = '%s.%d.fq' % (outprefix, a)
			p_outprefix = outprefix + '.' + str(a)
			try:
				fastq_primer_remover.main(infile, p_outprefix, trimfile)
			except IOError as err:
				print('[ERROR]: %s' % err)
	else:
		fastq_primer_remover.main(infile, outprefix, trimfile)
		infile = outprefix + '.trim.fq'

# 4) Perform QA, if specified
if perform_qa is True:
	if mid_list is not None:
		if graph_list is not None:
			try:
				graph_list = convert_csv_to_list(graph_list)
			except ValueError as err:
				print('[ERROR]: Unable to parse graph list: %s' % err)
				sys.exit(2)
			graph_list = [g for g in graph_list if g in mid_list]
		else:
			graph_list = mid_list
		print('[INFO]: Creating QA graphs for MIDs ' + str(graph_list))
		for a in graph_list:
			if remove_primers is True:
				infile = '%s.%d.trim.fq' % (outprefix, a)
				graphfile = '%s.%d.trim.jpg' % (outprefix, a)
			else:
				infile = '%s.%d.fq' % (outprefix, a)
				graphfile = '%s.%d.jpg' % (outprefix, a)
			try:
				qa.main(infile, graphfile, r_path, ascii_offset, end_length)
			except IOError as err:
				print('[ERROR]: %s' % err)
	else:
		print('[INFO]: Creating QA graphs of "%s"' % infile)
		if remove_primers is True:
			graphfile = outprefix + '.trim.jpg'
		else:
			graphfile = outprefix + '.jpg'
		qa.main(infile, graphfile, r_path, ascii_offset, end_length)