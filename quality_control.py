#! /software/bin/python3

'''This script takes in a FASTQ file, and for each sequence calculates the median quality
score. If this is lower than a user-specified cutoff, it trims from the 3' end until either
the read becomes too small or its median raises above the cutoff. Passed reads are written to
file. If paired-end, both reads muts pass to be written to their respective output files.
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
sys.path.append('/Users/sw10/Dropbox/Sanger/QUASR/QUASR_v6.09/modules/')
try:
	import qc, qa
except ImportError:
	print("[ERROR]: Path to QUASR modules not set in script. See 'docs/INSTALL' for more info")
	sys.exit(1)

if sys.version_info < (3,0):
	print("[ERROR]: QUASR requires Python3 to run. Please read 'docs/INSTALL' for more info")
	sys.exit(1)

prog = sys.argv[0]
examples = '''
[EG]: %s -f reads1.fastq -o example_output -m 20 -l 50
[EG]: %s -f reads1.fastq -r reads2.fastq -o example_output -m 30 -l 150 -g
[EG]: %s -f reads1.fastq -o ../example_output -m 25 -l 50 -g -p usr/bin/R -e 30''' % (prog, prog, prog)

usage = '''[USAGE]: %s <options>
-h   [None]\tDisplay help message with examples (--help)

--GENERAL--
-f * [File]\tFASTQ containing forward reads (--forward)
-r   [File]\tFASTQ containing reverse reads if paired-end (--reverse)
-o * [String]\tOutput directory and filename prefix (--outprefix)
-i   [None]\tIllumina ASCII offset (+64) used to encode quality (--illumina)

--QUALITY CONTROL--
-m * [Integer]\tMedian quality cutoff (--median)
-l * [Integer]\tRead length cutoff (--length)

--QUALITY-ASSURANCE--
-g   [None]\tGenerate QA graphs of output readset (--graphs)
-p   [File]\tPath to R binary for stats generation and graphing (--path)
-e   [Integer]\tNumber of bases to graph 3' mean quality dropoff [15] (--end)

[NOTE]: Options with * are mandatory. All others are optional.''' % prog

try:
	opts, args = getopt.getopt(sys.argv[1:], "hf:r:o:im:l:ge:p:", ["help", "forward=", "reverse=", "outprefix=", "illumina", "median=", "length=", "graphs", "end=", "path="])
except getopt.GetoptError as err:
	print(str(err))
	print(usage)
	sys.exit(2)

for_file = None
rev_file = None
outprefix = None
ascii_offset = 33
median_cutoff = None
length_cutoff = None
perform_qa = False
end_length = 15
r_path = None
paired = False

for o,a in opts:
	if o in ("-h", "--help"):
		print(usage)
		print(examples)
		sys.exit()
	elif o in ("-f", "--forward"):
		for_file = a
	elif o in ("-r", "--reverse"):
		rev_file = a
	elif o in ("-o", "--outprefix"):
		outprefix = a
	elif o in ("-i", "--illumina"):
		ascii_offset = 64
	elif o in ("-g", "--graphs"):
		perform_qa = True
	elif o in ("-e", "--end"):
		end_length = int(a)
	elif o in ("-p", "--path"):
		r_path = a
	elif o in ("-m", "--median"):
		median_cutoff = int(a)
	elif o in ("-l", "--length"):
		length_cutoff = int(a)
	else:
		assert False, "[ERROR]: Unhandled option"
	
if len(sys.argv) == 1:
		print(usage)
		sys.exit()
		
if for_file is None:
	print('[ERROR]: FASTQ file containing forward reads must be specified')
	sys.exit(2)
elif outprefix is None:
	print('[ERROR]: Output directory and filename prefix must be specified')
	sys.exit(2)
elif median_cutoff is None:
	print('[ERROR]: Median cutoff value must be specified')
	sys.exit(2)
elif length_cutoff is None:
	print('[ERROR]: Length cutoff value must be specified')
	sys.exit(2)

if rev_file is not None:
	paired = True
	outfile_r = outprefix + '.r.fq'
else:
	outfile_r = None
	
outfile_f = outprefix + '.f.fq'

print('[INFO]: Input parameters successfully parsed')
qc.main(for_file, outfile_f, rev_file, outfile_r, paired, ascii_offset, median_cutoff, length_cutoff)

if perform_qa is True:
	graphfile_f = outprefix + '.f.jpg'
	qa.main(outfile_f, graphfile_f, r_path, ascii_offset, end_length)
	
	if paired is True:
		graphfile_r = outprefix + '.r.jpg'
		qa.main(outfile_r, graphfile_r, r_path, ascii_offset, end_length)