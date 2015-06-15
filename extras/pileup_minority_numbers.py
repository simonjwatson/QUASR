#! /software/bin/python3

'''This script takes in a pileup file and writes to file all minority bases above 
a certain frequency for all positions. The minority bases are those that are not 
the most frequent base.
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
sys.path.append('/nfs/users/nfs_s/sw10/QUASR_v6.08/modules/')
try:
	from pileup import *
except ImportError:
	print("[ERROR]: Path to QUASR modules not set in script. See 'docs/INSTALL' for more info")
	sys.exit(1)

if sys.version_info < (3,0):
	print("[ERROR]: QUASR requires Python3 to run. Please read 'docs/INSTALL' for more info")
	sys.exit(1)

prog = sys.argv[0]

examples = '''
[EG]: %s -f input.pileup -o output_dir/out
[EG]: %s -f input.pileup -o ../output -c 15 -l 50
[EG]: %s -f input.pileup -o ~/out -c 20 -p /usr/bin/R -x 1000 -n 4''' % (prog, prog, prog)

usage = '''[USAGE]: %s <options>
-h   [None]\tDisplay this usage message with examples (--help)
-f * [File]\tInput pileup file (--infile)
-i   [None]\tIllumina ASCII offset (+64) used to encode quality (--illumina)
-c   [Integer]\tIgnore bases below this Phred score [0] (--cutoff)
-r   [File]\tReference file to confirm genome size if low coverage (--reference)
-t   [Integer]\tDisplay only those bases with a read depth above this value [0] (--depth)
-m   [Float]\tMinority value cutoff [0.99] (--minority)

[NOTE]: Options with * are mandatory. All others are optional.''' % prog

try:
	opts, args = getopt.getopt(sys.argv[1:], "hf:o:ic:dr:q:t:m:", ["help", "infile=", "outprefix=", "illumina", "cutoff=", "deletions", "reference=", "frequency=", "depth=", "minority="])
except getopt.GetoptError as err:
        # print help information and exit:
	print(str(err)) # will print something like "option -a not recognized"
	print(usage)
	sys.exit(2)

infile = None
outprefix = None
reference_file = None
phred_cutoff = 0
ascii_offset = 33
cutoff = 0.99
base_list = 'ACGT*'

for o, a in opts:
	if o in ("-h", "--help"):
		print(usage)
		print(examples)
		sys.exit()
	elif o in ("-f", "--infile"):
		infile = a
	elif o in ("-c", "--cutoff"):
		phred_cutoff = int(a)
	elif o in ("-i", "--illumina"):
		ascii_offset = 64
	elif o in ("-r", "--reference"):
		reference_file = a
	elif o in ("-d", "--deletions"):
		base_list = 'ACGT'
	elif o in ("-m", "--minority"):
		cutoff = float(a)

if len(sys.argv) == 1:
	print(usage)
	sys.exit()

if infile is None:
	print('[ERROR]: Input pileup file must be specified with the "-f" flag')
	sys.exit(2)

with open(infile, 'r') as pileup_fh:
	print('[INFO]: Parsing "%s"' % infile)
	if reference_file is not None:
		with open(reference_file, 'r') as ref_fh:
			pileup = PileupFile(pileup_fh, ref_fh)
	else:
		pileup = PileupFile(pileup_fh)
		
try:	
	bases = pileup.parse_read_bases(phred_cutoff, ascii_offset)
	segs = pileup.return_reference_names()
except AssertionError as err:
	print('[ERROR]: %s' % err)

count = 0
counter = 0
previous_seg = segs[0]
#outfile = outprefix + '.minority.txt'
#with open(outfile, 'w') as outfh:
#	outfh.write('SEGMENT\tPOS\tBASE\tFREQ\tDEPTH\tMAJORITY\n')	
for n in range(len(bases)):
	if 	segs[n] != previous_seg:
		counter = 1
		previous_seg = segs[n]
	else:
		counter += 1
	pos = bases[n]
	if pos == '':
		continue
	total = len(pos)
	max_base = None
	max_freq = 0.0
	freqs = {'A': 0.00, 'C': 0.00, 'G': 0.00, 'T': 0.00, '*': 0.00}
	# Calculate base frequencies	
	for base in base_list:
		f = pos.count(base) / total
		freqs[base] = f
		if f > max_freq:
			max_base = base
			max_freq = f
			
	if max_freq >= 0.99:
		count = count + 1
				
perc = float(count) / len(bases)
print('[INFO]: Genome percentage with major base above cutoff: %.2f' % perc)
