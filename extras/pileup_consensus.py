#! /software/bin/python3

'''Generates a consensus sequene from a pileup file. Works in 1 of 2 ways depending on command-line
flags used:
1) Quality-independent. Uses all bases - straight frequency at high depths, N's at low depths.
2) Quality dependent. Ignores all bases <cutoff. If depth < cutoff, lower-case sequence used. If
depth drops to 0 because of quality threshold, puts in an N, otherwise puts in a -
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
	from pileup import *
except ImportError:
	print("[ERROR]: Path to QUASR modules not set in script. See 'docs/INSTALL' for more info")
	sys.exit(1)

if sys.version_info < (3,0):
	print("[ERROR]: QUASR requires Python3 to run. Please read 'docs/INSTALL' for more info")
	sys.exit(1)

prog = sys.argv[0]

examples = '''
[EG]: %s -f input.pileup -o outDir/outPrefix -q -c 15
[EG]: %s -f input.pileup -o outDir/outPrefix -q -c 15 -l 50
[EG]: %s -f input.pileup -o outDir/outPrefix -d -l 50''' % (prog, prog, prog)

usage = '''[USAGE]: %s <options>
-h   [None]\tDisplay this usage message with examples (--help)

--GENERAL--
-f * [File]\tInput pileup file (--infile)
-o * [String]\tOutput directory and file prefix (--outprefix)
-i   [None]\tIllumina ASCII offset (+64) used to encode quality (--illumina) 
-a   [Float]\tMinority base frequency for inclusion as ambiguity code [0.3] (--ambiguity)
-r   [File]\tReference file to confirm genome size if low coverage (--reference)

--EITHER--
-q + [None]\tGenerate quality-dependent consensus sequence (--dependent)
-c   [Integer]\tIgnore bases below this Phred score [0] (--cutoff)
-l   [Integer]\tDepth below which bases are written in lowercase [10] (--lowcoverage)

--OR--
-d + [None]\tGenerate quality-independent consensus sequence (--independent)
-l   [Integer]\tDepth below which bases are written as 'N's [10] (--lowcoverage)

[NOTE]: Options with * are mandatory. Those with + are mandatory but mutually-exclusive.''' % prog

try:
	opts, args = getopt.getopt(sys.argv[1:], "hf:o:ic:r:a:ql:d", ["help", "infile=", "outprefix=", "illumina", "cutoff=", "reference=", "ambiguity=", "dependent", "lowcoverage=", "independent"])
except getopt.GetoptError as err:
        # print help information and exit:
	print(str(err)) # will print something like "option -a not recognized"
	print(usage)
	sys.exit(2)

pileup_file = None
outprefix = None
phred_cutoff = 0
ascii_offset = 33
ambiguity_threshold = 0.30
low_depth_cutoff = 10
quality_independent = None
reference_file = None
checker = 0 # This is just to check that both options haven't been specified

for o, a in opts:
	if o in ("-h", "--help"):
		print(usage)
		print(examples)
		sys.exit()
	elif o in ("-f", "--infile"):
		pileup_file = a
	elif o in ("-o", "--outprefix"):
		outprefix = a
	elif o in ("-c", "--cutoff"):
		phred_cutoff = int(a)
	elif o in ("-a", "--ambiguity"):
		ambiguity_threshold = float(a)
	elif o in ("-l", "--lowcoverage"):
		low_depth_cutoff = int(a)
	elif o in ("-d", "--independent"):
		quality_independent = True
		checker += 1
	elif o in ("-q", "--dependent"):
		quality_independent = False
		checker += 1
	elif o in ("-r", "reference="):
		reference_file = a

if len(sys.argv) == 1:
	print(usage)
	sys.exit()

if pileup_file is None:
	print('[ERROR]: Input pileup file must be specified with the "-f" flag')
	sys.exit(2)
elif outprefix is None:
	print('[ERROR]: Output directory and file prefix must be specified with the "-o" flag')
	sys.exit(2)
elif quality_independent is None:
	print('[ERROR]: Either quality-dependence or -independence must be chosen')
	sys.exit(2)
elif checker > 1:
	print('[ERROR]: Cannot specify both quality-dependent and -independent flags')
	sys.exit(2)

if quality_independent is True and phred_cutoff != 0:
	print('[WARNING]: Quality-independent consensus does not use Phred values. Ignoring "-c" flag')
	phred_cutoff = 0

AmbiguityCodes = { 	'AC': 'M', 'AG': 'R', 'AT': 'W',
					'CG': 'S', 'CT': 'Y', 'GT': 'K',
					'ACG': 'V', 'ACT': 'H', 'AGT': 'D',
					'CGT': 'B', 'ACGT': 'N' }

outfile = outprefix + '.consensus.fasta'
with open(pileup_file, 'r') as infh, open(outfile, 'w') as outfh:
	print('[INFO]: Parsing "%s"' % pileup_file)
	if reference_file is not None:
		with open(reference_file, 'r') as ref_fh:
			pileup = PileupFile(infh, ref_fh)
	else:
		pileup = PileupFile(infh)
	
	# Get the read bases at each position	
	if quality_independent is True:
		print('[INFO]: Generating quality-independent consensus')
		base_list = pileup.parse_read_bases(phred_cutoff, ascii_offset)
	else: # This removes all bases below the cutoff
		print('[INFO]: Generating quality-dependent consensus')
		base_list = pileup.parse_read_bases(phred_cutoff, ascii_offset)
	segment_list = pileup.return_reference_names()
	previous_segment = None
	genome_size = len(base_list)
	consensus = ''
	for i in range(genome_size):
		current_segment = segment_list[i]
		if current_segment != previous_segment:
			if previous_segment is not None:
				#outfh.write('>%s %s | amb:%.2f | phred:%d | depth:%d\n%s\n' % (outprefix, previous_segment, ambiguity_threshold, phred_cutoff, low_depth_cutoff, consensus))
				outfh.write('>%s\n%s\n' % (previous_segment, consensus))
				consensus = ''
				
		previous_segment = current_segment
		bases = base_list[i].upper()
		depth = len(bases)
		if depth == 0:
			if pileup._read_depth[i] > 0: # bad coding to access private variable, but easiest way to find out if there were bases covering!
				consensus += 'N'
			else:
				consensus += '-'
			continue
		elif depth < low_depth_cutoff and quality_independent is True:
			consensus += 'N'
			continue
		else:
			consensus_base = ''
			max_freq = 0.0
			max_base = ''
			if '*' in bases: # If more reads show a deletion than not, consider it as a deletion
				del_freq = bases.count('*') / depth
				if del_freq > 0.5:
					consensus += '-'
					continue
			for base in 'ACGT': # First calculate bases above threshold
				f = bases.count(base) / depth
				if f > max_freq:
					max_freq = f # Store the freqency of the highest frequency base to compare against deletion frequency
					max_base = base
				elif f == max_freq:
					max_base += base

				if f > ambiguity_threshold:
					consensus_base += base
			if len(consensus_base) == 0:
				consensus_base = max_base
			if len(consensus_base) > 1:
				consensus_base = AmbiguityCodes.get(consensus_base.upper(), 'N')
			if depth < low_depth_cutoff: # must be quality independent == False because the alternative is already dealt with above
				consensus += consensus_base.lower()
			else:
				consensus += consensus_base.upper()
				
	#outfh.write('>%s %s | amb:%.2f | phred:%d | depth:%d\n%s\n' % (outprefix, previous_segment, ambiguity_threshold, phred_cutoff, low_depth_cutoff, consensus))
	outfh.write('>%s\n%s\n' % (previous_segment, consensus))
	print('[INFO]: Consensus sequence written to "%s"' % outfile)
