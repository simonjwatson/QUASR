#! /software/bin/python3

'''This script takes in a pileup file and uses R to print out a stacked histogram of the
minority bases. The minority bases are those that are not the most frequent base.
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

import sys, os, math, tempfile, getopt
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
[EG]: %s -f input.pileup -o output_dir/out
[EG]: %s -f input.pileup -o ../output -c 15 -l 50
[EG]: %s -f input.pileup -o ~/out -c 20 -p /usr/bin/R -x 1000 -n 4''' % (prog, prog, prog)

usage = '''[USAGE]: %s <options>
-h   [None]\tDisplay this usage message with examples (--help)

--GENERAL--
-f * [File]\tInput pileup file (--infile)
-o * [String]\tOutput directory and file prefix (--outprefix)
-i   [None]\tIllumina ASCII offset (+64) used to encode quality (--illumina)
-c   [Integer]\tIgnore bases below this Phred score [0] (--cutoff)
-r   [File]\tReference file to confirm genome size if low coverage (--reference)

--GRAPH OPTIONS--
-p   [File]\tPath to R binary for stats generation and graphing (--path)
-x   [Integer]\tNumber of bases per plot (--xaxis)
-n   [Integer]\tMaximum number of plots per output JPG [4] (--num)
-l   [Integer]\tBackground shade bases with depth below cutoff [10] (--lowcoverage)

[NOTE]: Options with * are mandatory. All others are optional.''' % prog

try:
	opts, args = getopt.getopt(sys.argv[1:], "hf:o:ic:r:p:x:n:l:", ["help", "infile=", "outprefix=", "illumina", "cutoff=", "reference=", "path=", "xaxis=", "num=", "lowcoverage="])
except getopt.GetoptError as err:
        # print help information and exit:
	print(str(err)) # will print something like "option -a not recognized"
	print(usage)
	sys.exit(2)

pileup_file = None
outprefix = None
phred_cutoff = 0
depth_cutoff = 10
ascii_offset=33
r_binary_path = None
max_x_axis = None
num_graphs = 4
reference_file = None

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
	elif o in ("-x", "--xaxis"):
		max_x_axis = int(a)
	elif o in ("-n", "--num"):
		num_graphs = int(a)
	elif o in ("-i", "--illumina"):
		ascii_offset = 64
	elif o in ("-p", "path="):
		r_binary_path = a
	elif o in ("-l", "--lowcoverage"):
		depth_cutoff = int(a)
	elif o in ("-r", "--reference"):
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

def determine_R_path():
	if 'R_PATH' in os.environ:
		return os.environ['R_PATH']
	elif os.system('which R > /dev/null') == 0:
		return 'R'
	else:
		return None

# First determine R binary path
if r_binary_path is not None:
	if not os.path.exists(r_binary_path):
		print('[ERROR]: Specified path to R binary is invalid')
		sys.exit(2)
else:
	r_binary_path = determine_R_path()
	if r_binary_path is None:
		print('[ERROR]: Unable to determine R path. Please use flag to specify.')
		sys.exit(2)
	else:
		print('[INFO]: R binary found')
	
# now create temp files to store R commands and R input data
try:
	r_datafile = tempfile.NamedTemporaryFile(delete=False, mode='w')
	r_commandsfile = tempfile.NamedTemporaryFile(delete=False, mode='w')
except IOError as err:
	print('[ERROR]: Unable to open temporary files to write R commands: %s' % err)
	sys.exit(2)

print('[INFO]: Parsing "%s"' % pileup_file)
with open(pileup_file, 'r') as pileup_fh:
	print('[INFO]: Parsing "%s"' % pileup_file)
	if reference_file is not None:
		with open(reference_file, 'r') as ref_fh:
			pileup = PileupFile(pileup_fh, ref_fh)
	else:
		pileup = PileupFile(pileup_fh)
try:	
	bases = pileup.parse_read_bases(phred_cutoff, ascii_offset)
except AssertionError as err:
	print('[ERROR]: %s' % err)
	r_datafile.close()
	r_commandsfile.close()
	os.unlink(r_datafile.name)
	os.unlink(r_commandsfile.name)
	sys.exit()
	
# Bases is a list of the read bases for each reference position
# IE ['ATCGAT', 'CGATCGACTA', 'TACGATCGA']
max_sum = 0.0
for pos in bases:
	if pos == '':
		r_datafile.write('0.00\t0.00\t0.00\t0.00\t1\n') # A/C/G/T/binary depth (1: <cutoff, 0: >cutoff)
		continue
	
	freqs = {'A': 0.00, 'C': 0.00, 'G': 0.00, 'T': 0.00}
	total = len(pos)
	max_freq = 0.0
	max_base = ''
	sum = 0.0
	# First calculate base frequencies
	for base in 'ACGT':
		f = pos.count(base) / total
		if f > max_freq:
			max_freq = f
			max_base = base
		freqs[base] = f
	# If more reads show a deletion than not, consider it as a deletion
	if '*' in pos:
		del_freq = pos.count('*') / total
		if del_freq > max_freq:
			r_datafile.write('0.00\t0.00\t0.00\t0.00\t1\n')
			continue
			
	for base in 'ACGT':
		if base == max_base:
			r_datafile.write('0.00\t') # doesn't print out majority base for readability
		else:
			b = freqs[base]
			r_datafile.write('%.3f\t' % b)
			sum += b
	if total > depth_cutoff:
		r_datafile.write('0\n') # for plotting grey lines at low coverage
	else:
		r_datafile.write('1\n')
	if sum > max_sum:
		max_sum = sum
		
r_datafile.close()
			
print('[INFO]: Generating minority plots')

# First, work out how many plots will be produced
num_pages = None
# if a max num of nucleotides on X axis is not specified, each segment will be put on a new plot
if max_x_axis is None:
	seg_counter = pileup.calc_segment_sizes()
	num_segs = len(seg_counter)
	# If fewer segments than plots per page, make plots per page = number of segments
	if num_segs < num_graphs:
		num_graphs = num_segs 
	num_pages = math.ceil(num_segs/num_graphs)
else:
	genome_size = len(bases)
	graph_boundaries = [(b+1, b+max_x_axis) for b in range(0, genome_size, max_x_axis) if b+max_x_axis <= genome_size]
	remainder = genome_size % max_x_axis
	if remainder != 0:
		last_point = graph_boundaries[-1][1]
		graph_boundaries.append((last_point+1, (last_point + genome_size%max_x_axis)))
	total_graphs = len(graph_boundaries)
	if total_graphs < num_graphs:
		num_graphs = total_graphs
	num_pages = math.ceil(total_graphs/num_graphs)

r_commands = '''raw.data <- read.table	('%s', header=F, sep="\\t")
binary.depth <- raw.data[5]
raw.data <- raw.data[1:4]''' % r_datafile.name

graphfile = outprefix + '.minority%02d.jpg'
r_commands += '''
jpeg(file='%s', height=7016, width=4960, res=600)''' % graphfile

counter = 0
min = 0
max = 0
remainder = 0
for page in range(num_pages):
	r_commands += '''
par(oma=c(1,1,0,1))
par(mar=c(3,4,3,0))
par(mfrow=c(%d,1))
par(xpd=NA)
par(xaxs="i")
par(yaxs="i")
par(cex=0.8)''' % num_graphs

	if max_x_axis is None:
		if counter+num_graphs > num_segs:
			remainder = (counter+num_graphs) - num_segs 
		for s in range(num_graphs*page, (counter+num_graphs)-remainder):
			counter += 1
			name = seg_counter[s][0]
			size = seg_counter[s][1]
			min = max + 1
			max = max + size
			r_commands += '''
barplot(t(as.matrix(binary.depth[%d:%d,])), col="grey80", border=NA, xlab="", ylab="", space=0, axes=F)
par(new=T)
barplot(t(as.matrix(raw.data[%d:%d,])), ylim=c(0, %.3f), col=c("darkgreen", "darkblue", "goldenrod4", "red"), border=NA, ylab="Read frequency", names.arg=seq(1,%d), main="%s")
tmp.u <- par('usr')
legend(tmp.u[1], tmp.u[4], c("A", "C", "G", "T"), col=c("darkgreen", "darkblue", "goldenrod4", "red"), lty=1, cex=0.5, horiz=T, bty="n")''' % (min, max, min, max, max_sum, size, name)
		
	else:
		if counter+num_graphs > total_graphs:
			remainder = (counter+num_graphs) - total_graphs
		for s in range(num_graphs*page, (counter+num_graphs)-remainder):
			counter += 1
			lower = graph_boundaries[s][0]
			upper = graph_boundaries[s][1]
			if s is not 0:
				min = lower
			r_commands += '''
barplot(t(as.matrix(binary.depth[%d:%d,])), col="grey80", border=NA, xlab="", ylab="", space=0, axes=F)
par(new=T)
barplot(t(as.matrix(raw.data[%d:%d,])), ylim=c(0, %.3f), col=c("darkgreen", "darkblue", "goldenrod4", "red"), border=NA, ylab="Read frequency", names.arg=seq(%d,%d), main='Nucleotides %d-%d')
tmp.u <- par('usr')
legend(tmp.u[1], tmp.u[4], c("A", "C", "G", "T"), col=c("darkgreen", "darkblue", "goldenrod4", "red"), lty=1, cex=0.5, horiz=T, bty="n")''' % (min, upper, min, upper, max_sum, lower, upper, lower, upper)

r_commands += '''
dev.off()'''
r_commandsfile.write(r_commands)
r_commandsfile.close()
os.system('%s CMD BATCH %s' % (r_binary_path, r_commandsfile.name))
try:
	os.unlink(os.path.basename(r_commandsfile.name) + '.Rout')
	os.unlink(r_commandsfile.name)
	os.unlink(r_datafile.name)
except OSError as e:
	print('[ERROR]: Unable to complete clean up: %s' % e)
else:
	print('[INFO]: Depth plots written to:')
	for i in range(num_pages):
		print('[INFO]: %s.minority%.2d.jpg' % (outprefix, (i+1)))