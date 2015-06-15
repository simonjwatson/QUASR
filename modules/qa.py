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

import tempfile
import os
from fastq import *

def determine_R_path():
	if 'R_PATH' in os.environ:
		return os.environ['R_PATH']
	elif os.system('which R > /dev/null') == 0:
		return 'R'
	else:
		return None

def main(fastq_file, output_file, r_binary_path=None, ascii_offset=33, window_size=15):
	# First determine R binary path
	if r_binary_path is not None:
		if not os.path.exists(r_binary_path):
			print('[ERROR]: Specified path to R binary is invalid')
			return 1
	else:
		r_binary_path = determine_R_path()
		if r_binary_path is None:
			print('[ERROR]: Unable to determine R path. Please use flag to specify.')
			return 1
		else:
			print('[INFO]: R binary found')
		
	# now create temp files to store R commands and R input data
	try:
		r_datafile = tempfile.NamedTemporaryFile(delete=False, mode='w')
		r_commandsfile = tempfile.NamedTemporaryFile(delete=False, mode='w')
	except IOError as err:
		print('[ERROR]: Unable to open temporary files to write R commands: %s' % err)
		return 1
	
	# Next, go through the FASTQ file record-by-record, calculate metrics, and write to R datafile
	with open(fastq_file, 'r') as infh:
		print('[INFO]: Calculating QA metrics for "%s"' % fastq_file)
		for header, sequence, quality in fastq_iterator(infh):
			try:
				read = FastqRecord(header, sequence, quality)
			except IOError as err:
				print('[ERROR]: Unable to handle "%s": %s' % (header, err))
				continue
			
			# Want to calculate read length, GC%, median percentage
			read_length = read.get_sequence_length()
			gc = read.calculate_gc_percentage()
			median = read.calculate_median_quality(ascii_offset=ascii_offset)
			
			r_datafile.write('%.2f\t%d\t%d' % (gc, median, read_length))
			
			# Now handle 3' end cross-sectional window
			phreds = read.return_phred_scores(start=-window_size, ascii_offset=ascii_offset)
			phred_size = len(phreds)
			assert phred_size <= window_size, 'Window size incorrectly parsed'
			if phred_size < window_size:
				diff = window_size - phred_size
				r_datafile.write('%s' % '\tNA' * diff)
			for p in phreds:
				r_datafile.write('\t%d' % p)
			r_datafile.write('\n')
			
		r_datafile.close()
		
	# Finally generate the string to write into the R commands file
	r_commands = '''raw.data <- read.table('%s', header=F, sep='\\t')
jpeg(file='%s', height=7016, width=4960, res=600)
par(oma=c(0,0,2,0))
par(mar=c(4,4,4,2))
par(font.main=2)
par(xaxs="i")
par(yaxs="i")
par(cex.axis=0.9)
layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE))

mean.length <- round(mean(raw.data[,3]), digits=2)
sd.length <- round(sd(raw.data[,3]), digits=2)
hist(raw.data[,3], breaks=20, xlab="Read length", col="skyblue", xlim=c(0, max(raw.data[,3])), main=paste("Mean length:", mean.length, "+/-", sd.length))

mean.gc <- round(mean(raw.data[,1]), digits=2)
sd.gc <- round(sd(raw.data[,1]), digits=2)
hist(raw.data[,1], breaks=20, xlab="GC %%", col="lemonchiffon1", main=paste("Mean GC%%:", mean.gc, "+/-", sd.gc))

mean.median <- round(mean(raw.data[,2]), digits=2)
plot(raw.data[,2]~raw.data[,3], xlim=c(0,max(raw.data[,3])), ylim=c(0,max(raw.data[,2])), pch=18, col="gray70", xlab="Read length", ylab="Median quality", main="Read median quality as a function of length")
abline(h=mean.median, col="black", lty=2)
abline(v=mean.length, col="black", lty=2)

sd.median <- round(sd(raw.data[,2]), digits=2)
hist(raw.data[,2], breaks=20, xlab="Read median quality", xlim=c(0, max(raw.data[,2])), col="mistyrose2", main=paste("Mean median-quality:", mean.median, "+/-", sd.median))

par(mar=c(5,7,5,5))
par(xaxs="r")
means <- colMeans(raw.data[4:%d], na.rm=T)
y.max <- max(means, na.rm=T)
remainder <- y.max%%%%5
y.max <- y.max + (5-remainder)
plot(means, ylim=c(0,y.max), xlab="Position from end of read", ylab="Mean quality", axes=F, col="pink", pch=20, main="3'-end cross-sectional mean quality")
points(means, ylim=c(0,y.max), col="red", type='l')
axis(1, at=seq(0, %d, 10), lab=seq(-%d, 0, 10))
axis(2, at=seq(0, y.max, 5))

title(main=paste("%s total sequences:", length(raw.data[,1])), outer=T)
dev.off()''' % (r_datafile.name, output_file, window_size+3, window_size, window_size, os.path.basename(fastq_file))

	r_commandsfile.write(r_commands)
	r_commandsfile.close()
	if (os.system('%s CMD BATCH %s' % (r_binary_path, r_commandsfile.name))) == 0:
		print('[INFO]: QA graphs for "%s" written to "%s"' % (fastq_file, output_file))
		try:
			os.unlink(os.path.basename(r_commandsfile.name) + '.Rout')
			os.unlink(r_commandsfile.name)
			os.unlink(r_datafile.name)
		except OSError as err:
			print('[ERROR]: Unable to complete clean up: %s' % err)
	else:
		print('[ERROR]: Execution of "%s" failed' % r_binary_path)