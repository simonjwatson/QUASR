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

from fastq import *

def main(forward_file, outfile_f, reverse_file=None, outfile_r=None, paired=False, ascii_offset=33, median_cutoff=20, length_cutoff=30):
	with open(forward_file, 'r') as infh_f, open(outfile_f, 'w') as outfh_f:
		if paired is True:
			infh_r = open(reverse_file, 'r')
			outfh_r = open(outfile_r, 'w')
			reverse_reads = fastq_iterator(infh_r)
		
		total_reads = 0
		for_passed = 0
		for_failed = 0
		for_recovered = 0
		if paired is True:
			both_passed = 0
			rev_failed = 0
			rev_recovered = 0
			both_failed = 0
			print('[INFO]: Performing PE QC on "%s" and "%s"' % (forward_file, reverse_file))
		else:
			print('[INFO]: Performing SE QC on "%s"' % forward_file)
		# loop through the forward reads' generator, calling next() on the reverse reads one to keep in step.
		for forward_head, forward_seq, forward_qual in fastq_iterator(infh_f):
			total_reads += 1
			try:
				forward_read = FastqRecord(forward_head, forward_seq, forward_qual)
			except IOError as err:
				print('[ERROR]: Unable to handle "%s": %s' % (forward_head, err))
				continue
			else:
				for_length = forward_read.get_sequence_length()
			if paired is True:
				try:
					reverse_head, reverse_seq, reverse_qual = next(reverse_reads)
					reverse_read = FastqRecord(reverse_head, reverse_seq, reverse_qual)
				except IOError as e:
					print('[ERROR]: Unable to handle "%s": %s' % (reverse_head, err))
					continue
				except StopIteration as err:
					print('[ERROR]: More reverse reads than forward reads!')
					return 1
				else:
					rev_length = reverse_read.get_sequence_length()
					# First check read length > cutoff
					if rev_length < length_cutoff:
						rev_failed += 1
						continue
			
			if for_length < length_cutoff:
				for_failed += 1
				continue
			
			# Secondly, check if the read passes the median cutoff
			if forward_read.calculate_median_quality(ascii_offset=ascii_offset) >= median_cutoff:
				if paired is True:
					if reverse_read.calculate_median_quality(ascii_offset=ascii_offset) >= median_cutoff:
						both_passed += 1
						# PE foward and reverse passed
						forward_read.write_to_file(outfh_f)
						reverse_read.write_to_file(outfh_r)
					else:
						# PE forward passed, reverse failed
						# rev_length = reverse_read.get_sequence_length()
						while True:
							reverse_read.remove_nth_base(rev_length)
							rev_length -= 1
							if rev_length < length_cutoff:
								break
							elif reverse_read.calculate_median_quality(ascii_offset=ascii_offset) >= median_cutoff:
								break
						# out of loop - either length failed or median passed.
						if rev_length < length_cutoff:
							rev_failed += 1
							continue
						else:
							rev_recovered += 1
							# PE forward passed reverse passed
							forward_read.write_to_file(outfh_f)
							reverse_read.write_to_file(outfh_r)
				else:
					# SE forward passed
					for_passed += 1
					forward_read.write_to_file(outfh_f)
			
			else:
				# if both forward and reverse fail initial check, probably a bad spot
				if paired is True:
					if reverse_read.calculate_median_quality(ascii_offset=ascii_offset) >= median_cutoff:
						# for_length = forward_read.get_sequence_length()
						while True:
							forward_read.remove_nth_base(for_length)
							for_length -= 1
							if for_length < length_cutoff:
								break
							elif forward_read.calculate_median_quality(ascii_offset=ascii_offset) >= median_cutoff:
								break
						# out of loop - either length failed or median passed.
						if for_length < length_cutoff:
							for_failed += 1
							continue
						else:
							for_recovered += 1
							# PE forward passed reverse passed
							forward_read.write_to_file(outfh_f)
							reverse_read.write_to_file(outfh_r)
					else:
						both_failed += 1
						continue
				else:
					# SE forward failed. Try to recover
					for_length = forward_read.get_sequence_length()
					while True:
						forward_read.remove_nth_base(for_length)
						for_length -= 1
						if for_length < length_cutoff:
							break
						elif forward_read.calculate_median_quality(ascii_offset=ascii_offset) >= median_cutoff:
							break
					# out of loop - either length failed or median passed.
					if for_length < length_cutoff:
						for_failed += 1
						continue
					else:
						for_recovered += 1
						# PE forward passed reverse passed
						forward_read.write_to_file(outfh_f)

		# Summary data
		if paired is True:
			total_passed = both_passed + for_recovered + rev_recovered
			total_failed = both_failed + for_failed + rev_failed
			print('[TOTAL]: Read pairs in input file: %d' % total_reads)
			print('[TOTAL]: Read pairs passed: %d (%.2f%%)' % (total_passed, (total_passed/total_reads)*100))
			print('[STATS]: Read pairs passed on both forward and reverse reads: %d (%.2f%%)' % (both_passed, (both_passed/total_reads)*100))
			print('[STATS]: Read pairs passed after trimming forward read: %d (%.2f%%)' % (for_recovered, (for_recovered/total_reads)*100))
			print('[STATS]: Read pairs passed after trimming reverse read: %d (%.2f%%)' % (rev_recovered, (rev_recovered/total_reads)*100))
			print('[TOTAL]: Read pairs failed: %d (%.2f%%)' % (total_failed, (total_failed/total_reads)*100))
			print('[STATS]: Read pairs failed on both forward and reverse reads: %d (%.2f%%)' % (both_failed, (both_failed/total_reads)*100))
			print('[STATS]: Read pairs failed on forward read: %d (%.2f%%)' % (for_failed, (for_failed/total_reads)*100))
			print('[STATS]: Read pairs failed on reverse read: %d (%.2f%%)' % (rev_failed, (rev_failed/total_reads)*100))
		else:
			total_passed = for_passed + for_recovered
			print('[TOTAL]: Reads in input file: %d' % total_reads)
			print('[TOTAL]: Reads passed: %d (%.2f%%)' % (total_passed, (total_passed/total_reads)*100))
			print('[STATS]: Reads passed without trimming: %d (%.2f%%)' % (for_passed, (for_passed/total_reads)*100))
			print('[STATS]: Reads passed after trimming: %d (%.2f%%)' % (for_recovered, (for_recovered/total_reads)*100))
			print('[TOTAL]: Reads failed: %d (%.2f%%)' % (for_failed, (for_failed/total_reads)*100))