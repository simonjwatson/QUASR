#! /bin/bash

# This script takes in either a SAM or BAM file and outputs a pileup file
# using SAMtools

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

samtools_path='/nfs/users/nfs_s/sw10/QUASR_v6.07/samtools-0.1.8'

if [ "$#" = "0" ]; then
	echo "[USAGE]: convert_sam_to_pileup.sh <options>
	-h   [None]		Print this help message
	-i   [File]		Input SAM/BAM file
	-o   [String]	Output directory and file prefix
	-r   [File]		Reference FASTA file
"
	exit 1
fi

while getopts ":i:o:r:h" flag; do
	case $flag in
		h) echo "[USAGE]: convert_sam_to_pileup.sh <options>
	-h   [None]		Print this help message
	-i   [File]		Input SAM/BAM file
	-o   [String]	Output directory and file prefix
	-r   [File]		Reference FASTA file
	"
	    exit 1
		;;
		i) input_file=$OPTARG
		;;
		o) output_dir_prefix=$OPTARG
		;;
		r) ref_seq_fasta=$OPTARG
		;;
		:) echo "[ERROR]: Flag -$OPTARG requires an argument. See -h for details."
		   exit 1
		;;
		?) echo "[INFO]: Ignoring invalid option: -$OPTARG."
	esac
done

if [ "$input_file" = '' ] || echo "$input_file" | egrep -q -v '(.sam$|.bam$)'; then
	echo '[ERROR]: Input file must be given and must end in .sam or .bam'
	exit 1
elif [ "$output_dir_prefix" = '' ]; then
	echo '[ERROR]: Output directory and file prefix must be given'
	exit 1
fi

if echo "$input_file" | grep -q '.sam$'; then
	if [ "$ref_seq_fasta" = '' ]; then
		echo '[ERROR]: Reference FASTA file must be given if input file is in SAM format'
		exit 1
	fi
	#Convert sam to bam
	echo '[INFO]: Converting SAM file to BAM file.'
	bam_filename=${output_dir_prefix}.bam
	$samtools_path/samtools view -S -b -T $ref_seq_fasta -o $bam_filename $input_file > /dev/null
else
	bam_filename=${input_file}
fi

sorted_bam_prefix=${output_dir_prefix}.sorted
sorted_bam_filename=${output_dir_prefix}.sorted.bam
pileup_filename=${output_dir_prefix}.pileup

#Sort bam
echo '[INFO]: Creating sorted BAM file.'
$samtools_path/samtools sort $bam_filename $sorted_bam_prefix

#Index sorted bam
echo '[INFO]: Creating indexed sorted BAM file.'
$samtools_path/samtools index $sorted_bam_filename

#Generate pileup (for consensus/SNP detection)
echo '[INFO]: Creating SAMTools pileup of mapped reads.'
if [ "$ref_seq_fasta" = '' ]; then
	$samtools_path/samtools pileup -c $sorted_bam_filename > $pileup_filename
else
	$samtools_path/samtools pileup -c -f $ref_seq_fasta $sorted_bam_filename > $pileup_filename
fi