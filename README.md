# QUASR
## Installation
QUASR is intended to be a self-contained pipeline that has minimal requirements:

1. The Python3 interpreter must be installed.

2. Graphing requires R to be installed.

3. To use any of the `pileup_*` scripts in the `extras/` folder, SAMtools is recommended to create the pileup file.

If you have these installed, follow these steps to get QUASR running:

1.  
	* Change the top line of `readset_parser.py`, `quality_control.py`, and all scripts in the `extras/` folder to the location of your Python3 interpreter.
	* Or just call the script with python3 (or full path to python3 interpreter). _IE_: `python3 readset_parser.py` or `/path_to_python3_binary/python3 readset_parser.py`

2. In the same files, change the line near the top that reads `sys.path.append('/nfs/users/nfs_s/sw10/QUASR6/modules/')` to `sys.path.append('full_path_to_QUASR_folder/modules/')`

## Composition
QUASR is comprised of two main pre-assembly scripts that parse and process read files:

1. `readset_parser.py`

2. `quality_control.py`

and a number of post-assembly scripts that parse pileup or mpileup files:
* `pileup_consensus.py`
* `pileup_depth_graph.py`
* `pileup_minority_graph.py`
* `pileup_minority_list.py`
* `pileup_minority_numbers.py`

In addition, a number of useful additional scripts have been included to help handle read data:
* `convert_sam_to_pileup.sh`
* `fastq_duplicate_remover.py`
* `fastq_remove_Ns.py`

To find out what inputs any script takes, invoke it without any arguments, or with `-h/--help`.

## Usage
