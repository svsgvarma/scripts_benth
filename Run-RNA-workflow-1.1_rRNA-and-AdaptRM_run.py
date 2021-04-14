#!/usr/bin/python

"""
#Script to Pars vcf...
#####------Inputs-------
# python Run-RNA-workflow-1.1_rRNA-and-AdaptRM_run.py /home/gala0002/proj/proj_Ramesh/RNA-seq_2020/
"""

import sys
import re
import fnmatch, os
import tempfile
import commands
import subprocess
#import subprocess32
from subprocess import *
from subprocess import call

class fileHandler:
	def __init__(self):
		self.data = []
		#print "Calling fileHandler constructor"
	def open_file(self,readfl):
		self.rfile = open(readfl,'r').readlines()
		return self.rfile
	def write_file(self,writefl):
		self.wfile = open(writefl,'w')
		return self.wfile

class SearchDB(fileHandler):
	def __init__(self):
		self.data = []
		from collections import defaultdict
		self.ident_ranges_HMBM = defaultdict(list)
	def Search_CADD(self,readdir):
		"""
		Calling Search localsearch method
		"""
		SCR="/home/gala0002/proj/RNAseq-analysis-run/scripts_benth/"
		INDIR=readdir+"200130_A00605_0102_AHJ75GDRXX"
		OUTDIR=readdir+"SL-2400_1.1_sort-trim"
		if not os.path.exists(OUTDIR):
			os.makedirs(OUTDIR)
		#cmdFls1 = SCR+"Run-RNA-workflow-1.1_rRNA-and-AdaptRM.sh NG-14833_Pas_M_1_lib233714_5767_1"
		#grepout1 =  subprocess.check_output(cmdFls1, shell=True)
		#print grepout1
		file2 = self.write_file(SCR+"Run-RNA-workflow-1.1_rRNA-and-AdaptRM.sh.log")
		for file in os.listdir(INDIR):
			if fnmatch.fnmatch(file, 'SL-2400*.fastq.gz'):
				fl_spid = file.split("_")
				fl_samid = "_".join(fl_spid[:-2])
				fl_samlib = "_".join(fl_spid[:-1])
				cmdFls1 = SCR+"Run-RNA-workflow-1.1_rRNA-and-AdaptRM.sh "+fl_samid
				grepout1 = subprocess.check_output(cmdFls1, shell=True)
				file2.write(str(grepout1)+"\n")
		print "Done parsing...."
		return None

# file1: Input positions SNP
# write1: Output file

clF1 = SearchDB().Search_CADD(sys.argv[1])



