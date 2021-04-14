#!/usr/bin/python


"""
#Script to Pars vcf...

#####------Inputs-------
# ./Sript.py --Workdir
# python Run-RNA-workflow-4.1_express_pars_count.py /home/gala0002/proj/proj_Ramesh/RNA-seq_2020/

"""

import sys
import re
import os
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

	def Search_CADD(self,readfl1):
		"""
		Calling Search 
		"""
		INDIR=readfl1+"SL-2400_4.0_express_trained/"
		#OUTDIR=readfl1+"NG-14833_4.1_express_pars-tpm/"
		OUTDIR=readfl1+"SL-2400_4.1_express_pars-count/"
		
		if not os.path.exists(OUTDIR):
			os.makedirs(OUTDIR)
		
		for filename in os.listdir(INDIR):
			flout = filename.split("-")[2]
			#with open(INDIR+filename+"/results.xprs",'r') as f1, open(OUTDIR+flout+"_express.tpm",'w') as output:
			with open(INDIR+filename+"/results.xprs",'r') as f1, open(OUTDIR+flout+"_express.count",'w') as output:
				first_line = f1.readline().strip()
				for g in f1:
					g1 = g.strip()
					gg = g1.split("\t")
					#gg_count = gg[1]+"\t"+gg[-1]
					gg_count = gg[1]+"\t"+gg[4]
					output.write(str(gg_count+"\n"))
			print "Done parsing...."
		return None

# file1: Input positions SNP
# write1: Output file
clF1 = SearchDB().Search_CADD(sys.argv[1])







