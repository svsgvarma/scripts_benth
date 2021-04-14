#!/usr/bin/python


"""
#Script to Pars vcf...

#####------Inputs-------
# Run-Pars_KEGG_output.py INPUT=Inputfile OUTPUT=Outputfile DIR=Directory-fullpath

python Run-Pars_KEGG_output.py result_KEGG_PMTVWT.txt result_KEGG_PMTVWT_pars.txt /Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/RESULTS_2020_v3/KEGG_KofamKOALA/

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

"""
#Query_ID, Subject_ID, %Identity, Alignment_length, Mismatches, Gap_opens, Q.start, Q.end, S.start, S.end, E-value, Bit-score
#Query_ID	Subject_ID	%Identity	Alignment_length	Mismatches	Gap_opens	Q.start	Q.end	S.start	S.end	E-value	Bit-score
"""


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

	def Search_ReMM(self,readfl1,outfl1,workdir):
		"""
		Calling Search local DB
		"""
		def srchdb1(CDSName):
			DIR="/home/gala0002/proj/proj_Asa/REF_Genome/Ref_Arabidopsis-thaliana/GCF_000001735.4_TAIR10.1_genomic.gff"
			try:
				True
				cmdFls1 = "LANG=C grep -m 1 -w '"+str(CDSName)+"' "+str(DIR)+""
				cmdFls2 = subprocess.check_output(cmdFls1, shell=True)
				grepout = cmdFls2.strip()
			except:
				False
			return grepout

		with open(workdir+readfl1,'r') as f1, open(workdir+outfl1,'w') as output:
			first_line = f1.readlines()
			for lns in first_line:
				ln_all = lns.strip().split()
				Anno_out = str(ln_all[2]+"\t"+ln_all[1]+"\t"+str("_".join(ln_all[6:]))) #
				print Anno_out
				output.write(Anno_out+"\n")

		print "Done seach for ..."
		return None


clF1 = SearchDB().Search_ReMM(sys.argv[1],sys.argv[2],sys.argv[3])





