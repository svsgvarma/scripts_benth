#!/usr/bin/python


"""
#Script to Pars vcf...

#####------Inputs-------
# Run-Pars_KEGG_output.py INPUT=Inputfile OUTPUT=Outputfile DIR=Directory-fullpath

python Run-Pars_KEGG_output_pars_BLAST.py all_DESeq2.DE_results.P5e-2_C1.0_Mock_vs_PMTVWT_anno_pars_head.tsv all_DESeq2.DE_results.P5e-2_C1.0_Mock_vs_PMTVWT_anno_pars_KEGG-In.tsv /home/gala0002/proj/proj_Ramesh/RNA-seq_2020/RESULT_2020_v3/

##

python Run-Pars_KEGG_output_pars_BLAST.py all_DESeq2.DE_results.P5e-2_C1.0_Mock_vs_PMTVWT_anno_pars.tsv all_DESeq2.DE_results.P5e-2_C1.0_Mock_vs_PMTVWT_anno_pars_KEGG-InFC.tsv /home/gala0002/proj/proj_Ramesh/RNA-seq_2020/RESULTS_2020_v3/
python Run-Pars_KEGG_output_pars_BLAST.py all_DESeq2.DE_results.P5e-2_C1.0_Mock_vs_delta8K_anno_pars.tsv all_DESeq2.DE_results.P5e-2_C1.0_Mock_vs_delta8K_anno_pars_KEGG-InFC.tsv /home/gala0002/proj/proj_Ramesh/RNA-seq_2020/RESULTS_2020_v3/

python2.7 Run-Pars_KEGG_output_pars_BLAST_PMTVWT.py all_DESeq2.DE_results.PMTVWT-vs-Mock_anno.tsv all_DESeq2.DE_results.PMTVWT-vs-Mock_anno_KEGG-In.tsv /home/gala0002/proj/proj_Ramesh/RNA-seq_2020/RESULTS_2020_v4/
python2.7 Run-Pars_KEGG_output_pars_BLAST_PMTVWT.py all_DESeq2.DE_results.P5e-2_C1.0_PMTVWT-vs-Mock_anno.tsv all_DESeq2.DE_results.P5e-2_C1.0_PMTVWT-vs-Mock_anno_KEGG-In.tsv /home/gala0002/proj/proj_Ramesh/RNA-seq_2020/RESULTS_2020_v4/

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
		def srchdb1(GName):
			DIR="/home/gala0002/proj/proj_Ramesh/RNA-seq_2020/KEGG_KofamKOALA/Ref_Ntab-TN90/N.benQ-PMTVWT_vs_N.tab_blastp-max1_proteins.tsv"
			#DIR="/home/gala0002/proj/proj_Ramesh/RNA-seq_2020/KEGG_KofamKOALA/Ref_Ntab-TN90/N.benQ-Delta8K_vs_N.tab_blastp-max1_proteins.tsv"
			try:
				True
				cmdFls1 = "LANG=C grep -m 1 -w '"+str(GName)+"' "+str(DIR)+""
				cmdFls2 = subprocess.check_output(cmdFls1, shell=True)
				grepout = cmdFls2.strip()
			except:
				False
				grepout = "."
			return grepout

		def srchdb2(GName):
			DIR="/home/gala0002/proj/proj_Ramesh/RNA-seq_2020/KEGG_KofamKOALA/Ref_Ntab-TN90/N.benQ-PMTVWT_vs_N.tab_blastp-max1_proteins_XP-GID.tsv"
			#DIR="/home/gala0002/proj/proj_Ramesh/RNA-seq_2020/KEGG_KofamKOALA/Ref_Ntab-TN90/N.benQ-Delta8K_vs_N.tab_blastp-max1_proteins_XP-GID.tsv"
			try:
				True
				cmdFls1 = "LANG=C grep -m 1 -w '"+str(GName)+"' "+str(DIR)+""
				cmdFls2 = subprocess.check_output(cmdFls1, shell=True)
				grepout = cmdFls2.strip()
			except:
				False
				grepout = "."
			return grepout

		with open(workdir+readfl1,'r') as f1, open(workdir+outfl1,'w') as output:
			output.write(str("Gene"+"\t"+"FC"+"\t"+"P.Value"+"\t"+"FDR"+"\n"))

			first_line = f1.readlines()[1:]
			for lns in first_line:
				ln_all = lns.strip().split()
				ln_gid = ln_all[0]
				Anno_out1 = srchdb1(ln_gid)
				if Anno_out1 != ".":
					True
					Anno_XP = Anno_out1.split("\t")[1]
					Anno_out2 = srchdb2(Anno_XP)
					if Anno_out2 != ".":
						True
						Anno_out2_GID = Anno_out2.split(",")[0].split(":")[1]
						Anno_out_all = str(Anno_out2_GID+"\t"+ln_all[6]+"\t"+ln_all[9]+"\t"+ln_all[10])
						output.write(Anno_out_all+"\n")

		print "Done seach for ..."
		return None


clF1 = SearchDB().Search_ReMM(sys.argv[1],sys.argv[2],sys.argv[3])





