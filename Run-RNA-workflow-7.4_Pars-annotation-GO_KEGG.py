#!/usr/bin/python


"""
#Script to Pars vcf...

#####------Inputs-------
# Run-RNA-workflow-7.0_Pars-annotation-DEG.py INPUT=Inputfile OUTPUT=Outputfile DIR=Directory-fullpath

#python Run-RNA-workflow-7.4_Pars-annotation-GO_KEGG.py all_DESeq2.DE_results.P5e-2_C1.0_Mock_vs_delta8K_anno_pars_head.tsv all_DESeq2.DE_results.P5e-2_C1.0_Mock_vs_delta8K_anno_pars_GO.tsv /Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/RESULTS_2020/


#############################################------
#cut -d";" --output-delimiter=$'\t' -f3,4,5,7,8,11 | 

truncate -s 0 Thaliana_PIDs_anno

for i in $(awk '{print $2}' Quinoa_vs_Thaliana_blastp-max1_proteins.tsv)
do

#echo $i
LANG=C grep -w -m 1 $i /home/gala0002/proj/proj_Asa/REF_Genome/Ref_Arabidopsis-thaliana/GCF_000001735.4_TAIR10.1_genomic.gff | cut -f9 | awk -F';' '{ split($4,a,"=");print a[2],"\t",$3,"\t",$(NF-1)}' >> Thaliana_PIDs_anno
done

awk 'NR==FNR{a[$1]=$0; next} ($2 in a){print $0,'\t',a[$2]}' Thaliana_PIDs_anno Quinoa_vs_Thaliana_blastp-max1_proteins.tsv > Quinoa_vs_Thaliana_blastp-max1_proteins_anno.tsv

####------


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
		def srchdb1(GName1):
			DIR="/Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/RESULTS_2020/result_KEGG_Delta8K.txt"
			#DIR="/Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/RESULTS_2020/result_KEGG_PMTVWT.txt"
			try:
				True
				cmdFls1 = "LANG=C grep -m 1 -w '"+str(GName1)+"' "+str(DIR)+""
				cmdFls2 =  subprocess.check_output(cmdFls1, shell=True)
				grepout1 = str("\t".join(cmdFls2.split()[2:6]))+"\t"+str(" ".join(cmdFls2.split()[6:]))
			except:
				False
				grepout1 = str("."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+".")
			return grepout1

		def srchdb2(GName2):
			DIR="/Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/RESULTS_2020/result_Sign_GOSLIM_Delta8K_101559409.1.txt"
			#DIR="/Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/RESULTS_2020/result_Sign_GOSLIM_PMTVWT_278461476.1.txt"
			try:
				True
				cmdFls1 = "LANG=C grep -w '"+str(GName2)+"' "+str(DIR)+""
				cmdFls2 =  subprocess.check_output(cmdFls1, shell=True)
				#print re.split('\t|\n',cmdFls2)
				grepout2 = cmdFls2.split("\n") #re.split('\t|\n',cmdFls2) # #str("\t".join(cmdFls2.strip().split("\t")[0:9]))
			except:
				False
				grepout2 = str("."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+".")
			return grepout2


		with open(workdir+readfl1,'r') as f1, open(workdir+outfl1,'w') as output:
			first_line0 = f1.readline().strip().split("\t")
			first_lines = f1.readlines()[0:]
			output.write(str(str("\t".join(first_line0))+"\t"+"#KO_ID	Thrshld	Score	E-value	KO-definition	#GOSLIM_GO(P)ID	GOSLIM_GO(P)NAME	GOSLIM_GO(F)ID	GOSLIM_GO(F)NAME	GOSLIM_GO(C)ID	GOSLIM_GO(C)NAME"+"\n"))
			for lns in first_lines:
				lns_sp =  lns.strip().split("\t")
				lns_sp1 =  lns_sp[0]
				#print lns_sp1
				Anno_out1 = srchdb1(lns_sp1)
				Anno_out2 = srchdb2(lns_sp1)
				#Anno_out2_pars = re.split('\t|\n',Anno_out2)#Anno_out2
				list_P=[]; list_PN=[];list_F=[]; list_FN=[];list_C=[]; list_CN=[]
				for lns1 in Anno_out2:
					True
					lns_pars= lns1.split("\t")[:-1]
					try:
						if lns_pars != []:
							if lns_pars[1] == "P":
								list_P.append(lns_pars[0])
								list_PN.append(lns_pars[2])
						else:
							list_P.append(".")
							list_PN.append(".")
						if lns_pars != []:
							if lns_pars[1] == "F":
								list_F.append(lns_pars[0])
								list_FN.append(lns_pars[2])
						else:
							list_F.append(".")
							list_FN.append(".")
						if lns_pars != []:
							if lns_pars[1] == "C":
								list_C.append(lns_pars[0])
								list_CN.append(lns_pars[2])
						else:
							list_C.append(".")
							list_CN.append(".")
					except IndexError:
						#print lns_pars
						True

				#print str(";".join(list_F))
				Anno_out2_pars = str(";".join(list_P))+"\t"+str(";".join(list_PN))+"\t"+str(";".join(list_F))+"\t"+str(";".join(list_FN))+"\t"+str(";".join(list_C))+"\t"+str(";".join(list_CN))
				#print Anno_out1, Anno_out2_pars
				Anno_outall = str("\t".join(lns_sp)+"\t"+Anno_out1+"\t"+Anno_out2_pars+"\n")
				output.write(str(Anno_outall))

		print "Done seach for ..."
		return None

clF1 = SearchDB().Search_ReMM(sys.argv[1],sys.argv[2],sys.argv[3])





