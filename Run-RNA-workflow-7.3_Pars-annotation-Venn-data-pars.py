#!/usr/bin/python


"""
#Script to Pars vcf...

#####------Inputs-------
# Run-RNA-workflow-7.0_Pars-annotation-DEG.py INPUT=Inputfile OUTPUT=Outputfile DIR=Directory-fullpath

#python Run-RNA-workflow-7.3_Pars-annotation-Venn-data-pars.py venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_PasvTit_Inter_RegvTit_145_TIT_anno_uniq.tsv venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_PasvTit_Inter_RegvTit_145_TIT_anno_uniq_pars.tsv /Volumes/Seagate/Backup-MAC_HD2/proj_Asa/NG-14833_6.0_Venn-diagram/
#python Run-RNA-workflow-7.3_Pars-annotation-Venn-data-pars.py venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_PasvTit_Inter_PasvReg_2143_PAS_anno_uniq.tsv venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_PasvTit_Inter_PasvReg_2143_PAS_anno_uniq_pars.tsv /Volumes/Seagate/Backup-MAC_HD2/proj_Asa/NG-14833_6.0_Venn-diagram/
#python Run-RNA-workflow-7.3_Pars-annotation-Venn-data-pars.py venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_PasvReg_Inter_RegvTit_165_REG_anno_uniq.tsv venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_PasvReg_Inter_RegvTit_165_REG_anno_uniq_pars.tsv /Volumes/Seagate/Backup-MAC_HD2/proj_Asa/NG-14833_6.0_Venn-diagram/

#python Run-RNA-workflow-7.3_Pars-annotation-Venn-data-pars.py venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_PasvTit_Inter_PasvReg_Inter_RegvTit_PAS-TIT-REG_anno_uniq.tsv venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_PasvTit_Inter_PasvReg_Inter_RegvTit_PAS-TIT-REG_anno_uniq_pars.tsv /Volumes/Seagate/Backup-MAC_HD2/proj_Asa/NG-14833_6.0_Venn-diagram/

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

	def Search_all(self,readfl1,outfl1,workdir):
		"""
		Calling Search local DB
		"""
		with open(workdir+readfl1,'r') as f1, open(workdir+outfl1,'w') as output:
			first_line0 = f1.readline().strip().split("\t")
			first_lines = f1.readlines()[0:]
			output.write(str("\t".join(first_line0)+"\t"+"Gene-Anno"+"\n"))

			for lns in first_lines:
				TF_fam=[]
				lns_sp =  lns.strip().split("\t")
				#lns_sp1 =  lns_sp[0]
				lns_sp = [x.strip(' ') for x in lns_sp]
				#print lns_sp# lns_sp[46], lns_sp[49]
				if ((lns_sp[41] != "." ) and (lns_sp[44] == "." )):
					#Anno_out3 = srchdb3(PIDhit1)
					TF_fam.append(str("_".join(lns_sp[41].split("_")[1:])))
					#print str(lns_sp[41].split("_")[1]+"_family[PlantTFDB/ITAK]")
				elif ((lns_sp[44] != "." ) and (lns_sp[41] == "." )):
					#Anno_out3 = srchdb3(PIDhit1)
					TF_fam.append(str(lns_sp[44].split("|")[1]))
					#print str(lns_sp[44].split("|")[1]+"_family[PlantTFDB/ITAK]")
				elif ((lns_sp[41] != "." ) and (lns_sp[44] != "." )):
					#Anno_out3 = srchdb3(PIDhit1)
					db1=str("_".join(lns_sp[41].split("_")[1:]))
					db2=str(lns_sp[44].split("|")[1])
					#print str(db1+"/"+db2)
					TF_fam.append(str(db1+"/"+db2))
					#print str(lns_sp[44].split("|")[1]+"_family[PlantTFDB/ITAK]")

				if ((lns_sp[29] == "." )):
					GI_Desc = lns_sp[29].replace(".","NA")+"("+lns_sp[28]+")"
				else:
					GI_Desc = lns_sp[29]+"("+lns_sp[28]+")"
				
				if (TF_fam != []):
				 	True
					#print GI_Desc+"@"+str(TF_fam[0])+str("_family[PlantTFDB/ITAK]") 
					Anno_outall = str("\t".join(lns_sp))+"\t"+str(GI_Desc)+"@"+str(TF_fam[0])+str("_family[PlantTFDB/ITAK]")+"\n"
				else:
					True
					#print GI_Desc
					Anno_outall = str("\t".join(lns_sp))+"\t"+str(GI_Desc)+"\n"
				output.write(str(Anno_outall))
		print "Done seach for ..."
		return None

clF1 = SearchDB().Search_all(sys.argv[1],sys.argv[2],sys.argv[3])





