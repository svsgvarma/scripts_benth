#!/usr/bin/python


"""
#Script to Pars vcf...

#####------Inputs-------
# Run-RNA-workflow-7.0_Pars-annotation-DEG.py INPUT=Inputfile OUTPUT=Outputfile DIR=Directory-fullpath

#python Run-RNA-workflow-7.0_Pars-annotation-DEG_CDS-DEseq2.py Express_counts_all_Mock-PMTVWT-delta8k.matrix.Mock_vs_delta8K.DESeq2.DE_results.P1e-2_C1.0.DE.subset_head Express_counts_all_Mock-PMTVWT-delta8k.matrix.Mock_vs_delta8K.DESeq2.DE_results.P1e-2_C1.0.DE.subset_head_anno.tsv /home/gala0002/proj/proj_Ramesh/RNA-seq_2020/DESeq2_genes/


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
			DIR="/home/gala0002/proj/proj_Ramesh/REF/Ref_Nicotiana_benthamiana/Supplemental_dataset_5_NbC_Anno_noseq.txt"
			try:
				True
				cmdFls1 = "LANG=C grep -m 1 -w '"+str(GName1)+"' "+str(DIR)+""
				cmdFls2 =  subprocess.check_output(cmdFls1, shell=True)
				grepout1 = str("\t".join(cmdFls2.strip().split("\t")))
			except:
				False
				grepout1 = str("."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+".")
			return grepout1

		def srchdb2(GName2):
			DIR="/home/gala0002/proj/proj_Ramesh/RNA-seq_2020/SL-2400_5.0_Annotation/NbC-prot_blastp-max1_PlantTFDBv5.0_Ext-NbC-pep.tsv"
			try:
				True
				cmdFls1 = "LANG=C grep -m 1 -w '"+str(GName2)+"' "+str(DIR)+""
				cmdFls2 =  subprocess.check_output(cmdFls1, shell=True)
				ann = cmdFls2.strip().split("\t")
				#print AT_ann
				grepout2 = str(ann[1]+"\t"+ann[-2]+"\t"+ann[2])
			except:
				False
				grepout2 = str("."+"\t"+"."+"\t"+".")
			return grepout2

		def srchdb3(GName3):
			DIR="/home/gala0002/proj/proj_Ramesh/RNA-seq_2020/SL-2400_5.0_Annotation/NbC-cds_blastn-max1_ITAKv18.12_Ext-NbC-cds.tsv"
			try:
				True
				cmdFls1 = "LANG=C grep -m 1 -w '"+str(GName3)+"' "+str(DIR)+""
				cmdFls2 =  subprocess.check_output(cmdFls1, shell=True)
				ann = cmdFls2.strip().split("\t")
				#print AT_ann
				grepout3 = str(ann[1]+"\t"+ann[-2]+"\t"+ann[2])
			except:
				False
				grepout3 = str("."+"\t"+"."+"\t"+".")
			return grepout3

		with open(workdir+readfl1,'r') as f1, open(workdir+outfl1,'w') as output:
			first_line0 = f1.readline().strip().split("\t")
			first_lines = f1.readlines()[1:]
			output.write(str("NB_GeneID"+"\t"+str("\t".join(first_line0))+"\t"+"Name	ApoplastP.1.0.1_Prediction	SignalP	ApoplastP_SignalP	PFAM	GENENAME	DESCRIPTION	ENZYME	GO(P)ID	GO(P)NAME	GO(F)ID	GO(F)NAME	GO(C)ID	GO(C)NAME	KEYWORD	PATHWAY	GOSLIM	ANNOTATOR	PlantTFDBv5.0_Ext-NbC-pep_TFfamID	PlantTFDBv5.0_Evalue	PlantTFDBv5.0_Identity_Score	ITAKv18.12_Ext-NbC-cds_TFfamID	ITAKv18.12_Evalue	ITAKv18.12_Identity_Score"+"\n"))
			for lns in first_lines:
				lns_sp =  lns.strip().split("\t")
				lns_sp1 =  lns_sp[0]
				Anno_out1 = srchdb1(lns_sp1)
				Anno_out2 = srchdb2(lns_sp1)
				Anno_out3 = srchdb3(lns_sp1)
				#print Anno_out2, Anno_out3
				Anno_outall = str("\t".join(lns_sp)+"\t"+Anno_out1+"\t"+Anno_out2+"\t"+Anno_out3+"\n")
				output.write(str(Anno_outall))

		print "Done seach for ..."
		return None

clF1 = SearchDB().Search_ReMM(sys.argv[1],sys.argv[2],sys.argv[3])





