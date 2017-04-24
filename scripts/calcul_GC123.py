#!/usr/bin/python2.7.6
# -*- coding: utf-8 -*-
import os
import sys
import csv
from Bio import SeqIO
from Bio.SeqUtils import GC123
from Bio.SeqUtils import GC

def parsing_fasta(fichier): #cette fonction return une liste de seq génomiques
	li = []
	if os.path.getsize(fichier) != 0:
		handle = open(fichier, "rU")
		for record in SeqIO.parse(handle, "fasta") :
			li.append(str(record.seq))
		handle.close()
	else:
		pass
	return li


def CDS_Conformity(CDS): #Elle vérifie si le CDS est valide
	Start_Codon = 'ATG'
	Stop_Codon = ['TAA', 'TAG', 'TGA']

	if CDS[len(CDS)-3 : len(CDS)] in Stop_Codon and len(CDS)  % 3 == 0 :
		return(True)
	
	else :
	
		return (False)

#print(parsing_fasta("exfast.txt"))
def parsing_genom_CDS():
	chain = ""
	liste_gen = parsing_fasta("tmp.txt.genome")
	if liste_gen:
		##recuperer les TaxID
		with open('tmp.txt') as f:
			line = f.readline()
			line = line.rstrip()
			line = int(float(line))
		len_gen = 0
		N = 0
		gc = 0
		for i in liste_gen:
			len_gen = len_gen + len(str(i)) #lengueur des chromosomes
			for j in str(i):
				if j == 'N':
					N = N + 1 #nbr de N dans la seq génomique
	
				if j in ['G', 'C']:
					gc = gc + 1
		
		len_gen_valid = len_gen - N
	
		gc = gc*100/float(len_gen_valid) #taux de GC genome
		N = N/ float(len_gen) #%NA
		#print(gc, N, len_gen, len_gen_valid)

	####CDS
		li = parsing_fasta("tmp.txt.cds")
		if li:		
			lis = []
			a = 0 #taux GC CDSvalid
			b = 0 #taux GC1 CDSvalid
			c = 0 #taux GC2 CDSvalid
			d = 0 #taux GC3 CDSvalid
			LenCumCDS = 0
			NbrCDS_valid = 0
			LenCumCDS_valid = 0
			j = 0
			#NbrCDS = len(li)
			for i in li:
				i = str(i)
				j = j +1
				LenCumCDS = LenCumCDS + len(i)
				if CDS_Conformity(i):
					NbrCDS_valid += 1
					LenCumCDS_valid +=len(i)
					lis = GC123(i)
					a +=lis[0]
					b +=lis[1]
					c +=lis[2]
					d +=lis[3]
				else:
					pass
			NbrCDS = j
		#print(j, NbrCDS_valid, LenCumCDS, LenCumCDS_valid)
		#print(a/len(li),b/len(li),c/len(li),d/len(li))
	
			lst = [line, len_gen, round(gc), round(N), NbrCDS, LenCumCDS, round(a/len(li)), NbrCDS_valid, LenCumCDS_valid, round(b/len(li)), round(c/len(li)), round(d/len(li))]
			lst = map(str, lst)
			chain = "\t".join(lst)
			chain = chain + '\n'
	return chain


#print(parsing_genom_CDS())


#writer = csv.writer(sys.stdout, delimiter='\t', quotechar='"', quoting=csv.QUOTE_ALL)
sys.stdout.write(parsing_genom_CDS())
