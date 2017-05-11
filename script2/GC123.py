# -*- coding: utf-8 -*-

import os
import sys
from Bio import SeqIO
import re



def parsing_fasta(fichier): 
	''' Prends un fichier fasta et return une liste des séquences'''
	li = []

	handle = open(fichier, "rU")
	for record in SeqIO.parse(handle, "fasta") :
		li.append(str(record.seq))
	handle.close()

	return li


def CDS_Conformity(CDS): 
	'''Vérifier la validité d'une séquence codante, return Vrai ou Faux'''
	Start_Codon = 'ATG'
	Stop_Codon = ['TAA', 'TAG', 'TGA']

	if CDS[len(CDS)-3 : len(CDS)] in Stop_Codon and len(CDS)  % 3 == 0 :
		return(True)
	
	else :
	
		return (False)


#print(parsing_fasta("exfast.txt"))

def parsing_genom_CDS(tax_id, genome, cds):
	'''calcul le taux de GC123 des séquences'''
	
	chain = ""
	liste_gen = parsing_fasta(genome)
	if liste_gen:
		##recuperer les TaxID
		tax_id = str(tax_id)
		tax_id = tax_id.replace('\n', '')
		tax_id = tax_id[:-7]
		
			
		len_gen = 0 
		N = 0
		gc = 0.0
		for i in liste_gen:
			len_gen = len_gen + len(i) #lengueur des chromosomes
			for j in str(i):
				if j == 'N':
					N = N + 1 #nbr de N dans la seq génomique
	
				if j in ['G', 'C', 'g', 'c']:
					gc = gc + 1
		
		len_gen_valid = len_gen - N
		taux_GC_genome = gc*100/len_gen_valid
	

	####CDS
		li = parsing_fasta(cds)
		
		if li:		
			gc_CDS_tot = 0.0
			LenCumCDS = 0 #fait
			NbrCDS_valid = 0
			LenCumCDS_valid = 0
			gc_CDS_tot_valid = 0.0
			gc1 = 0.0
			gc2 = 0.0
			gc3 = 0.0
			NbrCDS = len(li) #nbre de CDS
			for i in li:
				i = str(i)
				
				LenCumCDS = LenCumCDS + len(i)
				gc_CDS_tot += (i.count('G') + i.count('C') + i.count('g') + i.count('c'))
				if CDS_Conformity(i):
					NbrCDS_valid += 1
					gc_CDS_tot_valid += i.count('G') + i.count('C')  + i.count('g') + i.count('c')
					LenCumCDS_valid += len(i)
					
					for j in range(2,len(i),3):
						if i[j] in ['G', 'C', 'g', 'c']:
							gc3 = gc3 + 1
						if i[j-1] in ['G', 'C', 'g', 'c']:
							gc2 = gc2 + 1
						if i[j-2] in ['G', 'C', 'g', 'c']:
							gc1 = gc1 + 1

			taux_gc_CDS_tot = gc_CDS_tot*100/LenCumCDS
			taux_gc_CDS_tot_valid = gc_CDS_tot_valid*100/LenCumCDS_valid

	
			lst = [tax_id, len_gen, len_gen_valid, gc, taux_GC_genome, N, NbrCDS, NbrCDS_valid,LenCumCDS, LenCumCDS_valid , gc_CDS_tot, taux_gc_CDS_tot, gc_CDS_tot_valid, taux_gc_CDS_tot_valid, gc1, gc1*100/(LenCumCDS_valid/3), gc2 , gc2*100/(LenCumCDS_valid/3), gc3 , gc3*100/(LenCumCDS_valid/3)]
			lst = map(str, lst)
			chain = "\t".join(lst)
			chain = chain + '\n'
	return chain


def main():
	tax_id = sys.argv[1]
	genome = sys.argv[2]
	cds = sys.argv[3]
	sys.stdout.write(parsing_genom_CDS(tax_id, genome, cds))


main()

