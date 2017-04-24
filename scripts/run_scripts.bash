#!/usr/bin/env bash
#shuf -n 100 ../data_raw/prokaryotes.txt > exemple100.txt

#rm tmp.*
touch resultats.csv
printf '%s\n' TaxId LenGenome GCgenome %NA NbrCDS LenCumCDS GC_CDS NbrCDS_valid LenCumCDS_valid GC1 GC2 GC3 | paste -sd '\t' >> resultats.csv
#génerer le fichier prokaryotes-complete-genome.txt
#bash extract-complete-genome.bash
#optimiser le scripts calcul_CG123.py 
cython calcul_GC123.py
gcc calcul_GC123.c -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing $(pkg-config python --cflags) -o calcul_GC123.so
while IFS='' read -r line || [[ -n "$line" ]]
do
	#extrait les numéros d'accession avec ton super script
	echo "$line" | python parsing_prokaryotes.py > tmp.txt
	# exécute la commande raa_query sur tmp.txt
	bash run_raa_query.sh tmp.txt
	#exécuter le script qui calcul le taux GC...
	python calcul_GC123.py >> resultats.csv


done < ../data/prokaryotes-complete-genome.txt

##bash run_scripts.bash

