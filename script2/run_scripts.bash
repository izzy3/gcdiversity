#!/usr/bin/env bash

touch resultats.csv
printf '%s\n' TaxId LenGenome len_gen_valid GCgenome taux_GC_genome N NbrCDS NbrCDS_valid LenCumCDS LenCumCDS_valid gc_CDS_tot taux_gc_CDS_tot gc_CDS_tot_valid  taux_gc_CDS_tot_valid GC1 tauc_gc1 GC2 taux_gc2 GC3 taux_gc3  | paste -sd '\t' >> resultats.csv

while IFS='' read -r LINE
do
    python GC123.py $LINE $LINE.genome $LINE.cds >> resultats.csv
    #printf "$LINE\t$LINE.genome\t$LINE.cds\n"

done < <(ls -1 | grep 'access$')
