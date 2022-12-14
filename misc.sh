awk -v OFS='\t' '{print $1,$2,$2,$3,$4}' Unrevised.none.txt > Unrevised.none.bed
awk -v OFS='\t' 'NR>1&&$23 {split($8,a,"_"); print a[1],a[2],a[3],a[4]}' /cluster/work/pausch/xena/eQTL/cis_all/testis/maf01/conditional/FDR_05/all.txt > eVariants.best.txt
