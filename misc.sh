awk -v OFS='\t' '{print $1,$2,$2,$3,$4}' Unrevised.none.txt > Unrevised.none.bed
awk -v OFS='\t' 'NR>1&&$23 {split($8,a,"_"); print a[1],a[2],a[3],a[4]}' /cluster/work/pausch/xena/eQTL/cis_all/testis/maf01/conditional/FDR_05/all.txt > eVariants.best.txt



zcat Bos_taurus.ARS-UCD1.2.108.chr.gtf.gz | awk -v OFS="\t" '$1~/^[0-9]/&&$3=="gene" { printf "%s\t%s\t%s\t",$1,$4,$5; for (i=9; i<=NF; ++i) { if ($i=="gene_id" || $i=="gene_biotype") printf "%s\t",$(i+1) } printf "+\n"}' |  sed 's/"//g; s/;//g'

bedtools complement -L -i /cluster/work/pausch/alex/REF_DATA/ARS.108.gtf.bed -g /cluster/work/pausch/alex/REF_DATA/ARS-UCD1.2_Btau5.0.1Y.fa.fai | awk '{printf "%s\tintergenic\tintergenic\t+\n",$0}'


cat /cluster/work/pausch/xena/eQTL/cis_all/epi_h/epi_h_samples_good.txt /cluster/work/pausch/xena/eQTL/cis_all/testis/testis_samples_good.txt /cluster/work/pausch/xena/eQTL/cis_all/vas_d/vas_d_samples_good.txt | awk '{print $1}' | sort | uniq -c | awk '$1==3 {print $2}'
