

#a
while read G P  
do
  echo $G $P
  for T in Testis Epididymis_head Vas_deferens
  do
    paste <(grep -P "(chr|${G})" /cluster/work/pausch/alex/RNA_call_test/Unnormalized_TPM/UNNORMALIZED_${T}_TPM.tsv | awk '{for(i=1;i<=NF;i++)a[i][NR]=$i}END{for(i in a)for(j in a[i])printf"%s"(j==NR?RS:FS),a[i][j]}' "${1+FS=$1}" | grep -f /cluster/work/pausch/alex/RNA_call_test/sample_order.txt | awk -v OFS='\t' '{print $1,$2}' | sort -k1,1)  <(bcftools query -r $P -f '[%GT\n]' /cluster/work/pausch/alex/RNA_call_test/eQTL/${T}_full/variants.normed.vcf.gz) <(bcftools query -r $P -f '[%GT\n]' /cluster/work/pausch/alex/RNA_call_test/eQTL/WGS_full/variants.normed.vcf.gz) > Editing_QTL/${G}/${T}.TPM
  done
done < <(echo -e "ENSBTAG00000053969 19:42546571\nENSBTAG00000027962 2:111893837\nENSBTAG00000000597 1:18578966")


