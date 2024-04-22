

# 
paste <(grep -P "(chr|ENSBTAG00000053969)" /cluster/work/pausch/xena/eQTL/gene_counts/testis/UNNORMALIZED_testis_TPM.tsv | awk '{for(i=1;i<=NF;i++)a[i][NR]=$i}END{for(i in a)for(j in a[i])printf"%s"(j==NR?RS:FS),a[i][j]}' "${1+FS=$1}" | grep -f /cluster/work/pausch/alex/RNA_call_test/sample_order.txt | awk -v OFS='\t' '{print $1,$2}')  <(bcftools query -r 19:42546571 -f '[%GT\n]' /cluster/work/pausch/alex/RNA_call_test/eQTL/Testis_full/variants.normed.vcf.gz) > Editing_QTL/ENSBTAG00000053969/Testis.TPM

paste <(grep -P "(chr|ENSBTAG00000053969)" /cluster/work/pausch/xena/eQTL/gene_counts/epi_h/UNNORMALIZED_epi_h_TPM.tsv | awk '{for(i=1;i<=NF;i++)a[i][NR]=$i}END{for(i in a)for(j in a[i])printf"%s"(j==NR?RS:FS),a[i][j]}' "${1+FS=$1}" | grep -f /cluster/work/pausch/alex/RNA_call_test/sample_order.txt | sed 's/EH//' | awk -v OFS='\t' '{print $1,$2}')  <(bcftools query -r 19:42546571 -f '[%GT\n]' /cluster/work/pausch/alex/RNA_call_test/eQTL/Epididymis_head_full/variants.normed.vcf.gz) > Editing_QTL/ENSBTAG00000053969/Epididymis_head.TPM

paste <(grep -P "(chr|ENSBTAG00000053969)" /cluster/work/pausch/xena/eQTL/gene_counts/vas_d/UNNORMALIZED_vas_d_TPM.tsv | awk '{for(i=1;i<=NF;i++)a[i][NR]=$i}END{for(i in a)for(j in a[i])printf"%s"(j==NR?RS:FS),a[i][j]}' "${1+FS=$1}" | grep -f /cluster/work/pausch/alex/RNA_call_test/sample_order.txt | sed 's/V//' | awk -v OFS='\t' '{print $1,$2}')  <(bcftools query -r 19:42546571 -f '[%GT\n]' /cluster/work/pausch/alex/RNA_call_test/eQTL/Vas_deferens_full/variants.normed.vcf.gz) > Editing_QTL/ENSBTAG00000053969/Vas_deferens.TPM
