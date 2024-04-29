if [ $# -ne 2 ]
  then
    echo -e "usage: ./TPM.sh <gene ID> <top variant location>\nFor example\n\t./TPM.sh ENSBTAG00000053969 19:42546571"
    exit 
fi

G=$1
P=$2

mkdir -p Editing_QTL/$G
for T in Testis Epididymis_head Vas_deferens
do
  paste <(grep -P "(chr|${G})" Unnormalized_TPM/UNNORMALIZED_${T}_TPM.tsv | awk '{for(i=1;i<=NF;i++)a[i][NR]=$i}END{for(i in a)for(j in a[i])printf"%s"(j==NR?RS:FS),a[i][j]}' "${1+FS=$1}" | grep -f sample_order.txt | awk -v OFS='\t' '{print $1,$2}' | sort -k1,1)  <(bcftools query -r $P -f '[%GT\n]' eQTL/${T}_full/variants.normed.vcf.gz) <(bcftools query -r $P -f '[%GT\n]' eQTL/WGS_full/variants.normed.vcf.gz) > Editing_QTL/${G}/${T}.TPM
done


