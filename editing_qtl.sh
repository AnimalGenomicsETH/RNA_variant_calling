G="ENSBTAG00000020116"
T="Epididymis_head"
V="23:27375827-29375827"

mkdir $G
zgrep -P "(chr|$G)" ../aligned_genes/${T}.bed.gz | bgzip -c > ${G}/${T}.bed.gz
tabix -p bed ${G}/${T}.bed.gz
 
QTLtools cis --vcf ../eQTL/${T}_full/variants.normed.vcf.gz --bed ${G}/${T}.bed.gz --cov ../covariates/${T}.full.${T}.full.txt.gz --normal --nominal 1 --window 1000000 --log ${G}/${T}.nominal.log --silent --out ${G}/${T}.nominal.txt
bcftools query -f '%ID\t%INFO/DR2\n' -r $V ../eQTL/${T}_full/variants.normed.vcf.gz > ${G}/${T}.DR2

QTLtools cis --vcf ../eQTL/WGS_full/variants.normed.vcf.gz --bed ${G}/${T}.bed.gz --cov ../covariates/WGS.full.${T}.full.txt.gz --normal --nominal 1 --window 1000000 --log ${G}/WGS.nominal.log --silent --out ${G}/WGS.nominal.txt
bcftools query -f '%ID\t%INFO/DR2\n' -r $V ../eQTL/WGS_full/variants.normed.vcf.gz > ${G}/WGS.DR2
