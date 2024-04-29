if [ $# -ne 3 ]
  then
    echo -e "usage: ./editing_qtl.sh <gene ID> <Tissue> <variant window>\nFor example\n\t./editing_qtl.sh ENSBTAG00000020116 Testis 23:27375827-29375827"
    exit 
fi

G=$1
T=$2
V=$3
mkdir -p Editing_QTL/$G
zgrep -P "(chr|$G)" aligned_genes/${T}.bed.gz | bgzip -c > Editing_QTL/${G}/${T}.bed.gz
tabix -p bed ${G}/${T}.bed.gz
 
QTLtools cis --vcf eQTL/${T}_full/variants.normed.vcf.gz --bed Editing_QTL/${G}/${T}.bed.gz --cov covariates/${T}.full.${T}.full.txt.gz --normal --nominal 1 --window 1000000 --log Editing_QTL/${G}/${T}.nominal.log --silent --out Editing_QTL/${G}/${T}.nominal.txt
bcftools query -f '%ID\t%INFO/DR2\n' -r $V eQTL/${T}_full/variants.normed.vcf.gz > Editing_QTL/${G}/${T}.DR2

QTLtools cis --vcf eQTL/WGS_full/variants.normed.vcf.gz --bed Editing_QTL/${G}/${T}.bed.gz --cov covariates/WGS.full.${T}.full.txt.gz --normal --nominal 1 --window 1000000 --log Editing_QTL/${G}/WGS.nominal.log --silent --out Editing_QTL/${G}/WGS.nominal.txt
bcftools query -f '%ID\t%INFO/DR2\n' -r $V eQTL/WGS_full/variants.normed.vcf.gz > Editing_QTL/${G}/WGS.DR2
