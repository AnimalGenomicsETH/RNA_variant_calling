#!/bin/bash
tp(){ awk '{for(i=1;i<=NF;i++)a[i][NR]=$i}END{for(i in a)for(j in a[i])printf"%s"(j==NR?RS:FS),a[i][j]}' "${1+FS=$1}";}
{ echo "ID age RIN" ; grep -P "(ID|age|RIN)" /cluster/work/pausch/xena/eQTL/cis_all/epi_h/epi_h_covariates.txt | tp | grep -Ff good_samples.txt - ;} | tp > covariates/Epididymis_head.fixed.txt

{ echo "ID age RIN" ; grep -P "(ID|age|RIN)" /cluster/work/pausch/xena/eQTL/cis_all/vas_d/vas_d_covariates.txt | tp | grep -Ff good_samples.txt - ;} | tp > covariates/Vas_deferens.fixed.txt

{ echo "ID age RIN" ; grep -P "(ID|age|RIN)" /cluster/work/pausch/xena/eQTL/cis_all/testis/testis_covariates.txt | tp | grep -Ff good_samples.txt - ;} | tp > covariates/Testis.fixed.txt
