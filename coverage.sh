for T in Epididymis_head Testis Vas_deferens
do
  parallel -j 8 "samtools flagstat --threads 2 ${T}/{}.full.bam | awk -v T=$T -v S={} 'NR==2 {print T,S,\$1}'" ::: $(ls Testis/*.full.bam | cut -d'/' -f 2 | cut -d'.' -f 1)
done
