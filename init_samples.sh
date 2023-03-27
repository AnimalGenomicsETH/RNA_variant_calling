
for ext in bam bam.bai
do
  for i in $(cat good_samples.txt); do ln -s /nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/BTA/eQTL/testis/RNA_alignments/normal/dedup_alignment/$i/${i}.${ext} subsampled_bams/Testis/${i}.full.${ext}; done
  for i in $(cat good_samples.txt); do ln -s /nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/BTA/eQTL/epi_h/RNA_alignments/normal/dedup_alignment/${i}EH/${i}EH.${ext} subsampled_bams/Epididymis_head/${i}.full.${ext}; done
  for i in $(cat good_samples.txt); do ln -s /nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/BTA/eQTL/vas_d/RNA_alignments/normal/dedup_alignment/${i}V/${i}V.${ext} subsampled_bams/Vas_deferens/${i}.full.${ext}; done
  for i in $(cat good_samples.txt); do ln -s /nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/BTA/eQTL/testis/RNA_alignments/normal/dedup_alignment/${i}/${i}.${ext} subsampled_bams/WGS/${i}.full.${ext}; done
done
