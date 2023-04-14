rule qtltools_quan:
    input:
        gtf = '',
        bam = ''
    output:
        stat = '',
        tpm = ''
    shell:
        '''
         QTLtools quan --bam {input.bam} --gtf {input.gtf} --sample {wildcards.sample} --out-prefix {params.prefix} --filter-mapping-quality 60 --tpm --check-proper-pairing --filter-failed-qc
         '''

rule combine_quan:
    input:
        ''
    output:
        ''
    localrule: True
    shell:
        '''
        input_path=/cluster/work/pausch/xena/eQTL/RNA_alignment/vas_d/star_aligned
 59 output_path=/cluster/work/pausch/xena/eQTL/gene_counts/all_samples/vas_d
 60 
 61 cat ${output_path}/BSWCHEM120145055733.LmWzrHf25Mj.gene.tpm.bed  | cut -f1-6 > ${output_path}/anno_info
 62 for bam in `ls ${input_path}`
 63 do  
 64     cat ${output_path}/${bam}.R0BRiR1nnu.gene.tpm.bed | cut -f7 > ${output_path}/${bam}_temp_quant
 65 done
 66 
 67 paste -d "\t" ${output_path}/anno_info ${output_path}/*_temp_quant > ${output_path}/FINAL_gene_TPM.tsv
        '''

rule featurecounts:
    input:
        gtf = '',
        bams = ''
    output:
        ''
    threads: 2
    resources:
        mem_mb = 2500,
        walltime = '4h'
    shell:
        '''
        /cluster/home/xmapel/miniconda3/envs/featurecounts/bin/featureCounts -O -M -p -T {threads} -t exon -g gene_id --fraction -Q 60 -s 2 --primary --tmpDir $TMPDIR -a {input.gft} -o {output[0]} {input.bams}
        '''

rule filter_TPM:
    input:
        quan = rules.combine_quan.output,
        FC = rules.featurecounts.output
    output:
        ''
    shell:
        '''
        testis_tpm <- read.table("/cluster/work/pausch/xena/eQTL/gene_counts/testis/UNFILTERED_gene_TPM.tsv", header=T)
  5 fc_testis <- read.table("/cluster/work/pausch/xena/eQTL/gene_counts/testis/featurecounts/featurecounts_testis", header =T)
  6
  7 #Raw TPM
  8 ##Testis filtering
  9 remove <- c("BSWCHEM120119981594", "BSWCHEM120134766497", "BSWCHEM120135927712", "BSWCHEM120144952301", "BSWCHEM120146480749", "BSWCHEM120149770120", "BSWCHEM120150918900", "BSWCHEM120153327877", "BSWCHEM120153643403", "BSWCHEM120155307242", "BSWCHEM120151536851", "BSWDEUM00955    2586445")
 10
 11 testis_tpm = testis_tpm[,!(names(testis_tpm) %in% remove)]
 12 testis_tpm <- testis_tpm[rowSums(testis_tpm[,c(7:123)] >= 0.1) / (ncol(testis_tpm)-6) >= 0.2,]
 13 fc_testis_good <- fc_testis[rowSums(fc_testis[,c(7:123)] >= 6) / (ncol(fc_testis)-6) >= 0.2,]
 14 testis_tpm = testis_tpm[testis_tpm$gene %in% fc_testis_good$Geneid,]
        '''

rule qtltools_pca:
    input:
        lambda wildcards: 'bed' if wildcards.mode == 'bed' else 'vcf'
    output:
        multiext('covariates/{tissue}.{coverage}.{mode}','.pca','.pca_stats')
    params:
        prefix = lambda wildcards, output: PurePath(output[0]).with_suffix('')
    shell:
        '''
        QTLtools pca --{wildcards.mode} {input} --out {params.prefix} --center --scale
        '''

rule make_covariates:
    input:
        ''
    output:
        ''
    shell:
        '''
        #ADD age and RIN (fixed)
        '''


# STEP 4: Filtering (R)
##### let me remake this so itll work in snakemake #####


# STEP 5: PEER Factors 


# STEP 6: 
