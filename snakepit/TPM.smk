from pathlib import PurePath
import pandas as pd
import subprocess
from tempfile import NamedTemporaryFile

rule qtltools_quan:
    input:
        gtf = config['GTF'],
        bam = 'subsampled_bams/{tissue}/{sample}.{coverage}.bam'
    output:
        temp(multiext('aligned_genes/{tissue}/{sample}.{coverage}','.gene.tpm.bed','.gene.count.bed','.exon.tpm.bed','.exon.count.bed','.stats'))
    params:
        prefix = lambda wildcards, output: PurePath(output[4]).with_suffix('')
    threads: 1
    resources:
        mem_mb = 2500,
        waltime = '30m'
    shell:
        '''
         QTLtools quan --bam {input.bam} --gtf {input.gtf} --sample {wildcards.sample} --out-prefix {params.prefix} --filter-mapping-quality 60 --tpm --check-proper-pairing --no-hash --filter-failed-qc
         '''

rule combine_quan:
    input:
        expand(rules.qtltools_quan.output[0],sample=config['samples'],allow_missing=True)
    output:
        temp('aligned_genes/{tissue}/gene_TPM.{coverage}.tsv.gz')
    params:
        slices = lambda wildcards, input:  '1-7,' + ','.join(map(str,range(14,len(input)*7+1,7)))
    localrule: True
    shell:
        '''
        paste {input} | cut -f {params.slices} | pigz -p 2 -c > {output}
        '''

rule featurecounts:
    input:
        gtf = config['GTF'],
        bams = expand('subsampled_bams/{tissue}/{sample}.{coverage}.bam',sample=config['samples'],allow_missing=True)
    output:
        temp(multiext('aligned_genes/{tissue}/genes.{coverage}.FC','','.summary'))
    threads: 16
    resources:
        mem_mb = 100,
        walltime = '1h'
    shell:
        '''
        featureCounts -O -M -p -T {threads} -t exon -g gene_id --fraction -Q 60 -s 2 --primary --tmpDir $TMPDIR -a {input.gtf} -o {output[0]} {input.bams}
        '''

rule filter_TPM:
    input:
        quan = rules.combine_quan.output[0],
        FC = rules.featurecounts.output[0]
    output:
        temp(multiext('aligned_genes/{tissue}/gene_counts.{coverage}.bed.gz','','.tbi'))
    localrule: True
    run:
        TPM = pd.read_csv(input.quan,delimiter='\t',low_memory=False)
        good_TPM = TPM[((TPM.iloc[:,7:]>= 0.1).sum(axis=1) / (TPM.shape[1]-6)) >= 0.2]
        FC = pd.read_csv(input.FC,delimiter='\t',low_memory=False,skiprows=1)
        good_FC = FC[((FC.iloc[:,7:] >= 6).sum(axis=1) / (FC.shape[1]-6)) >= 0.2]
        with NamedTemporaryFile(mode='w') as fout:
            good_TPM[good_TPM['gene'].isin(good_FC['Geneid'])].sort_values(['#chr', 'start']).to_csv(fout,index=False,sep='\t')
            subprocess.run(f'bgzip -c {fout.name} > {output[0]}; tabix -p bed {output[0]}', shell=True)

rule bcftools_prune:
    input:
        '{tissue}_{coverage}/autosomes.{imputed}.vcf.gz'
    output:
        temp(multiext('{tissue}_{coverage}/autosomes.{imputed}.pruned.vcf.gz','','.tbi'))
    resources:
        mem_mb = 2500
    shell:
        '''
        bcftools +prune -o {output[0]} -w 1000kb -n 5 -m r2=0.2 {input[0]}
        tabix -fp vcf {output[0]}
        '''

rule qtltools_pca:
    input:
        lambda wildcards: rules.filter_TPM.output if wildcards.mode == 'bed' else rules.bcftools_prune.output
    output:
        temp(multiext('covariates/{tissue}.{coverage}.{imputed}.{mode}','.pca','.pca_stats'))
    params:
        prefix = lambda wildcards, output: PurePath(output[0]).with_suffix('')
    resources:
        mem_mb = lambda wildcards: 1024 if wildcards.mode == 'bed' else 15000,
        walltime = '30m'
    shell:
        '''
        QTLtools pca --{wildcards.mode} {input[0]} --out {params.prefix} --center --scale
        '''

rule make_covariates:
    input:
        pca_expression = lambda wildcards: expand(rules.qtltools_pca.output,mode='bed',tissue=wildcards.expression,coverage=wildcards.exp_coverage,allow_missing=True),
        pca_variants = expand(rules.qtltools_pca.output,mode='vcf',allow_missing=True),
        fixed = lambda wildcards: config['fixed_covariates'][wildcards.expression]
    output:
        'covariates/{tissue}.{coverage}.{expression}.{exp_coverage}.{imputed}.txt.gz'
    localrule: True
    shell:
        '''
        EXP=10#$(grep -oP "#Thresh995 criterion: PC #\K\d+" {input.pca_expression[1]})
        VAR=3#$(grep -oP "#Thresh995 criterion: PC #\K\d+" {input.pca_variants[1]})
        {{ cat {input.fixed} ; awk -v T=$EXP 'NR>1 && NR<=(T+1)' {input.pca_expression[0]} ; awk -v T=$VAR 'NR>1 && NR<=(T+1)' {input.pca_variants[0]} ; }} |\
        sed 's/ /\\t/g' |\
        pigz -p 2 -c > {output[0]}
        '''
