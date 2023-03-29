from pathlib import PurePath

def format_MAF(maf):
    return str(maf)[2:].zfill(2)

rule qtltools_ase:
    input:
        vcf = '{tissue}_{coverage}/autosomes.imputed.vcf.gz',
        bam = 'subsampled_bams/{tissue}/{sample}.{coverage}.bam',
        reference = config['reference'],
        annotation = '/cluster/work/pausch/alex/RNA_call_test/Bos_taurus.ARS-UCD1.2.108.chr.gtf'
    output:
        'ase/{sample}.{tissue}.{coverage}.ase',
        #'ase/{sample}.{tissue}.metric',
        #temp('ase/{sample}.{tissue}.ref_bias')
    params:
        out = lambda wildcards, output: PurePath(output[0]).with_suffix('')
    threads: 1
    resources:
        mem_mb = 10000,
        walltime = '24h'
    shell:
        '''
        QTLtools ase --bam {input.bam} --vcf {input.vcf} --ind {wildcards.sample} --mapq 10 -f {input.reference} --gtf {input.annotation} --suppress-warnings --pvalue 0.001 --out {params.out}
        '''

rule calc_ase:
    input:
        rules.qtltools_ase.output[0]
    output:
        'ase/{sample}.{tissue}.{coverage}.metric'
    localrule: True
    shell:
        '''
        gawk '$21<0.001 {{ ++c; if ($23!="NA") {{ split($23,G,":"); Y[G[1]]; X[G[1]+G[3]] }} }} END {{ print "{wildcards.sample}","{wildcards.tissue}",c,length(Y),length(X) }}' {input} > {output}
        '''

rule gather_ase:
    input:
        expand(rules.calc_ase.output,sample=config['samples'],tissue=('Testis','Epididymis_head','Vas_deferens'),allow_missing=True)
    output:
        'ase/metrics.{coverage}.csv'
    localrule: True
    shell:
        '''
        cat <(echo -e "sample tissue count genic transcript") {input} > {output}
        '''

rule concat_genes:
    input:
        lambda wildcards: expand('/cluster/work/pausch/xena/eQTL/gene_counts/{tissue_code}/QTLtools/{chromosome}_gene_counts.gz',chromosome=range(1,30),tissue_code={'Testis':'testis','Epididymis_head':'epi_h','Vas_deferens':'vas_d'}[wildcards.expression])
    output:
        'aligned_genes/{expression}.bed.gz',
        'aligned_genes/{expression}.bed.gz.tbi' 
    localrule: True
    shell:
        '''
        zcat {input} | sort -u -k1,1n -k2,2n | bgzip -@ 2 -c > {output[0]}
        tabix -p bed {output[0]}
        '''

rule normalise_vcf:
    input:
        '{tissue}_{coverage}/autosomes.imputed.vcf.gz'
    output:
        'eQTL/{tissue}_{coverage}/variants.normed.vcf.gz'
    threads: 4
    resources:
        mem_mb = 2500,
        disk_scratch = 50
    shell:
        '''
        bcftools norm --threads {threads} -f {config[reference]} -m -any {input} -Ou | \
        bcftools sort -T $TMPDIR -Ou - | \
        bcftools annotate --threads {threads} --set-id '%CHROM\_%POS\_%TYPE\_%REF\_%ALT' -o {output} -
        tabix -fp vcf {output}
        '''

rule exclude_MAF:
    input:
        rules.normalise_vcf.output
    output:
        'eQTL/{tissue}_{coverage}/exclude_sites.{MAF}.txt'
    shell:
        '''
        bcftools query -f '%ID\n' -i 'MAF<0.{wildcards.MAF}' {input} > {output}
        '''

def get_pass(_pass,input):
    if _pass == 'permutations':
        return f'--permute {config["permutations"]}'
    elif _pass == 'conditionals':
        return f'--mapping {input.mapping}'
    elif _pass == 'nominals':
        return f'--nominal {config.get("nominal",0.05)}'

rule qtltools_parallel:
    input:
        vcf = rules.normalise_vcf.output,
        exclude = rules.exclude_MAF.output,
        bed = rules.concat_genes.output,
        cov = lambda wildcards: config['covariates'][wildcards.expression],
        mapping = lambda wildcards: rules.qtltools_FDR.output if wildcards._pass == 'conditionals' else []
    output:
        merged = temp('eQTL/{tissue}_{coverage}_{expression}/{_pass}.{chunk}.{MAF}.txt.gz')
    params:
        _pass = lambda wildcards,input: get_pass(wildcards._pass,input)
    threads: 1
    resources:
        mem_mb = 5500,
        walltime = '4h'
    shell:
        '''
        QTLtools cis --vcf {input.vcf} --bed {input.bed} --cov {input.cov} --std-err {params._pass} --window {config[window]} --normal --chunk {wildcards.chunk} {config[chunks]} --silent --log /dev/stderr --out /dev/stdout | pigz -p 2 -c > {output} 
        '''

rule qtltools_gather:
    input:
        expand(rules.qtltools_parallel.output,chunk=range(0,config['chunks']+1),allow_missing=True)
    output:
        'eQTL/{tissue}_{coverage}_{expression}/{_pass}.{MAF}.txt.gz'
    localrule: True
    shell:
        '''
        LC_ALL=C; pigz -p 2 -dc {input} | sort --parallel=2 -k9,9n -k10,10n | pigz -p 2 -c > {output}
        sort {params.sort_key} {input} > {output}
        '''

rule qtltools_FDR:
    input:
        expand(rules.qtltools_gather.output,_pass='permutations',allow_missing=True)
    output:
        'eQTL/{tissue}_{coverage}_{expression}/permutations_all.{MAF}.thresholds.txt'
    params:
        out = lambda wildcards, output: PurePath(output[0]).with_suffix('').with_suffix('')
    envmodules:
        'gcc/8.2.0',
        'r/4.2.2'
    shell:
        '''
        Rscript /cluster/work/pausch/alex/software/qtltools/scripts/qtltools_runFDR_cis.R {input} 0.05 {params.out}
        '''

### REPLICATION RESULTS

rule prepare_qtl:
    input:
        wgs = expand(rules.qtltools_gather.output,tissue='WGS',_pass='conditionals',MAF=format_MAF(config['MAF']),allow_missing=True),
        tissue = lambda wildcards: expand(rules.qtltools_gather.output,tissue=wildcards.expression,_pass='conditionals',MAF=format_MAF(config['MAF']),allow_missing=True)
    output:
        'replication/{expression}.{coverage}.qtl'
    #params:
    #    expression = lambda wildcards: '$A&&$B' if wildcards.mode == 'best' else '$B'
    localrule: True
    shell:
        '''
        mawk '$21&&$22 {{print $1,$8 }}' {input} > {output}
        '''

rule qtltools_replicate:
    input:
        phenotypes = rules.concat_genes.output[0],
        covariates = lambda wildcards: config['covariates'][wildcards.expression],
        vcf = lambda wildcards: expand(rules.normalise_vcf.output,tissue=wildcards.mode,coverage=wildcards.coverage),
        qtl = rules.prepare_qtl.output
    output:
        'replication/{expression}.{mode}.{coverage}.replicated'
    threads: 1
    resources:
        mem_mb = 2500
    shell:
        '''
        QTLtools rep --bed {input.phenotypes} --cov {input.covariates} --vcf {input.vcf} --qtl {input.qtl} --normal --out {output}
        '''

rule collate_replicate:
    input:
        replication = rules.qtltools_replicate.output,
        qtl = expand(rules.qtltools_gather.output,tissue='WGS',_pass='nominals',MAF=format_MAF(config['MAF']),allow_missing=True)
    output:
        'replication/{expression}.{mode}.{coverage}.csv'
    shell:
        '''
        while read a b c d e f g h i j k
        do
          mawk -v C=$a -v D=$f -v p=$j -v s=$k '{{ if (NR==1) {{ for (i = 1; i<=NF;i++) if ($i=="nom_pval") {{ A=i }}  else {{ if ($i=="slope") {{ B=i }} }} }} else {{ if ($0~C&&$0~D) {{ pri  nt C,D,$A,$B,p,s }} }} }}' {input.qtl} >> {output}
        done < {input.replication}
        '''

