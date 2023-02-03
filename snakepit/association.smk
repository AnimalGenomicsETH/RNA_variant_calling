from pathlib import PurePath

wildcard_constraints:
    _pass = r'permutations|conditionals|nominals',
    chunk = r'\d+',
    chrom = r'\d+',
    MAF = r'\d+',
    vcf = r'(eQTL|gwas)/\S+'

rule all:
    input:
        expand('eQTL/{tissue}/{_pass}.{MAF}.txt',tissue=config['vcf'],_pass=('conditionals',),MAF=(0.01))

localrules: concat_genes
rule concat_genes:
    input:
        lambda wildcards: expand('/cluster/work/pausch/xena/eQTL/gene_counts/{tissue_code}/QTLtools/{chromosome}_gene_counts.gz',chromosome=range(1,30),tissue_code={'Testis':'testis','Epididymis_head':'epi_h','Vas_deferens':'vas_d'}[wildcards.tissue])
    output:
        'aligned_genes/{tissue}.bed.gz',
        'aligned_genes/{tissue}.bed.gz.tbi' 
    shell:
        '''
        zcat {input} | sort -k1,1n -k2,2n | pigz -p 2 -c > {output[0]}
        tabix -p bed {output[0]}
        '''

rule normalise_vcf:
    input:
        '{tissue}/all.Unrevised.vcf.gz'
    output:
        'eQTL/{tissue}/variants.normed.vcf.gz'
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
        'eQTL/{tissue}/variants.normed.vcf.gz'
    output:
        'eQTL/{tissue}/exclude_sites.{MAF}.txt'
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
        return '--nominal 1.0'

rule qtltools_parallel:
    input:
        vcf = rules.normalise_vcf.output,
        exclude = rules.exclude_MAF.output,
        bed = rules.concat_genes.output,
        cov = lambda wildcards: config['covariates'][wildcards.tissue],
        mapping = lambda wildcards: 'eQTL/{tissue}/permutations_all.{MAF}.thresholds.txt' if wildcards._pass == 'conditionals' else []
    output:
        merged = temp('eQTL/{tissue}/{_pass}.{chunk}.{MAF}.txt')
    params:
        _pass = lambda wildcards,input: get_pass(wildcards._pass,input),#f'--permute {config["permutations"]}' if wildcards._pass == 'permutations' else f'--mapping {input.mapping}',
        debug = '--silent' if 'debug' in config else '',
        grp = lambda wildcards: '--grp-best' if wildcards.tissue == 'sQTL' else ''
    threads: 1
    resources:
        mem_mb = 12500,
        walltime = lambda wildcards: '24:00' if wildcards._pass == 'permutations' else '4:00'
    shell:
        '''
        QTLtools cis --vcf {input.vcf} --bed {input.bed} --cov {input.cov} {params._pass} {params.grp} --window {config[window]} --normal --chunk {wildcards.chunk} {config[chunks]} --out {output} {params.debug}
        '''

localrules: qtltools_gather, qtltools_postprocess

rule qtltools_gather:
    input:
        expand('eQTL/{tissue}/{_pass}.{chunk}.{MAF}.txt',chunk=range(1,config['chunks']+1),allow_missing=True)
    output:
        'eQTL/{tissue}/{_pass}.{MAF}.txt'
    resources:
        mem_mb = 3000,
        walltime = '20'
    params:
        sort_key = lambda wildcards: '-k9,9n -k10,10n' if wildcards.tissue != 'eQTL' else '-k11,11n -k12,12n'
    shell:
        '''
        sort {params.sort_key} {input} > {output}
        '''

rule qtltools_FDR:
    input:
        'eQTL/{tissue}/permutations.{MAF}.txt'
    output:
        'eQTL/{tissue}/permutations_all.{MAF}.thresholds.txt'
    params:
        out = lambda wildcards, output: PurePath(output[0]).with_suffix('').with_suffix('')
    shell:
        '''
        Rscript /cluster/work/pausch/alex/software/qtltools/scripts/qtltools_runFDR_cis.R {input} 0.05 {params.out}
        '''

rule qtltools_postprocess:
    input:
        'eQTL/{tissue}/conditionals.{MAF}.txt'
    output:
        'eQTL/{tissue}/significant_hits.{minS}.{MAF}.fastman'
    params:
        print_key = lambda wildcards: '$9"\t"($11-$10)"\t"$10"\t"$8"\tN\teQTL\t2\t"$20"\t"$19"\t"$18' if wildcards.tissue != 'eQTL' else '$11"\t"($13-$12)"\t"$12"\t"$10"\tN\teQTL\t2\t"$22"\t"$21"\t"$20'
    shell:
        '''
        echo "CHR\tsize\tBP\tSNP\tA1\tTEST\tNMISS\tBETA\tSTAT\tP" > {output}
        awk -v L={wildcards.minS} '($11-$10)>=L {{print {params.print_key}}}' {input} >>  {output}
        '''

