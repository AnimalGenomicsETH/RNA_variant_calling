from pathlib import PurePath

wildcard_constraints:
    _pass = r'permutations|conditionals|nominals',
    chunk = r'\d+',
    chrom = r'\d+',
    MAF = r'\d+',
    vcf = r'(eQTL|gwas)/\S+',
    tissue = r'WGS|Epididymis_head|Vas_deferens|Testis'

def format_MAF(maf):
    return str(maf)[2:].zfill(2)

rule allerst:
    input:
        expand('eQTL/WGS_{tissue}/{_pass}.{MAF}.txt',tissue=config['vcf'],_pass=('conditionals',),MAF=format_MAF(config['MAF'])),
        expand('eQTL/{tissue}_{tissue}/{_pass}.{MAF}.txt',tissue=config['vcf'],_pass=('conditionals',),MAF=format_MAF(config['MAF'])),
        expand('replication/{tissue}.{mode}.csv',tissue=config['vcf'],mode=('best',)),
        'ase_metrics.csv'

rule qtltools_ase:
    input:
        vcf = 'WGS/autosomes.Unrevised.imputed.vcf.gz',
        bam = 'bam_symlinks/{sample}.{tissue}.bam',
        reference = config['reference'],
        annotation = '/cluster/work/pausch/alex/RNA_call_test/Bos_taurus.ARS-UCD1.2.108.chr.gtf'
    output:
        'ase/{sample}.{tissue}.ase',
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
        awk 'NR>1 {{ print {wildcards.sample},{wildcards.tissue},$21 }}' {output[0]} > {output}
        '''

localrules: calc_ase
rule calc_ase:
    input:
        'ase/{sample}.{tissue}.ase'
    output:
        'ase/{sample}.{tissue}.metric'
    shell:
        '''
        gawk '$21<0.001 {{ ++c; if ($23!="NA") {{ split($23,G,":"); Y[G[1]]; X[G[1]+G[3]] }} }} END {{ print "{wildcards.sample}","{wildcards.tissue}",c,length(Y),length(X) }}' {input} > {output}
        '''

localrules: gather_ase
rule gather_ase:
    input:
        expand('ase/{sample}.{tissue}.metric',sample=config['samples'],tissue=('Testis','Epididymis_head','Vas_deferens'))
    output:
        'ase_metrics.csv'
    shell:
        '''
        cat <(echo -e "sample tissue count genic transcript") {input} > {output}
        '''

localrules: concat_genes
rule concat_genes:
    input:
        lambda wildcards: expand('/cluster/work/pausch/xena/eQTL/gene_counts/{tissue_code}/QTLtools/{chromosome}_gene_counts.gz',chromosome=range(1,30),tissue_code={'Testis':'testis','Epididymis_head':'epi_h','Vas_deferens':'vas_d'}[wildcards.expression])
    output:
        'aligned_genes/{expression}.bed.gz',
        'aligned_genes/{expression}.bed.gz.tbi' 
    shell:
        '''
        zcat {input} | sort -u -k1,1n -k2,2n | bgzip -@ 2 -c > {output[0]}
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
        rules.normalise_vcf.output
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
        return f'--nominal {config.get("nominal",0.05)}'

rule qtltools_parallel:
    input:
        vcf = rules.normalise_vcf.output,
        exclude = rules.exclude_MAF.output,
        bed = rules.concat_genes.output,
        cov = lambda wildcards: config['covariates'][wildcards.expression],
        mapping = lambda wildcards: 'eQTL/{tissue}_{expression}/permutations_all.{MAF}.thresholds.txt' if wildcards._pass == 'conditionals' else []
    output:
        merged = temp('eQTL/{tissue}_{expression}/{_pass}.{chunk}.{MAF}.txt')
    params:
        _pass = lambda wildcards,input: get_pass(wildcards._pass,input),#f'--permute {config["permutations"]}' if wildcards._pass == 'permutations' else f'--mapping {input.mapping}',
        debug = '--silent' if 'debug' in config else '',
        grp = lambda wildcards: '--grp-best' if wildcards.tissue == 'sQTL' else ''
    threads: 1
    resources:
        mem_mb = 12500,
        walltime = lambda wildcards: '24h' if wildcards._pass == 'permutationsX' else '4h'
    shell:
        '''
        QTLtools cis --vcf {input.vcf} --bed {input.bed} --cov {input.cov} {params._pass} {params.grp} --window {config[window]} --normal --chunk {wildcards.chunk} {config[chunks]} --out {output} {params.debug}
        '''

localrules: qtltools_gather, qtltools_postprocess

rule qtltools_gather:
    input:
        expand(rules.qtltools_parallel.output,chunk=range(0,config['chunks']+1),allow_missing=True)
        #'eQTL/{tissue}_{expression}/{_pass}.{chunk}.{MAF}.txt',chunk=range(1,config['chunks']+1),allow_missing=True)
    output:
        'eQTL/{tissue}_{expression}/{_pass}.{MAF}.txt'
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
        expand(rules.qtltools_gather.output,_pass='permutations',allow_missing=True)
    output:
        'eQTL/{tissue}_{expression}/permutations_all.{MAF}.thresholds.txt'
    params:
        out = lambda wildcards, output: PurePath(output[0]).with_suffix('').with_suffix('')
    envmodules:
        'gcc/8.2.0',
        'r/4.2.2'
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

### REPLICATION RESULTS

localrules: prepare_qtl
rule prepare_qtl:
    input:
        wgs = expand(rules.qtltools_gather.output,tissue='WGS',_pass='conditionals',MAF=format_MAF(config['MAF']),allow_missing=True),
        tissue = lambda wildcards: expand(rules.qtltools_gather.output,tissue=wildcards.expression,_pass='conditionals',MAF=format_MAF(config['MAF']),allow_missing=True)
    output:
        'replication/{expression}.qtl'
    #params:
    #    expression = lambda wildcards: '$A&&$B' if wildcards.mode == 'best' else '$B'
    shell:
        '''
        mawk '$21&&$22 {{print $1,$8 }}' {input} > {output}
        '''

rule qtltools_replicate:
    input:
        phenotypes = rules.concat_genes.output[0],
        covariates = lambda wildcards: config['covariates'][wildcards.expression],
        vcf = lambda wildcards: expand(rules.normalise_vcf.output,tissue=wildcards.mode),
        qtl = rules.prepare_qtl.output
    output:
        'replication/{expression}.{mode}.replicated'
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
        'replication/{expression}.{mode}.csv'
    shell:
        '''
        while read a b c d e f g h i j k
        do
          mawk -v C=$a -v D=$f -v p=$j -v s=$k '{{ if (NR==1) {{ for (i = 1; i<=NF;i++) if ($i=="nom_pval") {{ A=i }}  else {{ if ($i=="slope") {{ B=i }} }} }} else {{ if ($0~C&&$0~D) {{ pri  nt C,D,$A,$B,p,s }} }} }}' {input.qtl} >> {output}
        done < {input.replication}
        '''

rule qtltools_rtc:
    input:
        WGS = expand(rules.qtltools_gather.output,tissue='WGS',_pass='conditionals',MAF=format_MAF(config['MAF']),allow_missing=True),
        tissue = expand(rules.qtltools_gather.output,_pass='conditionals',MAF=format_MAF(config['MAF']),allow_missing=True),
        vcf = '',
        bed = '',
        cov = '',
        hotspot = '',
        IDs = '',
        qtl = ''
    output:
        'replication/{tissue}.rtc.out'
    shell:
        '''
        QTLtools rtc --vcf {input.vcf} --bed {input.bed} --cov {input.cov} --hotspot {input.hotspot} --mergeQTL-cis --window {config[window]} --normal --conditional --mergeQTL-cis {input.IDs} {input.qtl} --chunk {wildcards.chunk} {config[chunks]} --out {output}
        '''

