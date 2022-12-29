rule all:
    input:
        expand('replication/{tissue}.{mode}.csv',tissue=config['qtl'],mode=('best',))

localrules: prepare_qtl
rule prepare_qtl:
    input:
        qtl = lambda wildcards: config['qtl'][wildcards.tissue]
        #'/cluster/work/pausch/xena/eQTL/cis_all/testis/maf01/conditional/FDR_05/all.txt'
    output:
        'replication/{tissue}.{mode}.qtl'
    params:
        expression = lambda wildcards: '$A&&$B' if wildcards.mode == 'best' else '$B'
    shell:
        '''
        mawk '{{ if (NR==1) {{ for (i = 1; i<=NF;i++) if ($i=="bwd_best_hit") {{ A=i }}  else {{ if ($i=="bwd_sig") {{ B=i }} }} }} else {{ if ({params.expression}) {{ print $1" "$8 }} }} }}' {input} > {output}
        '''

rule qtltools_replicate:
    input:
        phenotypes = lambda wildcards: config['phenotypes'][wildcards.tissue],
        covariates = lambda wildcards: config['covariates'][wildcards.tissue],
        vcf = lambda wildcards: config['variants'][wildcards.tissue],
        qtl = rules.prepare_qtl.output
    output:
        'replication/{tissue}.{mode}.replicated'
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
        qtl = lambda wildcards: config['qtl'][wildcards.tissue]
    output:
        'replication/{tissue}.{mode}.csv'
    shell:
        '''
        while read a b c d e f g h i j k
        do
          mawk -v C=$a -v D=$f -v p=$j -v s=$k '{{ if (NR==1) {{ for (i = 1; i<=NF;i++) if ($i=="bwd_pval") {{ A=i }}  else {{ if ($i=="bwd_slope") {{ B=i }} }} }} else {{ if ($0~C&&$0~D) {{ print C,D,$A,$B,p,s }} }} }}' {input.qtl} >> {output}
        done < {input.replication}
        '''

