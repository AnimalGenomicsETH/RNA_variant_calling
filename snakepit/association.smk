from pathlib import PurePath

def format_MAF(maf):
    return str(maf)[2:].zfill(2)

rule qtltools_ase:
    input:
        vcf = 'WGS_{coverage}/autosomes.imputed.vcf.gz',
        bam = 'subsampled_bams/{tissue}/{sample}.{coverage}.bam',
        reference = config['reference'],
        annotation = rules.format_annotation.output['gtf']
    output:
        multiext('ase/{sample}.{tissue}.{coverage}','.ase','.metric','.ref_bias')
    params:
        out = lambda wildcards, output: PurePath(output[0]).with_suffix('')
    threads: 1
    resources:
        mem_mb = 10000,
        walltime = '24h'
    shell:
        '''
        QTLtools ase --bam {input.bam} --vcf {input.vcf} --ind {wildcards.sample} --both-alleles-seen --mapq 10 -f {input.reference} --gtf {input.annotation} --suppress-warnings --pvalue 0.001 --out {params.out}
        '''

rule calc_ase:
    input:
        rules.qtltools_ase.output[0]
    output:
        'ase/{sample}.{tissue}.{coverage}.ase_metric'
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

rule normalise_vcf:
    input:
        '{tissue}_{coverage}/autosomes.{imputed}.vcf.gz'
    output:
        'eQTL/{tissue}_{coverage}_{imputed}/variants.normed.vcf.gz'
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
        'eQTL/{tissue}_{coverage}_{imputed}/exclude_sites.{MAF}.txt'
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
        bed = lambda wildcards: expand(rules.filter_TPM.output[0],tissue=wildcards.expression,coverage=wildcards.exp_coverage,allow_missing=True),
        cov = rules.make_covariates.output,
        mapping = lambda wildcards: rules.qtltools_FDR.output if wildcards._pass == 'conditionals' else []
    output:
        merged = temp('eQTL/{tissue}_{coverage}_{expression}_{exp_coverage}_{imputed}/{_pass}.{chunk}.{MAF}.txt.gz')
    params:
        _pass = lambda wildcards,input: get_pass(wildcards._pass,input)
    localrule: lambda wildcards: wildcards._pass == 'conditionals'
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
        'eQTL/{tissue}_{coverage}_{expression}_{exp_coverage}_{imputed}/{_pass}.{MAF}.txt.gz'
    localrule: True
    shell:
        '''
        LC_ALL=C; pigz -p 2 -dc {input} | sort --parallel=2 -k9,9n -k10,10n | pigz -p 2 -c > {output}
        '''

rule qtltools_FDR:
    input:
        expand(rules.qtltools_gather.output,_pass='permutations',allow_missing=True)
    output:
        'eQTL/{tissue}_{coverage}_{expression}_{exp_coverage}_{imputed}/permutations_all.{MAF}.thresholds.txt'
    params:
        out = lambda wildcards, output: PurePath(output[0]).with_suffix('').with_suffix('')
    envmodules:
        'gcc/8.2.0',
        'r/4.2.2'
    shell:
        '''
        Rscript qtltools_runFDR_cis.R {input} 0.05 {params.out}
        '''

rule mismatched_QTL:
    input:
        expand('eQTL/WGS_full_{tissue}_{coverage}/conditionals.{MAF}.txt.gz',tissue=config['tissues'],MAF=format_MAF(config['MAF']),coverage=config['coverages']),
        expand('eQTL/{tissue}_{coverage}_{tissue}_{coverage}/conditionals.{MAF}.txt.gz',tissue=config['tissues'],MAF=format_MAF(config['MAF']),coverage=config['coverages'])
    output:
        'eQTL/mismatched_genes.csv'
    localrule: True
    shell:
        '''

        for i in $(find eQTL -type d -name '*_*_*_*'); do cp ${{i}}/conditionals.01.txt.gz ${{i}}.conditionals.csv.gz; done

        echo -e "tissue\\teGene\\tWGS threshold\\tWGS pval\\tRNA threshold\\tRNA pval" > {output}

coverage="full"
for tissue in Testis Vas_deferens Epididymis_head
do
  for variants in WGS $tissue
  do
    zcat eQTL/${{variants}}_${{coverage}}_${{tissue}}_${{coverage}}/conditionals.01.txt.gz | awk 'NR>1 {{print $1}}' | sort -u > eQTL/${{variants}}_${{tissue}}.uniq
  done
    comm -23 eQTL/WGS_${{tissue}}.uniq eQTL/${{tissue}}_${{tissue}}.uniq > eQTL/${{tissue}}.WGS.only
    comm -13 eQTL/WGS_${{tissue}}.uniq eQTL/${{tissue}}_${{tissue}}.uniq > eQTL/${{tissue}}.${{tissue}}.only

  while read G
  do
    echo -ne "${{tissue}}\\t$G\\t" >> {output}
    for variants in WGS ${{tissue}}
    do
      grep -hF $G eQTL/${{variants}}_${{coverage}}_${{tissue}}_${{coverage}}/permutations_all.01.thresholds.txt | awk '{{printf $2"\\t"}}' >> {output}
      zgrep -hF $G eQTL/${{variants}}_${{coverage}}_${{tissue}}_${{coverage}}/permutations.01.txt.gz | awk '{{printf $16"\\t"}}' >> {output}
    done
  echo >> {output}
  done < <(cat  eQTL/${{tissue}}.${{tissue}}.only  eQTL/${{tissue}}.WGS.only)

done
        '''
