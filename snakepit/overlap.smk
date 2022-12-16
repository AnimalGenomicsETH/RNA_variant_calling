

rule bcftools_isec:
    input:
        lambda wildcards: expand('{vcf}/all.{filter}.vcf.gz',filter=wildcards.filter,vcf=config['vcf'])
    output:
        'overlaps/{filter}.{mode}.txt'
    threads: 4
    resources:
        mem_mb = 1500
    shell:
        '''
        bcftools isec -n +1 --threads {threads} -c {wildcards.mode} -o {output} {input}
        '''

rule count_isec:
    input:
        'overlaps/{filter}.{mode}.txt'
    output:
        'overlaps/{filter}.{mode}.isec'
    shell:
        '''
        mawk 'length($3)==1&&length($4)==1 {{SNP[$5]+=1;next}} {{INDEL[$5]+=1}} END {{for (key in SNP) {{ print "SNP",key,SNP[key]}} for (key in INDEL) {{ print "INDEL",key,INDEL[key] }} }}' {input} > {output}
        '''

rule make_bed:
    input:
        'overlaps/{filter}.{mode}.txt'
    output:
        expand('overlaps/{filter}.{mode}.{sample}.bed',sample=config['vcf'],allow_missing=True)
    shell:
        '''
        awk '{{print $1,$2,$2+1}}' {input} > {output}
        '''

