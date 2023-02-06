from pathlib import PurePath

rule all:
    input:
        'happy_metrics.csv'

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

workflow.singularity_args = f'-B $TMPDIR -B {PurePath(config["reference"]).parent}'

rule bcftools_split:
    input:
        '{tissue}/all.Unrevised.vcf.gz'
    output:
        expand('split_vcfs/{tissue}/{sample}.vcf.gz',sample=config['samples'],allow_missing=True)
    params:
        _dir = lambda wildcards, output: PurePath(output[0]).parent
    resources:
        mem_mb = 2500
    shell:
        '''
        bcftools +split -i 'GT[*]="alt"' -Oz -o {params._dir} {input}
        '''

rule happy:
    input:
        vcfs = lambda wildcards: expand('split_vcfs/{tissues}/{sample}.vcf.gz',tissues=('WGS',wildcards.tissue),allow_missing=True),
        reference = config['reference']
    output:
        csv = 'happy/{sample}.{tissue}.summary.csv',
        others = temp(multiext('happy/{sample}.{tissue}','.bcf','.bcf.csi','.extended.csv','.roc.all.csv.gz','.runinfo.json'))
    params:
        _dir = lambda wildcards, output: PurePath(output.csv).with_suffix('').with_suffix('')
    container: '/cluster/work/pausch/alex/software/images/hap.py_latest.sif'
    threads: 4
    resources:
        mem_mb = 500,
        sratch = '10G'
    shell:
        '''
        /opt/hap.py/bin/hap.py -r {input.reference} --bcf --usefiltered-truth --no-roc --no-json -L --pass-only --scratch-prefix $TMPDIR -X --threads {threads} -o {params._dir} {input.vcfs}
        '''

localrules: gather_happy
rule gather_happy:
    input:
        expand('happy/{sample}.{tissue}.summary.csv',sample=config['samples'],tissue=('Testis','Epididymis_head','Vas_deferens'))
    output:
        'happy_metrics.csv'
    shell:
        '''
        echo -e "variant truth query recall precision truth_TiTv query_TiTv sample tissue" > {output}
        for i in {input}
        do
          awk -v I=$(basename $i) -F',' '$2=="PASS" {{ split(I,a,"."); print $1,$3,$6,$10,$11,$14,$15,a[1],a[2] }}' $i >> {output}
        done
        '''
