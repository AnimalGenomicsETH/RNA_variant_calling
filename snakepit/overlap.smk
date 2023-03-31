from pathlib import PurePath

rule bcftools_isec:
    input:
        tissues = lambda wildcards: expand('{tissue}_{coverage}/autosomes.{imputed}.vcf.gz',tissue=config['tissues'],allow_missing=True),
        WGS = lambda wildcards: expand('WGS_full/autosomes.{imputed}.vcf.gz',allow_missing=True)
    output:
        'overlaps/{imputed}.{mode}.{coverage}.txt'
    threads: 4
    resources:
        mem_mb = 1500
    shell:
        '''
        bcftools isec -n +1 --threads {threads} -c {wildcards.mode} -o {output} {input}
        '''

rule count_isec:
    input:
        isec = rules.bcftools_isec.output,
        annotation = 'annotation.bed'
    output:
        'overlaps/{imputed}.{mode}.{coverage}.isec'
    shell:
        '''
        awk -v OFS="\\t" '{{ print $1,$2,$2+1,$3,$4,$5 }}' {input.isec} |\
        bedtools intersect -wo -a - -b {input.annotation} |\
        awk ' {{ if($4~/,/||$5~/,/) {{ ++MA[$10=="intergenic"][$6] }} else {{ if (length($4)==1&&length($5)==1) {{++SNP[$10=="intergenic"][$6]}} else {{++INDEL[$10=="intergenic"][$6]}} }} }} END {{for (key in SNP[0]) {{ print "genic","SNP",key,SNP[0][key] }} for (key in SNP[1]) {{ print "intergenic","SNP",key,SNP[1][key] }} for (key in INDEL[0]) {{ print "genic","INDEL",key,INDEL[0][key] }} for (key in INDEL[1]) {{ print "intergenic","INDEL",key,INDEL[1][key] }} for (key in MA[0]) {{ print "genic","MA",key,MA[0][key]}} for (key in MA[1]) {{ print "intergenic","MA",key,MA[1][key]}} }}' > {output}
        '''

rule make_bed:
    input:
        rules.bcftools_isec.output
    output:
        expand('overlaps/{imputed}.{mode}.{coverage}.bed',sample=config['vcf'],allow_missing=True)
    shell:
        '''
        awk '{{print $1,$2,$2+1}}' {input} > {output}
        '''

workflow.singularity_args = f'-B $TMPDIR -B {PurePath(config["reference"]).parent}'

rule bcftools_split:
    input:
        '{tissue}_{coverage}/autosomes.{imputed}.vcf.gz'
    output:
        expand('split_vcfs/{tissue}_{coverage}_{imputed}/{sample}.vcf.gz',sample=config['samples'],allow_missing=True)
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
        vcf_tissue = lambda wildcards: expand('split_vcfs/{tissue}_{coverage}_{imputed}/{sample}.vcf.gz',tissues=wildcards.tissue,allow_missing=True),
        vcf_WGS = lambda wildcards: expand('split_vcfs/WGS_full_{imputed}/{sample}.vcf.gz',allow_missing=True),
        reference = config['reference']
    output:
        csv = 'F1/{sample}.{tissue}.{coverage}.{imputed}.summary.csv',
        others = temp(multiext('F1/{sample}.{tissue}.{coverage}.{imputed}','.bcf','.bcf.csi','.extended.csv','.roc.all.csv.gz','.runinfo.json'))
    params:
        _dir = lambda wildcards, output: PurePath(output.csv).with_suffix('').with_suffix('')
    container: '/cluster/work/pausch/alex/software/images/hap.py_latest.sif'
    threads: 2
    resources:
        mem_mb = 1000,
        scratch = '10G'
    shell:
        '''
        /opt/hap.py/bin/hap.py -r {input.reference} --bcf --usefiltered-truth --no-roc --no-json -L --pass-only --scratch-prefix $TMPDIR -X --threads {threads} -o {params._dir} {input.vcf_WGS} {input.vcf_tissue}
        '''

rule gather_happy:
    input:
        expand(rules.happy.output[0],sample=config['samples'],tissue=('Testis','Epididymis_head','Vas_deferens'),allow_missing=True)
    output:
        'F1/happy.{coverage}.{imputed}.csv'
    localrule: True
    shell:
        '''
        echo -e "variant truth query recall precision truth_TiTv query_TiTv sample tissue" > {output}
        for i in {input}
        do
          awk -v I=$(basename $i) -F',' '$2=="PASS" {{ split(I,a,"."); print $1,$3,$6,$10,$11,$14,$15,a[1],a[2] }}' $i >> {output}
        done
        '''
