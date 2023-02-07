from pathlib import PurePath

rule all:
    input:
        expand('{tissue}/autosomes.Unrevised.imputed.vcf.gz',tissue=config['tissues'])

rule bcftools_view:
    input:
        '{tissue}/all.Unrevised.vcf.gz'
    output:
        '{tissue}/autosomes.Unrevised.vcf.gz'
    params:
        autosomes = ','.join(map(str,range(1,30)))
    threads: 4
    resources:
        mem_mb = 1500
    shell:
        '''
        bcftools view --threads {threads} --regions {params.autosomes} -o {output} {input}
        tabix -p vcf {output}
        '''

rule beagle5_imputation:
    input:
        '{tissue}/autosomes.Unrevised.vcf.gz'
    output:
        multiext('{tissue}/autosomes.Unrevised.imputed.vcf.gz','','.tbi')
    params:
        lambda wildcards, output: PurePath(output[0]).with_suffix('').with_suffix('')
    threads: 8
    resources:
        mem_mb = 6500,
        walltime = '4:00'
    shell:
        '''
        java -jar -Xss25m -Xmx40G /cluster/work/pausch/alex/software/beagle.22Jul22.46e.jar gt={input} nthreads={threads} out={params}
        mv {output} $TMPDIR/temporary.vcf
        bcftools reheader -f {config[reference]} -o {output} $TMPDIR/temporary.vcf
        tabix -fp vcf {output}
        '''
