from pathlib import PurePath

rule wallaq:
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
        #'{tissue}/autosomes.Unrevised.imputed.vcf.gz'
        multiext('{tissue}/autosomes.Unrevised.imputed.vcf.gz','','.tbi')
    params:
        prefix = lambda wildcards, output: PurePath(output[0]).with_suffix('').with_suffix(''),
        name = lambda wildcards, output: PurePath(output[0]).name
    threads: 24
    resources:
        mem_mb = 3000,
        walltime = '4:00'
    shell:
        '''
        java -jar -Xss25m -Xmx65G /cluster/work/pausch/alex/software/beagle.22Jul22.46e.jar gt={input} nthreads={threads} out={params.prefix}
        mv {output[0]} $TMPDIR/{params.name}
        bcftools reheader -f {config[reference]}.fai -o {output[0]} $TMPDIR/{params.name}
        tabix -fp vcf {output[0]}
        '''
