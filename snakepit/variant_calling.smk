from pathlib import PurePath

rule samtools_subsample:
    input:
        bam = 'subsampled_bams/{tissue}/{sample}.full.bam'
    output:
        multiext('subsampled_bams/{tissue}/{sample}.{coverage}.bam','','.csi')
    params:
        sample_rate = lambda wildcards: str(config['coverages'][wildcards.coverage])[2:]
    threads: 2
    resources:
        mem_mb = 2500
    shell:
        '''
        samtools view -@ {threads} -s $(od -N 4 -t uL -An /dev/urandom | tr -d " ").{params.sample_rate} -o {output[0]} --write-index {input}
        '''

rule DeepVariant_call:
    input:
        config = 'config/analysis.yaml',
        bams = expand(rules.samtools_subsample.output[0],sample=config['samples'],allow_missing=True)
    output:
        '{tissue}_{coverage}/all.Unrevised.vcf.gz'
    params:
        model = lambda wildcards: 'WGS' if wildcards.tissue == 'WGS' else '/cluster/work/pausch/alex/REF_DATA/RNA_DV_models/model.ckpt',
        _index = lambda wildcards: '.bai' if wildcards.coverage == 'full' else '.csi'
    threads: 1
    resources:
        mem_mb = 2500,
        walltime = '24h'
    shell:
        '''
        snakemake -s /cluster/work/pausch/alex/BSW_analysis/snakepit/deepvariant.smk --configfile {input.config} \
                --config Run_name="{wildcards.tissue}_{wildcards.coverage}" bam_path="subsampled_bams/{wildcards.tissue}/" bam_index="{params._index}" bam_name="{{sample}}.{wildcards.coverage}.bam" model="{params.model}" \
                --profile "slurm/fullNT" --nolock
        '''

rule bcftools_view:
    input:
        rules.DeepVariant_call.output
    output:
        '{tissue}_{coverage}/autosomes.Unrevised.vcf.gz'
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
        rules.bcftools_view.output
    output:
        multiext('{tissue}_{coverage}/autosomes.imputed.vcf.gz','','.tbi')
    params:
        prefix = lambda wildcards, output: PurePath(output[0]).with_suffix('').with_suffix(''),
        name = lambda wildcards, output: PurePath(output[0]).name
    threads: 8
    resources:
        mem_mb = 6000,
        walltime = '4h'
    shell:
        '''
        java -jar -Xss25m -Xmx65G /cluster/work/pausch/alex/software/beagle.22Jul22.46e.jar gt={input} nthreads={threads} out={params.prefix}
        mv {output[0]} $TMPDIR/{params.name}
        bcftools reheader -f {config[reference]}.fai -o {output[0]} $TMPDIR/{params.name}
        tabix -fp vcf {output[0]}
        '''
