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
        multiext('{tissue}_{coverage}/all.Unrevised.vcf.gz','','.tbi')
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
        bcftools view --threads {threads} --regions {params.autosomes} -o {output} {input[0]}
        tabix -p vcf {output}
        '''

rule beagle4_imputation:
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
        walltime = lambda wildcards: '24h' if wildcards.tissue != 'WGS' else '48h'
    shell:
        '''
        java -jar -Xss25m -Xmx40G /cluster/work/pausch/alex/software/beagle.27Jan18.7e1.jar gl={input} nthreads={threads} out={params.prefix}
        mv {output[0]} $TMPDIR/{params.name}
        bcftools reheader -f {config[reference]}.fai -o {output[0]} $TMPDIR/{params.name}
        tabix -fp vcf {output[0]}
        '''

rule beagle5_panel_imputation:
    input:
        vcf = rules.beagle4_imputation.output,
        panel = 'panel.phased.vcf.gz'
    output:
        multiext('{tissue}_{coverage}/autosomes.panel.vcf.gz','','.tbi')
    params:
        prefix = lambda wildcards, output: PurePath(output[0]).with_suffix('').with_suffix(''),
        name = lambda wildcards, output: PurePath(output[0]).name
    threads: 24
    resources:
        mem_mb = 4000,
        waltime = '4h'
    shell:
        '''
        java -jar -Xss25m -Xmx80G /cluster/work/pausch/alex/software/beagle.22Jul22.46e.jar ref={input.panel} gt={input.vcf[0]} nthreads={threads} out={params.prefix}
        mv {output[0]} $TMPDIR/{params.name}
        bcftools reheader -f {config[reference]}.fai -o {output[0]} $TMPDIR/{params.name}
        tabix -fp vcf {output[0]}
        '''

rule plink_PCA:
    input:
        rules.beagle4_imputation.output[0]
        #'{tissue}_{coverage}/autosomes.pruned.vcf.gz'
    output:
        multiext('PCA/{tissue}.{coverage}','.eigenval','.eigenvec','.log')
    params:
        prefix = lambda wildcards, output: PurePath(output[0]).with_suffix(''),
        maf = 0.05
    threads: 1
    resources:
        mem_mb = 10000,
        walltime = '30m'
    shell:
        '''
        plink2 --pca --vcf {input} --out {params.prefix} --cow --maf {params.maf} --threads {threads} --vcf-half-call m
        '''

rule vep_download:
    output:
        directory('vep/bos_taurus')
    params:
        lambda wildcards, output: PurePath(output[0]).parent
    localrule: True
    shell:
        '''
        export TMPDIR=/cluster/scratch/alleonard
        curl -o $TMPDIR/bos_taurus.tar.gz https://ftp.ensembl.org/pub/release-110/variation/indexed_vep_cache/bos_taurus_vep_110_ARS-UCD1.2.tar.gz
        tar xzf $TMPDIR/bos_taurus.tar.gz -C {parmas}
        '''

rule vep_annotate:
    input:
        vcf = rules.beagle4_imputation.output,
        vep = rules.vep_download.output[0]
    output:
        'vep/{tissue}.{coverage}.vep'
    params:
        lambda wildcards, input: PurePath(input.vep).parent
    threads: 4
    resources:
        mem_mb = 1250
    container: '/cluster/work/pausch/alex/software/images/ensembl-vep_latest.sif'
    shell:
        '''
        export LC_ALL=C
        vep -i {input.vcf[0]} --tab --fields "Consequence,IMPACT" --species bos_taurus --offline --no_headers --no_stats --dir {params} --output_file STDOUT --fork {threads} |\
        mawk '{{ split($1,c,","); for (k in c) {{ ++CSQ[c[k]"_"$2] }} }} END {{ for (k in CSQ) {{ print k,CSQ[k] }} }}' |\
        sed 's/\(.*\)_/\\1 /' > {output}
        '''
