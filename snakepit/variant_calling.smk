
rule all:
    input:
        'subsampled_bams/A.50.bam'

rule samtools_subsample:
    input:
        bam = 'A.bam'
    output:
        multiext('subsampled_bams/{sample}.{sample_rate}.bam','','.csi')
    params:
        sample_rate = lambda wildcards: float(wildcards.sample_rate)/(10**len(wildcards.sample_rate))
    threads: 2
    resources:
        mem_mb = 2500
    shell:
        '''
        samtools view -@ {threads} -s $(od -N 4 -t uL -An /dev/urandom | tr -d " ").{wildcards.sample_rate} -o {output[0]} --write-index {input}
        '''
