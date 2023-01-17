rule all:
    input:
        'coverage/test.test.bed.gz'

rule perbase_depth:
    input:
        '/nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/BTA/eQTL/epi_h/RNA_alignments/normal/dedup_alignment/BSWCHEM120156182022EH/BSWCHEM120156182022EH.bam'
    output:
        'coverage/{sample}.{tissue}.bed.gz'
    threads: 4
    resources:
        mem_mb = 2500
    shell:
        '''
        perbase only-depth {input} --bed-format --bgzip -L 6 -T {threads} -t {threads} -o {output}
        '''
