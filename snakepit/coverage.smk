rule all:
    input:
        expand('coverage/{sample}.{tissue}.bed.gz',sample=config['samples'],tissue=config['bams'])

def get_ID(sample,tissue):
    match tissue:
        case 'Epididymis_head':
            return sample + 'EH'
        case 'Vas_deferens':
            return sample + 'V'
        case 'Testis':
            return sample
        case 'WGS':
            return sample

rule perbase_depth:
    input:
        'bam_symlinks/{sample}.{tissue}.bam'
        #'/cluster/work/pausch/vcf_UCD/2022_03/bam_repository/BSWCHEF120023224572.bam'
        #'/nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/BTA/eQTL/epi_h/RNA_alignments/normal/dedup_alignment/BSWCHEM120156182022EH/BSWCHEM120156182022EH.bam'
    output:
        'coverage/{sample}.{tissue}.bed.gz'
    params:
        window = 100000
    threads: 2
    resources:
        mem_mb = 5000
    shell:
        '''
        perbase only-depth {input} --bed-format -z -F 3848 -t {threads} --bgzip -L 6 -T {threads} -o {output}
        '''
        #'''
        #awk -v W={params.window} '$1~/^[0-9]/' {{ A[$1][int($2/W)]+=($3-$2)*$4 }} END {{ for (chr in A) {{ for (window in A[chr]) {{ print chr,window*W,window*W+W,A[chr][window] }} }} }}' > {output}
        #'''
