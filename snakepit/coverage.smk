rule all:
    input:
        'annotated_coverage.bed'
        #expand('coverage/{sample}.{tissue}.bed.gz',sample=config['samples'],tissue=config['bams'])

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

localrules: format_annotation
rule format_annotation:
    input:
        fai = '/cluster/work/pausch/alex/REF_DATA/ARS-UCD1.2_Btau5.0.1Y.fa.fai'
    output:
        'annotation.bed'
    params:
        temp_file = '$TMPDIR/annotation.gtf.gz'
    envmodules:
        'eth_proxy'
    shell:
        '''
        curl -s https://ftp.ensembl.org/pub/release-108/gtf/bos_taurus/Bos_taurus.ARS-UCD1.2.108.chr.gtf.gz |\
        gunzip |\
        awk -v OFS="\\t" '$1~/^[0-9]/&&$3=="gene" {{ printf "%s\\t%s\\t%s\\t",$1,$4,$5; for (i=9; i<=NF; ++i) {{ if ($i=="gene_id" || $i=="gene_biotype") printf "%s\\t",$(i+1) }} printf "+\\n" }}' |\
        sed 's/"//g; s/;//g' | sort -k1,1n -k2,2n |\
        tee {output} |\
        stdbuf -o1G bedtools complement -L -i - -g {input.fai} | awk ' {{ printf "%s\\tintergenic\\tintergenic\\t+\\n",$0 }} ' >> {output}
        sort -k1,1n -k2,2n -o {output} {output}
        '''

rule bedtools_annotate:
    input:
        beds = expand('coverage/{sample}.{tissue}.bed.gz',sample=config['samples'],tissue=config['bams']),
        annotation = 'annotation.bed'
    output:
        'annotated_coverage.bed'
    params:
        names = [f'{sample}_{tissue}' for sample,tissue in zip(config['samples'],config['bams'])]
    threads: 1
    resources:
        mem_mb = 15000
    shell:
        '''
        bedtools annotate -both -names {params.names} -i {input.annotation} -files {input.beds}
        '''
