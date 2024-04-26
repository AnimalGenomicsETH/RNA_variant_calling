rule perbase_depth:
    input:
        'subsampled_bams/{tissue}/{sample}.{coverage}.bam'
    output:
        'coverage/{tissue}/{sample}.{coverage}.bed.gz'
    params:
        window = 100000
    threads: 2
    resources:
        mem_mb = lambda wildcards: 15000 if wildcards.tissue == 'WGS' else 10000
    shell:
        '''
        perbase only-depth {input} --bed-format -z -F 3848 -t {threads} --bgzip -L 6 -T {threads} -o {output}
        '''

rule missing_genes:
    input:
        bed = rules.perbase_depth.output,
        annotation = 'annotation.bed'
    output:
        'coverage/{tissue}/{sample}.{coverage}.missing_genes.bed'
    shell:
        '''
        bedtools intersect -sorted -a {input.bed} -b {input.annotation} | awk '$4>100' | bedtools merge -i - -d 100 > {output}
        bedtools slop -b 10000 -i protein_coding.bed -g ../REF_DATA/ARS-UCD1.2_Btau5.0.1Y.fa.fai | bedtools subtract -sorted -b - -a coverage/BSWCHEM120154151075.Vas_deferens.bed.gz | awk '$4>200' | bedtools merge -i - -d 100 | awk '($3-$2)>2000 {print $1":"$2"-"$3}' > candidates.txt
        '''

rule estimate_coverage:
    input:
       rules.perbase_depth.output
    output:
        'coverage/{tissue}/{sample}.{coverage}.genome_coverage.csv'
    resources:
        walltime = '20m'
    shell:
        '''
        zcat {input} | mawk -v S={wildcards.sample} -v T={wildcards.tissue} -v OFS="\\t" '$1~/^[0-9]/ {{ c+=($3-$2); d+=($3-$2)*($4>1) }} END {{ print S,T,c,d }}' > {output}
        '''

rule collate_coverage:
    input:
        expand(rules.estimate_coverage.output,sample=config['samples'],tissue=config['bams'],allow_missing=True)
    output:
        'coverage/genome.{coverage}.csv'
    localrule: True
    shell:
        '''
        cat <(echo -e "sample\\ttissue\\tcoverage\\tcoverage >1") {input} > {output}
        '''

rule bedtools_annotate:
    input:
        beds = rules.perbase_depth.output,
        annotation = rules.format_annotation.output['bed']
    output:
        'coverage/{sample}.{tissue}.{coverage}.annotated.bed'
    threads: 1
    resources:
        mem_mb = lambda wildcards: 95000 if wildcards.tissue == 'WGS' else 25000
    shell:
        '''
        bedtools annotate -both -i {input.annotation} -files {input.beds} | awk 'NR>1 {{ print $0"\\t{wildcards.sample}\\t{wildcards.tissue}" }}' | sort -k1,1n -k2,2n >> {output}
        '''

rule collate_annotation:
    input:
        beds = expand(rules.bedtools_annotate.output,sample=config['samples'],tissue=config['bams'],allow_missing=True)
    output:
        'coverage/annotated.{coverage}.bed.gz'
    localrule: True
    shell:
        '''
        cat <(echo -e "chr\\tstart\\tstop\\tgene\\tfeature\\tstrand\\tdepth\\tcoverage\\tsample\\ttissue") {input} | pigz -c > {output}
        '''
