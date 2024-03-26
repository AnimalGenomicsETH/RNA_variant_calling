from pathlib import PurePath

rule bcftools_filter:
    input:
        '{tissue}_{coverage}/autosomes.imputed.vcf.gz'
    output:
        '{tissue}_{coverage}/autosomes.filtered.vcf.gz'
    threads: 2
    shell:
        '''
        bcftools filter --threads {threads} -i 'AF>0&&INFO/DR2>0' -o {output} {input}
        tabix -p vcf {output}
        '''

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
        annotation = 'annotation_detailed.regions.bed'
        #annotation = 'annotation.bed'
    output:
        'overlaps/{imputed}.{mode}.{coverage}.isec'
    shell:
        '''
        awk -v OFS="\\t" '{{ print $1,$2,$2+1,$3,$4,$5 }}' {input.isec} |\
        bedtools intersect -wo -a - -b {input.annotation} |\
        awk ' {{ if($4~/,/||$5~/,/) {{ ++MA[$10][$6] }} else {{ if (length($4)==1&&length($5)==1) {{++SNP[$10][$6]}} else {{++INDEL[$10][$6]}} }} }} END {{ for (R in SNP) {{ for (V in SNP[R]) {{ print R,"SNP",V,SNP[R][V] }}; for (V in INDEL[R]) {{ print R,"INDEL",V,INDEL[R][V] }}; for (V in MA[R]) {{ print R,"MA",V,MA[R][V] }} }} }} ' |\
        sed 's/ /_/' > {output}
        
        '''

rule make_bed:
    input:
        rules.bcftools_isec.output
    output:
        expand('overlaps/{imputed}.{mode}.{coverage}.bed',sample=config['tissues'],allow_missing=True)
    shell:
        '''
        awk '{{print $1,$2,$2+1}}' {input} > {output}
        '''

workflow._singularity_args = f'-B $TMPDIR -B {PurePath(config["reference"]).parent}'

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

rule bin_TPM_genes:
    input:
        #rules.combine_quan.output
        'aligned_genes/{tissue}/gene_TPM.full.tsv.gz'
    output:
        'happy/{tissue}.TPM_bins'
    localrule: True
    run:
        import polars as pl
        import numpy as np
        df = (pl.read_csv(input[0],separator='\t',ignore_errors=True)
                .drop_nulls()
             )
        df2 = (df.with_columns(df.select(pl.selectors.starts_with('BSWCHEM')).map_rows(np.median)).select(['gene','map'])
                .with_columns(pl.when(pl.col("map")<2).then(0)
                                .when(pl.col('map').is_between(2, 10, closed='left')).then(2)
                                .when(pl.col('map').is_between(10, 50, closed='left')).then(10)
                                .when(pl.col('map').is_between(50, 200, closed='left')).then(50)
                                .otherwise(200).alias('TPM_level'))
                .with_columns(TPM_bin=pl.col('TPM_level').replace({0:'TPM_zero',2:'TPM_low',10:'TPM_moderate',50:'TPM_moderate_high',200:'TPM_high'},default=None))
              )
        df2.select(['gene','TPM_bin']).write_csv(output[0],separator=' ',include_header=False)

rule make_CDS_regions:
    input:
        gtf = 'Bos_taurus.ARS-UCD1.2.108.gtf.gz',
        bins = rules.bin_TPM_genes.output
    output:
        CDS = 'happy/{tissue}.CDS.bed',
        TPMs = multiext('happy/{tissue}','.TPM_zero.bed','.TPM_low.bed','.TPM_moderate.bed','.TPM_moderate_high.bed','.TPM_high.bed'),
        noncoding_exons = 'happy/{tissue}.noncoding_exons.bed',
        non_exons = 'happy/{tissue}.others.bed',
        regions = 'happy/{tissue}.regions.tsv'
    localrule: True
    shell:
        '''
        zgrep -P "^\d" {input.gtf} | awk -v OFS='\\t' '$3=="CDS" {{(match($10,/[A-Z0-9]+/,m)); print $1,$4,$5,m[0]}}' | sort -k1,1n -k2,2n | uniq > {output.CDS}
        zgrep -P "^\d" {input.gtf} | awk -v OFS='\\t' '$3=="exon" {{(match($10,/[A-Z0-9]+/,m)); print $1,$4,$5,m[0]}}' | sort -k1,1n -k2,2n | uniq | bedtools subtract -A -a /dev/stdin -b {output.CDS} > {output.noncoding_exons}
        cat {output.CDS} {output.noncoding_exons} | sort -k1,1n -k2,2n |\
        bedtools complement -g {config[reference]}.fai -i /dev/stdin > {output.non_exons}

        awk -v OFS='\\t' -v T={wildcards.tissue} 'NR==FNR{{a[$1]=$2;next}}$4 in a{{print $1,$2,$3 > "happy/"T"."a[$4]".bed"}}' {input.bins} {output.CDS} 

        echo -e "CDS_zero\\thappy/{wildcards.tissue}.TPM_zero.bed\\nCDS_low\\thappy/{wildcards.tissue}.TPM_low.bed\\nCDS_moderate\\thappy/{wildcards.tissue}.TPM_moderate.bed\\nCDS_moderate_high\\thappy/{wildcards.tissue}.TPM_moderate_high.bed\\nCDS_high\\thappy/{wildcards.tissue}.TPM_high.bed\\nNCE\\t{output.noncoding_exons}\\nintergenic\\t{output.non_exons}" > {output.regions}
        '''

#x={0:'TPM<2',2:'2<TPM<=10',10:'10<TPM<=100',100:'100<TPM<=1000',1000:'1000<TPM'}
#zcat aligned_genes/Testis/gene_TPM.full.tsv.gz | awk 'NR>1 {T=0;n=0;for (i=7; i<=NF;++i) {T+=$i;++n}; print $4,T/n}'

rule happy:
    input:
        vcf_tissue = lambda wildcards: expand('split_vcfs/{tissue}_{coverage}_{imputed}/{sample}.vcf.gz',tissues=wildcards.tissue,allow_missing=True),
        vcf_WGS = lambda wildcards: expand('split_vcfs/WGS_full_{imputed}/{sample}.vcf.gz',allow_missing=True),
        reference = config['reference'],
        regions = rules.make_CDS_regions.output['regions']
    output:
        csv = 'F1/{sample}.{tissue}.{coverage}.{imputed}.summary.csv',
        others = temp(multiext('F1/{sample}.{tissue}.{coverage}.{imputed}','.bcf','.bcf.csi','.extended.csv','.roc.all.csv.gz','.runinfo.json'))
    params:
        _dir = lambda wildcards, output: PurePath(output.csv).with_suffix('').with_suffix('')
    container: '/cluster/work/pausch/alex/software/images/hap.py_latest.sif'
    threads: 1
    resources:
        mem_mb = 5000,
        scratch = '10G',
        #walltime = ''
    shell:
        '''
        /opt/hap.py/bin/hap.py -r {input.reference} --bcf --usefiltered-truth --no-roc --no-json -L --pass-only --scratch-prefix $TMPDIR -X --threads 2 -o {params._dir} --stratification {input.regions} {input.vcf_WGS} {input.vcf_tissue}
        '''

rule gather_happy:
    input:
        #need to consider extended
        #expand(rules.happy.output['others'][2],sample=config['samples'],tissue=config['tissues'],allow_missing=True)
        lambda wildcards: expand(rules.happy.output['others'][2],sample=config['samples'],tissue=(config['tissues'] if wildcards.coverage not in ['hundred','thirty'] else config['bams']),allow_missing=True)
    output:
        'F1/happy.{coverage}.{imputed}.csv'
    localrule: True
    shell:
        '''
        echo -e "variant region truth query recall precision truth_TiTv query_TiTv F1_score sample tissue" > {output}
        for i in {input}
        do
          awk -v I=$(basename $i) -F',' '$2=="*"&&($3~/CDS/||$3=="NCE"||$3=="intergenic")&&$4=="PASS" {{ split(I,a,"."); print $1,$3,$17,$38,$8,$9,$22,$43,$11,a[1],a[2] }}' $i >> {output}
        done
        '''

rule covered_by_variants:
    input:
        vcf = '{tissue}_{coverage}/autosomes.{imputed}.vcf.gz',
        genome = config['reference'] + '.fai'
    output:
        variants = 'coverage/{tissue}.{coverage}.{imputed}.variant.bed',
        covered = 'coverage/{tissue}.{coverage}.{imputed}.variant.csv'
    params:
        window = 1000
    threads: 1
    resources:
        mem_mb = 7500
    shell:
        '''
        bcftools query -i 'GT[*]="alt"' -f '%CHROM\\t%POS\\n' {input.vcf} |\
        awk -v OFS='\\t' '{{print $1,$2,$2+1}}' |\
        bedtools merge -i - -d {params.window} |\
        bedtools slop -g <(head -n 29 {input.genome} | cut -f -2)  -b {params.window} -i - > {output.variants}

        bedtools genomecov -g <(head -n 29 {input.genome} | cut -f -2) -i {output.variants} | awk '$2==1 {{print "{wildcards.tissue}","{wildcards.coverage}","{wildcards.imputed}",$1,$5}}' > {output.covered}
        '''

rule gather_covered:
    input:
        #expand(rules.covered_by_variants.output['covered'],tissue=config['bams'],coverage=config['coverages'],imputed=('imputed','panel'))
        expand(rules.covered_by_variants.output['covered'],tissue=config['bams'],coverage=('full','hundred','thirty'),imputed=('imputed','panel')),
        #expand(rules.covered_by_variants.output['covered'],tissue=config['bams'],coverage=('hundred','thirty'),imputed=('imputed',)),
        expand(rules.covered_by_variants.output['covered'],tissue=config['tissues'],coverage=('five',),imputed=('imputed',))
    output:
        'coverage/variants.csv'
    localrule: True
    shell:
        '''
        echo "tissue coverage imputed chromosome covered" > {output}
        cat {input} >> {output}
        '''

## ASE analysis
rule ASE_analysis:
    input:
        ase = 'ase/{sample}.{tissue}.{coverage}.ase',
        happy = 'F1/{sample}.{tissue}.{coverage}.imputed.bcf',
        exons = 'happy/exons.bed'
    output:
        'ase/{sample}.{tissue}.{coverage}.csv'
    resources:
        mem_mb = 5000
    shell:
        '''
        bedtools intersect -wao \
        -a <(bcftools view -H -R {input.exons} -i 'FORMAT/BVT[0]=="SNP"&FORMAT/BLT[0]=="het"&FORMAT/BD!="TP"' {input.happy} | awk -v OFS='\\t' '{{ print $1,$2,$2+1,substr($11,0,3)}}') \
        -b <(awk -v OFS='\\t' 'NR>1 {{ print $3,$4,$4+1,$21,$13/($11+$12),($8-$9)/$10 }}' {input.ase}) |\
        awk -v S={wildcards.sample} -v T={wildcards.tissue} -v C={wildcards.coverage} '$6!=-1 {{ print S,T,$4,$8,$9,$10,"0"}} END {{print S,T,"0","0","0","0",NR}}' > {output}


        bedtools intersect -wo \
        -a <(bcftools view -H -R {input.exons} -i 'FORMAT/BVT[0]=="SNP"&FORMAT/BLT[0]=="het"&FORMAT/BD=="TP"' {input.happy} | awk -v OFS='\\t' '{{ print $1,$2,$2+1,substr($11,0,3)}}') \
        -b <(awk -v OFS='\\t' 'NR>1 {{ print $3,$4,$4+1,$21,$13/($11+$12),($8-$9)/$10 }}' {input.ase}) |\
        awk -v S={wildcards.sample} -v T={wildcards.tissue} -v C={wildcards.coverage} '$6!=-1 {{ print S,T,$4,$8,$9,$10,"0"}}' >> {output}
        '''

rule gather_ASE:
    input:
        expand(rules.ASE_analysis.output,sample=config['samples'],tissue=config['tissues'],coverage=('full',))
    output:
        'ase/ASE.csv.gz'
    localrule: True
    shell:
        '''
        {{ echo "sample tissue GT pval ASE_adjusted ASE_raw het_errors" ; cat {input} ; }} | pigz -p 2 -c > {output}
        '''
