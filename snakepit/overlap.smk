from pathlib import PurePath

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

#awk '$5=="1110" {print $3,$4}' overlaps/imputed.none.full.txt | grep -v "," | awk 'length($0)==3' | sort | uniq -c
#awk '$5=="1110" {print $3,$4}' overlaps/imputed.none.full.txt | grep -v "," | awk 'length($0)>3 {print length($1)-length($2)}' | sort | uniq -c | sort -k2,2nr | less

rule count_isec:
    input:
        isec = rules.bcftools_isec.output,
        annotation = 'annotation.bed'
    output:
        'overlaps/{imputed}.{mode}.{coverage}.isec'
    shell:
        '''
        awk -v OFS="\\t" '{{ print $1,$2,$2+1,$3,$4,$5 }}' {input.isec} |\
        bedtools intersect -wo -a - -b {input.annotation} |\
        awk ' {{ if($4~/,/||$5~/,/) {{ ++MA[$10=="intergenic"][$6] }} else {{ if (length($4)==1&&length($5)==1) {{++SNP[$10=="intergenic"][$6]}} else {{++INDEL[$10=="intergenic"][$6]}} }} }} END {{for (key in SNP[0]) {{ print "genic","SNP",key,SNP[0][key] }} for (key in SNP[1]) {{ print "intergenic","SNP",key,SNP[1][key] }} for (key in INDEL[0]) {{ print "genic","INDEL",key,INDEL[0][key] }} for (key in INDEL[1]) {{ print "intergenic","INDEL",key,INDEL[1][key] }} for (key in MA[0]) {{ print "genic","MA",key,MA[0][key]}} for (key in MA[1]) {{ print "intergenic","MA",key,MA[1][key]}} }}' > {output}
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

rule happy:
    input:
        vcf_tissue = lambda wildcards: expand('split_vcfs/{tissue}_{coverage}_{imputed}/{sample}.vcf.gz',tissues=wildcards.tissue,allow_missing=True),
        vcf_WGS = lambda wildcards: expand('split_vcfs/WGS_full_{imputed}/{sample}.vcf.gz',allow_missing=True),
        reference = config['reference'],
        regions = 'happy/regions.tsv'
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
        walltime = '24h'
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
          awk -v I=$(basename $i) -F',' '$2=="*"&&($3=="others"||$3=="exons")&&$4=="PASS" {{ split(I,a,"."); print $1,$3,$17,$38,$8,$9,$22,$43,$11,a[1],a[2] }}' $i >> {output}
        done
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

# bcftools query -f '%CHROM\t%POS\n' Vas_deferens_full/autosomes.imputed.vcf.gz| awk -v OFS='\t' '{print $1,$2,$2+1}' | bedtools merge -i - -d 10000 > Vas_deferens_variant_cover.bed
#for i in WGS Testis Epididymis_head Vas_deferens; do bedtools genomecov -g <(head -n 29 ../REF_DATA/ARS-UCD1.2_Btau5.0.1Y.fa.fai | cut -f -2) -i ${i}_variant_cover.bed  | awk -v T=${i} '$2==1 {print T,$1,$5}'; done
