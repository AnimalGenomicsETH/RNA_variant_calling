
include: 'snakepit/association.smk'
include: 'snakepit/overlap.smk'
include: 'snakepit/coverage.smk'
include: 'snakepit/imputation.smk'

rule all:
    input:
        expand('eQTL/WGS_{tissue}/{_pass}.{MAF}.txt',tissue=config['vcf'],_pass=('conditionals',),MAF=format_MAF(config['MAF'])),
        expand('eQTL/{tissue}_{tissue}/{_pass}.{MAF}.txt',tissue=config['vcf'],_pass=('conditionals',),MAF=format_MAF(config['MAF'])),
        expand('replication/{tissue}.{tissue}.csv',tissue=config['vcf']),
        expand('replication/{tissue}.WGS.csv',tissue=config['vcf']),
        expand('{tissue}/autosomes.Unrevised.imputed.vcf.gz',tissue=config['vcf']),
        'ase_metrics.csv',
        'overlaps/Unrevised.none.isec',
        'happy_metrics.csv',
        'annotated_coverage.bed.gz',
        'genome_coverage.csv'
