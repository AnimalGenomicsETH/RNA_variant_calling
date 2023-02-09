
include: 'snakepit/association.smk'
include: 'snakepit/overlap.smk'
include: 'snakepit/coverage.smk'

rule all:
    input:
        expand('eQTL/WGS_{tissue}/{_pass}.{MAF}.txt',tissue=config['vcf'],_pass=('conditionals',),MAF=format_MAF(config['MAF'])),
        expand('eQTL/{tissue}_{tissue}/{_pass}.{MAF}.txt',tissue=config['vcf'],_pass=('conditionals',),MAF=format_MAF(config['MAF'])),
        expand('replication/{tissue}.{mode}.csv',tissue=config['vcf'],mode=('best',)),
        'ase_metrics.csv',
        'overlaps/Unrevised.none.isec',
        'happy_metrics.csv',
        'annotated_coverage.bed.gz',
        'genome_coverage.csv'
