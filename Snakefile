


include: 'snakepit/variant_calling.smk'
include: 'snakepit/association.smk'
include: 'snakepit/overlap.smk'
include: 'snakepit/coverage.smk'

rule all:
    input:
        expand('{tissue}_{coverage}/autosomes.imputed.vcf.gz',tissue=('Vas_deferens','Testis','Epididymis_head'),coverage=config['coverages']),
        'WGS_full/autosomes.imputed.vcf.gz',
        expand('eQTL/WGS_{tissue}.{MAF}.genevars.txt',tissue=config['vcf'],MAF=format_MAF(config['MAF'])),
        expand('eQTL/{tissue}_{tissue}.{MAF}.genevars.txt',tissue=config['vcf'],MAF=format_MAF(config['MAF'])),
        expand('replication/{tissue}.{tissue}.replicated',tissue=config['vcf']),
        expand('replication/{tissue}.WGS.replicated',tissue=config['vcf']),
        expand('{tissue}/autosomes.Unrevised.imputed.vcf.gz',tissue=config['vcf']),
        'ase_metrics.csv',
        'overlaps/Unrevised.none.isec',
        'happy_metrics.csv',
        'annotated_coverage.bed.gz',
        'genome_coverage.csv'
