wildcard_constraints:
    _pass = r'permutations|conditionals|nominals',
    chunk = r'\d+',
    MAF = r'\d+',
    tissue = r'WGS|Epididymis_head|Vas_deferens|Testis',
    expression = r'WGS|Epididymis_head|Vas_deferens|Testis',
    coverage = r'\w+'

config['tissues'] = list(config['bams'].keys())[1:]

include: 'snakepit/variant_calling.smk'
include: 'snakepit/TPM.smk'
include: 'snakepit/association.smk'
include: 'snakepit/overlap.smk'
include: 'snakepit/coverage.smk'

rule all:
    input:
        ## Variant calling
        expand('{tissue}_{coverage}/autosomes.imputed.vcf.gz',tissue=('Vas_deferens','Testis','Epididymis_head'),coverage=config['coverages']),
        expand('WGS_{coverage}/autosomes.imputed.vcf.gz',coverage=('full','hundred','thirty')),
        
        ## Variant analysis
        expand('overlaps/imputed.none.{coverage}.isec',coverage=config['coverages']),
        expand('F1/happy.{coverage}.{imputed}.csv',coverage=config['coverages'],imputed=('imputed','Unrevised')),
        
        ## Bam coverage/ASE
        expand('coverage/annotated.{coverage}.bed.gz',coverage=config['coverages']),
        expand('coverage/genome.{coverage}.csv',coverage=config['coverages']),
        expand('ase/metrics.{coverage}.csv',coverage=config['coverages']),
        
        ## Variant comparison
        'coverage/variants.csv',
        expand('F1/happy.{coverage}.{imputed}.csv',coverage=('full','hundred','thirty','five'),imputed=('imputed',)),
        
        ## eQTL analysis
        expand('eQTL/WGS_full_{tissue}_{coverage}_filtered/conditionals.{MAF}.txt.gz',tissue=config['tissues'],MAF=format_MAF(config['MAF']),coverage=config['coverages']),
        expand('eQTL/{tissue}_{coverage}_{tissue}_{coverage}_filtered/conditionals.{MAF}.txt.gz',tissue=config['tissues'],MAF=format_MAF(config['MAF']),coverage=config['coverages'])
