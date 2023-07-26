

wildcard_constraints:
    _pass = r'permutations|conditionals|nominals',
    chunk = r'\d+',
    MAF = r'\d+',
    tissue = r'WGS|Epididymis_head|Vas_deferens|Testis',
    expression = r'WGS|Epididymis_head|Vas_deferens|Testis'

config['tissues'] = list(config['bams'].keys())[1:]

include: 'snakepit/variant_calling.smk'
include: 'snakepit/TPM.smk'
include: 'snakepit/association.smk'
include: 'snakepit/overlap.smk'
include: 'snakepit/coverage.smk'


rule all:
    input:
        ## Variant calling
        #expand('{tissue}_{coverage}/autosomes.imputed.vcf.gz',tissue=('Vas_deferens','Testis','Epididymis_head'),coverage=config['coverages']),
        #'WGS_full/autosomes.imputed.vcf.gz',
        ## Variant analysis
        #expand('overlaps/imputed.none.{coverage}.isec',coverage=config['coverages']),
        #expand('F1/happy.{coverage}.{imputed}.csv',coverage=config['coverages'],imputed=('Unrevised','imputed')),
        ## Bam coverage/ASE
        #expand('coverage/annotated.{coverage}.bed.gz',coverage=config['coverages']),
        #expand('coverage/genome.{coverage}.csv',coverage=config['coverages']),
        #expand('ase/metrics.{coverage}.csv',coverage=config['coverages']),
        ## eQTL analysis
        expand('eQTL/WGS_full_{tissue}_{coverage}/nominals.{MAF}.txt.gz',tissue=config['tissues'],MAF=format_MAF(config['MAF']),coverage=config['coverages']),
        expand('eQTL/{tissue}_{coverage}_{tissue}_{coverage}/nominals.{MAF}.txt.gz',tissue=config['tissues'],MAF=format_MAF(config['MAF']),coverage=config['coverages']),
        #expand('replication/{tissue}.{tissue}.{coverage}.replicated',tissue=config['vcf'],coverage=config['coverages']),
        #expand('replication/{tissue}.WGS.{coverage}.replicated',tissue=config['vcf'],coverage=config['coverages'])
