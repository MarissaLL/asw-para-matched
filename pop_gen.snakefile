#!/usr/bin/env python3


###########
# GLOBALS #
###########

bayescan_container = 'shub://MarissaLL/singularity-containers:bayescan_2.1@e035d571b5e888e98e5b79f902c30388'
fastsimcoal_container = 'shub://MarissaLL/singularity-containers:fastsimcoal_2.6@37ca431784b209574f517ee09263fca2'
pgdspider_container = 'shub://MarissaLL/singularity-containers:pgdspider_2.1.1.5@e546f843e2b84401284745a766546c90'
stacks2beta_container = 'shub://TomHarrop/singularity-containers:stacks_2.0beta9@bb2f9183318871f6228b51104056a2d0'


#########
# RULES #
#########

subworkflow process_reads:
    snakefile: 'process_reads.snakefile'

subworkflow stacks:
    snakefile: 'stacks.snakefile'

rule target:
    input:
        'output/060_pop_genet/populations.snps.vcf',
        'output/070_pop_tests/compared_pops_prop.txt',
        'output/070_pop_tests/compared_para_prop.txt'








# Make a new SPID. CHANGE THIS TO GENEPOP?
rule make_arlequin_input:
     input:
        vcf = 'output/060_pop_genet/populations.snps.vcf',
        spid = 'data/convert_arlequin.spid'
    output:
        'output/071_fastsimcoal/.arp'
    params:
        in_format = 'VCF',
        out_format = 'ARLEQUIN'
    singularity:
        pgdspider_container
    log:
        'output/logs/071_fastsimcoal/pgdspider.log'
    shell:
        'java -jar /opt/pgdspider/PGDSpider2-cli.jar '
        '-inputfile {input.ped} '
        '-inputformat {params.in_format} '
        '-outputfile {output}'
        '-outputformat {params.out_format} '
        '-spid {input.spid} '
        '&> {log}'


rule bayescan_pops:
    input:
        genotypes = 'output/070_pop_tests/compare_pops.geste-outputformat'
    output:
        'output/070_pop_tests/compared_pops_prop.txt'
    params:
        outdir = 'output/070_pop_tests',
        outname = 'compared_pops'
    singularity:
        bayescan_container
    threads:
        50
    log:
        'output/logs/070_pop_tests/bayescan_pops.log'
    shell:
        'bayescan_2.1 '
        '{input.genotypes} '
        '-od {params.outdir} '
        '-o {params.outname} '
        '-pilot 5000 '
        '-burn 5000 '
        '-n 10000 '
        '-pr_jump 0.1 '
        '-pr_pref 0.5 '
        '-pr_odds 500 '
        '-out_pilot '
        '-out_freq '
        '&> {log}'

rule bayescan_para:
    input:
        genotypes = 'output/070_pop_tests/compare_para.geste-outputformat'
    output:
        'output/070_pop_tests/compared_para_prop.txt'
    params:
        outdir = 'output/070_pop_tests',
        outname = 'compared_para'
    singularity:
        bayescan_container
    log:
        'output/logs/070_pop_tests/bayescan_para.log'
    shell:
        'bayescan_2.1 '
        '{input.genotypes} '
        '-od {params.outdir} '
        '-o {params.outname} '
        '-pilot 15000 '
        '-burn 15000 '
        '-n 30000 '
        '-pr_jump 0.1 '
        '-pr_pref 0.5 '
        '-pr_odds 500 '
        '-out_pilot '
        '-out_freq '
        '&> {log}'

rule make_bayescan_input_para:
    input:
        ped = 'output/070_pop_tests/pop_para_snps.ped',
        spid = 'data/convert_para.spid'
    output:
        'output/070_pop_tests/compare_para.geste-outputformat'
    params:
        in_format = 'PED',
        out_format = 'GESTE_BAYE_SCAN'
    singularity:
        pgdspider_container
    log:
        'output/logs/070_pop_tests/pgdspider_2pops.log'
    shell:
        'java -jar /opt/pgdspider/PGDSpider2-cli.jar '
        '-inputfile {input.ped} '
        '-inputformat {params.in_format} '
        '-outputfile {output}'
        '-outputformat {params.out_format} '
        '-spid {input.spid} '
        '&> {log}'

rule make_bayescan_input_pops:
    input:
        ped = 'output/070_pop_tests/pop_para_snps.ped',
        spid = 'data/convert_pops.spid'
    output:
        'output/070_pop_tests/compare_pops.geste-outputformat'
    params:
        in_format = 'PED',
        out_format = 'GESTE_BAYE_SCAN'
    singularity:
        pgdspider_container
    log:
        'output/logs/070_pop_tests/pgdspider_4pops.log'
    shell:
        'java -jar /opt/pgdspider/PGDSpider2-cli.jar '
        '-inputfile {input.ped} '
        '-inputformat {params.in_format} '
        '-outputfile {output}'
        '-outputformat {params.out_format} '
        '-spid {input.spid} '
        '&> {log}'

rule add_popdata_to_ped:
    input:
        ped_file = 'output/060_pop_genet/snps.ped',
        para_info = process_reads('output/010_config/tidy_sample_info.tsv')
    output:
        'output/070_pop_tests/pop_para_snps.ped'
    log:
        'output/logs/070_pop_tests/add_popdata_to_ped.log'
    script:
        'src/add_popdata_to_ped.R'

# Run populations again on filtered data to get population summary statistics
rule populations_stats:
    input:
        stacks('output/040_stacks/catalog.fa.gz'),
        stacks('output/040_stacks/catalog.calls'),
        popmap = 'output/060_pop_genet/r0.8_filtered_popmap.txt',
        whitelist = 'output/060_pop_genet/whitelist.txt'
    output:
        'output/060_pop_genet/populations.snps.vcf',
        'output/060_pop_genet/populations.plink.ped'
    params:
        stacks_dir = 'output/040_stacks',
        outdir = 'output/060_pop_genet'
    singularity:
        stacks2beta_container
    threads:
        50
    log:
        'output/logs/060_pop_genet/pops_stats.log'
    shell:
        'populations '
        '-P {params.stacks_dir} '
        '-M {input.popmap} '
        '-O {params.outdir} '
        '-W {input.whitelist} '
        '-t {threads} '
        '-r 0 '
        '--genepop '
        '--plink '
        '--vcf '
        '--hwe '
        '--fstats '
        '&> {log}'

# Write a new popmap and a new whitelist with only the samples and SNPs that passed all the filtering
rule make_whitelist_popmap:
    input:
        plink_file = 'output/060_pop_genet/plink.raw'
    output:
        whitelist = 'output/060_pop_genet/whitelist.txt',
        popmap = 'output/060_pop_genet/r0.8_filtered_popmap.txt'
    log:
        'output/logs/060_pop_genet/make_whitelist_popmap.log'
    script:
        'src/make_whitelist_popmap.R'

# Convert ped to plink format, with unknown chromosomes
rule convert_to_plinkraw:
    input:
        ped = 'output/060_pop_genet/snps.ped',
        map = 'output/060_pop_genet/snps.map'
    output:
        'output/060_pop_genet/plink.raw'
    params:
        workdir = 'output/060_pop_genet'
    threads:
        25
    log:
        'output/logs/060_pop_genet/convert_plinkraw.log'
    shell:
        'cd {params.workdir} || exit 1 ; '
        'plink '
        '--ped snps.ped '
        '--map snps.map '
        '--recode A '
        '--aec'

# Filter SNPs on MAF and missing rate, also filter by sample missing rate
rule filter_snps_indivs:
    input:
        'output/060_pop_genet/snps.gds'
    output:
        'output/060_pop_genet/snps.ped',
        'output/060_pop_genet/snps.map'
    params:
        maf = 0.05,
        missing_rate = 0.2,
        sample_missing_quantile = 0.8,
        ped_file = 'output/060_pop_genet/snps'
    threads:
        25
    log:
        'output/logs/060_pop_genet/filter_snps_indivs.log'
    script:
        'src/filter_snps.R'

# Convert VCF to GDS format, keeping biallelic SNPs only
rule convert_to_gds:
    input:
        stacks('output/050_stacks_pops/r0/populations.snps.vcf')
    output:
        'output/060_pop_genet/snps.gds'
    threads:
        25
    log:
        'output/logs/060_pop_genet/gds_convert.log'
    script:
        'src/convert_gds.R'















