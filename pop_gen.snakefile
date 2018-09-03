#!/usr/bin/env python3


###########
# GLOBALS #
###########

bayescan_container = 'shub://MarissaLL/singularity-containers:bayescan_2.1@e035d571b5e888e98e5b79f902c30388'
fastsimcoal_container = 'shub://MarissaLL/singularity-containers:fastsimcoal_2.6@37ca431784b209574f517ee09263fca2'
pgdspider_container = 'shub://MarissaLL/singularity-containers:pgdspider_2.1.1.5@e546f843e2b84401284745a766546c90'
r_container = 'shub://TomHarrop/singularity-containers:r_3.5.0@490e801d406497fa461377d17b3b339b'
stacks2beta_container = 'shub://TomHarrop/singularity-containers:stacks_2.0beta9@bb2f9183318871f6228b51104056a2d0'

bayescan_runs = ['compared_4pops', 'compared_para', 'compared_2pops']

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
        'output/090_pop_genet/populations.snps.vcf',
        'output/060_pop_genet/filtered_snps.gds'
        #expand('output/070_bayescan/{run}.sel',
        #       run = bayescan_runs)






rule populations_three:
    input:
        stacks('output/040_stacks/catalog.fa.gz'),
        stacks('output/040_stacks/catalog.calls'),
        popmap = 'test/r0.2_filtered_popmap.txt',
        whitelist = 'test/m0_whitelist.txt'
    output:
        'output/090_pop_genet/populations.snps.vcf'
    params:
        stacks_dir = 'output/040_stacks',
        outdir = 'output/090_pop_genet'
    singularity:
        stacks2beta_container
    threads:
        50
    log:
        'output/logs/090_pop_genet/pops_stats.log'
    shell:
        'populations '
        '-P {params.stacks_dir} '
        '-M {input.popmap} '
        '-O {params.outdir} '
        '-W {input.whitelist} '
        '-t {threads} '
        '-r 0 '
        '--vcf '
        '&> {log}'

rule bayescan:
    input:
        genotypes = expand('output/070_bayescan/{run}.geste-outputformat',
                            run = bayescan_runs)
    output:
        expand('output/070_bayescan/{run}.sel',
               run = bayescan_runs),
        expand('output/070_bayescan/{run}_fst.txt',
               run = bayescan_runs)
    params:
        outdir = 'output/070_bayescan',
        outname = expand('{run}', 
                         run = bayescan_runs)
    singularity:
        bayescan_container
    threads:
        50
    log:
        expand('output/logs/070_bayescan/{run}.log',
               run = bayescan_runs)
    shell:
        'bayescan_2.1 '
        '{input.genotypes} '
        '-od {params.outdir} '
        '-o {params.outname} '
        '-pilot 5000 '
        '-nbp 20 '
        '-burn 15000 '
        '-n 30000 '
        '-thin 10 '
        '-pr_odds 500 '
        '-out_pilot '
        '-out_freq '
        '&> {log}'


#Convert vcf data into geste format for bayescan. 
# Note that the popmap is actually specified within the spid file
# This has only produced one of the three output files. FIX...
rule convert_bayescan_inputs: 
    input:
        vcf = 'output/060_pop_genet/populations.snps.vcf',
        spid = expand('data/{run}.spid',
                      run = bayescan_runs),
        popmap = expand('output/070_bayescan/popmap_{run}.txt',
                        run = bayescan_runs)
    output:
        expand('output/070_bayescan/{run}.geste-outputformat',
               run = bayescan_runs)
    params:
        in_format = 'VCF',
        out_format = 'GESTE_BAYE_SCAN',
        out_path = expand('output/070_bayescan/{run}.geste',
                          run = bayescan_runs)
    singularity:
        pgdspider_container
    threads:
        50
    log:
        expand('output/logs/070_bayescan/pgdspider_{run}.log',
               run = bayescan_runs)
    shell:
        'java -jar /opt/pgdspider/PGDSpider2-cli.jar '
        '-inputfile {input.vcf} '
        '-inputformat {params.in_format} '
        '-outputfile {params.out_path}'
        '-outputformat {params.out_format} '
        '-spid {input.spid} '
        '&> {log}'

rule make_bayescan_popmaps:
    input: 
        popmap = 'output/060_pop_genet/r0.8_filtered_popmap.txt',
        para_data = process_reads('output/010_config/tidy_sample_info.tsv')
    output:
        popmap_4pops = 'output/070_bayescan/popmap_compared_4pops.txt',
        popmap_para = 'output/070_bayescan/popmap_compared_para.txt',
        popmap_2pops = 'output/070_bayescan/popmap_compared_2pops.txt'
    singularity:
        r_container
    threads:
        25
    log:
        'output/logs/070_bayescan/make_bayescan_popmaps.log'
    script:
        'src/make_bayescan_popmaps.R'

rule convert_gds_again:
    input:
        'output/060_pop_genet/populations.snps.vcf'
    output:
        'output/060_pop_genet/filtered_snps.gds'
    threads:
        25
    log:
        'output/logs/060_pop_genet/filtered_gds_convert.log'
    script:
        'src/convert_gds.R'


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

# Used for DAPC here

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

# Filter SNPs on MAF and missing rate, also filter samples by missing rate
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















