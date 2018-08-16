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
        expand('output/070_bayescan/{run}.sel',
               run = bayescan_runs)
        #'output/080_pop_simulations/fsc_sfs.txt',
        #'fsc_test/test.arp'





# rule subset_snps:
#     input:
#         stacks('output/040_stacks/catalog.fa.gz'),
#         stacks('output/040_stacks/catalog.calls'),
#         popmap = 'output/080_test/shortened_popmap.txt',
#         whitelist = 'output/080_test/whitelist.txt'
#     output:
#         'output/080_test/populations.snps.genepop'
#     params:
#         stacks_dir = 'output/040_stacks',
#         outdir = 'output/080_test'
#     singularity:
#         stacks2beta_container
#     threads:
#         50
#     log:
#         'output/logs/080_test/pops_stats.log'
#     shell:
#         'populations '
#         '-P {params.stacks_dir} '
#         '-M {input.popmap} '
#         '-O {params.outdir} '
#         '-W {input.whitelist} '
#         '-t {threads} '
#         '-r 0 '
#         '--genepop '
#         '--write_random_snp '
#         '&> {log}'  

# rule simulate_sfs:
#     input:
#         par = 'fsc_test/test.par',
#         obs = 'fsc_test/test_MAFpop0.obs'
#     output:
#         'fsc_test/test.arp'
#     singularity:
#         fastsimcoal_container
#     shell:
#         'fsc26 '
#         '--ifile {input.par} '
#         '--noarloutput '
#         '--dnatosnp 0 '
#         '--msfs '
#         '--numsims 1000 '
#         '--seed 1234 '


# rule calculate_obs_sfs:
#     input:
#         vcf = 'output/060_pop_genet/populations.snps.vcf',
#         popmap = 'output/060_pop_genet/r0.8_filtered_popmap.txt'
#     output:
#         sfs_popmap = 'output/080_pop_simulations/sfs_popmap.txt',
#         dadi = 'output/080_pop_simulations/dadi.txt',
#         fsc = 'output/080_pop_simulations/fsc_sfs.txt'
#     params:
#         vcf2sfs = 'vcf2sfs.R'
#     singularity:
#         r_container
#     threads:
#         20
#     log:
#         'output/logs/080_pop_simulations/sfs_calcs.log'
#     script:
#         'src/sfs_calcs.R'







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

rule make_bayescan_input_2pops:
    input: 
        vcf = 'output/060_pop_genet/populations.snps.vcf',
        spid = 'data/convert_2pops.spid',
        popmap = 'output/070_bayescan/popmap_2pops.txt'
    output:
        'output/070_bayescan/compared_2pops.geste-outputformat'
    params:
        in_format = 'VCF',
        out_format = 'GESTE_BAYE_SCAN',
        out_path = 'output/070_bayescan/compared_2pops.geste'
    singularity:
        pgdspider_container
    threads:
        50
    log:
        'output/logs/070_bayescan/pgdspider_para.log'
    shell:
        'java -jar /opt/pgdspider/PGDSpider2-cli.jar '
        '-inputfile {input.vcf} '
        '-inputformat {params.in_format} '
        '-outputfile {params.out_path}'
        '-outputformat {params.out_format} '
        '-spid {input.spid} '
        '&> {log}'

rule make_bayescan_input_para:
    input: 
        vcf = 'output/060_pop_genet/populations.snps.vcf',
        spid = 'data/convert_para.spid',
        popmap = 'output/070_bayescan/popmap_para.txt'
    output:
        'output/070_bayescan/compared_para.geste-outputformat'
    params:
        in_format = 'VCF',
        out_format = 'GESTE_BAYE_SCAN',
        out_path = 'output/070_bayescan/compared_para.geste'
    singularity:
        pgdspider_container
    threads:
        50
    log:
        'output/logs/070_bayescan/pgdspider_para.log'
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
        popmap_para = 'output/070_bayescan/popmap_para.txt',
        popmap_2pops = 'output/070_bayescan/popmap_2pops.txt'
    singularity:
        r_container
    log:
        'output/logs/070_bayescan/make_bayescan_popmaps.log'
    script:
        'src/make_bayescan_popmaps.R'


#Convert vcf data into geste format (with 4 populations)for bayescan. 
# Note that the popmap is actually specified within the spid file
rule make_bayescan_input_4pops: 
    input:
        vcf = 'output/060_pop_genet/populations.snps.vcf',
        spid = 'data/convert_4pops.spid',
        popmap = 'output/060_pop_genet/r0.8_filtered_popmap.txt'
    output:
        'output/070_bayescan/compared_4pops.geste-outputformat'
    params:
        in_format = 'VCF',
        out_format = 'GESTE_BAYE_SCAN',
        out_path = 'output/070_bayescan/compared_4pops.geste'
    singularity:
        pgdspider_container
    threads:
        50
    log:
        'output/logs/070_bayescan/pgdspider_4pops.log'
    shell:
        'java -jar /opt/pgdspider/PGDSpider2-cli.jar '
        '-inputfile {input.vcf} '
        '-inputformat {params.in_format} '
        '-outputfile {params.out_path}'
        '-outputformat {params.out_format} '
        '-spid {input.spid} '
        '&> {log}'


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















