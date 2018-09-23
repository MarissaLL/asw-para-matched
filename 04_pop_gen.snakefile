#!/usr/bin/env python3


###########
# GLOBALS #
###########

bayescan_container = 'shub://MarissaLL/singularity-containers:bayescan_2.1@e035d571b5e888e98e5b79f902c30388'
fastsimcoal_container = 'shub://MarissaLL/singularity-containers:fastsimcoal_2.6@37ca431784b209574f517ee09263fca2'
pgdspider_container = 'shub://MarissaLL/singularity-containers:pgdspider_2.1.1.5@e546f843e2b84401284745a766546c90'
plink_container = 'shub://TomHarrop/singularity-containers:plink_1.90beta5@43e5ecf38b3490b64a5d7e1f5ead046d'
r_container = 'shub://MarissaLL/singularity-containers:r_3.5.0@1078bd77b7e550e72486881defed9bad'
stacks2beta_container = 'shub://TomHarrop/singularity-containers:stacks_2.0beta9@bb2f9183318871f6228b51104056a2d0'
vcftools_container = 'shub://MarissaLL/singularity-containers:vcftools_0.1.17@230db32b3097775cd51432092f9cbcb1'


bayescan_runs = ['compared_island', 'compared_para', 'compared_2pops', 
                 'compared_invermay', 'compared_lincoln', 'compared_ruakura',
                 'compared_ruakura_poa', 'compared_4pops']

bayescan_subsets = ['compared_lincoln', 'compared_ruakura', 'compared_ruakura_poa', 
                    'compared_invermay', 'compared_2pops']

bayescan_full = ['compared_island', 'compared_4pops', 'compared_para']

#########
# RULES #
#########

subworkflow process_reads:
    snakefile: 'process_reads.snakefile'

subworkflow stacks:
    snakefile: 'stacks.snakefile'

rule target:
    input:
        # 'output/060_pop_genet/populations.snps.vcf',
        # 'output/070_bayescan/popmap_compared_invermay.txt',
        expand('output/070_bayescan/{bayescan_run}.sel',
              bayescan_run = bayescan_runs),
        'output/070_bayescan/compared_island_prhi.sel'
        # 'output/071_DAPC/dapc_para_results.tsv',
        # expand('output/070_bayescan/{bayescan_subset}.geste-outputformat',
        #         bayescan_subset = bayescan_subsets),
        # expand('output/070_bayescan/{bayescan_full}.geste-outputformat',
        #         bayescan_full = bayescan_full)

# Run bayescan with higher prior odds to detect outlying SNPs
rule bayescan_prhi:
    input:
        genotypes = 'output/070_bayescan/compared_island.geste-outputformat'
    output:
        'output/070_bayescan/compared_island_prhi.sel',
        'output/070_bayescan/compared_island_prhi.txt'
    params:
        outdir = 'output/070_bayescan',
        outname = 'compared_island_prhi'
    singularity:
        bayescan_container
    threads:
        50
    log:
        'output/logs/070_bayescan/compared_island_prhi.log'
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
        '-pr_odds 10000 '
        '-out_pilot '
        '-out_freq '
        '&> {log}'



# Run bayescan to detect outlying SNPs
rule bayescan:
    input:
        genotypes = 'output/070_bayescan/{bayescan_run}.geste-outputformat'
    output:
        'output/070_bayescan/{bayescan_run}.sel',
        'output/070_bayescan/{bayescan_run}_fst.txt'
    params:
        outdir = 'output/070_bayescan',
        outname = '{bayescan_run}'
    singularity:
        bayescan_container
    threads:
        50
    log:
        'output/logs/070_bayescan/{bayescan_run}.log'
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


#Convert SNP data in VCF format into geste format for bayescan. For the full vcf
# Note that the location of the popmap is actually specified within the spid file
# It is only specified here to link the dependencies of the rules
# rule convert_full_pop_bayescan_inputs: 
#     input:
#         vcf = 'output/060_pop_genet/populations.snps.vcf',
#         spid = 'data/{bayescan_full}.spid',
#         popmap = 'output/070_bayescan/popmap_{bayescan_full}.txt'
#     output:
#         'output/070_bayescan/{bayescan_full}.geste-outputformat'
#     params:
#         in_format = 'VCF',
#         out_format = 'GESTE_BAYE_SCAN',
#         out_path = 'output/070_bayescan/{bayescan_full}.geste'
#     singularity:
#         pgdspider_container
#     threads:
#         50
#     log:
#         'output/logs/070_bayescan/pgdspider_{bayescan_full}.log'
#     shell:
#         'java -jar /opt/pgdspider/PGDSpider2-cli.jar '
#         '-inputfile {input.vcf} '
#         '-inputformat {params.in_format} '
#         '-outputfile {params.out_path}'
#         '-outputformat {params.out_format} '
#         '-spid {input.spid} '
#         '&> {log}'


## Convert data to geste format for bayescan, for the subset vcfs
## Note that the location of the popmap is actually specified within the spid file
## It is only specified here to link the dependencies of the rules
# rule convert_subset_bayescan_inputs:
#     input:
#         vcf = 'output/070_bayescan/{bayescan_subset}.recode.vcf',
#         spid = 'data/{bayescan_subset}.spid',
#         popmap = 'output/070_bayescan/popmap_{bayescan_subset}.txt'
#     output:
#         'output/070_bayescan/{bayescan_subset}.geste-outputformat'
#     params:
#         in_format = 'VCF',
#         out_format = 'GESTE_BAYE_SCAN',
#         out_path = 'output/070_bayescan/{bayescan_subset}.geste'
#     singularity:
#         pgdspider_container
#     threads:
#         50
#     log:
#         'output/logs/070_bayescan/pgdspider_{bayescan_subset}.log'
#     shell:
#         'java -jar /opt/pgdspider/PGDSpider2-cli.jar '
#         '-inputfile {input.vcf} '
#         '-inputformat {params.in_format} '
#         '-outputfile {params.out_path}'
#         '-outputformat {params.out_format} '
#         '-spid {input.spid} '
#         '&> {log}'

# Subset vcfs to only include individuals needed in each bayescan run
rule subset_vcfs:
    input:
        individual_list = 'output/070_bayescan/{bayescan_subset}_indivs.txt',
        full_vcf = 'output/060_pop_genet/populations.snps.vcf'
    output: 
        subset_vcf = 'output/070_bayescan/{bayescan_subset}.recode.vcf'
    params:
        outname = 'output/070_bayescan/{bayescan_subset}'
    singularity:
        vcftools_container
    threads:
        25
    log:
        'output/logs/070_bayescan/{bayescan_subset}_subset_vcf.log'
    shell:
        'vcftools '
        '--vcf {input.full_vcf} '
        '--out {params.outname} '
        '--keep {input.individual_list} '
        '--recode '
        '&> {log}'

#  Make lists of individuals to keep when subsetting vcfs
rule make_indiv_lists:
    input:
        popmap_2pops = 'output/070_bayescan/popmap_compared_2pops.txt',
        popmap_ruakura = 'output/070_bayescan/popmap_compared_ruakura.txt',
        popmap_ruakura_poa = 'output/070_bayescan/popmap_compared_ruakura_poa.txt',
        popmap_lincoln = 'output/070_bayescan/popmap_compared_lincoln.txt',
        popmap_invermay = 'output/070_bayescan/popmap_compared_invermay.txt'
    output:
        indivs_2pops = 'output/070_bayescan/compared_2pops_indivs.txt',
        indivs_ruakura = 'output/070_bayescan/compared_ruakura_indivs.txt',
        indivs_ruakura_poa = 'output/070_bayescan/compared_ruakura_poa_indivs.txt',
        indivs_lincoln = 'output/070_bayescan/compared_lincoln_indivs.txt',
        indivs_invermay = 'output/070_bayescan/compared_invermay_indivs.txt'
    singularity:
        r_container
    log:
        'output/logs/070_bayescan/make_indiv_lists.log'
    script:
        'src/define_vcf_individuals.R'


# Make popmap files specifying which individuals to use for each bayescan run
rule make_bayescan_popmaps:
    input: 
        popmap = 'output/060_pop_genet/r0.8_filtered_popmap.txt',
        para_data = process_reads('output/010_config/tidy_sample_info.tsv')
    output:
        popmap_4pops = 'output/070_bayescan/popmap_compared_4pops.txt',
        popmap_para = 'output/070_bayescan/popmap_compared_para.txt',
        popmap_2pops = 'output/070_bayescan/popmap_compared_2pops.txt',
        popmap_island = 'output/070_bayescan/popmap_compared_island.txt',
        popmap_ruakura = 'output/070_bayescan/popmap_compared_ruakura.txt',
        popmap_ruakura_poa = 'output/070_bayescan/popmap_compared_ruakura_poa.txt',
        popmap_lincoln = 'output/070_bayescan/popmap_compared_lincoln.txt',
        popmap_invermay = 'output/070_bayescan/popmap_compared_invermay.txt'
    singularity:
        r_container
    threads:
        25
    log:
        'output/logs/070_bayescan/make_bayescan_popmaps.log' 
        # The more informative logs actually get output automatically as Verif.txt files
    script:
        'src/make_bayescan_popmaps.R'


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
    singularity:
        r_container
    log:
        'output/logs/060_pop_genet/make_whitelist_popmap.log'
    script:
        'src/make_whitelist_popmap.R'

# Do DAPC on the SNP data, with various subsets of the data and groupings
rule run_DAPC:
    input:
        plink_file = 'output/060_pop_genet/plink.raw',
        para_file =  'output/010_config/tidy_sample_info.tsv'
    output:
        pca_results = 'output/071_DAPC/pca_results.tsv',
        dapc_para_results = 'output/071_DAPC/dapc_para_results.tsv',
        dapc_pop_results = 'output/071_DAPC/dapc_pop_results.tsv',
        dapc_invermay_results = 'output/071_DAPC/dapc_invermay_results.tsv',
        dapc_ruakura_results = 'output/071_DAPC/dapc_ruakura_results.tsv',
        dapc_rpoa_results = 'output/071_DAPC/dapc_rpoa_results.tsv',
        dapc_lincoln_results = 'output/071_DAPC/dapc_lincoln_results.tsv'
    singularity:
        r_container
    log:
        'output/logs/071_DAPC/run_DAPC.log'
    script:
        'src/DAPC.R'


# Convert ped to plink format, specify that chromosomes are unknown so that it behaves
rule convert_to_plinkraw:
    input:
        ped = 'output/060_pop_genet/snps.ped',
        map = 'output/060_pop_genet/snps.map'
    output:
        'output/060_pop_genet/plink.raw'
    params:
        workdir = 'output/060_pop_genet'
    singularity:
        plink_container
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
    singularity:
        r_container
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
    singularity:
        r_container
    threads:
        25
    log:
        'output/logs/060_pop_genet/gds_convert.log'
    script:
        'src/convert_gds.R'















