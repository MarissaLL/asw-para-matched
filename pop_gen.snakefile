#!/usr/bin/env python3


###########
# GLOBALS #
###########

stacks2beta_container = 'shub://TomHarrop/singularity-containers:stacks_2.0beta9@bb2f9183318871f6228b51104056a2d0'
bayescan_container = 'shub://MarissaLL/singularity-containers:bayescan_2.1@e035d571b5e888e98e5b79f902c30388'

#########
# RULES #
#########

subworkflow stacks:
    snakefile: 'stacks.snakefile'

rule target:
    input:
        'output/060_pop_genet/populations.snps.vcf',
        'output/070_pop_tests/compare_para.geste',
        'output/070_pop_tests/compared_para_prop.txt'


rule bayescan_para:
    input:
        genotypes = 'output/070_pop_tests/compare_para.geste'
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
        '-pr_odds 1000 '
        '-out_pilot '
        '-out_freq '
        '&> {log}'

# rule make_bayescan_input_para:
#     input:
#         ped = 'output/070_pop_tests/pop_para_snps.ped',
#         spid = 'data/convert_para.spid'
#     output:
#         'output/070_pop_tests/compare_para.geste'
#     params:
#         in_format = 'PED'
#         out_format = 'GESTE_BAYE_SCAN'
#     singularity:
#         pgdspider_container
#     log:
#         'output/logs/070_pop_tests/pgdspider_2pops.log'
#     shell:
#         'somehow run PGDSpider2-cli.jar' # FIX THIS
#         '-inputfile {input.ped} '
#         '-inputformat {params.in_format} '
#         '-outputfile {output}'
#         '-outputformat {params.out_format} '
#         '-spid {input.spid} '
#         '&> {log}'

# rule make_bayescan_input_pop:
#     input:
#         ped = 'output/070_pop_tests/pop_para_snps.ped',
#         spid = 'data/convert_pops.spid'
#     output:
#         'output/070_pop_tests/compare_pops.geste'
#     params:
#         in_format = 'PED'
#         out_format = 'GESTE_BAYE_SCAN'
#     singularity:
#         pgdspider_container
#     log:
#         'output/logs/070_pop_tests/pgdspider_4pops.log'
#     shell:
#         'somehow run PGDSpider2-cli.jar' # FIX THIS
#         '-inputfile {input.ped} '
#         '-inputformat {params.in_format} '
#         '-outputfile {output}'
#         '-outputformat {params.out_format} '
#         '-spid {input.spid} '
#         '&> {log}'

# rule add_popdata_to_ped:
#     input:
#         ped_file = 'output/060_pop_genet/populations.plink.ped'
#         para_info = process_reads('output/010_config/tidy_sample_info.tsv')
#     output:
#         'output/070_pop_tests/pop_para_snps.ped'
#     log:
#         'output/logs/070_pop_tests/add_popdata_to_ped.log'
#     script:
#         'src/add_popdata_to_ped.R'

# Run populations again on filtered data to get Fst etc.
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















