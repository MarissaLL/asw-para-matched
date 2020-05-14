#!/usr/bin/env python3

###########
# GLOBALS #
###########

bayescan_container = 'shub://MarissaLL/singularity-containers:bayescan_2.1'
pgdspider_container = 'shub://MarissaLL/singularity-containers:pgdspider_2.1.1.5'
plink_container = 'shub://MarissaLL/singularity-containers:plink_1.9'


dataset = ['geo_final', 'geo_ns', 'para']

# , 'geo_ns'
#########
# RULES #
#########


rule target:
    input:
       expand('output/090_geo_survey/{dataset_dir}/geo_survey.geste-outputformat',
                dataset_dir = dataset)
       

# GO according to approach of Laporte et al 2016 (Lamprey paper)
# GO for SNPs within genes exon or interior intron
# Go for SNPs within 5000bp of a gene



#### Unfinished
# rule list_sig_snps:
#     input:
#         fst_file = 'output/090_geo_survey/{dataset}/geo_survey_fst.txt'
#         vcf = 'data/{dataset}/pop.vcf'
#     output:
#            'sig_snp_list.txt'
#     params:
#         info = '{dataset}'
#     script:
#         'plot_and_export_sig.R'


rule bayescan:
    input:
        genotypes = 'output/090_geo_survey/{dataset}/geo_survey.geste-outputformat'
    output:
        'output/090_geo_survey/{dataset}/geo_survey.sel',
        'output/090_geo_survey/{dataset}/geo_survey_fst.txt'
    params:
        outdir = 'output/090_geo_survey/{dataset}/',
        outname = 'geo_survey'
    singularity:
        bayescan_container
    threads:
        20
    log:
        'output/logs/090_geo_survey/{dataset}_bayescan.log'
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

#Convert SNP data in VCF format into geste format for bayescan.
# Note that the location of the popmap is actually specified within the spid file
# It is only specified here to link the dependencies of the rules
rule convert_bayescan_input: 
    input:
        vcf = 'data/{dataset}/pop.vcf',
        # vcf = 'output/090_geo_survey/{dataset}/subset_ldpruned.vcf',
        spid = 'data/{dataset}/geo_survey.spid',
        popmap = 'data/{dataset}/popmap_subset.txt'
    output:
        'output/090_geo_survey/{dataset}/geo_survey.geste-outputformat'
    params:
        in_format = 'VCF',
        out_format = 'GESTE_BAYE_SCAN',
        out_path = 'output/090_geo_survey/{dataset}/geo_survey.geste'
    singularity:
        pgdspider_container
    log:
        'output/logs/090_geo_survey/pgdspider_{dataset}.log'
    shell:
        'java -jar /opt/pgdspider/PGDSpider2-cli.jar '
        '-inputfile {input.vcf} '
        '-inputformat {params.in_format} '
        '-outputfile {params.out_path}'
        '-outputformat {params.out_format} '
        '-spid {input.spid} '
        '&> {log}'


rule filter_vcf:
    input:
        vcf = 'data/{dataset}/populations.snps.vcf',
        to_filter = 'output/090_geo_survey/{dataset}/geo_survey_sorted.prune.in'
    output:
         'output/090_geo_survey/{dataset}/subset_ldpruned.vcf'
    shell:
         "bcftools filter -e 'ID!=@{input.to_filter}' {input.vcf} > {output}"



rule thin_snps_12:
    input:
        sorted_bed = 'output/090_geo_survey/{dataset}/geo_survey_sorted.bed'
    output:
        'output/090_geo_survey/{dataset}/geo_survey_sorted.prune.in'
    params:
        filename_prefix = 'output/090_geo_survey/{dataset}/geo_survey_sorted',
        outfile = 'output/090_geo_survey/{dataset}/geo_survey_sorted'
    log:
        'output/logs/090_geo_survey/plink_prune_{dataset}.log'
    shell:
        'plink --bfile {params.filename_prefix} --aec --indep 50 5 1.2 '
        '--out {params.outfile}'

rule sort_bed:
    input:
        bed = 'output/090_geo_survey/{dataset}/geo_survey.bed'
    output:
        'output/090_geo_survey/{dataset}/geo_survey_sorted.bed'
    params:
        filename_prefix = 'output/090_geo_survey/{dataset}/geo_survey',
        outfile = 'output/090_geo_survey/{dataset}/geo_survey_sorted'
    shell:
        'plink --bfile {params.filename_prefix} --make-bed --aec --out {params.outfile}'
        


rule make_bed_file:
    input:
        vcf = 'data/{dataset}/pop.vcf'
    output:
        'output/090_geo_survey/{dataset}/geo_survey.bed'
    params:
        outfile = 'output/090_geo_survey/{dataset}/geo_survey'
    shell:
        'plink --vcf {input.vcf} --make-bed --aec --out {params.outfile}'


        
