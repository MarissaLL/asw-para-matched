#!/usr/bin/env python3

import numpy
import os
import pandas
import re

#############
# FUNCTIONS #
#############


def get_ustacks_individuals(counts_file):
	counts_data = pandas.read_csv(counts_file)
	indivs = sorted(set(counts_data.loc[counts_data['#Kept'] >1e6]['#Individual']))
	passed_read_filter = counts_data.loc[counts_data['#Individual'].isin(indivs)]
	q90 = numpy.percentile(passed_read_filter['mean_gc'], 90)
	passed_all_filters = sorted(set(passed_read_filter.loc[passed_read_filter['mean_gc'] < q90]['#Individual']))
	return(passed_all_filters)

###########
# GLOBALS #
###########

counts_file = 'output/010_config/filtering_stats.csv'
stacks_container = 'shub://TomHarrop/singularity-containers:stacks_2.0b@099f0c7d8c8ff2baf7ad763ad7bcd17b'

#########
# SETUP #
#########

ustacks_individuals = get_ustacks_individuals(counts_file)

#########
# RULES #
#########

subworkflow process_reads:
    snakefile: 'process_reads.snakefile'


rule target:
    input:
        'output/060_pop_genet/plink.raw'

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
        'output/050_stacks_pops/r0/populations.snps.vcf'
    output:
        'output/060_pop_genet/snps.gds'
    threads:
        25
    log:
        'output/logs/060_pop_genet/gds_convert.log'
    script:
        'src/convert_gds.R'


rule populations:
    input:
        'output/040_stacks/catalog.fa.gz',
        'output/040_stacks/catalog.calls',
        popmap = 'output/010_config/filtered_popmap.txt'
    output:
        'output/050_stacks_pops/r0/populations.snps.vcf'        
    params:
        stacks_dir = 'output/040_stacks',
        outdir = 'output/050_stacks_pops/r0'
    singularity:
        'shub://TomHarrop/singularity-containers:stacks_2.0beta9@bb2f9183318871f6228b51104056a2d0'
    threads:
        50
    log:
        'output/logs/050_stacks_pops/pops.log'
    shell:
        'populations '
        '-P {params.stacks_dir} '
        '-M {input.popmap} '
        '-O {params.outdir} '
        '-t {threads} '
        '-r 0 '
        '--genepop --vcf '
        '&> {log}'

rule gstacks:
    input:
        expand('output/040_stacks/{individual}.matches.bam',
                individual=ustacks_individuals),
        popmap = 'output/010_config/filtered_popmap.txt'
    output:
        'output/040_stacks/catalog.fa.gz',
        'output/040_stacks/catalog.calls'
    params:
        stacks_dir = 'output/040_stacks'
    singularity:
        stacks_container
    threads:
        75
    log:
        'output/logs/040_stacks/gstacks.log'
    shell:
        'gstacks '
        '-P {params.stacks_dir} '
        '-M {input.popmap} '
        '-t {threads} '
        '&> {log}'

rule tsv2bam:
    input:
        expand('output/040_stacks/{individual}.matches.tsv.gz',
                individual=ustacks_individuals),
        popmap = 'output/010_config/filtered_popmap.txt'
    output:
        expand('output/040_stacks/{individual}.matches.bam',
                individual=ustacks_individuals)
    params:
        stacks_dir = 'output/040_stacks'
    singularity:
        stacks_container
    threads:
        75
    log:
        'output/logs/040_stacks/tsv2bam.log'
    shell:
        'tsv2bam '
        '-P {params.stacks_dir} '
        '-M {input.popmap} '
        '-t {threads} '
        '&> {log}'

rule sstacks:
    input:
        catalog = 'output/040_stacks/catalog.tags.tsv.gz',
        popmap = 'output/010_config/filtered_popmap.txt'
    output:
        expand('output/040_stacks/{individual}.matches.tsv.gz',
                individual=ustacks_individuals)
    singularity:
        stacks_container
    threads:
        75
    params:
        stacks_dir = 'output/040_stacks'
    log:
        'output/logs/040_stacks/sstacks.log'
    shell:
        'sstacks '
        '-P {params.stacks_dir} '
        '-M {input.popmap} '
        '-p {threads} '
        '&> {log}'

rule cstacks:
    input:
        expand('output/040_stacks/{individual}.alleles.tsv.gz',
               individual=ustacks_individuals),
        popmap = 'output/010_config/filtered_popmap.txt'
    output:
        'output/040_stacks/catalog.tags.tsv.gz',
        'output/040_stacks/catalog.snps.tsv.gz',
        'output/040_stacks/catalog.alleles.tsv.gz'
    threads:
        75
    params:
        stacks_dir = 'output/040_stacks'
    log:
        'output/logs/040_stacks/cstacks.log'
    shell:
        'cstacks '
        '-p {threads} '
        '-P {params.stacks_dir} '
        '-M {input.popmap} '
        '-n 5 '
        '&> {log}'

rule combine_coverage:
	input:
		coverage_file = expand('output/041_ustacks_coverage/{individual}.csv',
			individual=ustacks_individuals),
		full_popmap = process_reads('output/010_config/full_popmap.txt')
	output:
		coverage_file = 'output/010_config/combined_coverage_ustacks.csv',
		filtered_popmap = 'output/010_config/filtered_popmap.txt'
	log:
		'output/logs/041_ustacks_coverage/combine_coverage.log'
	script:
		'src/combine_coverage.R'

rule calculate_coverage:
	input:
		tags = 'output/040_stacks/{individual}.tags.tsv.gz'
	output:
		csv = 'output/041_ustacks_coverage/{individual}.csv'
	log:
		'output/logs/041_ustacks_coverage/{individual}.log'
	script:
		'src/calculate_coverage.R'

rule ustacks:
    input:
        fq = process_reads('output/021_filtered/{individual}.fq.gz'),
        pickle = 'output/010_config/individual_i.p'
    output:
        'output/040_stacks/{individual}.alleles.tsv.gz',
        'output/040_stacks/{individual}.snps.tsv.gz',
        'output/040_stacks/{individual}.tags.tsv.gz'
    threads:
        1
    params:
        wd = 'output/040_stacks',
        m = '3',
        M = '4'
    log:
        'output/logs/040_stacks/{individual}_ustacks.log'
    benchmark:
        'output/benchmarks/040_stacks/{individual}_ustacks.log'
    run:
        with open(input.pickle, 'rb') as f:
            individual_i = pickle.load(f)
        sample_i = individual_i[wildcards.individual]
        shell('ustacks '
              '-p {threads} '
              '-t gzfastq '
              '-f {input.fq} '
              '-o {params.wd} '
              '-i {sample_i} '
              '-m {params.m} '
              '-M {params.M} '
              '&> {log} ',
              bench_record=bench_record)

rule enumerate_samples:
    input:
        process_reads(expand('output/021_filtered/{individual}.fq.gz',
               individual=ustacks_individuals))
    output:
        pickle = 'output/010_config/individual_i.p'
    run:
        my_files = [re.sub('\.fq\.gz$', '', os.path.basename(x))
                    for x in input]
        my_individuals = enumerate(sorted(set(my_files)))
        individual_i = {y: x for x, y in my_individuals}
        # pickle the individual_i dict for other rules to use
        with open(output.pickle, 'wb+') as f:
            pickle.dump(individual_i, f)