#!/usr/bin/env python3

import os
import pandas
import re

#############
# FUNCTIONS #
#############


def get_ustacks_individuals(counts_file): 
	counts_data = pandas.read_csv(counts_file)
	indivs = sorted(set(counts_data.loc[counts_data['#Kept'] >1e6]['#Individual']))
	return(indivs)

###########
# GLOBALS #
###########

counts_file = 'output/010_config/individual_counts.csv'

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
        'output/010_config/combined_coverage_ustacks.csv'

rule combine_coverage:
	input:
		coverage_file = expand('output/041_ustacks_coverage/{individual}.csv',
			individual=ustacks_individuals)
	output:
		coverage_file = 'output/010_config/combined_coverage_ustacks.csv'
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
    params:
        wd = 'output/040_stacks',
        m = '3',
        M = '4'
    output:
        'output/040_stacks/{individual}.alleles.tsv.gz',
        'output/040_stacks/{individual}.snps.tsv.gz',
        'output/040_stacks/{individual}.tags.tsv.gz'
    threads:
        1
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