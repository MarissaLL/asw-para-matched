#!/usr/bin/env python3

import numpy
import os
import pandas
import pickle
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

def lookup_indiv(pickle_file, individual):
    with open(pickle_file, 'rb') as f:
        individual_i = pickle.load(f)
        sample_i = individual_i[individual]
        return(sample_i)

###########
# GLOBALS #
###########

counts_file = 'output/010_config/filtering_stats.csv'
r_container = 'shub://MarissaLL/singularity-containers:r_3.5.0@1f713cf93d67765d2d37eb87339495de5'
stacks2b_container = 'shub://TomHarrop/singularity-containers:stacks_2.0b@099f0c7d8c8ff2baf7ad763ad7bcd17b'
stacks2beta_container = 'shub://TomHarrop/singularity-containers:stacks_2.0beta9@bb2f9183318871f6228b51104056a2d0'

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
        'output/050_stacks_pops/r0/populations.snps.vcf'


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
        stacks2beta_container
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
        stacks2b_container
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
        stacks2b_container
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
        stacks2b_container
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
        stacks2b_dir = 'output/040_stacks',
        n = '5'
    singularity:
        stacks2b_container
    log:
        'output/logs/040_stacks/cstacks.log'
    shell:
        'cstacks '
        '-p {threads} '
        '-P {params.stacks2b_dir} '
        '-M {input.popmap} '
        '-n {params.n} '
        '&> {log}'

rule combine_coverage:
    input:
        coverage_file = expand('output/041_ustacks_coverage/{individual}.csv',
            individual=ustacks_individuals),
        full_popmap = process_reads('output/010_config/full_popmap.txt')
    output:
        coverage_file = 'output/010_config/combined_coverage_ustacks.csv',
        filtered_popmap = 'output/010_config/filtered_popmap.txt'
    singularity:
        r_container
    log:
        'output/logs/041_ustacks_coverage/combine_coverage.log'
    script:
        'src/combine_coverage.R'

# Calculate the number of reads and coverage within individuals
rule calculate_coverage:
    input:
        tags = 'output/040_stacks/{individual}.tags.tsv.gz'
    output:
        csv = 'output/041_ustacks_coverage/{individual}.csv'
    singularity:
        r_container
    log:
        'output/logs/041_ustacks_coverage/{individual}.log'
    script:
        'src/calculate_coverage.R'

# Form loci from reads within individuals
rule ustacks:
    input:
        fq = process_reads('output/021_filtered/{individual}.fq.gz'),
        pickle = 'output/010_config/individual_i.p'
    output:
        'output/040_stacks/{individual}.alleles.tsv.gz',
        'output/040_stacks/{individual}.snps.tsv.gz',
        'output/040_stacks/{individual}.tags.tsv.gz'
    singularity:
        stacks2b_container
    threads:
        1
    params:
        wd = 'output/040_stacks',
        m = '3',
        M = '4',
        i = lambda wildcards, input: lookup_indiv(input.pickle, wildcards.individual)
    log:
        'output/logs/040_stacks/{individual}_ustacks.log'
    benchmark:
        'output/benchmarks/040_stacks/{individual}_ustacks.log'
    shell:
        'ustacks '
        '-p {threads} '
        '-t gzfastq '
        '-f {input.fq} '
        '-o {params.wd} '
        '-i {params.i} '
        '-m {params.m} '
        '-M {params.M} '
        '&> {log}'

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