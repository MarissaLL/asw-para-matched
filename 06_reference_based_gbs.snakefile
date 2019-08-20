#!/usr/bin/env python3

import numpy
import os
import pandas

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

bbmap_container = 'shub://MarissaLL/singularity-containers:bbmap_37.92@fa973c37883055c243c5e37f82f68f4d'
bwa_container = 'shub://TomHarrop/singularity-containers:bwa_0.7.17'

#########
# SETUP #
#########

ustacks_individuals = get_ustacks_individuals(counts_file)

#########
# RULES #
#########

subworkflow process_reads:
    snakefile: '01_process_reads.snakefile'


rule target:
    input:
        expand('output/080_against_genome/{individual}_bwa.sam',
                individual=ustacks_individuals)


# Map GBS reads to the genome using BWA
rule map_tidied_reads_bwa:
    input:
        tidy_reads = 'output/021_filtered/{individual}.fq.gz',
        index = expand('output/080_against_genome/flye_denovo_full.racon.fasta.{suffix}',
               suffix=['amb', 'ann', 'bwt', 'pac', 'sa'])
    output:
        'output/080_against_genome/{individual}_bwa.sam'
    params:
        prefix = 'output/080_against_genome/flye_denovo_full.racon.fasta',
    threads:
        30
    log:
        'output/logs/080_against_genome/bwa_map_reads_{individual}.log'
    singularity:
        bwa_container
    shell:
        'bwa mem '
        '-t {threads} '
        '{params.prefix} '
        '-B 1 ' # with -B 1 I get 18 alignments, with 2 I get 4
        '{input.tidy_reads} '
        '1> {output} '
        '2> {log}'

# Map GBS reads to the genome using bbmap
rule map_tidied_reads_bbmap:
    input:
        genome = 'data/flye_denovo_full.racon.fasta',
        tidy_reads = process_reads('output/021_filtered/{individual}.fq.gz')
    output:
        'output/080_against_genome/{individual}_bbmap.sam'
    params:
        ref_path = 'output/080_against_genome/'
    log:
        'output/logs/080_against_genome/bbmap_map_reads_{individual}.log'
    singularity:
        bbmap_container
    shell:
        'bbmap.sh '
        'in={input.tidy_reads} '
        'out={output} '
        'ref={input.genome} '
        'path={params.ref_path} '
        'trimreaddescriptions=t '
        '&> {log}'
    


# Index the genome for BWA
rule bwa_index:
    input:
        'data/flye_denovo_full.racon.fasta'
    output:
        expand('output/080_against_genome/flye_denovo_full.racon.fasta.{suffix}',
               suffix=['amb', 'ann', 'bwt', 'pac', 'sa'])
    params:
        prefix = 'output/080_against_genome/flye_denovo_full.racon.fasta',
        algorithm = 'is'
    threads:
        1
    log:
        'output/logs/080_against_genome/bwa_index.log'
    singularity:
        bwa_container
    shell:
        'bwa index '
        '-p {params.prefix} '
        '-a {params.algorithm} '
        '{input} '
        '&> {log}'
