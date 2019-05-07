#!/usr/bin/env python3

import multiprocessing

###########
# GLOBALS #
###########

bwa_container = 'shub://TomHarrop/singularity-containers:bwa_0.7.17'

#########
# RULES #
#########

rule target:
    input:
        'output/080_against_genome/sig_regions.fa',
        expand('output/080_against_genome/flye_denovo_full.racon.fasta.{suffix}',
               suffix=['amb', 'ann', 'bwt', 'pac', 'sa']),
        temp('output/080_against_genome/aln.sam')


# Find loci that differed significantly between north and south island populations in the GBS SNPs 
# Done without aligning anything to the genome or LD-pruning
rule find_seqs:
    input:
        sig_loci = 'output/040_stacks/loc_num_file.txt',
        catalog = 'output/040_stacks/catalog.fa.gz'
    output:
        'output/080_against_genome/sig_regions.fa' 
    shell:
        'zcat {input.catalog} | '
        'grep '
        '-f {input.sig_loci} '
        '-A 1 '
        '-F '
        '--no-group-separator '
        '&> {output}'


# Map the full stacks catalog to the genome
rule map_full_catalog:
    input:
        catalog = 'output/040_stacks/catalog.fa.gz',
        index = expand('output/080_against_genome/flye_denovo_full.racon.fasta.{suffix}',
               suffix=['amb', 'ann', 'bwt', 'pac', 'sa'])
    output:
        temp('output/080_against_genome/aln.sam')
    params:
        prefix = 'output/080_against_genome/flye_denovo_full.racon.fasta',
    threads:
        multiprocessing.cpu_count()
    log:
        'output/logs/080_against_genome/map_full_catalog.log'
    singularity:
        bwa_container
    shell:
        'bwa mem '
        '-t {threads} '
        '{params.prefix} '
        '{input.catalog} '
        '> {output} '
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
