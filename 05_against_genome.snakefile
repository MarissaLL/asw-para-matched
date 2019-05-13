#!/usr/bin/env python3

import multiprocessing

###########
# GLOBALS #
###########

bwa_container = 'shub://TomHarrop/singularity-containers:bwa_0.7.17'
stacks_container = 'shub://TomHarrop/singularity-containers:stacks_2.3e'
samtools_container = 'shub://TomHarrop/singularity-containers:samtools_1.9'


#########
# RULES #
#########

rule target:
    input:
        # 'output/080_against_genome/locus_coordinates.tsv',
        expand('output/080_against_genome/flye_denovo_full.racon.fasta.{suffix}',
               suffix=['amb', 'ann', 'bwt', 'pac', 'sa']),
        'output/080_against_genome/aln.sam',
        'output/080_against_genome/aln_mapped_1.sam'


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



# # Integrate alignment info back into stacks workflow
# rule integrate_alignments:
#     input:
#         catalog = 'output/040_stacks/populations.loci.fa',
#         bam = 'output/080_against_genome/aln.bam'
#     output:
#         aln_catalog = 'output/080_against_genome/catalog.fa.gz',
#         tsv = 'output/080_against_genome/locus_coordinates.tsv',
#         calls = 'output/080_against_genome/catalog.calls'
#     params:
#         stacks_dir = 'output/040_stacks',
#         out_dir = 'output/080_against_genome'
#     threads:
#         1
#     log:
#         'output/logs/080_against_genome/integrate_alignments.log'
#     singularity:
#         stacks_container
#     shell:
#         'stacks-integrate-alignments '
#         '-P {params.stacks_dir} '
#         '-B {input.bam} '
#         '-O {params.out_dir} '
#         '&> {log}' 

# # Format conversion
# rule sam_to_bam:
#     input:
#         aln = 'output/080_against_genome/aln.sam'
#     output:
#         bam = 'output/080_against_genome/aln.bam'
#     threads:
#         1
#     log:
#         'output/logs/080_against_genome/sam_to_bam.log'
#     singularity:
#         stacks_container
#     shell:
#         'samtools view '
#         '--threads {threads} '
#         '-O BAM '
#         '-bh '
#         '{input.aln} '
#         '1> {output.bam} '
#         '2> {log}'


rule filter_unaligned:
    input:
        'output/080_against_genome/aln.sam'
    output:
        'output/080_against_genome/aln_mapped_1.sam'
    log:
        'output/logs/080_against_genome/filter_unaligned.log'
    singularity:
        samtools_container
    shell:
        'samtools view '
        '-F4 {input} '
        '1> {output} '
        '2> {log}'


# Map the full stacks catalog to the genome
rule map_filtered_catalog:
    input:
        catalog = 'output/060_pop_genet/populations.loci.fa',
        index = expand('output/080_against_genome/flye_denovo_full.racon.fasta.{suffix}',
               suffix=['amb', 'ann', 'bwt', 'pac', 'sa'])
    output:
        'output/080_against_genome/aln.sam'
    params:
        prefix = 'output/080_against_genome/flye_denovo_full.racon.fasta',
    threads:
        30
    log:
        'output/logs/080_against_genome/map_filtered_catalog.log'
    singularity:
        bwa_container
    shell:
        'bwa mem '
        '-t {threads} '
        '{params.prefix} '
        '-B 1 ' # with -B 1 I get 18 alignments, with 2 I get 4
        '{input.catalog} '
        '1> {output} '
        '2> {log}'

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