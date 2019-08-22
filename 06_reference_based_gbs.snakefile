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

def lookup_indiv(pickle_file, individual):
    with open(pickle_file, 'rb') as f:
        individual_i = pickle.load(f)
        sample_i = individual_i[individual]
        return(sample_i)

###########
# GLOBALS #
###########

counts_file = 'output/010_config/filtering_stats.csv'

bbmap_container = 'shub://MarissaLL/singularity-containers:bbmap_37.92@fa973c37883055c243c5e37f82f68f4d'
bwa_container = 'shub://TomHarrop/singularity-containers:bwa_0.7.17'
stacks2b_container = 'shub://TomHarrop/singularity-containers:stacks_2.0b@099f0c7d8c8ff2baf7ad763ad7bcd17b'

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
                individual=ustacks_individuals),
        expand('output/080_against_genome/{individual}_bbmap.sam',
                individual=ustacks_individuals),
        expand('output/081_genome_mapped_stacks/{individual}_bwa.snps.tsv.gz',
                individual=ustacks_individuals)

rule run_sstacks:
    input:
        catalog = 'output/081_genome_mapped_stacks/catalog.tags.tsv.gz',
        popmap = 'output/010_config/filtered_popmap.txt'
    output:
        expand('output/081_genome_mapped_stacks/{individual}.matches.tsv.gz',
                individual=ustacks_individuals)
    singularity:
        stacks2b_container
    threads:
        75
    params:
        stacks_dir = 'output/081_genome_mapped_stacks'
    log:
        'output/logs/081_genome_mapped_stacks/sstacks.log'
    shell:
        'sstacks '
        '-P {params.stacks_dir} '
        '-M {input.popmap} '
        '-p {threads} '
        '&> {log}'


# b — database/batch ID of the catalog to consider (default: guess).
# P — path to the directory containing Stacks files.
# M — path to a population map file from which to take sample names.
# s — filename prefix from which to load sample loci.
# c — path to the catalog.
# g,--aligned — base matching on alignment position, not sequence identity.
# p — enable parallel execution with num_threads threads.
# o — output path to write results.
# x — don't verify haplotype of matching locus.




## Do I need to worry about the sample prefix -s? "s — sample prefix from which to load loci into the catalog."

# Generate a catalog
rule run_cstacks;
    input:
        'output/081_genome_mapped_stacks/{individual}_bwa.tags.tsv.gz',
        'output/081_genome_mapped_stacks/{individual}_bwa.models.tsv.gz', 
        'output/081_genome_mapped_stacks/{individual}_bwa.snps.tsv.gz', 
        'output/081_genome_mapped_stacks/{individual}_bwa.alleles.tsv.gz',
        popmap = 'output/010_config/filtered_popmap.txt'
    output:
        'output/081_genome_mapped_stacks/catalog.fa.gz'
    params:
        input_dir = 'output/081_genome_mapped_stacks/',
        output_dir = 'output/081_genome_mapped_stacks/',
        n = '5'
    threads:
            12
    log:
        'output/logs/081_genome_mapped_stacks/cstacks.log'
    singularity:
            stacks2b_container
        shell:
            'cstacks '
            '--aligned '
            '-P {params.input_dir} '
            '-M {input.popmap} '
            '-n {params.n} '
            '-p {threads} '
            '-s '
            '-o {params.outdir} '
            '&> {log}'


# Assess read mapping rates based on the pstacks and process radtags logs at this point. (number of primary 
# alignments/initial num reads from process_radtags; proportion discarded due to soft-clipping; avg per locus 
# coverage


###### Might need to consider if mapping rate is low
# --max_clipped [float] — alignments with more than this fraction of soft-clipped bases are discarded (default 15%).
# --min_mapq [int] — minimum required quality (default 10).
# --keep_sec_alns — keep secondary alignments (default: false, only keep primary alignments).



# Extract reference-mapped stacks
rule run_pstacks:
    input:
        'output/080_against_genome/{individual}_bwa.sam',
        pickle = 'output/010_config/individual_i.p'
    output:
        'output/081_genome_mapped_stacks/{individual}_bwa.tags.tsv.gz',
        'output/081_genome_mapped_stacks/{individual}_bwa.models.tsv.gz', 
        'output/081_genome_mapped_stacks/{individual}_bwa.snps.tsv.gz', 
        'output/081_genome_mapped_stacks/{individual}_bwa.alleles.tsv.gz'
    params:
        input_filetype = 'sam',
        outdir = 'output/081_genome_mapped_stacks/',
        i = lambda wildcards, input: lookup_indiv(input.pickle, wildcards.individual),
        m = '3'
    threads:
        12
    log:
        'output/logs/081_genome_mapped_stacks/{individual}_pstacks.log'
    singularity:
        stacks2b_container
    shell:
        'pstacks '
        '-f {input} '
        '-i {params.i} '
        '-m 3 '
        '-p {threads} '
        '-t {params.input_filetype} '
        '-o {params.outdir} '
        '&> {log}'




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
        tidy_reads = 'output/021_filtered/{individual}.fq.gz'
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
