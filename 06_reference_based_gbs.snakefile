#!/usr/bin/env python3

import numpy
import pandas
import pickle

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

counts_file = 'output/010_config/old_filtering_stats.csv'

bbmap_container = 'shub://MarissaLL/singularity-containers:bbmap_37.92@fa973c37883055c243c5e37f82f68f4d'
bwa_container = 'shub://TomHarrop/singularity-containers:bwa_0.7.17'
samtools_container = 'shub://MarissaLL/singularity-containers:samtools_1.9'
stacks2b_container = 'shub://TomHarrop/singularity-containers:stacks_2.0b@099f0c7d8c8ff2baf7ad763ad7bcd17b'

#########
# SETUP #
#########

ustacks_individuals = get_ustacks_individuals(counts_file)
test_individuals = ['I16', 'I21']


#########
# RULES #
#########

subworkflow process_reads:
    snakefile: '01_process_reads.snakefile'

map_method = ['bbmap']


rule target:
    input:
     'output/081_genome_mapped_stacks/catalog.fa.gz',
        'output/081_genome_mapped_stacks/catalog.calls'
        # # expand('output/080_against_genome/{individual}_bwa.sam',
        # #         individual=ustacks_individuals),
        # expand('output/081_genome_mapped_stacks/{individual}_bbmap.bam',
        #         individual=ustacks_individuals),
        # # 'output/081_genome_mapped_stacks/catalog.fa.gz',
        # 'output/logs/081_genome_mapped_stacks/collated_read_bbmap_stats.log'
        #        # # expand('output/081_genome_mapped_stacks/{individual}_bwa.snps.tsv.gz',
        # #         individual=ustacks_individuals)





rule run_gstacks:
    input:
        bam_files = expand('output/081_genome_mapped_stacks/{individual}_bbmap.bam',
                            individual=ustacks_individuals),
       # popmap = 'output/081_genome_mapped_stacks/testing_popmap.txt'
        popmap = 'output/010_config/filtered_popmap.txt'
    output:
        'output/081_genome_mapped_stacks/catalog.fa.gz',
        'output/081_genome_mapped_stacks/catalog.calls'
    params:
        stacks_dir = 'output/081_genome_mapped_stacks/'
    threads:
        18
    log:
        'output/logs/081_genome_mapped_stacks/gstacks_bbmapped_reads.log'
    singularity:
        stacks2b_container
    shell:
        'gstacks '
        '-I {params.stacks_dir} '
        '-S _bbmap.bam '
        '-O {params.stacks_dir} '
        '-M {input.popmap} '
        '--max-clipped 0.5 '
        '--min-mapq 1 '
        '--details '
        '--phasing-dont-prune-hets '
        '--unpaired '
        '-t {threads} '
        '2> {log}'



# Sort the sam file by locus name. Can this output a bam?
rule sort_sam:
    input:
        sam = 'output/081_genome_mapped_stacks/{individual}_bbmap.sam',
        ref = 'data/flye_denovo_full.racon.fasta'
    output:
        'output/081_genome_mapped_stacks/{individual}_bbmap.bam'
    log:
        'output/logs/081_genome_mapped_stacks/samtools_sort_{individual}.log'
    singularity:
        samtools_container
    shell:
        'samtools sort '
        '{input.sam} '
        '-o {output} '
        '-O BAM '
        '--reference {input.ref} '
        '&> {log}'

# # Map GBS reads to the genome using BWA
# rule map_tidied_reads_bwa:
#     input:
#         tidy_reads = 'output/021_filtered/{individual}.fq.gz',
#         index = expand('output/080_against_genome/flye_denovo_full.racon.fasta.{suffix}',
#                suffix=['amb', 'ann', 'bwt', 'pac', 'sa'])
#     output:
#         'output/081_genome_mapped_stacks/{individual}_bwa.sam'
#     params:
#         prefix = 'output/080_against_genome/flye_denovo_full.racon.fasta',
#     log:
#         'output/logs/081_genome_mapped_stacks/bwa_map_reads_{individual}.log'
#     singularity:
#         bwa_container
#     shell:
#         'bwa mem '
#         '-t 8 '
#         '{params.prefix} '
#         '-B 4 ' # with -B 1 I get 18 alignments, with 2 I get 4. default is 4 I think
#         '{input.tidy_reads} '
#         '1> {output} '
#         '2> {log}'


rule get_bbmap_mapping_stats:
    input:
        expand('output/logs/081_genome_mapped_stacks/bbmap_map_reads_{individual}.log',
                individual=ustacks_individuals)
    output:
        'output/logs/081_genome_mapped_stacks/collated_read_bbmap_stats.log' 
    shell:
        'grep '
        '--with-filename '
        '"unambiguous" '
        '-C 1 '
        '{input} > {output}'       


# Map GBS reads to the genome using bbmap
rule map_tidied_reads_bbmap:
    input:
        genome = 'data/flye_denovo_full.racon.fasta',
        tidy_reads = 'output/021_filtered/{individual}.fq.gz'
    output:
        'output/081_genome_mapped_stacks/{individual}_bbmap.sam'
    params:
        ref_path = 'output/081_genome_mapped_stacks/'
    log:
        'output/logs/081_genome_mapped_stacks/bbmap_map_reads_{individual}.log'
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
    


# # Index the genome for BWA
# rule bwa_index:
#     input:
#         'data/flye_denovo_full.racon.fasta'
#     output:
#         expand('output/080_against_genome/flye_denovo_full.racon.fasta.{suffix}',
#                suffix=['amb', 'ann', 'bwt', 'pac', 'sa'])
#     params:
#         prefix = 'output/080_against_genome/flye_denovo_full.racon.fasta',
#         algorithm = 'is'
#     threads:
#         30
#     log:
#         'output/logs/080_against_genome/bwa_index.log'
#     singularity:
#         bwa_container
#     shell:
#         'bwa index '
#         '-p {params.prefix} '
#         '-a {params.algorithm} '
#         '{input} '
#         '&> {log}'
