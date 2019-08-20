#!/usr/bin/env python3

import multiprocessing

###########
# GLOBALS #
###########

bwa_container = 'shub://TomHarrop/singularity-containers:bwa_0.7.17'
stacks_container = 'shub://TomHarrop/singularity-containers:stacks_2.3e'
stacks2beta_container = 'shub://TomHarrop/singularity-containers:stacks_2.0beta9@bb2f9183318871f6228b51104056a2d0'
samtools_container = 'shub://TomHarrop/singularity-containers:samtools_1.9'
bbmap_container = 'shub://TomHarrop/singularity-containers:bbmap_38.45'


#########
# RULES #
#########

rule target:
    input:
        'output/080_against_genome/sig_regions.fa',

        # 'output/080_against_genome/locus_coordinates.tsv',
        expand('output/080_against_genome/flye_denovo_full.racon.fasta.{suffix}',
               suffix=['amb', 'ann', 'bwt', 'pac', 'sa']),
        'output/080_against_genome/aln.sam',
        'output/080_against_genome/aln_mapped_1.sam',
        # 'output/080_against_genome/ref/genome/1/summary.txt',
        # 'output/080_against_genome/bbmapped.sam',
        'output/080_against_genome/loci_noheader.fa',
        # 'output/080_against_genome/bbmapped_full.sam',
        # 'output/080_against_genome/bbmapped_full.bam',
        'output/080_against_genome/locus_coordinates.tsv',
        # 'output/080_against_genome/bbmapped_full_sorted.bam',
        # 'output/080_against_genome/bbmapped_full_sorted_sam3.sam',
        'output/080_against_genome/bbmapped_checked_sorted.bam',
        'output/080_against_genome/populations.snps.vcf'
        
# grep -E -o ".{0,5}scaffold.{0,5}" populations.snps.vcf > scaffold_names.txt


# 
rule genome_population_stats:
    input:
        aln_catalog = 'output/080_against_genome/catalog.fa.gz',
        calls = 'output/080_against_genome/catalog.calls',
    output:
        'output/080_against_genome/populations.snps.vcf'
    params:
        stacks_dir = 'output/080_against_genome',
        outdir = 'output/080_against_genome'
    singularity:
        stacks2beta_container
    threads:
        50
    log:
        'output/logs/080_against_genome/genome_pop_stats.log'
    shell:
        'populations '
        '-P {params.stacks_dir} '
        '-O {params.outdir} '
        '-t {threads} '
        '-r 0 '
        '--genepop '
        '--plink '
        '--vcf '
        '--hwe '
        '--fstats '
        '--fasta_loci '
        '&> {log}'







# Integrate alignment info back into stacks workflow
rule integrate_alignments:
    input:
        catalog = 'output/040_stacks/catalog.fa.gz',
        bam = 'output/080_against_genome/bbmapped_checked_sorted.bam'
    output:
        aln_catalog = 'output/080_against_genome/catalog.fa.gz',
        tsv = 'output/080_against_genome/locus_coordinates.tsv',
        calls = 'output/080_against_genome/catalog.calls'
    params:
        stacks_dir = 'output/040_stacks',
        out_dir = 'output/080_against_genome'
    threads:
        1
    log:
        'output/logs/080_against_genome/integrate_alignments.log'
    singularity:
        stacks_container
    shell:
        'stacks-integrate-alignments '
        '-P {params.stacks_dir} '
        '-B {input.bam} '
        '-O {params.out_dir} '
        '&> {log}' 

# Format conversion
rule sam_to_bam:
    input:
        aln = 'output/080_against_genome/bbmapped_checked_sorted.sam'
    output:
        bam = 'output/080_against_genome/bbmapped_checked_sorted.bam'
    threads:
        1
    log:
        'output/logs/080_against_genome/sam_to_bam.log'
    singularity:
        stacks_container
    shell:
        'samtools view '
        '--threads {threads} '
        '-O BAM '
        '-bh '
        '{input.aln} '
        '1> {output.bam} '
        '2> {log}'


# Sort the sam file by locus name. Can this output a bam?
rule sort_sam:
    input:
        sam = 'output/080_against_genome/bbmapped_checked_14.sam',
        ref = 'data/flye_denovo_full.racon.fasta'
    output:
        'output/080_against_genome/bbmapped_checked_sorted.sam'
    log:
        'output/logs/080_against_genome/samtools_sort.log'
    shell:
        'samtools sort '
        '{input.sam} '
        '-o {output} '
        '-n '
        '--reference {input.ref} '
        '&> {log}'


# Should be possible to combine this with the previous rule using outm and sam=1.3
rule check_sam:
    input:
        'output/080_against_genome/bbmapped_full.sam'
    output:
        'output/080_against_genome/bbmapped_checked_14.sam'
    log:
        'output/logs/080_against_genome/check_sam.log'
    singularity:
        bbmap_container
    shell:
        'reformat.sh '
        'in={input} '
        'out={output} '
        'sam=1.3 '
        'mappedonly=t '
        '&> {log}'


# Indexes the genome, then maps the provided loci to it. 
# ref_path defines where the reference created by indexing goes
rule bbmap_full:
    input:
        genome = 'data/flye_denovo_full.racon.fasta',
        loci = 'output/040_stacks/catalog.fa.gz'
    output:
        'output/080_against_genome/bbmapped_full.sam'
    params:
        ref_path = 'output/080_against_genome/'
    log:
        'output/logs/080_against_genome/bbmap_index_map.log'
    singularity:
        bbmap_container
    shell:
        'bbmap.sh '
        'in={input.loci} '
        'out={output} '
        'ref={input.genome} '
        'path={params.ref_path} '
        'trimreaddescriptions=t '
        '&> {log}'


# Remove the comment line from the filtered catalog file so that bbmap can read it later
rule remove_header:
    input:
        loci = 'output/060_pop_genet/populations.loci.fa'
    output: 
        loci_noheader = 'output/080_against_genome/loci_noheader.fa'
    shell:
        'cp {input.loci} {output} && '
        'sed -i \'1d\' {output}'
    

# Map the full stacks catalog to the genome with bwa mem. Seems to only work on a subset of the data.
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