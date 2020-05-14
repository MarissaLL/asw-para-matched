#!/usr/bin/env python3

import pandas


###########
# GLOBALS #
###########

bayescan_container = 'shub://MarissaLL/singularity-containers:bayescan_2.1'
bbmap_container = 'shub://TomHarrop/singularity-containers:bbmap_38.45'
bwa_container = 'shub://TomHarrop/singularity-containers:bwa_0.7.17'
pgdspider_container = 'shub://MarissaLL/singularity-containers:pgdspider_2.1.1.5'
plink_container = 'shub://TomHarrop/singularity-containers:plink_1.90beta5'
r_container = 'shub://MarissaLL/singularity-containers:r_3.5.0'
samtools_container = 'shub://TomHarrop/singularity-containers:samtools_1.9'
stacks_container = 'shub://TomHarrop/singularity-containers:stacks_2.3e'
stacks2beta_container = 'shub://TomHarrop/singularity-containers:stacks_2.0beta9'
vcftools_container = 'shub://MarissaLL/singularity-containers:vcftools_0.1.17@230db32b3097775cd51432092f9cbcb1'


filtered_popmap = 'data/geo_survey_filtered_popmap.txt'

indivs = pandas.read_csv(filtered_popmap, sep = '\t', header=None)
filtered_individuals = indivs[0]

catalog_version = ['ref_based']
# , 'catalog_mapped'

#########
# RULES #
#########


rule target:
	input:
		'output/090_geo_survey/catalog_mapped/pre_filter/populations.snps.vcf',
		'output/090_geo_survey/ref_based/pre_filter/populations.snps.vcf'
		


#### Rules 01,02,03 and 06 are to map GBS reads then run stacks
#### Rules 04, 045, 05 and 07 are to map the exiting catalog then integrate back into stacks

# 07
rule run_populations_mapped:
    input:
        aln_catalog = 'output/090_geo_survey/catalog_mapped/catalog.fa.gz',
        calls = 'output/090_geo_survey/catalog_mapped/catalog.calls',
        popmap = 'data/geo_survey_filtered_popmap.txt'
    output:
        'output/090_geo_survey/catalog_mapped/pre_filter/populations.snps.vcf'
    params:
        stacks_dir = 'output/090_geo_survey/catalog_mapped/',
        outdir = 'output/090_geo_survey/catalog_mapped/pre_filter/'
    singularity:
        stacks2beta_container
    threads:
        50
    log:
        'output/logs/090_geo_survey/mapped_populations_pre_filter.log'
    shell:
        'populations '
        '-P {params.stacks_dir} '
        '-O {params.outdir} '
        '-M {input.popmap} '
        '-t {threads} '
        '-r 0 '
        '--vcf '
        '&> {log}'


# 06
rule run_populations_refbased:
    input:
        aln_catalog = 'output/090_geo_survey/ref_based/catalog.fa.gz',
        calls ='output/090_geo_survey/ref_based/catalog.calls',
        popmap = 'data/geo_survey_filtered_popmap.txt'
    output:
        'output/090_geo_survey/ref_based/pre_filter/populations.snps.vcf'
    params:
        stacks_dir = 'output/090_geo_survey/ref_based/',
        outdir = 'output/090_geo_survey/ref_based/pre_filter/'
    singularity:
        stacks2beta_container
    threads:
        50
    log:
        'output/logs/090_geo_survey/refbased_populations_pre_filter.log'
    shell:
        'populations '
        '-P {params.stacks_dir} '
        '-O {params.outdir} '
        '-M {input.popmap} '
        '-t {threads} '
        '-r 0 '
        '--vcf '
        '&> {log}'



# 05 Integrate alignment info back into stacks workflow. Using edited file because the original only outputs
# loci without a t in their name
rule integrate_alignments:
    input:
        catalog =  'data/geo_survey_catalog/catalog.fa.gz',
        calls = 'data/geo_survey_catalog/catalog.calls',
        bam = 'output/090_geo_survey/bbmapped_full.bam'
    output:
        aln_catalog = 'output/090_geo_survey/catalog_mapped/catalog.fa.gz',
        tsv = 'output/090_geo_survey/catalog_mapped/locus_coordinates.tsv',
        calls = 'output/090_geo_survey/catalog_mapped/catalog.calls'
    params:
        stacks_dir = 'data/geo_survey_catalog/',
        out_dir = 'output/090_geo_survey/catalog_mapped/'
    threads:
        1
    log:
        'output/logs/090_geo_survey/integrate_bwa_catalog_alignments.log'
    
    shell:
        ' ./stacks-integrate-alignments-edited '
        '-P {params.stacks_dir} '
        '-B {input.bam} '
        '-O {params.out_dir} '
        '&> {log}' 

#045
rule sam_to_bam:
    input:
        'output/090_geo_survey/bbmapped_full.sam'
    output:
        'output/090_geo_survey/bbmapped_full.bam'
    singularity:
        samtools_container
    shell:
        'samtools view -S -b {input} > {output}'

# 04
rule bbmap_full_catalog:
    input:
        genome = 'data/flye_denovo_full.racon.fasta',
        catalog = 'data/geo_survey_catalog/catalog.fa.gz'
    output:
        'output/090_geo_survey/bbmapped_full.sam'
    params:
        ref_path = 'output/081_genome_mapped_stacks/' # Keeping the old path to re-use existing index
    log:
        'output/logs/090_geo_survey/bbmap_index_map.log'
    singularity:
        bbmap_container
    shell:
        'bbmap.sh '
        'in={input.catalog} '
        'out={output} '
        'ref={input.genome} '
        'path={params.ref_path} '
        'trimreaddescriptions=t '
        '&> {log}'


# 03
rule run_gstacks:
    input:
        bam_files = expand('output/090_geo_survey/{individual}_bwa.bam',
                            individual=filtered_individuals),
        popmap = 'data/geo_survey_filtered_popmap.txt'
    output:
        'output/090_geo_survey/ref_based/catalog.fa.gz',
        'output/090_geo_survey/ref_based/catalog.calls'
    params:
        stacks_dir = 'output/090_geo_survey/',
        output_dir = 'output/090_geo_survey/ref_based',
        file_suffix = '_bwa.bam'
    threads:
        18
    log:
        'output/logs/090_geo_survey/gstacks_bwamapped_reads.log'
    singularity:
        stacks2beta_container
    shell:
        'gstacks '
        '-I {params.stacks_dir} '
        '-S {params.file_suffix} '
        '-O {params.output_dir} '
        '-M {input.popmap} '
        '--max-clipped 0.5 '
        '--min-mapq 1 '
        '--details '
        '--phasing-dont-prune-hets '
        '--unpaired '
        '-t {threads} '
        '&> {log}'



# 02 Sort the sam file by locus name, convert to bam. 
rule sort_sam:
    input:
        sam = 'output/090_geo_survey/{individual}_bwa.sam',
        ref = 'data/flye_denovo_full.racon.fasta'
    output:
        'output/090_geo_survey/{individual}_bwa.bam'
    log:
        'output/logs/090_geo_survey/samtools_sort_{individual}_bwa.log'
    singularity:
        samtools_container
    shell:
        'samtools sort '
        '{input.sam} '
        '-o {output} '
        '-O BAM '
        '--reference {input.ref} '
        '&> {log}'


# # 01 Map GBS reads to the genome using BWA
# # NOTE this step was run separately to this script (bwa difficulty) on the 143 fq.gz files (filtered indivs)
# # and the logs for all indivs are combined into bwa_mapping_samples.log
# rule map_tidied_reads_bwa:
#     input:
#         tidy_reads = '/Volumes/archive/deardenlab/tomharrop/projects/asw-stacks-singleplate/output/combined/{individual}.fq.gz',
#         index = expand('output/080_against_genome/flye_denovo_full.racon.fasta.{suffix}',
#                suffix=['amb', 'ann', 'bwt', 'pac', 'sa'])
#     output:
#         'output/090_geo_survey/{individual}_bwa.sam'
#     params:
#         prefix = 'output/090_geo_survey/flye_denovo_full.racon.fasta',
#     log:
#         'output/logs/090_geo_survey/bwa_map_reads_{individual}.log'
#     singularity:
#         bwa_container
#     shell:
#         'bwa mem '
#         '-t 8 '
#         '{params.prefix} '
#         '{input.tidy_reads} '
#         '1> {output} '
#         '2> {log}'