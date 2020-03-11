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

bayescan_container = 'shub://MarissaLL/singularity-containers:bayescan_2.1'
bbmap_container = 'shub://MarissaLL/singularity-containers:bbmap_37.92@fa973c37883055c243c5e37f82f68f4d'
bwa_container = 'shub://TomHarrop/singularity-containers:bwa_0.7.17'
pgdspider_container = 'shub://MarissaLL/singularity-containers:pgdspider_2.1.1.5'
plink_container = 'shub://TomHarrop/singularity-containers:plink_1.90beta5'
r_container = 'shub://MarissaLL/singularity-containers:r_3.5.0'
samtools_container = 'shub://MarissaLL/singularity-containers:samtools_1.9'
stacks2b_container = 'shub://TomHarrop/singularity-containers:stacks_2.0b@099f0c7d8c8ff2baf7ad763ad7bcd17b'
stacks2beta_container = 'shub://TomHarrop/singularity-containers:stacks_2.0beta9'
vcftools_container = 'shub://MarissaLL/singularity-containers:vcftools_0.1.17'

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
bayescan_subsets = ['compared_island', 'compared_4pops', 'compared_para']


rule target:
    input:
        # expand('output/081_genome_mapped_stacks/{individual}_bwa.sam',
        #     individual=ustacks_individuals),
        'output/081_genome_mapped_stacks/catalog.fa.gz',
        'output/081_genome_mapped_stacks/populations.snps.vcf',
        'output/081_genome_mapped_stacks/whitelist.txt',
        'output/081_genome_mapped_stacks/final_popgen/populations.snps.vcf',
        # expand('trial/{bayescan_subset}.sel',
        #       bayescan_subset=bayescan_subsets)
        # 'output/081_genome_mapped_stacks/catalog.fa.gz',
        # 'output/081_genome_mapped_stacks/catalog.calls'
        # # expand('output/080_against_genome/{individual}_bwa.sam',
        # #         individual=ustacks_individuals),
        # expand('output/081_genome_mapped_stacks/{individual}_bbmap.bam',
        #         individual=ustacks_individuals),
        # # 'output/081_genome_mapped_stacks/catalog.fa.gz',
        # 'output/logs/081_genome_mapped_stacks/collated_read_bbmap_stats.log'
        #        # # expand('output/081_genome_mapped_stacks/{individual}_bwa.snps.tsv.gz',
        # #         individual=ustacks_individuals)




rule bayescan:
    input:
        genotypes = 'output/081_genome_mapped_stacks/final_popgen/{bayescan_subset}.geste-outputformat'
    output:
        'output/081_genome_mapped_stacks/final_popgen/{bayescan_subset}.sel',
        'output/081_genome_mapped_stacks/final_popgen/{bayescan_subset}_fst.txt'
    params:
        outdir = 'output/081_genome_mapped_stacks/final_popgen/',
        outname = '{bayescan_subset}'
    singularity:
        bayescan_container
    threads:
        6
    log:
        'output/logs/081_genome_mapped_stacks/{bayescan_subset}_bayescan_filtered_catalog.log'
    shell:
        'bayescan_2.1 '
        '{input.genotypes} '
        '-od {params.outdir} '
        '-o {params.outname} '
        '-pilot 5000 '
        '-nbp 20 '
        '-burn 15000 '
        '-n 30000 '
        '-thin 10 '
        '-pr_odds 500 '
        '-out_pilot '
        '-out_freq '
        '&> {log}'

        

# ## Note that the location of the popmap is actually specified within the spid file
# ## It is only specified here to link the dependencies of the rules
rule convert_subset_bayescan_inputs:
    input:
        vcf = 'output/081_genome_mapped_stacks/final_popgen/{bayescan_subset}.recode.vcf',
        spid = 'data/{bayescan_subset}.spid',
        popmap = 'output/070_bayescan/popmap_{bayescan_subset}.txt'
    output:
        'output/081_genome_mapped_stacks/final_popgen/{bayescan_subset}.geste-outputformat'
    params:
        in_format = 'VCF',
        out_format = 'GESTE_BAYE_SCAN',
        out_path = 'output/081_genome_mapped_stacks/final_popgen/{bayescan_subset}.geste'
    singularity:
        pgdspider_container
    threads:
        50
    log:
        'output/logs/081_genome_mapped_stacks/pgdspider_{bayescan_subset}_filtered_catalog.log'
    shell:
        'java -jar /opt/pgdspider/PGDSpider2-cli.jar '
        '-inputfile {input.vcf} '
        '-inputformat {params.in_format} '
        '-outputfile {params.out_path}'
        '-outputformat {params.out_format} '
        '-spid {input.spid} '
        '&> {log}'

# # Subset vcfs to only include individuals needed in each bayescan run
rule subset_vcfs:
    input:
        individual_list = 'output/070_bayescan/{bayescan_subset}_indivs.txt',
        full_vcf = 'output/081_genome_mapped_stacks/final_popgen/populations.snps.vcf'
    output: 
        subset_vcf = 'output/081_genome_mapped_stacks/final_popgen/{bayescan_subset}.recode.vcf'
    params:
        outname = 'output/081_genome_mapped_stacks/final_popgen/{bayescan_subset}'
    singularity:
        vcftools_container
    threads:
        25
    log:
        'output/logs/081_genome_mapped_stacks/{bayescan_subset}_subset_vcf.log'
    shell:
        'vcftools '
        '--vcf {input.full_vcf} '
        '--out {params.outname} '
        '--keep {input.individual_list} '
        '--recode '
        '&> {log}'



# Run populations again on filtered data to get population summary statistics. 
# Output a fasta file of the consensus sequences for the filtered loci
rule populations_stats:
    input:
        aln_catalog = 'output/081_genome_mapped_stacks/catalog.fa.gz',
        calls = 'output/081_genome_mapped_stacks/catalog.calls',
         whitelist = 'output/081_genome_mapped_stacks/whitelist.txt',
        popmap = 'output/081_genome_mapped_stacks/ref_based_popmap.txt'
    output:
        'output/081_genome_mapped_stacks/final_popgen/populations.snps.vcf',
        'output/081_genome_mapped_stacks/final_popgen/populations.plink.ped',
        'output/081_genome_mapped_stacks/final_popgen/populations.loci.fa'
    params:
        stacks_dir = 'output/081_genome_mapped_stacks',
        outdir = 'output/081_genome_mapped_stacks/final_popgen'
    singularity:
        stacks2beta_container
    threads:
        20
    log:
        'output/logs/081_genome_mapped_stacks/pops_stats.log'
    shell:
        'populations '
        '-P {params.stacks_dir} '
        '-M {input.popmap} '
        '-O {params.outdir} '
        '-W {input.whitelist} '
        '-t {threads} '
        '-r 0 '
        '--genepop '
        '--plink '
        '--vcf '
        '--hwe '
        '--fstats '
        '--fasta_loci '
        '&> {log}'



rule make_whitelist_popmap:
    input:
        plink_file = 'output/081_genome_mapped_stacks/plink.raw'
    output:
        whitelist = 'output/081_genome_mapped_stacks/whitelist.txt',
        popmap = 'output/081_genome_mapped_stacks/ref_based_popmap.txt'
    singularity:
        r_container
    log:
        'output/logs/081_genome_mapped_stacks/make_whitelist_popmap.log'
    script:
        'src/make_whitelist_popmap.R'


# Convert ped to plink format, specify that chromosomes are unknown so that it behaves
rule convert_to_plinkraw:
    input:
        ped = 'output/081_genome_mapped_stacks/snps.ped',
        map = 'output/081_genome_mapped_stacks/snps.map'
    output:
        'output/081_genome_mapped_stacks/plink.raw'
    params:
        workdir = 'output/081_genome_mapped_stacks'
    singularity:
        plink_container
    threads:
        25
    log:
        'output/logs/081_genome_mapped_stacks/convert_plinkraw.log'
    shell:
        'cd {params.workdir} || exit 1 ; '
        'plink '
        '--ped snps.ped '
        '--map snps.map '
        '--recode A '
        '--aec'


# Filter SNPs on MAF and missing rate, also filter samples by missing rate
rule filter_snps_indivs:
    input:
        'output/081_genome_mapped_stacks/snps.gds'
    output:
        'output/081_genome_mapped_stacks/snps.ped',
        'output/081_genome_mapped_stacks/snps.map'
    params:
        maf = 0.05,
        missing_rate = 0.2,
        sample_missing_quantile = 1,
        ped_file = 'output/081_genome_mapped_stacks/snps'
    singularity:
        r_container
    threads:
        25
    log:
        'output/logs/081_genome_mapped_stacks/filter_snps_indivs.log'
    script:
        'src/filter_snps.R'

# Convert VCF to GDS format, keeping biallelic SNPs only
rule convert_to_gds:
    input:
        'output/081_genome_mapped_stacks/populations.snps.vcf'
    output:
        'output/081_genome_mapped_stacks/snps.gds'
    singularity:
        r_container
    threads:
        25
    log:
        'output/logs/081_genome_mapped_stacks/gds_convert.log'
    script:
        'src/convert_gds.R'


rule genome_population_stats:
    input:
        aln_catalog = 'output/081_genome_mapped_stacks/catalog.fa.gz',
        calls = 'output/081_genome_mapped_stacks/catalog.calls',
        popmap = 'output/010_config/filtered_popmap.txt'
    output:
        'output/081_genome_mapped_stacks/populations.snps.vcf'
    params:
        stacks_dir = 'output/081_genome_mapped_stacks',
        outdir = 'output/081_genome_mapped_stacks'
    singularity:
        stacks2beta_container
    threads:
        50
    log:
        'output/logs/081_genome_mapped_stacks/genome_pop_stats.log'
    shell:
        'populations '
        '-P {params.stacks_dir} '
        '-O {params.outdir} '
        '-M {input.popmap} '
        '-t {threads} '
        '-r 0 '
        '--genepop '
        '--plink '
        '--vcf '
        '--hwe '
        '--fstats '
        '--fasta_loci '
        '&> {log}'


# Currently only works on bwa mapped files. No idea why
rule run_gstacks:
    input:
        bam_files = expand('output/081_genome_mapped_stacks/{individual}_bwa.bam',
                            individual=ustacks_individuals),
       # popmap = 'output/081_genome_mapped_stacks/testing_popmap.txt'
        popmap = 'output/010_config/filtered_popmap.txt'
    output:
        'output/081_genome_mapped_stacks/catalog.fa.gz',
        'output/081_genome_mapped_stacks/catalog.calls'
    params:
        stacks_dir = 'output/081_genome_mapped_stacks/',
      #  file_suffix = '_{map_method}.bam'
         file_suffix = '_bwa.bam'
    threads:
        18
    log:
        'output/logs/081_genome_mapped_stacks/gstacks_bwamapped_reads.log'
    singularity:
        stacks2b_container
    shell:
        'gstacks '
        '-I {params.stacks_dir} '
        '-S {params.file_suffix} '
        '-O {params.stacks_dir} '
        '-M {input.popmap} '
        '--max-clipped 0.5 '
        '--min-mapq 1 '
        '--details '
        '--phasing-dont-prune-hets '
        '--unpaired '
        '-t {threads} '
        '&> {log}'



# Sort the sam file by locus name. Only needed for bbmap
rule sort_sam:
    input:
        sam = 'output/081_genome_mapped_stacks/{individual}_bbmap.sam',
        ref = 'data/flye_denovo_full.racon.fasta'
    output:
        'output/081_genome_mapped_stacks/{individual}_bbmap.bam'
    log:
        'output/logs/081_genome_mapped_stacks/samtools_sort_{individual}_bbmap.log'
    singularity:
        samtools_container
    shell:
        'samtools sort '
        '{input.sam} '
        '-o {output} '
        '-O BAM '
        '--reference {input.ref} '
        '&> {log}'

# Map GBS reads to the genome using BWA
# NOTE this step was run separately and the logs for all indivs are combined into bwa_mapping_samples.log
rule map_tidied_reads_bwa:
    input:
        tidy_reads = 'output/021_filtered/{individual}.fq.gz',
        index = expand('output/080_against_genome/flye_denovo_full.racon.fasta.{suffix}',
               suffix=['amb', 'ann', 'bwt', 'pac', 'sa'])
    output:
        'output/081_genome_mapped_stacks/{individual}_bwa.sam'
    params:
        prefix = 'output/080_against_genome/flye_denovo_full.racon.fasta',
    log:
        'output/logs/081_genome_mapped_stacks/bwa_map_reads_{individual}.log'
    singularity:
        bwa_container
    shell:
        'bwa mem '
        '-t 8 '
        '{params.prefix} '
        '{input.tidy_reads} '
        '1> {output} '
        '2> {log}'

# Collate per sample statistics for bbmap mapping
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
    #    algorithm = 'is'
    threads:
        30
    log:
        'output/logs/080_against_genome/bwa_index.log'
    singularity:
        bwa_container
    shell:
        'bwa index '
        '-p {params.prefix} '
     #   '-a {params.algorithm} '
        '{input} '
        '2> {log}'
