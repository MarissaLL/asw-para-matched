#!/usr/bin/env python3


###########
# GLOBALS #
###########

bayescan_container = 'shub://MarissaLL/singularity-containers:bayescan_2.1'
bbmap_container = 'shub://TomHarrop/singularity-containers:bbmap_38.45'
bwa_container = 'shub://TomHarrop/singularity-containers:bwa_0.7.17'
pgdspider_container = 'shub://MarissaLL/singularity-containers:pgdspider_2.1.1.5'
samtools_container = 'shub://TomHarrop/singularity-containers:samtools_1.9'
stacks_container = 'shub://TomHarrop/singularity-containers:stacks_2.3e'
stacks2beta_container = 'shub://TomHarrop/singularity-containers:stacks_2.0beta9@bb2f9183318871f6228b51104056a2d0'
vcftools_container = 'shub://MarissaLL/singularity-containers:vcftools_0.1.17@230db32b3097775cd51432092f9cbcb1'

bayescan_subsets = ['compared_island', 'compared_4pops', 'compared_para']
catalog_version = ['filtered', 'full']

#########
# RULES #
#########
# grep -E -o ".{0,5}scaffold.{0,5}" populations.snps.vcf > scaffold_names.txt

rule target:
    input:
        expand('output/080_against_genome/{catalog_ver}_catalog/populations.snps.vcf',
                catalog_ver=catalog_version),
        'output/080_against_genome/bbmapped_full.sam',
        expand('output/080_against_genome/{catalog_ver}_catalog/maf_stats.INFO',
                catalog_ver=catalog_version),
        expand('output/080_against_genome/{catalog_ver}_catalog/missingness.lmiss',
                catalog_ver=catalog_version)
        # expand('trial/{bayescan_subset}.sel',
        #        bayescan_subset=bayescan_subsets)



# rule bayescan:
#     input:
#         genotypes = 'output/080_against_genome/filtered_catalog/{bayescan_subset}.geste-outputformat'
#     output:
#         'trial/{bayescan_subset}.sel',
#         'trial/{bayescan_subset}_fst.txt'
#     params:
#         outdir = 'trial/',
#         outname = '{bayescan_subset}'
#     singularity:
#         bayescan_container
#     threads:
#         6
#     log:
#         'trial/{bayescan_subset}_bayescan_filtered_catalog.log'
#     shell:
#         'bayescan_2.1 '
#         '{input.genotypes} '
#         '-od {params.outdir} '
#         '-o {params.outname} '
#         '-pilot 5000 '
#         '-nbp 20 '
#         '-burn 15000 '
#         '-n 30000 '
#         '-thin 10 '
#         '-pr_odds 500 '
#         '-out_pilot '
#         '-out_freq '
#         '&> {log}'

        

# ## Note that the location of the popmap is actually specified within the spid file
# ## It is only specified here to link the dependencies of the rules
# rule convert_subset_bayescan_inputs:
#     input:
#         vcf = 'output/080_against_genome/filtered_catalog/{bayescan_subset}.recode.vcf',
#         spid = 'data/{bayescan_subset}.spid',
#         popmap = 'output/070_bayescan/popmap_{bayescan_subset}.txt'
#     output:
#         'output/080_against_genome/filtered_catalog/{bayescan_subset}.geste-outputformat'
#     params:
#         in_format = 'VCF',
#         out_format = 'GESTE_BAYE_SCAN',
#         out_path = 'output/080_against_genome/filtered_catalog/{bayescan_subset}.geste'
#     singularity:
#         pgdspider_container
#     threads:
#         50
#     log:
#         'output/logs/080_against_genome/pgdspider_{bayescan_subset}_filtered_catalog.log'
#     shell:
#         'java -jar /opt/pgdspider/PGDSpider2-cli.jar '
#         '-inputfile {input.vcf} '
#         '-inputformat {params.in_format} '
#         '-outputfile {params.out_path}'
#         '-outputformat {params.out_format} '
#         '-spid {input.spid} '
#         '&> {log}'

# # Subset vcfs to only include individuals needed in each bayescan run
# rule subset_vcfs:
#     input:
#         individual_list = 'output/070_bayescan/{bayescan_subset}_indivs.txt',
#         full_vcf = 'output/080_against_genome/filtered_catalog/populations.snps.vcf'
#     output: 
#         subset_vcf = 'output/080_against_genome/filtered_catalog/{bayescan_subset}.recode.vcf'
#     params:
#         outname = 'output/080_against_genome/filtered_catalog/{bayescan_subset}'
#     singularity:
#         vcftools_container
#     threads:
#         25
#     log:
#         'output/logs/080_against_genome/{bayescan_subset}_subset_vcf.log'
#     shell:
#         'vcftools '
#         '--vcf {input.full_vcf} '
#         '--out {params.outname} '
#         '--keep {input.individual_list} '
#         '--recode '
#         '&> {log}'

## Extract loci that fulfil the MAF criteria. ALSO need to do missingness for the full catalog.
#rule extract_passing_loci:
    input:
        expand('output/080_against_genome/{catalog_ver}_catalog/maf_stats.INFO',
                catalog_ver=catalog_version),
        expand('output/080_against_genome/{catalog_ver}_catalog/missingness.lmiss',
                catalog_ver=catalog_version)
    output:
        'output/080_against_genome/pass_loci.txt'


## Extract missingness info
rule extract_missingness_rate:
    input:
        'output/080_against_genome/{catalog_ver}_catalog/populations.snps.vcf'
    output:
        'output/080_against_genome/{catalog_ver}_catalog/missingness.lmiss'
    params:
        output_name = 'output/080_against_genome/{catalog_ver}_catalog/missingness'
    log:
        'output/logs/080_against_genome/extract_missingness_{catalog_ver}.log'
    shell:
        'vcftools '
        '--vcf {input} '
        '--missing-site '
        '--out {params.output_name} '
        '2> {log}'

## Extract maf info
rule extract_allele_frequencies:
    input:
        'output/080_against_genome/{catalog_ver}_catalog/populations.snps.vcf'
    output:
        'output/080_against_genome/{catalog_ver}_catalog/maf_stats.INFO'
    params:
        output_name = 'output/080_against_genome/{catalog_ver}_catalog/maf_stats'
    log:
        'output/logs/080_against_genome/extract_AF_{catalog_ver}.log'
    shell:
        'vcftools '
        '--vcf {input} '
        '--get-INFO AF '
        '--out {params.output_name} '
        '2> {log}'

# Run stacks populations on the mapped filtered catalog loci
rule genome_population_stats:
    input:
        aln_catalog = 'output/080_against_genome/{catalog_ver}_catalog/catalog.fa.gz',
        calls = 'output/080_against_genome/{catalog_ver}_catalog/catalog.calls',
        popmap = 'output/010_config/filtered_popmap.txt'
    output:
        'output/080_against_genome/{catalog_ver}_catalog/populations.snps.vcf'
    params:
        stacks_dir = 'output/080_against_genome/{catalog_ver}_catalog',
        outdir = 'output/080_against_genome/{catalog_ver}_catalog'
    singularity:
        stacks2beta_container
    threads:
        50
    log:
        'output/logs/080_against_genome/genome_pop_stats_{catalog_ver}.log'
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


# Integrate alignment info back into stacks workflow. Using edited file because the original only outputs
# loci without a t in their name
rule integrate_alignments:
    input:
        catalog = 'output/040_stacks/catalog.fa.gz',
        calls = 'output/040_stacks/catalog.calls',
        bam = 'output/080_against_genome/bbmapped_{catalog_ver}.bam'
    output:
        aln_catalog = 'output/080_against_genome/{catalog_ver}_catalog/catalog.fa.gz',
        tsv = 'output/080_against_genome/{catalog_ver}_catalog/locus_coordinates.tsv',
        calls = 'output/080_against_genome/{catalog_ver}_catalog/catalog.calls'
    params:
        stacks_dir = 'output/040_stacks',
        out_dir = 'output/080_against_genome/{catalog_ver}_catalog/'
    threads:
        1
    log:
        'output/logs/080_against_genome/integrate_{catalog_ver}_alignments.log'
    
    shell:
        ' ./stacks-integrate-alignments-edited '
        '-P {params.stacks_dir} '
        '-B {input.bam} '
        '-O {params.out_dir} '
        '&> {log}' 




# Sort the sam file by locus name. Output as a bam file
rule sort_sam:
    input:
        sam = 'output/080_against_genome/bbmapped_{catalog_ver}.sam',
        ref = 'data/flye_denovo_full.racon.fasta'
    output:
        'output/080_against_genome/bbmapped_{catalog_ver}.bam'
    log:
        'output/logs/080_against_genome/samtools_sort_{catalog_ver}.log'
    shell:
        'samtools sort '
        '{input.sam} '
        '-o {output} '
        '-O BAM '
        '-n '
        '--reference {input.ref} '
        '&> {log}'




rule bbmap_filtered:
    input:
        genome = 'data/flye_denovo_full.racon.fasta',
        loci = 'output/080_against_genome/loci_noheader.fa'
    output:
        'output/080_against_genome/bbmapped_filtered.sam'
    params:
        ref_path = 'output/080_against_genome/'
    log:
        'output/logs/080_against_genome/bbmap_index_map_filtered.log'
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





# Indexes the genome, then maps the full de novo catalog to it. 
# Not currently used for anything further. 
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