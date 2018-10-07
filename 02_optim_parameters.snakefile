#!/usr/bin/env python3

import os
import pandas

#############
# FUNCTIONS #
#############

def all_read_files(data_dir):
    '''
    Return the fastq.gz files that matches wildcards.fc
    '''
    data_dir_files = list((dirpath, filenames)
                          for (dirpath, dirnames, filenames)
                          in os.walk(data_dir))
    all_read_files = []
    for dirpath, filenames in data_dir_files:
        for filename in filenames:
            if filename.endswith('.fastq.gz'):
                all_read_files.append(os.path.join(dirpath, filename))
    return(all_read_files)

def generate_fc_dicts(key_file, data_dir):
    '''
    Parse the key_file and return dicts
    '''
    # initialise the dictionaries
    fc_to_indiv = {}
    fc_to_readfile = {}
    # find files
    read_files = all_read_files(data_dir)
    # read the key data
    key_data = pandas.read_csv(key_file)
    grouped_key_data = key_data.groupby('key')
    # populate dicts
    for name, group in grouped_key_data:
        fc_to_indiv[name] = sorted(set(x for x in group['sample_name']))
        fc_to_readfile[name] = [x for x in read_files
                                if name in os.path.basename(x)][0]
    return(fc_to_indiv, fc_to_readfile)

###########
# GLOBALS #
###########

data_dir = 'data/asw_para_matched'
key_file = 'data/asw_para_matched/combined_key_data.csv'

#########
# SETUP #
#########

fc_to_indiv, fc_to_readfile = generate_fc_dicts(
    key_file,
    data_dir)
all_fcs = list(set(fc_to_indiv.keys()))
all_indivs = sorted(set(y for x in all_fcs for y in fc_to_indiv[x]))

#########
# RULES #
#########

subworkflow process_reads:
    snakefile: 'process_reads.snakefile'

# Compare optimisation results with results using default parameters
rule compare_defaults:
    input:
        'output/030_optim/stats_n/samplestats_combined.csv',
        expand('output/021_filtered/{individual}.fq.gz',
               individual=all_indivs),
        popmap = 'output/010_config/full_popmap.txt'
    output:
        'output/030_optim/compare_defaults/optimised_samplestats_combined.csv'
    threads:
        50
    params:
        outdir = 'output/030_optim',
        indir = 'output/021_filtered'
    log:
        'output/logs/030_optim/compare_defaults.log'
    shell:
        'stacks_parameters '
        '--mode compare_defaults '
        '-m 3 '
        '-M 4 '
        '-n 5 '
        '-o {params.outdir} '
        '--individuals 8 '
        '--replicates 3 '
        '--threads {threads} '
        '{input.popmap} '
        '{params.indir} '
        '&> {log} '

# Optimise n
rule optim_n:
    input:
        'output/030_optim/stats_Mm/samplestats_combined.csv',
        expand('output/021_filtered/{individual}.fq.gz',
               individual=all_indivs),
        popmap = 'output/010_config/full_popmap.txt'
    output:
        'output/030_optim/stats_n/samplestats_combined.csv'
    threads:
        50
    params:
        outdir = 'output/030_optim',
        indir = 'output/021_filtered'
    log:
        'output/logs/030_optim/optim_n.log'
    shell:
        'stacks_parameters '
        '--mode optim_n '
        '-m 3 '
        '-M 4 '
        '-o {params.outdir} '
        '--individuals 8 '
        '--replicates 3 '
        '--threads {threads} '
        '{input.popmap} '
        '{params.indir} '
        '&> {log} '

# Optimise m and M
rule optim_mM:
    input:
        'output/030_optim/filtering/replicate_1_popmap.txt',
        expand('output/021_filtered/{individual}.fq.gz',
               individual=all_indivs),
        popmap = 'output/010_config/full_popmap.txt'
    output:
        'output/030_optim/stats_Mm/samplestats_combined.csv'
    threads:
        50
    params:
        outdir = 'output/030_optim',
        indir = 'output/021_filtered'
    log:
        'output/logs/030_optim/optim_mM.log'
    shell:
        'stacks_parameters '
        '--mode optim_Mm '
        '-o {params.outdir} '
        '--individuals 8 '
        '--replicates 3 '
        '--threads {threads} '
        '{input.popmap} '
        '{params.indir} '
        '&> {log} '

# Select individuals to use for parameter optimisation
rule optim_setup:
    input:
        process_reads(expand('output/021_filtered/{individual}.fq.gz',
               individual=all_indivs)),
        popmap = process_reads('output/010_config/full_popmap.txt')
    output:
        'output/030_optim/filtering/replicate_1_popmap.txt'
    threads:
        50
    params:
        outdir = 'output/030_optim',
        indir = 'output/021_filtered'
    log:
        'output/logs/030_optim/optim_setup.log'
    shell:
        'stacks_parameters '
        '--mode setup '
        '-o {params.outdir} '
        '--individuals 8 '
        '--replicates 3 '
        '--threads {threads} '
        '{input.popmap} '
        '{params.indir} '
        '&> {log} '

