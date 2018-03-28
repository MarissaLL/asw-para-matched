#!/usr/bin/env python3

import os
import pandas
import pathlib
import re


#############
# FUNCTIONS #
#############

# def all_key_files(data_dir):
#     '''
#     Return a list of all .txt from data_dir
#     '''
#     data_dir_files = list((dirpath, filenames)
#                           for (dirpath, dirnames, filenames)
#                           in os.walk(data_dir))
#     all_key_files = []
#     for dirpath, filenames in data_dir_files:
#         for filename in filenames:
#             if filename.endswith('.txt'):
#                 all_key_files.append(os.path.join(dirpath, filename))
#     return(all_key_files)


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


def get_full_path(x):
    return str(pathlib.Path(x).resolve())


def generate_fc_dicts(key_file, data_dir):
    '''
    Parse the key_file and return dicts
    '''
    # initialise the dictionaries
    fc_to_indiv = {}
    fc_to_readfile = {}
    # fc_to_keyfile = {}
    # find files
    read_files = all_read_files(data_dir)
    #key_files = all_key_files(data_dir)
    # read the key data
    key_data = pandas.read_csv(key_file)
    grouped_key_data = key_data.groupby('key')
    # populate dicts
    for name, group in grouped_key_data:
        fc_to_indiv[name] = sorted(set(x for x in group['sample_name']))
        fc_to_readfile[name] = [x for x in read_files
                                if name in os.path.basename(x)][0]
        # fc_to_keyfile[name] = [x for x in key_files
        #                        if name in os.path.basename(x)][0]
    return(fc_to_indiv, fc_to_readfile)


def read_keydata_and_write_config(key_file, outdir):
    '''
    Parse the keyfile and generate config files in outdir
    '''
    # read the key data
    key_data = pandas.read_csv(key_file)
    grouped_key_data = key_data.groupby('key')
    for name, group in grouped_key_data:
            config_file = os.path.join(outdir, '{}_config'.format(name))
            subset = group[['barcode', 'sample_name']]
            if len(subset) > 0:
                subset.to_csv(config_file,
                              sep='\t',
                              header=False,
                              index=False)


###########
# GLOBALS #
###########

data_dir = 'data/test_data'
key_file = 'data/test_data/combined_key_data.csv'

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
rule target:
    input:
        expand('output/020_demux/{individual}.fq.gz',
               individual=all_indivs)

for fc in all_fcs:
    rule:
        input:
            config = 'output/010_config/{}_config'.format(fc),
            reads = fc_to_readfile[fc]
        output:
            expand('output/020_demux/{individual}.fq.gz',
                   individual=fc_to_indiv[fc])
        params:
            outdir = 'output/020_demux'
        threads:
            1
        log:
            'output/logs/020_demux/{}.log'.format(fc)
        shell:
            'process_radtags '
            '-f {input.reads} '
            '-i gzfastq -y gzfastq '
            '-b {input.config} '
            '-o {params.outdir} '
            '-c -q '
            '-t 91 '
            '--inline_null '
            '--renz_1 apeKI --renz_2 mspI '
            '&> {log} '

# 010 generate stacks config
rule generate_config_files:
    input:
        key_file = key_file
    threads:
            1
    output:
        expand('output/010_config/{fc_name}_config',
               fc_name=all_fcs)
    params:
        outdir = 'output/010_config'
    run:
        read_keydata_and_write_config(input.key_file, params.outdir)

