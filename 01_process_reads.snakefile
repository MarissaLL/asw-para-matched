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

data_dir = 'data/asw_para_matched'
key_file = 'data/asw_para_matched/combined_key_data.csv'

bbmap_container = 'shub://MarissaLL/singularity-containers:bbmap_37.92@fa973c37883055c243c5e37f82f68f4d'
fastqc_container = 'shub://MarissaLL/singularity-containers:fastqc_0.11.5@2fa47d3690d60343361e575ae6c54a33'
r_container = 'shub://MarissaLL/singularity-containers:r_3.5.0@1f713cf93d67765d2d37eb87339495de5'
stacks2b_container = 'shub://TomHarrop/singularity-containers:stacks_2.0b@099f0c7d8c8ff2baf7ad763ad7bcd17b'

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
        'output/010_config/filtering_stats.csv',
        'output/010_config/full_popmap.txt',
        'output/010_config/tidy_sample_info.tsv',
        expand('output/022_fastqc/before_filter/{individual}_fastqc.zip',
               individual=all_indivs),
        expand('output/022_fastqc/after_filter/{individual}_fastqc.zip',
               individual=all_indivs),
        expand('output/021_filtered/readlength_hist/{individual}_readlength_hist.txt',
              individual=all_indivs)
    
        
# Generate the list of individuals to use
rule generate_popmap:
    input:
        key_file = key_file
    output:
        popmap = 'output/010_config/full_popmap.txt'
    threads:
        1
    run:
        key_data = pandas.read_csv(input.key_file)
        subset = key_data.loc[key_data['population'] != 'GBSNEG',
                              ['sample_name', 'population']]
        subset.to_csv(output.popmap,
                      sep='\t',
                      header=False,
                      index=False)

# Run FastQC again after the filtering, to check improvement
rule fastqc_after_filter:
    input:
        expand('output/021_filtered/{individual}.fq.gz',
               individual=all_indivs)
    output:
        expand('output/022_fastqc/after_filter/{individual}_fastqc.zip',
               individual=all_indivs)
    params:
        outdir='output/022_fastqc/after_filter'
    singularity:
        fastqc_container
    threads:
        10
    log:
        'output/logs/022_fastqc/after_filter_fastqc.log'
    shell:
        'fastqc '
        '-o {params.outdir} '
        '-t {threads} '
        '{input} '
        '&> {log}'

rule combine_stats:
    input:
        stats_file = expand('output/021_filtered/adapter_stats/{individual}.txt',
            individual=all_indivs),
        gc_file = expand('output/021_filtered/gc_hist/{individual}.txt',
            individual=all_indivs)
    output:
        stats_file = 'output/010_config/filtering_stats.csv'
    singularity:
        r_container
    log:
        'output/logs/021_filtered/combined_stats.log'
    script:
        'src/combine_stats.R'

# Remove any sequence matching that of the adapter, and truncate reads to 80 bp 
rule filter_adapters:
    input:
        FQ = 'output/020_demux/{individual}.fq.gz'
    output:
        kept_FQ = 'output/021_filtered/{individual}.fq.gz',
        discarded_FQ = 'output/021_filtered/discarded_FQ/{individual}.fq.gz',
        adapter_stats = 'output/021_filtered/adapter_stats/{individual}.txt',
        truncation_stats = 'output/021_filtered/truncation_stats/{individual}.txt',
        gc_hist = 'output/021_filtered/gc_hist/{individual}.txt'
    params:
        adapters = 'data/bbduk_adapters.fa'
    singularity:
        bbmap_container
    threads:
        1
    log:
        adapter_log = 'output/logs/021_filtered/{individual}_adapters.log',
        truncation_log = 'output/logs/021_filtered/{individual}_truncation.log'
    benchmark: 
        'output/benchmarks/021_filtered/{individual}.txt'
    shell:
        'bbduk.sh '
        'threads={threads} '
        'in={input.FQ} '
        'interleaved=f '
        'out=stdout.fq '
        'stats={output.adapter_stats} '
        'ref={params.adapters} '
        'ktrim=r k=23 mink=11 hdist=1 '
        'findbestmatch=t '
        'gchist={output.gc_hist} '
        'gcbins=auto '
        '&> {log.adapter_log} '
        ' | '
        'bbduk.sh '
        'threads={threads} '
        'in=stdin.fq '
        'interleaved=f '
        'outnonmatch={output.kept_FQ} '
        'outmatch={output.discarded_FQ} '
        'stats={output.truncation_stats} '
        'overwrite=t '
        'forcetrimright=79 '
        'minlength=80 '
        'ziplevel=9 '
        '&> {log.truncation_log}'

# Calculate how much of the read is ASW sequence
rule readlength:
    input:
        FQ = 'output/020_demux/{individual}.fq.gz'
    output:
        lhist = 'output/021_filtered/readlength_hist/{individual}_readlength_hist.txt'
    params:
        adapters = 'data/bbduk_adapters.fa'
    singularity:
        bbmap_container
    threads:
        1
    log:
        match_log = 'output/logs/021_filtered/{individual}_match.log',
        lhist_log = 'output/logs/021_filtered/{individual}_readlength.log'
    shell:
        'bbduk.sh '
        'threads={threads} '
        'in={input.FQ} '
        'interleaved=f '
        'out=stdout.fq '
        'ref={params.adapters} '
        'ktrim=r k=23 mink=11 hdist=1 '
        'findbestmatch=t '
        '2> {log.match_log}'
        ' | '
        'reformat.sh '
        'in=stdin.fq '
        'int=f '
        'out=stdout.fq '
        'lhist={output.lhist} '
        '2> {log.lhist_log}'

# Run FastQC on the demultiplexed data
rule fastqc_before_filter:
    input:
        expand('output/020_demux/{individual}.fq.gz',
               individual=all_indivs)
    output:
        expand('output/022_fastqc/before_filter/{individual}_fastqc.zip',
               individual=all_indivs)
    params:
        outdir='output/022_fastqc/before_filter'
    singularity:
        fastqc_container
    threads:
        10
    log:
        'output/logs/022_fastqc/before_filter_fastqc.log'
    shell:
        'fastqc '
        '-o {params.outdir} '
        '-t {threads} '
        '{input} '
        '&> {log}'       

# Demultiplex
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
        singularity:
            stacks2b_container
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
            '--barcode_dist_1 0 '
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


rule tidy_sample_info:
    input:
        'data/sample_catalog.csv'
    output:
        'output/010_config/tidy_sample_info.tsv'
    singularity:
        r_container
    log:
        'output/logs/010_config/tidy_sample_info.log'
    script:
        'src/tidy_sample_info.R'
