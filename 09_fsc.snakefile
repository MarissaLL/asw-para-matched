#!/usr/bin/env python3

from pathlib import Path


def get_best_run(wildcards):
    co = checkpoints.find_best_run.get(
        model=wildcards.model,
        mig=wildcards.mig).output['outdir']
    best_file = Path(co, 'best_is_{n}.txt').as_posix() # might need .as_posix()
    best_n = glob_wildcards(best_file).n[0]
    lhood_file = ('output/0xx_fastsimcoal/best_'
                f'{wildcards.model}.{wildcards.mig}.{{n}}/lhood_dist.txt')
    formatted_lhood = lhood_file.format(n=best_n)
    return {'params': formatted_lhood}


###########
# GLOBALS #
###########


fastsimcoal_container = 'shub://MarissaLL/singularity-containers:fastsimcoal_2.6'
r_container = 'shub://MarissaLL/singularity-containers:r_3.5.0'


#########
# RULES #
#########

rule target:
    input:
        expand('output/0xx_fastsimcoal/final_res/{model}.{mig}_lhood_dist.txt',
                model=['model0','model1','model2','model3','model4'],
                mig=['no_mig', 'mig'])

# Maybe I should move the files in best_run (contain maxl params) into this final_results dir, 
# then get rid of all the other outputs and horrendous amount of junk created


# This isn't really necessary, just easier to find the files
rule move_results:
    input:
        unpack(get_best_run)
    output: 
        'output/0xx_fastsimcoal/final_res/{model}.{mig}_lhood_dist.txt'
    shell:
        'cp {input} {output}'


rule cat_results:
    input:
        lhood = expand('output/0xx_fastsimcoal/best_{{model}}.{{mig}}.{{n}}/{rep}/{{model}}.{{mig}}_maxL/{{model}}.{{mig}}_maxL.lhoods',
                        rep = range(1, 101, 1))
    output:
        'output/0xx_fastsimcoal/best_{model}.{mig}.{n}/lhood_dist.txt'
    shell:
        'grep -v "MAFLHood_1" {input} > {output}'


rule fsc_distribs:
    input:
        par_file = 'output/0xx_fastsimcoal/{model}.{mig}/run{n}/{model}.{mig}/{model}.{mig}_maxL.par', # this is 'output/0xx_fsc/{model}.{mig}/run{n}/..MaxL.par'
        sfs = 'populations_jointMAFpop1_0.obs'
    output:
        sfs = temp('output/0xx_fastsimcoal/best_{model}.{mig}.{n}/{rep}/{model}.{mig}_maxL_jointMAFpop1_0.obs'),
        par_file = temp('output/0xx_fastsimcoal/best_{model}.{mig}.{n}/{rep}/{model}.{mig}_maxL.par'),
        out = 'output/0xx_fastsimcoal/best_{model}.{mig}.{n}/{rep}/{model}.{mig}_maxL/{model}.{mig}_maxL.lhoods' 
    params:
        wd = 'output/0xx_fastsimcoal/best_{model}.{mig}.{n}/{rep}',
        model_prefix = '{model}.{mig}',
        numsims = 100
    singularity:
        fastsimcoal_container
    log:
        'fsc_distribs_{model}.{mig}.{n}.{rep}.log'
    shell:
        'cp {input.par_file} {output.par_file} ; '
        'cp {input.sfs} {output.sfs} ; '
        'cd {params.wd} ; '
        'fsc26 '
        '--ifile {params.model_prefix}_maxL.par ' 
        '--numsims {params.numsims} '
        '--msfs '
        '&> {log}'
      
      

# Find best run & calc AIC
checkpoint find_best_run:
    input:
        bestlhoods = expand('output/0xx_fastsimcoal/{{model}}.{{mig}}/run{run}/{{model}}.{{mig}}/{{model}}.{{mig}}.bestlhoods',
               run=range(1, 5, 1))
    output:
        outdir = directory('output/0xx_fastsimcoal/{model}.{mig}/best_run/')
    singularity:
        r_container
    log:
        'output/logs/0xx_fastsimcoal/find_best_run_{model}.{mig}.log'
    script: 
        'src/find_best_run.R'


### Run fsc. Why is this outputting the log in two places?
rule fsc:
    input:
        tpl = 'src/fsc_models/{model}.{mig}.tpl',
        est = 'src/fsc_models/{model}.{mig}.est',
        sfs = 'populations_jointMAFpop1_0.obs'
    output:
        bestl = 'output/0xx_fastsimcoal/{model}.{mig}/run{run}/{model}.{mig}/{model}.{mig}.bestlhoods',
        maxlpar = 'output/0xx_fastsimcoal/{model}.{mig}/run{run}/{model}.{mig}/{model}.{mig}_maxL.par',
        tpl = temp('output/0xx_fastsimcoal/{model}.{mig}/run{run}/{model}.{mig}.tpl'),
        est = temp('output/0xx_fastsimcoal/{model}.{mig}/run{run}/{model}.{mig}.est'),
        sfs = temp('output/0xx_fastsimcoal/{model}.{mig}/run{run}/{model}.{mig}_jointMAFpop1_0.obs')
    params:
        wd = 'output/0xx_fastsimcoal/{model}.{mig}/run{run}',
        model_prefix = '{model}.{mig}',
        numsims = 100
    singularity:
        fastsimcoal_container
    log:
        'fsc_{model}.{mig}.{run}.log'
    shell:
        'cp {input.tpl} {output.tpl} ; '
        'cp {input.est} {output.est} ; '
        'cp {input.sfs} {output.sfs} ; '
        'cd {params.wd} ; '
        'fsc26 '
        '--tplfile {params.model_prefix}.tpl '
        '--estfile {params.model_prefix}.est '
        '--numsims {params.numsims} '
        '--msfs '
        '--maxlhood '
        '--numloops 60 '
        '--cores 2 '
        '&> {log}'