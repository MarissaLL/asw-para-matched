#!/usr/bin/env python3


###########
# GLOBALS #
###########


fastsimcoal_container = 'shub://MarissaLL/singularity-containers:fastsimcoal_2.6'
r_container = 'shub://MarissaLL/singularity-containers:r_3.5.0'


#########
# RULES #
#########

models = ['model0_mig','model0_no_mig', 'model1_mig', 'model1_no_mig',
 'model2_mig','model2_no_mig']
 # , 'model3_mig', 'model3_no_mig'
runs = range(1,3,1)


rule target:
	input:
		expand('{model}_run_{run}/asw_{model}.tpl',
				model = models, run = runs),
		expand('{model}_run_{run}/asw_{model}/asw_{model}.brent_lhoods',
				model = models, run = runs),
		expand('{model}.AIC',
				model = models)



### Think this needs an obs file called model_maxL_joint...obs
# rule compare_likelihood_dists:
	# input:
	# output:
	# singularity:
	# 	fastsimcoal_container
	# shell:
	# 	'fsc26 '
	# 	'--ifile asw_{model}_maxL.par '
	# 	'--numsims 1000000 '
	# 	'-m -q -0'

########### This requires that bestlikelihood, est and tpl files are in the same dir, 
########### so needs some files moving around before it can just run
rule best_aic_calc:
	input:
		'{model}_run_{run}/asw_{model}/asw_{model}.bestlhoods'
	output:
		'{model}.AIC'
	params:
		model_prefix = '{model}_run_{run}/asw_{model}/asw_{model}'
	singularity:
		r_container
	script:
		'compareAIC.R {params.model_prefix}'


rule fastsimcoal:
	input:
		template = '{model}_run_{run}/asw_{model}.tpl',
		est_file = '{model}_run_{run}/asw_{model}.est',
		sfs = '{model}_run_{run}/asw_{model}_jointMAFpop1_0.obs'
	output:
		'{model}_run_{run}/asw_{model}/asw_{model}.brent_lhoods',
		'{model}_run_{run}/asw_{model}/asw_{model}.bestlhoods'
	params:
		numsims = 100,
		current_run = '{model}_run_{run}',
		model_prefix = 'asw_{model}'
	singularity:
		fastsimcoal_container
	log:
		'{model}_run_{run}/fsc.log'
	shell:
		'cd {params.current_run} && '
		'fsc26 '
		'--tplfile {params.model_prefix}.tpl '
		'--numsims {params.numsims} '
		'--msfs '
		'--estfile {params.model_prefix}.est '
		'--maxlhood '
		'--numloops 60 '
		'--cores 20 '
		'&> fsc.log && '
		'cd .. '


rule input_dir_struct:
	input:
		template = 'asw_{model}.tpl',
		est_file = 'asw_{model}.est',
		sfs = 'asw_{model}_jointMAFpop1_0.obs'
	output:
		tpl = '{model}_run_{run}/asw_{model}.tpl',
		est = '{model}_run_{run}/asw_{model}.est',
		sfs = '{model}_run_{run}/asw_{model}_jointMAFpop1_0.obs'
	shell:
		'cp {input.template} {output.tpl} && '
		'cp {input.est_file} {output.est} && '
		'cp {input.sfs} {output.sfs}'






