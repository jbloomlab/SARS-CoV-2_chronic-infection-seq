### ======= Plot figures from dedicated R scipts for manuscript ======= ###
#
# Author: Will Hannon 
# Email: wwh22@uw.edu
# Date: 11/19/2020
#

rule plot_rbd_escape_mutations:
    """ Plot the escape mutations in the RBD for figure 2C of the manuscript.
    """
    input: join(config['pileup_dir'], "pysam_variants.csv") 
    output: join(config['figure_dir'], "rbd_escape_mutations.svg") 
    conda: "../envs/r.yml"
    script: "../scripts/plot_rbd_escape_mutations.R"


rule plot_all_spike_mutations:
    """ Plot all mutations in spike for supplemental figure S4.
    """
    input: join(config['pileup_dir'], "pysam_variants.csv") 
    output: join(config['figure_dir'], "all_spike_mutations.svg") 
    conda: "../envs/r.yml"
    script: "../scripts/plot_all_spike_mutations.R"


rule plot_coverage_in_spike:
    """ Plot the read overage over spike (capped at 500X) for supplemental figure S4.
    """
    input: join(config['coverage_dir'], "merged.depth")
    output: join(config['figure_dir'], "coverage_over_spike.svg") 
    params: config['samples']['file']
    conda: "../envs/r.yml"
    script: "../scripts/plot_coverage_in_spike.R"

rule plot_phased_rbd_mutations:
    """ Plot the phased mutations in the RBD for timepoints 143,146,152.
    """
    input: join(config['pileup_dir'], "pysam_variants_phasing.csv")
    output: join(config['figure_dir'], "phased_rbd_mutations.svg") 
    params: config['samples']['file']
    conda: "../envs/r.yml"
    script: "../scripts/plot_phased_rbd_mutations.R"
