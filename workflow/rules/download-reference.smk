### ======= Download and format reference genomes and annotations ======= ###
#
# Author: Will Hannon 
# Email: wwh22@uw.edu
# Date: 10/30/2020
#

rule get_ref:
    """ Download the reference genome fasta.
    """
    output: join(config['ref_dir'], '{genome}.fa')
    params: ftp = lambda wildcards: config[wildcards.genome]['ref']
    wildcard_constraints: genome="[^.]+"
    shell: 'wget -O - {params.ftp} | gunzip -c > {output}'


rule get_gtf:
    """ Download the gtf format annotation. 
    """
    output: join(config['gtf_dir'], '{genome}.gtf')
    params: ftp = lambda wildcards: config[wildcards.genome]['gtf']
    wildcard_constraints: genome="[^.]+"
    shell: 'wget -O - {params.ftp} | gunzip -c > {output}'


rule cat_ref:
    """ Concatenate virus/host genomes if both genomes are present.
    """
    input: 
        virus_ref=join(config['ref_dir'], '{virus}.fa'),
        host_ref=join(config['ref_dir'], '{host}.fa')
    output: join(config['ref_dir'], '{virus}.{host}.fa')
    shell: "cat {input.virus_ref} {input.host_ref} > {output}"


rule cat_gtf:
    """ 
    Concatenate virus/host gtf files. 
    This is primarily for `STAR` which requires 
    custom GTF files, so this is not downstream of
    `get_ref`. See run instructions for details. 
    """
    input: 
        virus_gtf=join(config['gtf_dir'], '{virus}.gtf'),
        host_gtf=join(config['gtf_dir'], '{host}.gtf')
    output: join(config['gtf_dir'], '{virus}.{host}.gtf')
    shell: "cat {input.virus_gtf} {input.host_gtf} > {output}"

