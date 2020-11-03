### ======= Index reference genomes in specified format ======= ###
#
# Author: Will Hannon 
# Email: wwh22@uw.edu
# Date: 10/30/2020
#

rule bwa_index:
    """ Index the genome with `BWA` before mapping.
    """
    input: join(config['ref_dir'], '{genome}.fa')
    output: join(config['index_dir']['bwa'], '{genome}.fa')
    conda: '../envs/align.yml'
    params: algorithm="bwtsw"
    shell:
        """
        cp {input} {output}
        bwa index -a {params.algorithm} {output}
        """

rule samtools_index:
    """ Index genome with `samtools` for `BSQR`. 
    """
    input: join(config['ref_dir'], '{genome}.fa')
    output: 
        fa=join(config['index_dir']['samtools'], '{genome}.fa'),
        idx=join(config['index_dir']['samtools'], '{genome}.fa.fai'),
        idxdict=join(config['index_dir']['samtools'], '{genome}.fa.dict')
    conda: '../envs/samtools.yml'
    shell:
        """
        cp {input} {output.fa}
        samtools faidx {output.fa}
        samtools dict -o {output.idxdict} {output.fa}
        """

