### ======= Import runs as fastq formatted files ======= ###
#
# Author: Will Hannon 
# Email: wwh22@uw.edu
# Date: 10/30/2020
#

rule fetch_fastq_se:
    """ Get single-ended fastq files from either local or remote source.
    """
    output: join(config['fastq_dir'], "{accession}", "{accession}.fastq.gz")
    params: join(config['fastq_dir'], "{accession}")
    threads: config['threads']['fastq_download']
    conda: '../envs/fastq.yml'
    script: '../scripts/fetch_fastq.py'


rule fetch_fastq_pe:
    """ Get paired-end fastq files from either local or remote source.
    """
    output: join(config['fastq_dir'], "{accession}", "{accession}_1.fastq.gz"),
            join(config['fastq_dir'], "{accession}", "{accession}_2.fastq.gz")
    params: join(config['fastq_dir'], "{accession}")
    threads: config['threads']['fastq_download']
    conda: '../envs/fastq.yml'
    script: '../scripts/fetch_fastq.py'

