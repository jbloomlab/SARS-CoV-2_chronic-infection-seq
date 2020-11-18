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


rule deinterleave_fastq:
    """ 
    Fastq reads are paired-end, but they are interleaved into a single file. 
    Parse these into two separate files for downstream analysis. 
    """
    input: join(config['fastq_dir'], "{accession}", "{accession}.fastq.gz")
    output: join(config['fastq_dir'], "deinterleaved", "{accession}", "{accession}_1.fastq.gz"),
            join(config['fastq_dir'], "deinterleaved", "{accession}", "{accession}_2.fastq.gz")
    conda: '../envs/seqtk.yml'
    shell:
        """
        pigz -c -d {input} | tee >(seqtk seq -1 - | pigz > {output[0]}) | seqtk seq -2 - | pigz > {output[1]} 
        """
