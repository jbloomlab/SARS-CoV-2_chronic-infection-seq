### ======= Filter processed fastq files by kmers with BBduk ======= ###
#
# Author: Will Hannon 
# Email: wwh22@uw.edu
# Date: 10/30/2020
#

rule filter_reads_se:
    """ 
    Filter out reads using kmer filtering w/ BBduk Single-Ended implementation.
    This will only work for the viral genomes because larger genomes
    begin to require huge amounts of heap memory.

    Genome is determined from the `samples.csv` file configuration.
    The genomes are indexed with samtools. 
    """
    input: 
        reads=get_avaliable_trimmed_fastqs,
        genome=get_genome
    output: 
        matched=join(config['filter_dir'], "{accession}", "{accession}.filtered.fastq.gz"),
        unmatched=join(config['filter_dir'], "{accession}", "{accession}.unfiltered.fastq.gz"),
        stats=join(config['qc_dir'], "{accession}", "BBduk", "{accession}.filter.stats")
    threads: config['threads']['max_cpu']
    params: error=join(config['filter_dir'], "{accession}", "{accession}.error.log")
    conda: '../envs/filter.yml'
    shell:
        """
        bbduk.sh -Xmx80g \
            in={input.reads} \
            out={output.unmatched} \
            outm={output.matched} \
            ref={input.genome} \
            k=31 \
            hdist=2 \
            stats={output.stats} \
            overwrite=TRUE \
            t={threads} \
            &> {params.error}
        """

rule filter_reads_pe:
    """ 
    Filter out reads using kmer filtering w/ BBduk Paried-End implementation.
    This will only work for the viral genomes because larger genomes
    begin to require huge amounts of heap memory.

    Genome is determined from the `samples.csv` file configuration.
    The genomes are indexed with samtools. 
    """
    input: 
        reads=get_avaliable_trimmed_fastqs,
        genome=get_genome
    output: 
        matched=[join(config['filter_dir'], "{accession}", "{accession}_1.filtered.fastq.gz"), 
                 join(config['filter_dir'], "{accession}", "{accession}_2.filtered.fastq.gz")],
        unmatched=[join(config['filter_dir'], "{accession}", "{accession}_1.unfiltered.fastq.gz"), 
                   join(config['filter_dir'], "{accession}", "{accession}_2.unfiltered.fastq.gz")],
        stats=join(config['qc_dir'], "{accession}", "BBduk", "{accession}.filter.stats")
    threads: config['threads']['max_cpu']
    params: error=join(config['filter_dir'], "{accession}", "{accession}.error.log")
    conda: '../envs/filter.yml'
    shell:
        """
        bbduk.sh -Xmx80g \
            in1={input.reads[0]} \
            in2={input.reads[1]} \
            out1={output.unmatched[0]} \
            out2={output.unmatched[1]} \
            outm1={output.matched[0]} \
            outm2={output.matched[1]} \
            ref={input.genome} \
            k=31 \
            hdist=2 \
            stats={output.stats} \
            overwrite=TRUE \
            t={threads} \
            &> {params.error}
        """