### ======= Short read alignment of filtered reads to normal or hybrid genome ======= ###
#
# Author: Will Hannon 
# Email: wwh22@uw.edu
# Date: 11/02/2020
#

rule bwa_align:
    """ 
    Perform short read alignment with `bwa-mem`.
    Sort the aligned reads with samtools sort.
    Will align to virus or virus and host depending 
    on the status of `config[filter_reads]`.
    """
    input: 
        reads=lambda wildcards: get_avaliable_filtered_fastqs(wildcards) if config['filter_reads'] else get_avaliable_trimmed_fastqs(wildcards),
        genome=lambda wildcards: get_genome(wildcards, index = "BWA") if config['filter_reads'] else get_hybrid_genome(wildcards)
    output: join(config['align_dir'], "BWA", "{accession}", "{accession}.BWA.sorted.bam")
    threads: config['threads']['max_cpu']
    params: tmp=join(config['align_dir'], "BWA", "{accession}", "{accession}.final.tmp")
    conda: '../envs/align.yml'
    shell: 
        """
        bwa mem -t {threads} \
            {input.genome} \
            {input.reads} | \
            samtools sort -o {output} -T {params.tmp} - 
        """


rule split_bam:
    """ 
    Split bam file into reads that map to virus (as identified
    by the `contig` name) and reads that dont map to virus. 
    """
    input: join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.sorted.bam"),
    output: virus_bam=join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.sorted.bam"),
            other_bam=join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.other.sorted.bam"),
            virus_sam=temp(join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.sorted.sam")),
            other_sam=temp(join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.other.sorted.sam"))
    params: contig=lambda wildcards: config[pd.read_csv(config['samples']['file']).set_index('Run').at[wildcards.accession, 'Virus']]['contig']
    conda: '../envs/samtools.yml'
    shell:
        """
        # Get the header of the bam file
        samtools view -H {input} | tee {output.virus_sam} {output.other_sam}

        # Extract only the reads that do map to virus
        samtools view {input} \
            | awk '$3== "{params.contig}"' \
            >> {output.virus_sam} 

        # Sort the results and convert to bam
        samtools view -Sb {output.virus_sam} | \
            samtools sort - -o {output.virus_bam}

        # Extract only the reads that don't map to virus
        samtools view {input} \
            | awk '$3!= "{params.contig}"' \
            >> {output.other_sam} 

        # Sort the results and convert to bam
        samtools view -Sb {output.other_sam} | \
            samtools sort - -o {output.other_bam}
        """


rule mark_duplicates:
    """ 
    This rule uses `Picard` to mark/remove probable PCR duplicates.
    Only removes duplicates if `true` is specified in the config file, 
    otherwise this tool only marks them.
    """
    input: bam=join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.{split}.sorted.bam")
    output: bam=join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.{split}.sorted.marked.bam"),
            metrics=join(config['qc_dir'], "{accession}", "Picard", "{accession}.{aligner}.{split}.metrics.txt")
    params: remove_duplicates=config['remove_duplicates']
    conda: '../envs/picard.yml'
    shell:
        """
        picard MarkDuplicates REMOVE_DUPLICATES={params.remove_duplicates} INPUT={input.bam} \
            OUTPUT={output.bam} METRICS_FILE={output.metrics} 
        """


rule index_bam:
    """ Index mapped, sorted, and marked bam files. 
    """
    input: join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.{split}.sorted.marked.bam"),
    output: join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.{split}.sorted.marked.bam.bai")
    conda: '../envs/samtools.yml'
    shell: "samtools index {input}"


