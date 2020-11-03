### ======= Genome guided assembly of consensus ======= ###

rule trinity:
    """ Assembly transcriptome from RNA-seq alignments.
    """
    input: 
        bam=join(config['dedup_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.bam"),
        bam_index=join(config['dedup_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.bam.bai")
    output:
        directory(join(config['trinity_dir'], "{accession}.{aligner}.trinity"))
    threads: config['max_cpu']
    conda: "../envs/assembly.yml"
    shell:
        """
        Trinity --genome_guided_bam {input.bam} \
            --max_memory 50G --genome_guided_max_intron 1 --CPU {threads} --output {output}
        """
    

rule vicuna: 
    """ Assembly of genomes with mVicuna. 
    """
    input:
        join(config['sra_dir'], "{accession}", "{accession}_1.fastq"),
        join(config['sra_dir'], "{accession}", "{accession}_2.fastq")
    output:
        directory(join(config['vicuna_dir'], "{accession}.{aligner}.vicuna"))
    params: 
        outfq=[join(config['vicuna_dir'], "{accession}", "{accession}_1.out.fastq"),join(config['vicuna_dir'], "{accession}", "{accession}_2.out.fastq")],
        outsfq=join(config['vicuna_dir'], "{accession}", "{accession}.singleton.fastq")
    threads: config['max_cpu']
    conda: "../envs/assembly.yml"
    shell:
        """
        mvicuna -tasks DupRm,SFrqEst \
        -ipfq {input[0]},{input[1]} \
        -opfq {params.outfq[0]},{params.outfq[1]} \
        -osfq test_singleton.fastq \
        -drm_op duprm_fas_1.fastq,duprm_fas_2.fastq
        """
   








mvicuna -tasks DupRm,SFrqEst -ipfq /home/whannon/whannon/2020/viral-deepseq/results/sra/RCD_02_NoIndex_L001/RCD_02_NoIndex_L001_1.fastq,/home/whannon/whannon/2020/viral-deepseq/results/sra/RCD_02_NoIndex_L001/RCD_02_NoIndex_L001_2.fastq -opfq test_fas_1.fastq,test_fas_2.fastq -osfq test_singleton.fastq -drm_op duprm_fas_1.fastq,duprm_fas_2.fastq    

