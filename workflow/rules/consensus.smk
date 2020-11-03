### ======= Generate consensus sequence from BAM ======= ###

rule generate_consensus:
    """ Use `bcftools` to generate a consensus sequence.
    """
    input:
        genome=get_genome,
        bam=join(config['dedup_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.bam"),
        bam_index=join(config['dedup_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.bam.bai")
    output:
        vcf=join(config['consensus_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.bcftools.vcf.gz"),
        fa=join(config['consensus_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.consensus.fa")
    params:
        name=">{accession}.{aligner}"
    conda: '../envs/consensus.yml'    
    shell:
        """
        # Call variants with `bcftools`
        bcftools mpileup -Ou -f {input.genome} {input.bam} | \
            bcftools call -Ou -mv -A | \
            bcftools norm -f {input.genome} -Oz -o {output.vcf}

        # Index the variant calls
        tabix {output.vcf}

        # Generate the consensus
        bcftools consensus -f {input.genome} {output.vcf} > {output.fa}

        # Fix Fasta header
        sed -i "1s/.*/{params.name}/" {output.fa}
  
        """

rule join_consensus:
    """ Join the consensus sequences for MSA
    """
    input: 
        expand(join(config['consensus_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.consensus.fa"), accession=pd.read_csv(config['samples']['file'])['Run'], aligner=['BWA','STAR'])
    output: 
        join(config['consensus_dir'], "all_samples_consensus.fa")
    shell:
        """
        cat {input} > {output}
        """
    