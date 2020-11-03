rule lofreq_indelqual:
    """
    This should adjust indel qualities instead of BSQR from GATK.
    It works for most alignments, but seems to have trouble with random
    STAR-made BAMs. I think it has something to do with how STAR makes
    their cigar strings. 
    """
    input:
        genome=get_genome,
        bam=join(config['split_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.bam"),
        output=join(config['split_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.bam.bai")
    output:
        bam=join(config['realign_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.bam"),
        bam_index=join(config['realign_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.bam.bai")
    conda: '../envs/variant.yml'
    shell:
        """
        if [[ $(samtools view {input.bam} | head -n 5) ]]; then

            lofreq indelqual --dindel -f {input.genome} -o {output.bam} {input.bam}  
            cp {input.output} {output.bam_index}

        else
            touch {output.bam}; touch {output.bam_index}
        fi
        """

def get_virus(wildcards):
    """ Get the samtools input genomes for GATK BSQR
    """
    # Get the sample dataframe
    sample_df = pd.read_csv(config['samples']['file'])
    # Get the correct wildcard
    genome = sample_df.loc[sample_df['Run'] == wildcards.accession, 'Virus'].iloc[0]
    # Generate the ouput file
    return genome


rule bqsr:
    """ 
    Re-calibrate base-qualities with `BSQR` before SNP/indel calling. This also
    does not work. In this case it's because there is now input VCF of "known sites". 
    I'm not sure how important this is for small genomes, though. 
    """
    input:
    	fa=lambda wildcards: expand(join(config['index_dir']['samtools'], '{genome}.fa'), genome=get_virus(wildcards)),
    	idx=lambda wildcards: expand(join(config['index_dir']['samtools'], '{genome}.fa'), genome=get_virus(wildcards)),
        idxdict=lambda wildcards: expand(join(config['index_dir']['samtools'], '{genome}.fa'), genome=get_virus(wildcards)),
        bam=join(config['split_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.bam"),
        bai=join(config['split_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.bam.bai")
    output:
        outbam=join(config['recal_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.bam"),
        outbai=join(config['recal_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.bam.bai")
    params: recal=join(config['recal_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.recal")
    conda: '../envs/variant.yml'
    shell:
        """
        gatk BaseRecalibrator \
            -R {input.fa} \
            -I {input.bam} \
            --use-original-qualities \
            -O {params.recal}

        gatk ApplyBQSR \
            --add-output-sam-program-record \
            -R {input.fa} \
            -I {input.bam} \
            --use-original-qualities \
            -O {output.outbam} \
            --bqsr-recal-file {params.recal}
        """


rule remove_duplicates:
    """ This rule uses `Picard` to mark/remove probable PCR duplicates.
    """
    input: 
        bam=join(config['split_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.{organism}.bam")
    output:
        bam=join(config['dedup_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.{organism}.bam"),
        metrics=join(config['dedup_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.{organism}.metrics.txt")
    conda: '../envs/picard.yml'
    shell:
        """
        picard MarkDuplicates REMOVE_DUPLICATES=true INPUT={input.bam} \
            OUTPUT={output.bam} METRICS_FILE={output.metrics} 
        """

