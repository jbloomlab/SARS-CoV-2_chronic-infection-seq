### ======= Identify variants using vatiant calling/annotation programs ======= ###
#
# Author: Will Hannon 
# Email: wwh22@uw.edu
# Date: 10/19/2020
#
rule varscan_calling:
    """ SNP calling with Varscan. Parameters are controlled from the config file.  
    """
    input: 
        varscan=join(config['tools'], "VarScan.v2.4.0.jar"),
        pileup=join(config['pileup_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.mpileup.txt")
    output: 
        variants=join(config['variant_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.varscan.vcf"),
        tmp_snp=temp(join(config['variant_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.varscan.snp.vcf")),
        tmp_indel=temp(join(config['variant_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.varscan.indel.vcf"))
    params:
        minimum_coverage=config['min_coverage'],
        minumum_supporting_reads=config['min_reads_supporting'],
        minimum_base_quality=config['BQ'],
        minimum_variant_freq=config['min_allele_frequency'],
        strand_filter=config['strand_bias_filter']
    conda: '../envs/variant.yml'    
    shell:
        """
        # Call SNPs using the mpileup file
        java -jar {input.varscan} \
            mpileup2snp {input.pileup} \
            --output-vcf 1 \
            --min-coverage {params.minimum_coverage} \
            --min-reads2 {params.minumum_supporting_reads} \
            --min-avg-qual {params.minimum_base_quality} \
            --strand-filter {params.strand_filter} \
            --min-var-freq {params.minimum_variant_freq} > {output.tmp_snp}

        # Call InDels using mpileup file
        java -jar {input.varscan} \
            mpileup2indel {input.pileup}  \
            --output-vcf 1 \
            --min-coverage {params.minimum_coverage} \
            --min-reads2 {params.minumum_supporting_reads} \
            --min-avg-qual {params.minimum_base_quality} \
            --min-var-freq {params.minimum_variant_freq} > {output.tmp_indel}

        bcftools view {output.tmp_snp} -Oz -o {output.tmp_snp}.gz;tabix -f {output.tmp_snp}.gz

        bcftools view {output.tmp_indel} -Oz -o {output.tmp_indel}.gz;tabix -f {output.tmp_indel}.gz
        
        bcftools concat -a -o {output.variants} -O v {output.tmp_snp}.gz {output.tmp_indel}.gz

        """


rule lofreq_calling:
    """ Call variants (SNP and indels) with the reference. Using only defaut filtering.
    """
    input: bam=join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.sorted.marked.bam"),
           bai=join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.sorted.marked.bam.bai"),        
           genome=get_genome
    output: join(config['variant_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.lofreq.vcf")
    params: 
        maxdepth=1000000,
        minimum_coverage=config['min_coverage']
    conda: '../envs/variant.yml'
    threads: config['threads']['max_cpu']
    shell:
        """
        if [[ $(samtools view {input.bam} | head -n 5) ]]; then
            lofreq call-parallel --pp-threads {threads} \
                -d {params.maxdepth} \
                -f {input.genome} \
                --call-indels \
                {input.bam} \
                -o {output}
        else
            touch {output}
        fi
        """ 


def get_virus(wildcards):
    """ Get the virus name from accession
    """
    # Get the sample dataframe
    sample_df = pd.read_csv(config['samples']['file'])
    # Get the correct wildcard
    genome = sample_df.loc[sample_df['Run'] == wildcards.accession, 'Virus'].iloc[0]
    # Generate the ouput file
    return genome


rule annotate_vcf:
    """
    Use the program SnpEff to annotate the effect of 
    mutations using the custom virus genome.
    """
    input:
        virusdir=lambda wildcards: expand(join(config['tools'], 'snpEff/data/{genome}'), genome=get_virus(wildcards)),
        vcf=join(config['variant_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.{caller}.vcf")
    output: join(config['variant_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.{caller}.ann.vcf")
    conda: '../envs/java.yml'
    params:
        snpEff=join(config['tools'], 'snpEff/snpEff.jar'),
        config=join(config['tools'], "snpEff/snpEff.config"),
        genome=get_virus
    shell:
        """
        java -jar {params.snpEff} -c {params.config} -noStats \
            -v {params.genome} {input.vcf} > \
            {output}
        """


rule vcf_to_table:
    """
    Convert varscan VCF files to tables for easy data
    analysis in R or Python. This will only work with 
    Varscan2 due to the program specific data fields. 
    """
    input: join(config['variant_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.{caller}.ann.vcf")
    output: join(config['variant_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.{caller}.ann.txt")
    conda: '../envs/variant.yml'    
    shell:
        """
        if [ -s {input} ]; then
            gatk VariantsToTable \
                -V {input} \
                -F CHROM -F POS -F QUAL -F REF -F ALT \
                -F DP -F AF -F FILTER -GF DP \
                -GF RD -GF FREQ -GF SDP -GF AD -F ANN \
                -O {output}
        else
            touch {output}
        fi
        """


rule add_metadata:
    """
    This rule adds in the metadata from the csv file
    that is used to run the experiment. 
    """
    input: join(config['variant_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.{caller}.ann.txt")
    output: join(config['variant_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.{caller}.ann.csv")
    params: metadata=config['samples']['file']
    conda: '../envs/r.yml'
    script: "../scripts/add_metadata.R"


rule aggregate_variants:
    """
    This rule aggregates all of the variants. 
    """
    input: expand([join(config['variant_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.{caller}.ann.csv")], accession=pd.read_csv(config['samples']['file'])['Run'], aligner=['BWA'], caller=['varscan', 'lofreq'])
    output: join(config['variant_dir'], "variants.csv")
    run:
        paths = []
        for f in input:
            try:
                pd.read_csv(f)
                paths.append(f)
            except:
                pass
        df = pd.concat(map(pd.read_csv, paths))
        df.to_csv(output[0], index = False)


