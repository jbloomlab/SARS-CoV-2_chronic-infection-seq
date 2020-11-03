### ======= Variant calling and annotation ======= ###

rule varscan_calling:
    """ SNP calling with Varscan. 
    """
    input: 
        varscan=join(config['tools'], "VarScan.v2.4.0.jar"),
        genome=get_genome,
        bam=join(config['split_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.bam"),
        bam_index=join(config['split_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.bam.bai")
    output: join(config['variant_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.varscan.vcf")
    params:
        tmp_snp=temp(join(config['variant_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.varscan.snp.vcf")),
        tmp_indel=temp(join(config['variant_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.varscan.indel.vcf")),
        maxdepth=0
    conda: '../envs/variant.yml'    
    shell:
        """
        if [[ $(samtools view {input.bam} | head -n 5) ]]; then
            samtools mpileup -d {params.maxdepth} -f {input.genome} \
                {input.bam} | \
                java -jar {input.varscan} \
                    mpileup2snp \
                    --output-vcf 1 \
                    --min-coverage 1 \
                    --min-reads2 1 \
                    --min-avg-qual 30 \
                    --min-var-freq 0.0 > {params.tmp_snp}

            samtools mpileup -d {params.maxdepth} -f {input.genome} \
                {input.bam} | \
                java -jar {input.varscan} \
                    mpileup2indel \
                    --output-vcf 1 \
                    --min-coverage 100 \
                    --min-reads2 5 \
                    --min-avg-qual 30 \
                    --min-var-freq 0.01 > {params.tmp_indel}

            bcftools view {params.tmp_snp} -Oz -o {params.tmp_snp}.gz
            tabix -f {params.tmp_snp}.gz
            bcftools view {params.tmp_indel} -Oz -o {params.tmp_indel}.gz
            tabix -f {params.tmp_indel}.gz
            bcftools concat -a -o {output} -O v {params.tmp_snp}.gz {params.tmp_indel}.gz
        else
            touch {output}        
        fi
        """


rule lofreq_calling:
    """ Call variants (SNP and indels) with the reference. Using only defaut filtering.
    """
    input: 
        genome=get_genome,
        bam=join(config['split_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.bam"),
        bam_index=join(config['split_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.bam.bai")
    output: join(config['variant_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.lofreq.vcf")
    params: maxdepth=10000
    conda: '../envs/variant.yml'
    threads: config['max_cpu']
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
    input: expand([join(config['variant_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.{caller}.ann.csv")], accession=pd.read_csv(config['samples']['file'])['Run'], aligner=['BWA'], caller=['varscan'])
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


rule raw_variants:
    """
    Calculate the raw varaints using `pileup` from Rsamtools. 
    """
    input:         
        join(config['split_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.bam"),
        join(config['split_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.bam.bai"),
        get_genome
    output: join(config['variant_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.rawvariants.csv")
    conda: '../envs/r.yml'
    script: "../scripts/variant_calling.R"

rule aggregate_raw_variants:
    """
    This rule aggregates all of the raw variants. 
    """
    input: expand([join(config['variant_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.rawvariants.csv")], accession=pd.read_csv(config['samples']['file'])['Run'], aligner=['BWA'])
    output: join(config['variant_dir'], "raw-variants.csv")
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