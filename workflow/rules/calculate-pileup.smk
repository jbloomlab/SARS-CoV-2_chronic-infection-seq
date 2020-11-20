### ======= Identify variants by calculating and parsing a pileup file ======= ###
#
# Author: Will Hannon 
# Email: wwh22@uw.edu
# Date: 11/02/2020
#

rule calculate_samtools_pileup:
    """ 
    Calculate the pileup of bases at every position in virus genome.
    Only calculates bases with Phred scaled quality score higher than 30.
    This is prased by a custom python script and used by Varscan2 to vallidate
    observed SNPs with multiple methods. 
    """
    input: bam=join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.sorted.marked.bam"),
           bai=join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.sorted.marked.bam.bai"),        
           genome=get_genome
    output: join(config['pileup_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.mpileup.txt")
    params: score=config['BQ']
    conda: '../envs/samtools.yml'
    shell: "samtools mpileup -d 0 -E --excl-flags UNMAP,SECONDARY,QCFAIL -q {params.score} -Q {params.score} -f {input.genome} {input.bam} -O -s --reverse-del -a -o {output}"


rule parse_pileup:
    """ Parse the pileup file generated by samtools.
    """
    input: join(config['pileup_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.mpileup.txt")
    output: join(config['pileup_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.mpileup.csv")
    script: '../scripts/process_pileup.py'


rule process_pileup:
    """ Process the parsed variants from the above python script, aggregating variants for easy plotting in R.
    """
    input: join(config['pileup_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.mpileup.csv")
    output: join(config['pileup_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.mpileup.processed.csv")
    conda: '../envs/r.yml'
    script: '../scripts/process_pileup.R'


rule aggregate_samtools_pileup:
    """ Aggregate variants called from samtools mpileup + parsing + processing into a single table. 
    """
    input: expand(join(config['pileup_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.mpileup.processed.csv"), accession=pd.read_csv(config['samples']['file'])['Run'], aligner=['BWA'])
    output: join(config['pileup_dir'], "samtools_variants.csv")
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


rule calculate_pysam_pileup:
    """ Calculate pileup statistics and process with python/pysam.
    """
    input: bam=join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.sorted.marked.bam"),
           bai=join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.sorted.marked.bam.bai"),        
           genome=get_genome
    output: join(config['pileup_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.pysam.pileup.csv")
    params: score=config['BQ']
    conda: "../envs/pysam.yml"
    script: "../scripts/pysam_pileup.py"

    
rule aggregate_pysam_pileup:
    """ Aggregate variants called with python and pysam.. 
    """
    input: expand(join(config['pileup_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.pysam.pileup.csv"), accession=pd.read_csv(config['samples']['file'])['Run'], aligner=['BWA'])
    output: join(config['pileup_dir'], "pysam_variants.csv")
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
