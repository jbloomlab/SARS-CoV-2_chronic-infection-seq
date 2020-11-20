### ======= Identify variants by calculating and parsing a pileup file ======= ###
#
# Author: Will Hannon 
# Email: wwh22@uw.edu
# Date: 11/02/2020
#

rule calculate_pileup:
    """ 
    Calculate the pileup of bases at every position in virus genome.
    Only calculates bases with Phred scaled quality score higher than 30.
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
    """ Parse the pileup file generated by samtools.
    """
    input: join(config['pileup_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.mpileup.csv")
    output: join(config['pileup_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.mpileup.processed.csv")
    conda: '../envs/r.yml'
    script: '../scripts/process_pileup.R'

rule aggregate_pileup:
    """
    """
    input: expand(join(config['pileup_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.mpileup.processed.csv"), accession=pd.read_csv(config['samples']['file'])['Run'], aligner=['BWA'])
    output: join(config['pileup_dir'], "raw-variants.csv")
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


