### ======= Calculate the coverage in bins over the genome ======= ###
#
# Author: Will Hannon 
# Email: wwh22@uw.edu
# Date: 10/31/2020
#

rule bedtools_coverage:
    """ Calculate the average read coverage over bins in the genome.
    """
    input: bam=join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.sorted.marked.bam"),
           bai=join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.sorted.marked.bam.bai"),
           genome=get_genome
    output: join(config['coverage_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.bedgraph")
    params: binsize = config['bin_size']
    conda: '../envs/samtools.yml'
    shell:
        """
        awk 'FS=OFS="\t"{{print $1, 0, $2}}' {input.genome}.fai \
            | bedtools makewindows -b - -w {params.binsize} \
            | bedtools coverage -a stdin -b {input.bam} -mean \
            > {output}
            
        sed -i "s/$/\t{wildcards.accession}/" {output} 
        """


rule merge_coverage:
    """ Merge the read coverage tables for all of the accessions into a single file.
    """
    input: expand(join(config['coverage_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.bedgraph"), accession=pd.read_csv(config['samples']['file'])['Run'], aligner=['BWA'])
    output: bedgraph=join(config['coverage_dir'], "merged.bedgraph"),
            header=temp(join(config['coverage_dir'], "merged.bedgraph.tmp"))
    shell:
        """
        cat {input} > {output.header}

        awk 'BEGIN{{print "Name\tStart\tStop\tDepth\tAccession"}}1' {output.header} > {output.bedgraph}
        """


rule samtools_depth:
    """ Calculate the depth over each position filtering by the phred base score. 
    """
    input: bam=join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.sorted.marked.bam"),
           bai=join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.sorted.marked.bam.bai")
    output: join(config['coverage_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.depth")
    params: score=config['BQ'],
            binsize=config['bin_size']
    conda: '../envs/samtools.yml'
    shell: 
        """
        samtools depth -a -m 0 -q {params.score} -g DUP {input.bam} \
            | awk '{{sum+=$3}} NR%{params.binsize}==0 {{print $2-{params.binsize} "\t" $2 "\t" sum/{params.binsize} "\t"; sum=0}}' - \
            > {output}

        sed -i "s/$/\t{wildcards.accession}/" {output} 
        """


rule merge_depth:
    """ Merge the samtools depth tables for all of the accessions into a single file.
    """
    input: expand(join(config['coverage_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.depth"), accession=pd.read_csv(config['samples']['file'])['Run'], aligner=['BWA'])
    output: depth=join(config['coverage_dir'], "merged.depth"),
            header=temp(join(config['coverage_dir'], "merged.depth.tmp"))
    shell:
        """
        cat {input} > {output.header}

        awk 'BEGIN{{print "Start\tStop\tDepth\tAccession"}}1' {output.header} > {output.depth}
        """

rule average_depth:
    input: bam=join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.sorted.marked.bam"),
           bai=join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.sorted.marked.bam.bai")
    output: join(config['coverage_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.average.depth")
    params: score=config['BQ'], 
            min_coverage=config['min_coverage']
    conda: '../envs/samtools.yml'
    shell:
        """
        samtools depth -a -q {params.score} -g DUP {input.bam} \
            | awk '{{ if ($3 >= {params.min_coverage}) sum+=1}} END {{print sum/NR*100}}' - \
            > {output}

        sed -i "s/$/\t{wildcards.accession}/" {output} 
        """

rule merge_average_depth:
    """ Merge the samtools depth tables for all of the accessions into a single file.
    """
    input: expand(join(config['coverage_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.average.depth"), accession=pd.read_csv(config['samples']['file'])['Run'], aligner=['BWA'])
    output: depth=join(config['coverage_dir'], "merged.average.depth"),
            header=temp(join(config['coverage_dir'], "merged.average.depth.tmp"))
    shell:
        """
        cat {input} > {output.header}

        awk 'BEGIN{{print "Percent\tAccession"}}1' {output.header} > {output.depth}
        """

rule coverage_stats: 
    input: expand(join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.sorted.marked.bam"), accession=pd.read_csv(config['samples']['file'])['Run'], aligner=['BWA'])
    output: join(config['coverage_dir'], "coverage.stats")
    params: score=config['BQ']
    conda: '../envs/samtools.yml'
    shell:
        """
        echo "rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\tfilename" > {output}
        
        list="{input}"

        for f in $list; do
            echo -e "$(samtools coverage -H -Q {params.score} --excl-flags UNMAP,SECONDARY,QCFAIL $f)\t$f" >> {output}
        done
        """

