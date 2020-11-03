### ======= Trim adaptors, poly A tails, and low quality reads ======= ###
#
# Author: Will Hannon 
# Email: wwh22@uw.edu
# Date: 10/30/2020
#

rule trim_adapters_se:
    """
    Fast all-in-one processing of single `fastq` files. 
    Automatic adaptor trimming,
    low-qual base filtering (50% of bases w/ Phred >20), 
    and reporting by HTML and JSON.

    Exclude reads < 50bp and cut poly X enabled (--trim_poly_x).
    """
    input: get_avaliable_fastqs
    output: 
        trimmed=join(config['trim_dir'], "{accession}", "{accession}.trimmed.fastq.gz"),
        html=join(config['qc_dir'], "{accession}", "fastp", "{accession}.fastp.html"),
        json=join(config['qc_dir'], "{accession}", "fastp", "{accession}.fastp.json")
    conda: '../envs/trim.yml'
    shell:
        """ 
        fastp \
            -i {input} \
            -q 20 \
            -u 50 \
            -l 50 \
            --trim_poly_x \
            -o {output.trimmed} \
            --html {output.html} \
            --json {output.json}
        """

rule trim_adapters_pe:
    """
    Fast all-in-one processing of single `fastq` files. 
    Automatic adaptor trimming,
    low-qual base filtering (50% of bases w/ Phred >20), 
    and reporting by HTML and JSON.

    Exclude reads < 50bp and cut poly X enabled (--trim_poly_x).
    """
    input: get_avaliable_fastqs
    output:
        trimmed=[join(config['trim_dir'], "{accession}", "{accession}_1.trimmed.fastq.gz"), 
                 join(config['trim_dir'], "{accession}", "{accession}_2.trimmed.fastq.gz")],
        html=join(config['qc_dir'], "{accession}", "fastp", "{accession}.fastp.html"),
        json=join(config['qc_dir'], "{accession}", "fastp", "{accession}.fastp.json")
    conda: '../envs/trim.yml'
    shell:
        """ 
        fastp \
            -q 20 \
            -u 50 \
            -l 50 \
            --trim_poly_x \
            --detect_adapter_for_pe \
            -o {output.trimmed[0]} \
            -O {output.trimmed[1]} \
            --html {output.html} \
            --json {output.json} 
        """
