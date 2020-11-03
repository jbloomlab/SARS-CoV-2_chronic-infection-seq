### ======= Python utilities for running the pipeline ======= ###
#
# Author: Will Hannon 
# Email: wwh22@uw.edu
# Date: 10/30/2020
#

def single_ended(Run):
    """
    Returns `True` if the Library Layout is single-end reads.
    """
    # Get the Library Layout
    layout = pd.read_csv(config['samples']['file']).set_index('Run').at[Run, 'LibraryLayout']
    # Check the layout
    if layout == "SINGLE":
        return True
    else:
        return False


def get_avaliable_fastqs(wildcards):
    """
    This function fills in the avaliable fastqs depending on 
    the library layout for the rule `trim_adapters_(pe | se)`. 
    """
    # If the layout is single-ended.
    if single_ended(wildcards.accession):
        # Return the target files.
        return expand(join(config["fastq_dir"], "{accession}", "{accession}.fastq.gz"), accession=wildcards.accession)
    # Otherwise the layout is assumed to be paired-ended. 
    return expand([join(config["fastq_dir"], "{accession}", "{accession}_1.fastq.gz"),
                   join(config["fastq_dir"], "{accession}", "{accession}_2.fastq.gz")], accession=wildcards.accession)


def get_avaliable_trimmed_fastqs(wildcards):
    """
    This function fills in the avaliable trimmed fastqs 
    depending on the library layout for the rule `filter_reads_(pe | se)`. 
    """
    # If the layout is single-ended.
    if single_ended(wildcards.accession):
        # Return the target files.
        return expand(join(config["trim_dir"], "{accession}", "{accession}.trimmed.fastq.gz"), accession=wildcards.accession)
    # Otherwise the layout is assumed to be paired-ended. 
    return expand([join(config["trim_dir"], "{accession}", "{accession}_1.trimmed.fastq.gz"),
                   join(config["trim_dir"], "{accession}", "{accession}_2.trimmed.fastq.gz")], accession=wildcards.accession)
     

def get_avaliable_filtered_fastqs(wildcards):
    """
    This function fills in the avaliable filtered fastq's 
    depending on the library layout. 
    """
    # If the layout is single-ended.
    if single_ended(wildcards.accession):
        # Return the target files.
        return expand(join(config["filter_dir"], "{accession}", "{accession}.filtered.fastq.gz"), accession=wildcards.accession)
    # Otherwise the layout is assumed to be paired-ended. 
    return expand([join(config["filter_dir"], "{accession}", "{accession}_1.filtered.fastq.gz"),
                   join(config["filter_dir"], "{accession}", "{accession}_2.filtered.fastq.gz")], accession=wildcards.accession)
     


def get_hybrid_genome(wildcards):
    """
    This function determine which genomes to use or generate
    for alignment of a given run using the `samples.csv` file.
    """
    # Read the metadata into Pandas dataframe
    samples_df = pd.read_csv(config['samples']['file'])
    # Get the viral genome for a given accession
    virusname = samples_df.loc[samples_df.Run == wildcards.accession, ["Virus"]].values.flatten().tolist()[0]
    # Get the host genome for a given accession
    hostname = samples_df.loc[samples_df.Run == wildcards.accession, ["Host"]].values.flatten().tolist()[0]
    # Handle sequences that do not need to be mapped to host\
    if hostname == 'none':
        # Virus only genome - just downloaded 
        return expand(join(config['index_dir']['bwa'], '{virus}.fa'), virus=virusname)
    # Virus and host genome - downloaded and concatenated
    return expand(join(config['index_dir']['bwa'], '{virus}.{host}.fa'), virus=virusname, host=hostname)


def get_genome(wildcards, index = "samtools"):
    """ Function to get the correct genome for variant calling, indexed with samtools. 
    """
    # Read the metadata into Pandas dataframe
    samples_df = pd.read_csv(config['samples']['file'])
    # Get the viral genome for a given accession
    virusname = samples_df.loc[samples_df.Run == wildcards.accession, ["Virus"]].values.flatten().tolist()[0]
    # Return the viral genome file
    if index == "samtools": 
        return expand(join(config['index_dir']['samtools'], '{virus}.fa'), virus=virusname)
    elif index == "BWA":
        return expand(join(config['index_dir']['bwa'], '{virus}.fa'), virus=virusname)


def get_organism(wildcards): 
    """ Determine if organism is host or virus depending on wildcard.
    """
    # Read in the sample dataframe
    samples_df = pd.read_csv(config['samples']['file'])
    # Get the host and virus names as a list
    virus = samples_df.Virus.tolist()
    host = samples_df.Host.tolist()
    # Figure out if the organism is host or virus 
    if wildcards.organism in host:
        # Invert grep
        return "-v"
    # Don't invert grep
    return ""


def get_filtered_bams(wildcards):
    """ 
    This function determines which filtered bam files to generate
    based on the accession.
    """
    # Read the metadata into Pandas dataframe
    samples_df = pd.read_csv(config['samples']['file'])
    # Get the viral genome for a given accession
    virusname = samples_df.loc[samples_df.Run == wildcards.accession, ["Virus"]].values.flatten().tolist()[0]
    # Get the host genome for a given accession
    hostname = samples_df.loc[samples_df.Run == wildcards.accession, ["Host"]].values.flatten().tolist()[0]
    # Handle sequences that do not need to be mapped to host
    if hostname == 'none':

        return expand(join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.{organism}.bam"), organism=virusname)

    return expand(join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.{organism}.bam"), organism=[virusname, hostname])

