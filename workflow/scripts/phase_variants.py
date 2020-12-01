__author__ = "Will Hannon"
__copyright__ = "Copyright 2020 Will Hannon"
__email__ = "wwh22@uw.edu"
__license__ = "MIT"

## === Import Libraries === ##
import pysam #count variant alleles from BAM
import pandas as pd #data frames
import numpy as np #arrays
import os #interacting with files
from Bio import SeqIO #reading fasta format
import re #regular expressions

## === Functions === ## 

def check_read(read):
    """
    Helper function to decide what reads should
    be keep when parsing alignment file with `pysam`. 

    Parameters
    ----------
    read : AlignedSegment
        read from alignment file parsed with `pysam`.

    Returns
    -------
    bool
        True/False if read should be included
        
    """
    # Exclude Quality Failures
    if read.is_qcfail:
        return False
    # Exclude Secondary Mappings
    if read.is_secondary:
        return False
    # Exclude Unmapped Reads
    if read.is_unmapped:
        return False
    else:
        return True


def get_day(path):
    """
    If the samples are downloaded from the SRA, use the
    SRA accession to get the day the sample was collected on. 

    Parameters
    ----------
    path : str
        path to the sample

    Returns
    -------
    int
        the day the sample was collected

    """
    
    day_mapping = { 
    'SRR13160722':152, 
    'SRR13160723':146, 
    'SRR13160724':143, 
    'SRR13160725':130, 
    'SRR13160726':128, 
    'SRR13160727':81, 
    'SRR13160728':75, 
    'SRR13160729':25, 
    'SRR13160730':18
    } 
    
    basename = os.path.basename(path)

    regex = re.compile("(?P<run>SRR\d+)\.BWA")
    
    m = regex.match(basename)

    return day_mapping[str(m.group('run'))]
    

def phase_variants(aa_change_list, joined_spike_count_df, filepath, contig = "NC_045512.2"):
    """
    Uses a list of variant alleles, a data-frame of SNPs corresponding to BAM, and an alignment file to make a N (SNPs) X P (Reads)
    dataframe with the read support for every obseved phase of varaints. 

    Parameters
    ----------
    joined_spike_count_df : Pandas.DataFrame
        Dataframe with SNPs summarized over a region.

    filepath : str
        path to the bam file to be parsed
        
    aa_change_list : list
        a list of alleles to phase
   
    conitg : str
        name of the chromosome, in this case the whole reference.


    Returns
    -------
    Pandas.DataFrame
       A Pandas.DataFrame of 0, 1, or NaN, describing read support for a given phase.
        
    """
    
    ## ===== Format inputs and get SNP list ===== ##

    # Get the SNPs to query over specific day
    SNPs_df = joined_spike_count_df[joined_spike_count_df.AA_CHANGE.isin(aa_change_list)]
    SNPs_set = set(pair for pair in zip(SNPs_df.POS, SNPs_df.ALT))                   
    
    ## ===== Move through the BAM file and get haplotypes ===== ##
    
    # Save the haplotypes in a dictionary
    haplotype_dict = {}

    # Get the start and stop
    start = sorted([pos for pos, alt in SNPs_set])[0]
    stop = sorted([pos for pos, alt in SNPs_set])[-1]

    # Open alignment with pysam
    with pysam.AlignmentFile(filepath, "rb") as bamfile:

        # Get the pileup column for a specific region
        for pileupcolumn in bamfile.pileup(contig, start = start, stop = stop, stepper = 'nofilter'):

            # Check the if the position has a target SNP (converted to 0-indexed)
            if pileupcolumn.pos in [pos-1 for pos, alt in SNPs_set]:

                # Iterate over every alignment in the pileup column. 
                for pileupread in pileupcolumn.pileups:

                    # Check if the read is valid and can be parsed
                    if check_read(pileupread.alignment) and not pileupread.is_del and not pileupread.is_refskip:

                        # Save the query name
                        qname = pileupread.alignment.query_name
                        
                        # Save the 1-indexed position
                        pos = pileupcolumn.pos + 1

                        # Save the base at that position in the read
                        alt = pileupread.alignment.query_sequence[pileupread.query_position]

                        # Check if this read has the SNP or not
                        if (pos, alt) in SNPs_set:
                            phase = 1 # The read has the SNP

                        else:
                            phase = 0 # The read doesn't have the SNP

                        # Add the readname to the dictionary  
                        if qname in haplotype_dict.keys():
                            haplotype_dict[qname].append((pos, phase, alt)) 
                        else:
                            haplotype_dict[qname] = [(pos, phase, alt)]
                     
    ## ===== Convert the haplotype dictionary to a dataframe filling in missings ===== ##
    
    # Save the haplotypes with missing values (not covered by reads) filled
    completed_haplotype_dict = {}

    # Iterate over the haplotype dictionary created above
    for qname, alleles in haplotype_dict.items():
        
        # Make a dictionary to fill with observed SNPs for each read
        aa_dict = {tup:"-" for tup in SNPs_set}
        
        # Iterate over all observed alleles
        for allele in alleles:
            
            # If it's 1, then an allele was observed
            if allele[1] == 1:
                # Add this observation based on the key
                aa_dict[(allele[0], allele[2])] = 1
            # If it's a wild type allele 
            elif allele[1] == 0:
                # Check every possible wt allele
                for key in aa_dict.keys():
                    # Check by position
                    if key[0] == allele[0]:
                        # If it's the same positon add a 0
                        aa_dict[key] = 0

        # Add this populated dictionary to the haplotype dictionary
        completed_haplotype_dict[qname] = aa_dict

    
    # Convert this dicitonary into a dataframe and take the transpose
    haplotype_df = pd.DataFrame(completed_haplotype_dict).T
    
    ## ===== Format the dataframe ===== ##
    
    # * FIX THE BUG WITH ALLELES IN SAME BASE *
    
    # Dict to format the column names
    aa_conversion = {pair:aa for pair, aa in zip(zip(SNPs_df.POS, SNPs_df.ALT),SNPs_df.AA_CHANGE)}
    
    # Flatten the column indicies
    haplotype_df.columns = [aa_conversion[(pos,alt)] for pos,alt in haplotype_df.columns.values]

    # Replace the '-' with 'NaN'
    haplotype_df = haplotype_df.replace('-', np.nan)
    
    # Return the phased SNP df
    return haplotype_df


## === Main === ## 

# Input path to SNP df
inputpath_csv = str(snakemake.input.csv)
inputpath_bam = str(snakemake.input.bam)

# Output path to csv file
outpath = str(snakemake.output)

# Name of sample
accession = os.path.basename(inputpath_bam)

# Get the day of the sample 
if accession.startswith("Day"):

    regex = re.compile("Day(?P<Day>\d+)_")
    m = regex.match(accession)
    day = int(m.group('Day'))

else:

   day = get_day(inputpath)



# Read in the join_calling df 
input_snp_df = pd.read_csv(inputpath_csv)  

# List of SNPs to phase
snp_list = ["N440D", "E484A", "F486I", "F486L", "Y489H", "Q493K"]

# Phase the SNPs
haplotypes_df = phase_variants(snp_list, input_snp_df, inputpath_bam, contig = "NC_045512.2")

# Add accession and day to data frame
haplotypes_df.insert(0, 'ACCESSION', accession)
haplotypes_df.insert(1, 'DAY', day)

# Write out the results
haplotypes_df.to_csv(outpath, index=False)
