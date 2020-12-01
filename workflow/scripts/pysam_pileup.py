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

def translate(codon):
    """
    Translate a three letter DNA string into 
    a one letter amino acid code. 

    Parameters
    ----------
    codon : str
        three letter DNA sequence

    Returns
    -------
    str
        one letter amino acid code

    Raises
    ------
    AssertionError
        error if codon sequence is invalid
        
    """
    
    table = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W', 
    } 
    
    assert codon in table.keys(), "Not a valid codon sequence."
    
    return table[codon]

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


def mutate(codon, alt, index):
    """
    Replace (mutate) a base in a codon with an 
    alternate base. 
        
    Parameters
    ----------
    codon : str
        three letter DNA sequence
        
    alt : str
        alternative base
        
    index : int
        index of the alt base in codon (0|1|2). 

    Returns
    -------
    str
        codon with alternative base

    Raises
    ------
    AssertionError
        error if index is not valid (0|1|2)
        
    AssertionError
        error if base is not valid (A|T|C|G)
        
    """
    
    assert index in [0,1,2], "Not a valid index."
    
    assert alt in ["A", "T", "C", "G"], "Not a valid base."
    
    return "".join([alt if i == index else b for i,b in enumerate(codon)])


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
        
    
def build_af_df(filepath, 
                callback_function=check_read, 
                ref = "NC_045512.2", 
                ref_path = "../../config/ref/SARS2.fa", 
                minimum_AF = 0.01, 
                minimum_qual = 25):
    """
    Read in BAM file and convert to a dataframe containing the frequency of 
    of any bases present at a given position in the reference genome using the
    `pysam` command `count_coverage`. 

    Parameters
    ----------
    filepath : str
        path to the bam file to be parsed
        
    callback_function : function
        function that decides which reads to keep/exclude
        
    ref : str
        name of the contig to count coverage over
        
    ref_path : str
        path to the reference genome as fasta
        
    minimum_AF : float
        minimum allele frequency to include allele in dataframe
        
    minimum_qual : int
        minimum QUAL score to count read at a position.

    Returns
    -------
    Pandas.DataFrame
       Data Frame containing the bases represented at each positon in the genome
        
    """

    # Open alignment with pysam
    with pysam.AlignmentFile(filepath, "rb") as bamfile:
        
        # Get a dataframe of the counts
        count_df = pd.DataFrame.from_dict({base:counts for base, counts in zip("ACGT", bamfile.count_coverage(contig = ref, read_callback=callback_function, quality_threshold=minimum_qual))})
        
        # Add the depth at each position
        count_df['DP'] = count_df.sum(axis = 1)
        
        # Add the position 
        count_df['POS'] = count_df.index + 1
        
        # Add the reference allele
        count_df['REF'] = [base.upper() for base in list(SeqIO.parse(ref_path, "fasta"))[0].seq]
        
        # convert counts to frequency 
        count_df.iloc[:,0:4] = count_df.iloc[:,0:4].div(count_df.DP, axis = 0)
        
        # handle any NaNs created by dividing by 0 coverage
        count_df = count_df.fillna(0)
        
        # Melt the data frame to a longer ('tidy') form
        count_df = pd.melt(count_df, 
                           id_vars=['POS', 'DP', 'REF'],
                           value_vars=[base for base in 'ATGC'],
                           value_name='AF',
                           var_name='ALT')
        
        # Filter out anything less than minimum allele freq
        count_df = count_df[count_df['AF'] >= minimum_AF]
        
        # TRUE/FALSE if it's a SNP
        count_df['SNP'] = np.where(count_df['ALT'] != count_df['REF'], True, False)
        
        # Is a base consensus or not.
        count_df['CONS'] = count_df['AF'].map(lambda x: x >= 0.5)
    
        # Sort by position and reset the index
        return count_df.sort_values('POS').reset_index(drop=True)
     
    
def annotate_coding_change_in_spike(count_df, ref_genome):
    """
    Annotate the protein coding changes of mutations in the Spike gene 
    of SARS-CoV-2 (i.e. in the format `D614G`). 

    Parameters
    ----------
    count_df : Pandas.DataFrame
       Data Frame containing the bases represented at each positon in the genome
        
    ref_genome : str
        reference genome for all of SARS-CoV-2
        

    Returns
    -------
    Pandas.DataFrame
       Data Frame containing annotation for all position over Spike
        
    """
    # Lists to hold residue changes and positions
    residue_change_list = []
    residue_position_list = []
    resiude_wt_list = []
    residue_mut_list = []
    
    # Spike Sequence by coordinates of spike 0-indexed
    spike_sequence = "".join(ref_genome[21562:25384])

    # Pull out only the SNPs in Spike 1-indexed
    ## spike_SNP_df = count_df[(count_df.SNP) & (count_df.POS >= 21563) & (count_df.POS <= 25384)]
    
    spike_SNP_df = count_df[(count_df.POS >= 21563) & (count_df.POS <= 25384)]

    # Iterate over every row in the data frame by alternative allele and position
    for alt, pos in zip(spike_SNP_df.ALT.tolist(),spike_SNP_df.POS.tolist()):
        
        # Get the position in the genome relative to Spike
        position_in_spike = pos - 21562
        
        ## === Check the position of the SNP in the codon == ##

        # First base in codon
        if position_in_spike % 3 == 1:
            
            # Save codon position
            position_in_codon = 1

            # Get the sequence of the wt codon by index
            wt_codon = spike_sequence[position_in_spike-1:position_in_spike+2]

            # Translate wt codon to residue
            wt_aa = translate(wt_codon)

            # Get the sequence of the mut codon by index
            mut_codon = mutate(wt_codon, alt, position_in_codon-1)

            # Translate mut codon to residue
            mut_aa = translate(mut_codon)

            # Calculate the position of the residue in Spike
            residue_position_in_spike = int((position_in_spike + 2) / 3)
            
            # Build the list of mutations and positions
            residue_change_list.append(f"{wt_aa}{residue_position_in_spike}{mut_aa}")
            residue_position_list.append(residue_position_in_spike)
            resiude_wt_list.append(wt_aa)
            residue_mut_list.append(mut_aa)


        # Second base in codon
        elif position_in_spike % 3 == 2:
            
            # Save codon position
            position_in_codon = 2
            
            # Get the sequence of the wt codon by index
            wt_codon = spike_sequence[position_in_spike-2:position_in_spike+1]

            # Translate wt codon to residue            
            wt_aa = translate(wt_codon)

            # Get the sequence of the mut codon by index           
            mut_codon = mutate(wt_codon, alt, position_in_codon-1)
            
            # Translate mut codon to residue
            mut_aa = translate(mut_codon)

            # Calculate the position of the residue in Spike
            residue_position_in_spike = int((position_in_spike + 1) / 3)

            # Build the list of mutations and positions
            residue_change_list.append(f"{wt_aa}{residue_position_in_spike}{mut_aa}")
            residue_position_list.append(residue_position_in_spike)
            resiude_wt_list.append(wt_aa)
            residue_mut_list.append(mut_aa)

       
        # Third base in codon
        elif position_in_spike % 3 == 0:

            # Save codon position
            position_in_codon = 3
            
            # Get the sequence of the wt codon by index
            wt_codon = spike_sequence[position_in_spike-3:position_in_spike]

            # Translate wt codon to residue            
            wt_aa = translate(wt_codon)
            
            # Get the sequence of the mut codon by index  
            mut_codon = mutate(wt_codon, alt, position_in_codon-1)

            # Translate mut codon to residue
            mut_aa = translate(mut_codon)

            # Calculate the position of the residue in Spike
            residue_position_in_spike = int((position_in_spike) / 3)

            # Build the list of mutations and positions
            residue_change_list.append(f"{wt_aa}{residue_position_in_spike}{mut_aa}")
            residue_position_list.append(residue_position_in_spike)
            resiude_wt_list.append(wt_aa)
            residue_mut_list.append(mut_aa)

    return spike_SNP_df.assign(AA_CHANGE = residue_change_list,
                               PROT_POS = residue_position_list,
                               WT_AA = resiude_wt_list,
                               MUT_AA = residue_mut_list)

## ===== Main ===== ##

# Input path to BAM
inputpath = str(snakemake.input.bam)

# Output path to csv file
outpath = str(snakemake.output)

# Reference genome
ref_genome = [base.upper() for base in list(SeqIO.parse(str(snakemake.input.genome), "fasta"))[0].seq]

# Name of sample
accession = os.path.basename(inputpath)

if accession.startswith("Day"):

    # Get the day of the sample 
    regex = re.compile("Day(?P<Day>\d+)_")
    m = regex.match(accession)
    day = int(m.group('Day'))

else:

   day = get_day(inputpath)

# Build the count df
count_df = build_af_df(inputpath, 
            callback_function=check_read, 
            ref = "NC_045512.2", 
            ref_path = str(snakemake.input.genome), 
            minimum_AF = 0.01, 
            minimum_qual = int(snakemake.params.score))

# Add accession and day to data frame
count_df.insert(0, 'ACCESSION', accession)
count_df.insert(1, 'DAY', day)

# Annotate the count df   
joined_spike_count_df = annotate_coding_change_in_spike(count_df, ref_genome)

# TRUE/FALSE if there is a missense mutation
joined_spike_count_df['MISSENSE'] = np.where(joined_spike_count_df['WT_AA'] != joined_spike_count_df['MUT_AA'], True, False)

# TRUE/FALSE if the mutation is in the RBD
joined_spike_count_df['RBD'] = np.where((joined_spike_count_df['PROT_POS'] >= 331) & (joined_spike_count_df['PROT_POS'] <= 531), True, False)

# Write this dataframe into R for plotting and testing.
joined_spike_count_df.to_csv(outpath, index=False)

## ===== Main ===== ##
