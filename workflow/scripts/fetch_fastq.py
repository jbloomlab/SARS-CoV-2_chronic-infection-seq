__author__ = "Will Hannon"
__copyright__ = "Copyright 2020, Will Hannon"
__email__ = "wwh22@uw.edu"
__license__ = "MIT"

# Import necessary packages - should be installed in the main environment running the pipeline
import os
import pandas as pd
from snakemake.shell import shell

# Determine if the run is local or on deposited on the SRA from ./config/samples.csv
source = pd.read_csv(snakemake.config['samples']['file']).set_index('Run').at[snakemake.wildcards.accession, 'Source']

# If the sample needs to be downloaded from the NCBI SRA
if source == "SRA":

    # Attempt the download at least 3 times
    attempt = 0

    # Attempt at three times
    while attempt < 3:

    # Attempt the download
        try:
            shell("fasterq-dump {snakemake.wildcards.accession} --outdir {snakemake.params} --temp {snakemake.params} --threads {snakemake.threads} -f; gzip -c {snakemake.params}/*_1.fastq > {snakemake.output[0]}; gzip -c {snakemake.params}/*_2.fastq > {snakemake.output[1]}")  
            
            break 

        # If the download fails increment the attempt and try again
        except: 
            # Increment the number of attempts left 
            attempt += 1

# If the run exists as a local file
elif source == "local":

    # Get the path to the files
    path = pd.read_csv(snakemake.config['samples']['file']).set_index('Run').at[snakemake.wildcards.accession, 'Path'].strip('][').split(', ')

    # Check if the run is paired or single-ended 
    layout = pd.read_csv(snakemake.config['samples']['file']).set_index('Run').at[snakemake.wildcards.accession, 'LibraryLayout']

    # Two approaches depending on the layout
    if layout == "PAIRED":

        # Make sure there are two paths 
        assert len(path) == 2

        # Check if runs are zipped
        if path[0].endswith(".gz"):
            # Make directory and move into the correct direcotry with naming scheme 
            shell("cp {path[0]} {snakemake.output[0]}; cp {path[1]} {snakemake.output[1]}")  
        else:
            # Zip into a new directory with the correct naming scheme
            shell("gzip -c {path[0]} > {snakemake.output[0]}; gzip -c {path[1]} > {snakemake.output[1]}")

    # Otherwise the runs are assumed to be single ended
    else: 

        # Make sure there is only one path
        assert len(path) == 1

        # Check if runs are zipped
        if path[0].endswith(".gz"):
            # Make directory and move into the correct direcotry with naming scheme 
            shell("cp {path[0]} {snakemake.output[0]}")  
        else:
            # Zip into a new directory with the correct naming scheme
            shell("gzip -c {path[0]} > {snakemake.output[0]}")

