## ---------------------------
##
## Script name: `rule plot_rbd_all_mutations`
##
## Purpose of script: Plot the all of the mutations in RBD for every timepoint day 18 to day 152.
##
## Author: Will Hannon
##
## Date Created: 2020-11-19
##
## Copyright (c) Will Hannon, 2020
## Email: wwh22@uw.edu
##
## ---------------------------

require(tidyverse)
require(grid)

## ---------------------------

## ==== Path to pileup data from pysam script ==== ##
# snakemake@input[[1]] 
pysam.data = "../../config/pysam_pileup.csv" # Make the snakemake input file

## ==== Add matplotlib Tableau color scheme ==== ##
tab_colors = c( "#7f7f7f", "#1f77b4", "#ff7f0e")

## ==== Read in the pysam/pileup data ==== ##
pysam.df = read_csv(pysam.data)

## ==== Get a list of the amino acid positions to plot ==== ##
sites.to.plot = pysam.df %>% 
  # Must have at least one missense mutation
  filter(MISSENSE == "TRUE") %>% 
  # Must be in the RBD
  filter(RBD == "TRUE") %>% 
  # At least one mutation over 10%
  filter(AF >= 0.1) %>% 
  # Get a list of the relevant protein positions
  pull(PROT_POS) %>% 
  # only the unqiue sites are necessary
  unique()





