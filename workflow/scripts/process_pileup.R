## ---------------------------
##
## Script name: `rule process_pileup`
##
## Purpose of script: Process the parsed pileup for analysis
##
## Author: Will Hannon
##
## Date Created: 2020-11-02
##
## Copyright (c) Will Hannon, 2020
##
## Email: wwh22@uw.edu
##
## ---------------------------

require(tidyverse)

## ==== Get file paths from snakemake object ==== ##

pileup.filepath = snakemake@input[[1]] 

## ==== Process the pileup file ==== ##

pileup.df = read_csv(pileup.filepath) %>%
  # Remove insertions and deletions
  select(!starts_with("Ins") & !starts_with("Del")) %>%
  # Separate sample names into parts
  separate(col = Accession, into = c("Accession", "Replicate"), sep = "-") %>%
  # Fix the replicate names
  mutate(Replicate = ifelse(is.na(Replicate), 1, Replicate)) %>%
  # Rename the columns
  rename(REF = "Reference Base", CONS = "Consensus Base", POS = "Position") %>%
  # Determine the total support for each base
  mutate(`A Total` = `A For` + `A Rev`,
         `C Total` = `C For` + `C Rev`,
         `T Total` = `T For` + `T Rev`,
         `G Total` = `G For` + `G Rev`) %>%
  # Pivot the columns
  pivot_longer(cols = c("A For", "A Rev",
                        "C For", "C Rev",
                        "G For", "G Rev",
                        "T For", "T Rev",
                        "A Total", "C Total",
                        "T Total", "G Total"),
               names_to = "ALT", values_to = "Count") %>%
  # Determine the strand
  separate(col = ALT, into = c("ALT", "Strand")) %>%
  # Rename the strands
  mutate(Strand = case_when(Strand == "For" ~ "+",
                            Strand == "Rev" ~  "-",
                            Strand == "Total" ~ "+/-"),
         SNP = paste0(REF, POS, ALT)) %>%
  # Calculate Allele Freq and identify SNPs only
  mutate(AF = Count/Depth, Type = case_when(ALT != REF ~ "SNP",
                                            ALT == REF ~ "REF")) %>%
  # Filter only SNPs with non-zero AF
  filter(Type == "SNP", AF != 0, Strand == "+/-")
  
## ==== Export the results ==== ##

write_csv(pileup.df, snakemake@output[[1]])

## ==== End of script ==== ##









