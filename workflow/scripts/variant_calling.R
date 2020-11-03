## ---------------------------
##
## Script name: Variant Calling
##
## Purpose of script: Simple Variant Calling with the pileup command from Rsamtools
##
## Author: William Hannon 
##
## Date Created: 2020-10-03
##
## Copyright (c) Will Hannon, 2020
##
## Email: wwh22@uw.edu
##
## ---------------------------

require(tidyverse)
require(Rsamtools)
require(seqinr)

## ---------------------------

## ====================== Input Files ====================== ##

bamfile = snakemake@input[[1]] 
indexfile = snakemake@input[[2]] 
genome = snakemake@input[[3]]

# bamfile = "../../results/split/BWA/SRR11939953-2/SRR11939953-2.BWA.virus.bam"
# indexfile = "../../results/split/BWA/SRR11939953-2/SRR11939953-2.BWA.virus.bam.bai"
# genome = "../../config/ref/SARS2.fasta"

## ====================== Read Fasta ====================== ##

# Generate a dataframe of reference alleles from the fasta reference
ref.df = data.frame((strsplit(toupper(read.fasta(file = genome, as.string = TRUE, seqtype = "DNA")[[1]][1]), split = ""))) %>% 
  mutate(pos = 1:nrow(.))

# Rename the columns
names(ref.df ) = c("REF", "pos")

## ====================== Pileup Params ====================== ##

pileup.params = PileupParam(max_depth=1,
            min_base_quality=0,
            min_mapq=0,
            min_nucleotide_depth=1,
            min_minor_allele_depth=0, 
            distinguish_strands=TRUE,
            distinguish_nucleotides=TRUE,
            ignore_query_Ns=TRUE, 
            include_deletions=FALSE, 
            include_insertions=FALSE, 
            left_bins=NULL, 
            query_bins=NULL, 
            cycle_bins=NULL)

## ====================== Generate Pileup ====================== ##

pileup.df = pileup(bamfile, index = indexfile, pileupParam = pileup.params) %>% 
  group_by(pos, strand, nucleotide) %>%
  summarize(n = sum(count)) %>% 
  pivot_wider(names_from = nucleotide, values_from = n) %>%
  mutate_at(c("A", "T", "C", "G"), function(x) ifelse(is.na(x), 0, x)) %>%  
  mutate(sum = sum(A, `T`, G, C)) %>% 
  left_join(.,ref.df, by = "pos")


result.df = cbind(as.data.frame(pileup.df[,3:6]/pileup.df$sum), as.data.frame(pileup.df[,c(1,2,7,8)])) %>% 
  mutate_at(c("A", "T", "C", "G"), function(x) ifelse(is.na(x), 0, x)) %>% 
  pivot_longer(cols = c("A", "T", "C", "G"), names_to = "ALT", values_to = "AF") %>% 
  mutate(is_reference = ifelse(REF == ALT, "yes", "no")) %>% 
  filter(is_reference == "no", AF != 0) %>% 
  select(!is_reference, count = "sum") %>% 
  mutate(Run = snakemake@wildcards[['accession']], Aligner = snakemake@wildcards[['aligner']]) 
  
## ====================== Export Results ====================== ##            
  
write_csv(result.df, snakemake@output[[1]])






