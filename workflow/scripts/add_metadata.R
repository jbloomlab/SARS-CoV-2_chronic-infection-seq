## ---------------------------
##
## Script name: `rule add_metadata`
##
## Purpose of script: read in file with given accession and add in the requisite metadata
##
## Author: Will Hannon
##
## Date Created: 2020-06-26
##
## Copyright (c) Will Hannon, 2020
## Email: wwh22@uw.edu
##
## ---------------------------

require(tidyverse)

## ---------------------------

## ==== Get file paths from snakemake object ==== ##

vcf.filepath = snakemake@input[[1]] 
metadata.filepath = snakemake@params[[1]] 

## ==== Boolean to see if file is empty and there are no columns ==== ##
info = file.info(vcf.filepath)
empty = (info$size == 0)

## ==== Import Metadata, same table used to import samples from `Run Selector w/ minimal changes` ==== ##
#
# Add in any constrains on metadata here.
#
# Download an <-                            Effect  ==  "downstream_gene_variant"  ~ "Synonymous",
Effect  ==  "stop_retained_variant"  ~ "Synonymous",
Effect  ==  "stop_lost" ~ "Nonsense",
Effect  ==  "start_lost"  ~ "Nonsense",
                                Effect  ==  "stop_gained" ~ "Nonsense")) %>% 
      # add sample ID
      mutate(Accession = sample.name) %>%
      mutate(Caller = sample.caller) %>% 
      mutate(Aligner = sample.aligner)
    
    sample.df$AF = as.numeric(sample.df$AF) / 100
  }
  
  ## ==== Join the variant data with the metadata by 'Run' id ==== ##
  
  # Join with metadata
  sample.df = left_join(sample.df, metadata.df, by = "Accession")
  
  # Write out to a file
  write.csv(sample.df, file = snakemake@output[[1]], row.names=FALSE, sep="\t")

}else{
  # Write an empty output if file is empty
  cat(NULL,file=snakemake@output[[1]])
}


## ==== END ==== ##













