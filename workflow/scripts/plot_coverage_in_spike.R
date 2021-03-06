## ---------------------------
##
## Script name: `rule plot_coverage_in_spike`
##
## Purpose of script: Plot the coverage over all of Spike.
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

## ---------------------------

## ==== Path to depth data from samtools and sample data that runs analysis ==== ##
depth.data = snakemake@input[[1]]
sample.data = snakemake@params[[1]]

# == Define the regions of interest - Spike, RBD, and RBM == #
spike.coords = c(21563, 25384)
rbd.coords = c(22556, 23158)

# == Import the coverage data and process == #
# Calculate the number of base observavtions in 10bp bins
samtools.depth.bins = left_join(read.table(depth.data, header = T), read_csv(sample.data), by = c("Accession" = "Run"))

# Plot the depth over Spike - no cuttoff
samtools.depth.bins %>% 
  mutate(Depth = ifelse(Depth <= 500, Depth, 500)) %>%
  mutate(Day2 = fct_reorder(paste("Day", Day, sep = " "), Day, .desc = F)) %>% 
  ggplot(aes(x = (Start+Stop)/2, y = Depth)) + 
    geom_line(size = 1.5) +
    facet_wrap(~Day2, ncol = 6) +
    scale_y_continuous(breaks = c(0, 100, 200, 300, 400, 500), labels = c("0", "100", "200", "300", "400", expression(phantom(x) >=500))) +
    scale_x_continuous(breaks = c(22000, 23500, 25000), limits = c(spike.coords[1], spike.coords[2])) + 
    xlab("nucleotide position in the genome") + 
    ylab("depth of coverage") + 
    annotate("rect", xmin = spike.coords[1], xmax = spike.coords[2], ymin = -10, ymax = -50, fill = "darkgoldenrod1", col = "black", alpha = 0.75) + 
    annotate("text", x = (spike.coords[1] + spike.coords[2])/2, y = (-5 + -50)/2, label = "Spike") + 
    annotate("rect", xmin = rbd.coords[1], xmax = rbd.coords[2], ymin = -55, ymax = -100, fill = "coral2", col = "black", alpha = 0.75) + 
    annotate("text", x = (rbd.coords[1] + rbd.coords[2])/2, y = (-55 + -100)/2, label = "RBD") + 
    theme_classic() +
    theme(legend.position="none") +
    theme(legend.box.background = element_rect(colour = "black")) +
    theme(text=element_text(size=22,  family="Helvetica")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(panel.grid.major.y = element_line(colour="grey", linetype="dashed")) 

# Save the file to figures directory
ggsave(snakemake@output[[1]], width = 20, height = 7, dpi = 300)
