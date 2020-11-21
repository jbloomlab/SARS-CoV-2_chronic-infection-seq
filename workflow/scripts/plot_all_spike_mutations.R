## ---------------------------
##
## Script name: `rule plot_rbd_all_mutations`
##
## Purpose of script: Plot the all of the mutations in Spike for every timepoint day 18 to day 152.
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

## ==== Path to pileup data from pysam script ==== ##
pysam.data = snakemake@input[[1]]  # Make the snakemake input file

## ==== Add matplotlib Tableau color scheme ==== ##
tab_colors = c("#1f77b4","#ff7f0e")

## ==== Read in the pysam/pileup data ==== ##
pysam.df = read_csv(pysam.data)

## ==== Get a list of the amino acid positions to plot ==== ##
sites.to.plot = pysam.df %>% 
  # Must have at least one missense mutation
  filter(MISSENSE == "TRUE") %>% 
  # At least one mutation over 10%
  filter(AF >= 0.1) %>% 
  # Must have at least 100X coverage (not an issue for RBD)
  filter(DP >= 100) %>% 
  # Get a list of the relevant protein positions
  pull(PROT_POS) %>% 
  # only the unqiue sites are necessary
  unique()

## ==== Further filtering to remove sites that don't change over the interval ==== ##
figure.df = pysam.df[which(pysam.df$PROT_POS %in% sites.to.plot),] %>% 
  # Get only the missense mutations
  filter(MISSENSE == "TRUE") %>% 
  # Drop the extra columns 
  select(!c("ACCESSION", "MISSENSE", "DP", "CONS", "SNP")) %>% 
  # Make a new column for pivot 
  mutate(DAY = paste0("DAY_", DAY)) %>% 
  # Pivot wider to have information for every single time point
  pivot_wider(names_from = DAY, values_from = AF, values_fill = 0) %>%  
  # Pivot longer to tidy version. 
  pivot_longer(starts_with("DAY_"), names_to = "DAY", values_to = "AF") %>% 
  # Split the day back into a number 
  mutate(DAY = as.integer(str_replace_all(DAY, "DAY_", ""))) %>% 
  # Add the color information for plotting
  mutate(COLOR_SCHEME = case_when(AA_CHANGE == "N440D" ~ "ESCAPE",
                                  AA_CHANGE == "F486I" ~ "ESCAPE",
                                  AA_CHANGE == "Q493K" ~ "ESCAPE",
                                  AA_CHANGE == "E484A" ~ "ESCAPE",
                                  AA_CHANGE == "Y489H" ~ "ESCAPE",
                                  TRUE ~ "NON-ESCAPE"))


## ==== Plot Figure ==== ##
y.pos = 0
point.size = 1
adj = 0.05
tick.color = "#bcbd22"

figure.df %>% 
  # Remove D614G which is present in all timepoints.
  filter(AA_CHANGE != "D614G") %>%
  # Add "residue" to the name for plotting
  mutate(RESIDUE = fct_reorder(paste("Residue", PROT_POS, sep = " "), PROT_POS, .desc = F)) %>% 
  ## ========= Plot the figure ========= ##
  
  ggplot(aes(x = DAY, y = AF, fill = COLOR_SCHEME, group = MUT_AA)) +
    # Annotation to represent the WT allele frequency
    annotate(geom = "rect", xmin = 18, xmax = 152, ymin=0, ymax=1, fill = "#7f7f7f") +
    geom_area() +
    #geom_point(col ="black") +
    facet_wrap(~RESIDUE, ncol = 5) +
    ylab("mutation frequency") + 
    xlab("days after diagnosis") + 
    scale_fill_manual(values = tab_colors) +
    scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(-.05,1) , expand = c(0,0)) +
    scale_x_continuous(limits = c(18,152) , expand = c(0,0)) +
    ### ==== Segments/Point Annotations X-Axis ==== ###
    annotate(geom = "segment", x = 145, y = 0 ,  xend = 145, yend = 1, size = 1, col = "black", linetype = 4) +
    annotate(geom = "segment", x = 18, y = 0 ,  xend = 152, yend = 0, size = 1, col = "black", linetype = 1) +
  
    annotate(geom = "segment", x = 18, y = y.pos ,  xend = 18, yend = y.pos+adj, size = 1, col = tick.color) +
    annotate(geom = "segment", x = 25, y = y.pos ,  xend = 25, yend = y.pos+adj, size = 1, col = tick.color) +
    annotate(geom = "segment", x = 75, y = y.pos ,  xend = 75, yend = y.pos+adj, size = 1, col = tick.color) +
    annotate(geom = "segment", x = 81, y = y.pos ,  xend = 81, yend = y.pos+adj, size = 1, col = tick.color) +
    annotate(geom = "segment", x = 128, y = y.pos ,  xend = 128, yend = y.pos+adj, size = 1, col = tick.color) +
    annotate(geom = "segment", x = 130, y = y.pos ,  xend = 130, yend = y.pos+adj, size = 1, col = tick.color) +
    annotate(geom = "segment", x = 143, y = y.pos ,  xend = 143, yend = y.pos+adj, size = 1, col = tick.color) +
    annotate(geom = "segment", x = 146, y = y.pos ,  xend = 146, yend = y.pos+adj, size = 1, col = tick.color) +
    annotate(geom = "segment", x = 152, y = y.pos ,  xend = 152, yend = y.pos+adj, size = 1, col = tick.color) +
    ### ==== Segments/Point Annotations X-Axis ==== ###
    
    ### ==== Themes for plot aestetics ==== ###
  
  
    theme_classic() +
    theme(
      # Set the appropriate font size
      text=element_text(size=20,  family="Helvetica"),
      # Change x axis title position
      axis.title.x = element_text(vjust= -0.5),
      # Change y axis title position
      axis.title.y = element_text(vjust= 1),
      # Change the distace of lables from x axis
      axis.text.x = element_text(vjust= -1),
      # Change the distace of lables from y axis
      axis.text.y = element_text(hjust= 2),
      # Remove the legend
      legend.position = "none",
      # Change the margin to make sure x axis labels aren't cut-off
      plot.margin = margin(10, 30, 10, 10)
    )
  
# Save the file to figures directory - will be the Snakemake output path
ggsave(snakemake@output[[1]], width = 20, height = 7, dpi = 300)

# NOTES
# I need to change the plotting so that extra mutations are shown even though they aren't escape mutations.
