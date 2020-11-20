## ---------------------------
##
## Script name: `rule plot_rbd_escape_mutations`
##
## Purpose of script: Plot the escape mutations in RBD for timepoints 143, 146, and 152
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
pysam.data = snakemake@input[[1]] # Make the snakemake input file

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
  # Must be in the timeframe of interest (Day 143, 146, and 152)
  filter(DAY >= 143) %>% 
  # At least one mutation over 10%
  filter(AF >= 0.1) %>% 
  # Get a list of the relevant protein positions
  pull(PROT_POS) %>% 
  # only the unqiue sites are necessary
  unique()

## ==== Further filtering to remove sites that don't change over the interval ==== ##
figure.df = pysam.df[which(pysam.df$PROT_POS %in% sites.to.plot),] %>% 
  # Must be in the timeframe of interest (Day 143, 146, and 152)
  filter(DAY >= 143) %>% 
  # Get only the missense mutations
  filter(MISSENSE == "TRUE") %>% 
  # Drop the extra columns 
  select(!c("ACCESSION", "RBD", "MISSENSE", "DP", "CONS", "SNP")) %>% 
  # Make a new column for pivot 
  mutate(DAY = paste0("DAY_", DAY)) %>% 
  # Pivot wider to have information for every single time point
  pivot_wider(names_from = DAY, values_from = AF, values_fill = 0) %>%  
  # Filter out the mutations that change minimally
  mutate(MEAN_CHANGE = rowMeans(select(.,starts_with("DAY_")))) %>% 
  # Filter if the change is minimal i.e., between around 98% frequency and greater than 2%
  filter(MEAN_CHANGE <= 0.98 & MEAN_CHANGE >= 0.02) %>% 
  # Pivot longer to tidy version. 
  pivot_longer(starts_with("DAY_"), names_to = "DAY", values_to = "AF") %>% 
  # Split the day back into a number 
  mutate(DAY = as.integer(str_replace_all(DAY, "DAY_", ""))) %>% 
  # Remove the `MEAN_CHANGE` columns 
  select(!MEAN_CHANGE) 

## ==== Get the matching WT positions and their frequencies ==== ##
wt.df = pysam.df[which(pysam.df$POS %in% figure.df$POS),] %>% 
  # Must be in the timeframe of interest (Day 143, 146, and 152)
  filter(DAY >= 143) %>% 
  # Get only the missense mutations
  filter(MISSENSE == "FALSE") %>% 
  # Drop the extra columns 
  select(!c("ACCESSION", "RBD", "MISSENSE", "DP", "CONS", "SNP")) %>% 
  # Make a new column for pivot 
  mutate(DAY = paste0("DAY_", DAY)) %>% 
  # Pivot wider to have information for every single time point
  pivot_wider(names_from = DAY, values_from = AF, values_fill = 0) %>%  
  # Pivot longer to tidy version. 
  pivot_longer(starts_with("DAY_"), names_to = "DAY", values_to = "AF") %>% 
  # Split the day back into a number 
  mutate(DAY = as.integer(str_replace_all(DAY, "DAY_", ""))) 

## ==== Combine the WT and MUT data and plot figure ==== ##
rbind(figure.df, wt.df) %>% 
  # Is the row a WT or MUT amino acid
  mutate(IS_WT = if_else(WT_AA == MUT_AA, 'WT', 'MUT')) %>%
  # If there are more than one mutation at a site, manually alter the coloring
  mutate(IS_WT = case_when(AA_CHANGE == "F486L" ~ "MUT-2",
                           TRUE ~ IS_WT)) %>% 
  # Relevel the factors for plotting
  mutate(IS_WT = fct_relevel(IS_WT, "WT", "MUT", "MUT-2")) %>% 
  
  ## == Plot the labeled stacked area plot == ##
  
  ggplot(aes(x = DAY, y = AF, fill = IS_WT)) +
  geom_area() +
  facet_wrap(~PROT_POS, nrow = 1) +
  ylab("mutation frequency") + 
  xlab("days after diagnosis") + 
  scale_x_continuous(breaks = c(143, 146, 152), limits = c(143, 152), expand = c(0,0)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), expand = c(0,0)) +
  scale_fill_manual(values = tab_colors) +
  ### ==== Segments/Point Annotations X-Axis ==== ###
  annotate(geom = "segment", x = 145, y = 0 ,  xend = 145, yend = 1, size = 1.5, col = "black", linetype = 4) +
  ### ==== Themes for plot aestetics ==== ###
  theme_classic() +
  theme(
    # Set the appropriate font size
    text=element_text(size=19,  family="Helvetica"),
    # Change the spacing of the facets
    panel.spacing = unit(2, "lines"),
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
ggsave(snakemake@output[[1]], width = 18, height = 4, dpi = 300)

