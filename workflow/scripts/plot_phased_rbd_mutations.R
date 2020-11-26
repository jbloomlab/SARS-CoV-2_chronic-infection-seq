## ---------------------------
##
## Script name: `rule plot_phased_rbd_mutations`
##
## Purpose of script: Plot the haplotypes of the mutations in RBD for timepoints 143, 146, and 152
##
## Author: Will Hannon
##
## Date Created: 2020-11-26
##
## Copyright (c) Will Hannon, 2020
## Email: wwh22@uw.edu
##
## ---------------------------

require(tidyverse)

## ---------------------------

## ==== Path to haplotype data from pysam script ==== ##
haplotype.data = snakemake@input[[1]] # Make the snakemake input file

## ==== Read in the haplotype data ==== ##
haplotype.df = read_csv(haplotype.data) %>% 
  filter(DAY >= 143) %>% 
  # Remove the N440D column becuase the haplotype is know (100% at day 146 0% all other days.)
  select(!N440D)


## ==== Save the annotations ==== ##
# N440D, F486I, F486L, E484A, Y489H, Q493K
#   1      1      1      1      1      1 

HAP = c("000000",
        "000010",
        "001000",
        "010000",
        "000001",
        "000100", 
        "001100", 
        "010100", 
        "000110", 
        "100001", 
        "000101",
        "001110", 
        "010101", 
        "010110"
)

SNPs = c("Wildtype",
         "Y489H",
         "F486L",
         "F486I",
         "Q493K",
         "E484A", 
         "F486L-E484A", 
         "F486I-E484A", 
         "E484A-Y489H", 
         "N440D-Q493K", 
         "E484A-Q493K",
         "F486L-E484A-Y489H", 
         "F486I-E484A-Q493K", 
         "F486I-E484A-Y489H"
)

annotations = data.frame(HAP , SNPs)

## === Haplotypes at day 143  === ##
day143 = haplotype.df %>% 
  filter(DAY == 143) %>% 
  select(!Q493K) %>% 
  mutate(F486L = case_when(F486I == 1 ~ 0,
                           TRUE ~ F486L)) %>% 
  mutate(F486I = case_when(F486L == 1 ~ 0,
                           TRUE ~ F486I)) %>% 
  drop_na() %>% 
  mutate(HAP = paste0(0, F486I, F486L, E484A, Y489H, 0)) %>% 
  group_by(HAP) %>% 
  summarise(Population=n()) %>% 
  mutate(Frequency = Population/sum(Population)) %>% 
  arrange(-Frequency) %>% 
  mutate(DAY = 143)

## === Haplotypes at day 146  === ##
day146 = haplotype.df %>% 
  filter(DAY == 146) %>% 
  mutate(F486L = case_when(F486I == 1 ~ 0,
                           TRUE ~ F486L)) %>% 
  mutate(F486I = case_when(F486L == 1 ~ 0,
                           TRUE ~ F486I)) %>% 
  drop_na() %>% 
  mutate(HAP = paste0(1, F486I, F486L, E484A, Y489H, Q493K)) %>% 
  group_by(HAP) %>% 
  summarise(Population=n()) %>% 
  mutate(Frequency = Population/sum(Population)) %>% 
  arrange(-Frequency) %>% 
  mutate(DAY = 146)

## === Haplotypes at day 143  === ##
day152 = haplotype.df %>% 
  filter(DAY == 152) %>% 
  #select(!Y489H) %>% 
  mutate(F486L = case_when(F486I == 1 ~ 0,
                           TRUE ~ F486L)) %>% 
  mutate(F486I = case_when(F486L == 1 ~ 0,
                           TRUE ~ F486I)) %>% 
  drop_na() %>% 
  mutate(HAP = paste0(0, F486I, F486L, E484A, Y489H, Q493K)) %>% 
  group_by(HAP) %>% 
  summarise(Population=n()) %>% 
  mutate(Frequency = Population/sum(Population)) %>% 
  arrange(-Frequency) %>% 
  mutate(DAY = 152)



## === Main Plot === ##

y.pos = 0
point.size = 1
adj = 0.05
tick.color = "#bcbd22"

left_join(rbind(rbind(day143, day146), day152) , annotations) %>% 
  filter(Frequency > 0.01) %>% 
  select(!Population) %>% 
  # Pivot wider to have information for every single time point
  pivot_wider(names_from = DAY, values_from = Frequency, values_fill = 0) %>%  
  # Pivot longer to tidy version. 
  pivot_longer(!c(HAP, SNPs), names_to = "DAY", values_to = "Frequency") %>% 
  ggplot(aes(x = as.integer(DAY), y = Frequency, fill = SNPs)) + 
  geom_area() + 
  ylab("mutation frequency") + 
  xlab("days after diagnosis") + 
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1) , expand = c(0,0)) +
  scale_x_continuous(limits = c(143,152) , breaks = c(143, 146, 152), expand = c(0,0)) +
  ### ==== Segments/Point Annotations X-Axis ==== ###
  annotate(geom = "segment", x = 145, y = 0 ,  xend = 145, yend = 1, size = 1, col = "black", linetype = 4) +
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
    plot.margin = margin(50, 50, 50, 50)
  ) + 
  theme(plot.title = element_text(hjust = 0.5))

# Save the file to figures directory - will be the Snakemake output path
ggsave(snakemake@output[[1]], width = 10, height = 8, dpi = 300)






