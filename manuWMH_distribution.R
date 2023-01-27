library(tidyverse)
library(ggpubr)
library(glue)

# Set wd to where 'data' with necessary data are located
setwd("/PATH-to-working-dir/")

##############################################################################
#  This script can be used to recreate Figure 2, which plots the distribution of 
# cluster size and total lesion volume in subjects with manual tracing of WMH 
# (train + test of MRi-Share, UKB and MWC).
##############################################################################

# Folder to store outputs of this script
dir.create('manuWMH_distribution', showWarnings = FALSE)

##### 1) Load and prepare data
cohorts <- c('MRi-Share', 'MWC', 'UKB')

dat <- read_csv('data/manuWMH_cluster_distribution_data.csv') %>%
  mutate(
    cohort = factor(cohort, levels = cohorts),
    log_cluster_vol = log10(cluster_volume),
    log_total_vol = log10(total_raw_truth_vol)
  )


##### 2) Plot and save distribution
p <- dat %>%
  ggscatterhist(
    x = "log_total_vol", y = "log_cluster_vol",
    color = "cohort", size = 1.5,
    palette = c(adjustcolor("#00AFBB", 0.3), adjustcolor("#FC4E07", 0.3), adjustcolor("#E7B800", 0.65)),
    margin.params = list(fill = "cohort", color = "black", size = 0.2),
    title = "WMH cluster size and total volume distribution",
    xlab = expression(paste("log10(lesion volume in mm"^3, ")")),
    ylab = expression(paste("log10(lesion cluster size in mm"^3, ")"))
  ) 

p$sp <- p$sp + 
  geom_vline(xintercept = 3, linetype = "dashed", color = "lightslategray") +
  geom_vline(xintercept = log10(5000), linetype = "dashed", color = "lightslategray") +
  geom_vline(xintercept = log10(15000), linetype = "dashed", color = "lightslategray") +
  annotate("text", x = 2.95, y = 4.25, hjust = 1,
           label = "Total WMH volume", color = 'lightslategray') +
  annotate("text", x = c(2.95, log10(5000)-0.05, log10(15000)-0.05), y = 4, hjust = 1,
           label = c("< 1mL", "< 5mL", "<15mL"), color = 'lightslategray')

p
ggsave(
  filename = glue("manuWMH_distribution/Figure2.png"),
  device="png",
  width=6,
  height=5.5,
  units="in",
  dpi=300)
