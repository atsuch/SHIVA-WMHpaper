library(tidyverse)
library(ggpubr)
library(glue)

# Set wd to where 'data' with necessary data are located
setwd("/PATH-to-working-dir/")

##############################################################################
#  This script is used to create the plot comparing the predicted and ground
# truth (i.e. manually-traced WMH) volumes, and can be used to generate Figure 4
# of the SHIVA-WMH paper main text.
##############################################################################

# Folder to store outputs of this script
dir.create('volume_comparison', showWarnings = FALSE)

##### 1) Load and prepare data
main_methods <- c('LST-LPA', 'WMH_pgs', 'HPM', 'enh270')
main_labels <- c('LST-LPA', 'PGS', 'HPM', 'SHIVA')

dat <- read_csv('data/method_comparison_data.csv') %>%
  filter(methods %in% main_methods, output_type == 'Voxel') %>%
  mutate(
    log_pred_vol = log10(total_pred_vol),
    log_truth_vol = log10(total_raw_truth_vol)
    ) %>%
  mutate(methods = factor(methods, levels = main_methods, labels = main_labels))

##### 2) Plot and save
p <- dat %>%
  ggscatter(x = 'log_truth_vol', y = 'log_pred_vol', 
            color = 'methods', palette = c('magenta2', 'cyan3', 'darkorange1', 'green'),
            add = 'reg.line', corr.coef = TRUE, conf.int = TRUE, 
            title = 'Estimated against true WMH volume for each method',
            xlab = expression(paste("log10(lesion volume in mm"^3, ")")),
            ylab = expression(paste("log10(predicted lesion volume in mm"^3, ")")),
            xlim = c(0.9, 5), ylim=c(0.9, 5))  +
  geom_abline(slope = 1, intercept = 0) +
  stat_cor(aes(color = methods), label.x = 1, size = 5, show.legend = FALSE) +
  font("title", size = 14) +
  font("xy", size = 12) +
  font("xy.text", size = 10)

p

ggsave(glue("volume_comparison/Figure4.png"),
       p,
       device="png",
       width=5.5,
       height=6,
       units="in",
       dpi=300)