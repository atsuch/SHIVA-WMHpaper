library(tidyverse)
library(rstatix)
library(ggpubr)
library(glue)

# Set wd to where 'data' with necessary data are located
setwd("/PATH-to-working-dir/")

##############################################################################
#  This script can be used to recreate components of Figure 3, which plots
# comparisons for DWI-based WM properties in WMH vs NAWM in 50 MRi-Share 
# subjects with manual WMH segmentation.
##############################################################################

# Folder to store outputs of this script
dir.create('manuWMH_diff_profiles', showWarnings = FALSE)

##### 1) Load data
dat <- read_csv('data/MRiSHARE_manuWMH_diff_profile_data.csv') %>%
  convert_as_factor(SUB_ID)

##### 2) Comparison of WMH vs NAWM

# Function to perform paired t-test and plot
plot_wmh_profile <- function(metric = 'NDI') {
  if (metric == 'MD') {
    metric_col = glue('scaled{metric}')
    ylabel = expression("mean MD"~ ( x ~ 10^{-4} ~ mm^{2}/sec))
  } else{
    metric_col = glue('mean{metric}')
    ylabel = glue('mean {metric}')
  }
  metric_dat <- dat %>%
    select('SUB_ID', 'mask_type', metric_col)
  
  f = reformulate(termlabels = 'mask_type', response = metric_col)
  
  # stat test with mean and conf int info
  res <- t.test(formula = f, data = metric_dat, paired = TRUE, conf.level=0.95)
  
  stat.test <- metric_dat  %>%
    t_test(formula = f, paired = TRUE) %>%
    mutate(estimate = res$estimate,
           conf_l = res$conf.int[[1]],
           conf_h = res$conf.int[[2]]) %>%
    add_significance("p") %>%
    add_xy_position(x = "mask_type")
  
  scaleFUN <- function(x) sprintf("%.1f", x)
  
  p <- ggpaired(
    dat, x = 'mask_type', y = metric_col, id = 'SUB_ID',
    color = 'slategray', line.color = "gray", line.size = 0.2,
    point.size = 2, fill = 'mask_type', palette = c("turquoise3", "violetred3"),
    ylab = ylabel, xlab = '') +
    stat_pvalue_manual(
      stat.test, hide.ns = TRUE, step.increase = 0.1
    ) +
    labs(
      subtitle = get_test_label(stat.test, detailed = TRUE)
    ) +
    scale_y_continuous(labels = scaleFUN) +
    font("subtitle", size = 18) +
    font("ylab", size = 24) +
    font("xy.text", size = 16) +
    theme(legend.position = "none")

  ggsave(glue("manuWMH_diff_profiles/{metric}_WMH_vs_NAWM.png"),
         p,
         device="png",
         width=6,
         height=5,
         units="in",
         dpi=300)
  
  return(stat.test)
  
}

# Now plot and show comparisons (Figure 3) 
ndi_comparison <- plot_wmh_profile(metric = 'NDI')
fa_comparison <- plot_wmh_profile(metric = 'FA')
md_comparison <- plot_wmh_profile(metric = 'MD')
