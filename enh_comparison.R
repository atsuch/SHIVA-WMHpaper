library(tidyverse)
library(ggpubr)
library(glue)

# Set wd to where 'data' with necessary data are located
setwd("/PATH-to-working-dir/")

##############################################################################
#  This script can be used to recreate Supp Figure 1, which plots performance 
# improvements from training data enhancements with outputs of a model trained  
# only with 40 MRiShare training data (i.e. WMH segmentations based on the model
# specialised for MRiShare data are used as surrogate reference segmentations.) 
##############################################################################

# Folder to store outputs of this script
dir.create('enh_comparison', showWarnings = FALSE)

##### 1) Load and prepare data
enh_models <- c('ishare_only', 'base', 'enh90', 'enh270', 'enh300', 'enh400')
enh_labels <- c('MRi-Share-specific', 'base', 'enh90', 'enh270', 'enh300', 'enh400')

dat <- read_csv('data/enh_comparison_data.csv') %>%
  mutate(model_name = factor(model_name, levels = enh_models, labels = enh_labels)) %>%
  pivot_wider(id_cols = c(SUB_ID, cohort, model_name), names_from = output_type, values_from = F1) %>%
  mutate(cv_mean_F1 = (Voxel + Cluster)/2)

# data without MRi-Share-specific model, which is not comparable to other models
comp_dat <- dat %>%
  filter(model_name != 'MRi-Share-specific')

##### 2) Plot and save
pal = c(
  'MRi-Share-specific' = 'gold', 
  'base' = 'green4',
  'enh90' = 'green3',
  'enh270' = 'green2',
  'enh300' = 'chartreuse1',
  'enh400' = 'greenyellow'
)

cohort.labs <- c("MRi-Share (n = 40)", "MWC (n = 50)")
names(cohort.labs) <- c("MRi-Share", "MWC")
p <- ggplot( data = dat, aes(x = model_name, y = cv_mean_F1, fill = model_name)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual("model name", values = pal) +
  geom_jitter(color='lightsteelblue3', size=0.5, alpha=0.9) +
  scale_x_discrete(guide = guide_axis(n.dodge=2)) +
  stat_summary(data = comp_dat, fun = median, geom = 'line', color = 'black', aes(group = 1)) +
  xlab("") +
  ylab("Average of VL- and CL-Dice scores") +
  ylim(0, 1.0) +
  facet_wrap(
    ~cohort,
    labeller = labeller(cohort = cohort.labs)
    ) +
  theme_bw() +
  theme(
    strip.text.x = element_text(size = 14),
    strip.background = element_rect(color="white", fill="white")
  )

p

ggsave("enh_comparison/SuppFig1.png",
       p,
       device="png",
       width=10,
       height=5.5,
       units="in",
       dpi=300)
