library(tidyverse)
library(rstatix)
library(ggpubr)
library(gt)
library(gtsummary)
library(glue)

# Set wd to where 'data' with necessary data are located
setwd("/PATH-to-working-dir/")

##############################################################################
#  This script performs statistical comparisons and visualization of similarity
# metrics for comparing following methods for WMH detection:
#
#  - SHIVA-WMH
#  - LST-LPA
#  - PGS
#  - HyperMapper (HPM)
#
# It can be used to generate:
# 
#  - Table 2&3
#  - Figure 5
#  - Supp Table 1
#  - Supp Figure 2
# 
# of SHIVA-WMH paper.
#
##############################################################################

# Folder to store outputs of this script
dir.create('metric_comparison', showWarnings = FALSE)

##### 1) Load and prepare data
#####
# Extract metrics used for comparison and put in long form: Cluster-level HD95 
# is discarded since there is only single set of HD95 scores and 'Voxel' vs 
# 'Cluster' level distinction does not apply for this metric

dat <- read_csv('data/method_comparison_data.csv') %>%
  select(SUB_ID, cohort, methods, output_type, precision, sensitivity, F1, HD95) %>%
  pivot_longer(cols = precision:HD95, names_to = 'metric', values_to = 'score') %>%
  mutate(keep = if_else((output_type == "Cluster" & metric == 'HD95'), 0, 1)) %>%
  filter(keep == 1) %>%
  select(-keep)

# Add combined data across cohort in the nested data so that separate
# comparisons can be performed collapsed across the cohort.

all_dat <- dat %>%
  mutate(cohort = 'All')

dat <- dat %>%
  bind_rows(all_dat) 

# Summary of each performance metric
summary <- dat %>%
  group_by(cohort, methods, output_type, metric) %>%
  summarise(n = n(), mean = mean(score, na.rm = TRUE), sd = sd(score, na.rm = TRUE))

#####
  
##### 2) Pairwise comparison of each performance metric
#####
# Create nested data to perform pairwise comparisons for each cohort, output_type
# and metric. 
nested_dat <- dat %>%
  nest_by(cohort, output_type, metric)

# Function to perform pairwise comparison of reference models against SHIVA 
# (enh270 OR enh270_mono) using paired t-test in each cohort 
get_pwc <- function(df, cohort, ref = 'enh270') {
  if (cohort == 'MWC') {
    comp_models = c('LST-LPA', 'HPM')
  }
  else {
    comp_models = c('LST-LPA', 'WMH_pgs', 'HPM')
  }
  pwc = NULL
  
  formula = reformulate(termlabels = 'methods', response = 'score')
  
  for (comp_model in comp_models) {
    comparison = c(ref, comp_model)
    
    if (cohort == 'All' & comp_model == 'WMH_pgs') {
      pwc_df <- df %>%
        filter(!grepl('MWC', SUB_ID))
        
    } else {
        pwc_df <- df
    }
    
    pwc_df <- pwc_df %>%
      filter(methods == ref | methods == comp_model) %>%
      convert_as_factor(SUB_ID) %>%        
      mutate(methods = factor(methods, levels = comparison))
        
    pwc_to_bind <- pwc_df %>%
        t_test(
          formula, paired = TRUE, p.adjust.method = "bonf"
        ) %>%
        mutate(test_type = 'paired_t')

    pwc <- pwc %>%
      bind_rows(pwc_to_bind)
    }

  pwc <- pwc %>%
    mutate(p.adj = ifelse(p < 0.05, p*length(comp_models), p)) %>%
    add_significance("p.adj") %>%
    add_xy_position(x = "methods") %>%
    mutate(xmax = seq(from = 2, by = 1, length.out = length(comp_models)))
    
  return(pwc)
}

# Summary of pairwise comparisons for SHIVA model 
main_pwc = NULL
for (i in 1:nrow(nested_dat)) {
  pwc_out <- get_pwc(nested_dat$data[[i]], cohort = nested_dat$cohort[[i]]) %>%
    mutate(
      cohort = nested_dat$cohort[[i]],
      output_type = nested_dat$output_type[[i]],
      metric = nested_dat$metric[[i]]
      ) %>%
    select(cohort, output_type, test_type, metric, everything())
  
  main_pwc <- main_pwc %>%
    bind_rows(pwc_out)
}

# Summary of pairwise comparisons for SHIVA (FLAIR-only) model
supp_pwc = NULL
for (i in 1:nrow(nested_dat)) {
  pwc_out <- get_pwc(nested_dat$data[[i]], cohort = nested_dat$cohort[[i]], ref = 'enh270_mono') %>%
    mutate(
      cohort = nested_dat$cohort[[i]],
      output_type = nested_dat$output_type[[i]],
      metric = nested_dat$metric[[i]]
    ) %>%
    select(cohort, output_type, test_type, metric, everything())
  
  supp_pwc <- supp_pwc %>%
    bind_rows(pwc_out)
}

#####

##### 3) Summary Table containing mean (SD) for each metric/method
#####
# Make a summary table of Mean (SD) of each metric for each method, with asterisk
# indicating any significant difference between SHIVA and other methods

# Function to create summary table for a given cohort
make_summary_table = function(cohort ='All', ref='enh270') {
  if (ref == 'enh270') {
    fname_prefix = 'SHIVA'
    method_level = c('enh270', 'lstlpa', 'WMH_pgs', 'HPM')
    method_labels = c('SHIVA', 'LST-LPA', 'PGS', 'HPM')
    summary_dat <- summary %>%
      filter(methods != 'enh270_mono')
    pwc_dat <- main_pwc
  } else if (ref == 'enh270_mono') {
    fname_prefix = 'SHIVA-FLonly'
    method_level = c('enh270_mono', 'lstlpa', 'WMH_pgs', 'HPM')
    method_labels = c('SHIVA (FLAIR only)', 'LST-LPA', 'PGS', 'HPM')
    summary_dat <- summary %>%
      filter(methods != 'enh270')
    pwc_dat <- supp_pwc
  }
  
  pwc_pvals <- pwc_dat %>%
    select(cohort, output_type, metric, group2, p.adj) %>%
    rename(methods = group2)
  
  # Prepare data for table
  table_dat <- summary_dat %>%
    filter(cohort == {{cohort}}) %>%
    left_join(pwc_pvals, by = c('cohort', 'output_type', 'metric', 'methods')) %>%
    mutate(
      metric_label = case_when(
        (output_type == 'Voxel' & metric == 'sensitivity') ~ 'VLTPR',
        (output_type == 'Voxel' & metric == 'precision') ~ 'VLPPV',
        (output_type == 'Voxel' & metric == 'F1') ~ 'VLDice',
        (output_type == 'Cluster' & metric == 'sensitivity') ~ 'CLTPR',
        (output_type == 'Cluster' & metric == 'precision') ~ 'CLPPV',
        (output_type == 'Cluster' & metric == 'F1') ~ 'CLDice',
        TRUE ~ metric
      ),
      methods = ifelse(methods == 'LST-LPA', 'lstlpa', methods)
    ) %>%
    mutate(
      metric_label = factor(
        metric_label,
        levels = c('VLTPR', 'VLPPV', 'VLDice', 'CLTPR', 'CLPPV', 'CLDice', 'HD95')
        ),
      methods = factor(
        methods, 
        levels = method_level,
        labels = method_labels
        ),
      fmt_mean =
        format(mean, digit = 1, nsmall = 2, scientific = FALSE),
      fmt_sd =
        format(sd, digit = 1, nsmall = 2, scientific = FALSE)
    ) %>%
    mutate(fmt_mean = case_when(
      p.adj < 0.0001 ~ glue("{fmt_mean}****"),
      p.adj < 0.001 & p.adj >= 0.0001 ~ glue("{fmt_mean}***"),
      p.adj < 0.01 & p.adj >= 0.001 ~ glue("{fmt_mean}**"),
      p.adj < 0.05 & p.adj >= 0.01 ~ glue("{fmt_mean}*"),
      TRUE ~ glue("{fmt_mean}"))
    ) %>%
    arrange(metric_label, methods) %>%
    pivot_wider(id_cols = c(cohort, methods),
                names_from = metric_label,
                values_from = starts_with("fmt_"),
                names_sep = "_") %>%
    group_by(cohort) 
  
  # Create a gt!
  gt_table <- table_dat %>%
    gt(rowname_col = "methods") %>%
    tab_spanner(
      label = "Mean (SD)",
      columns = starts_with("fmt_")
    ) %>%
    cols_merge(
      columns = c(fmt_mean_VLTPR, fmt_sd_VLTPR),
      hide_columns = c(fmt_sd_VLTPR),
      pattern = "{1} ({2})"
    ) %>%
    cols_merge(
      columns = c(fmt_mean_VLPPV, fmt_sd_VLPPV),
      hide_columns = c(fmt_sd_VLPPV),
      pattern = "{1} ({2})"
    ) %>%
    cols_merge(
      columns = c(fmt_mean_VLDice, fmt_sd_VLDice),
      hide_columns = c(fmt_sd_VLDice),
      pattern = "{1} ({2})"
    ) %>%
    cols_merge(
      columns = c(fmt_mean_CLTPR, fmt_sd_CLTPR),
      hide_columns = c(fmt_sd_CLTPR),
      pattern = "{1} ({2})"
    ) %>%
    cols_merge(
      columns = c(fmt_mean_CLPPV, fmt_sd_CLPPV),
      hide_columns = c(fmt_sd_CLPPV),
      pattern = "{1} ({2})"
    ) %>%
    cols_merge(
      columns = c(fmt_mean_CLDice, fmt_sd_CLDice),
      hide_columns = c(fmt_sd_CLDice),
      pattern = "{1} ({2})"
    ) %>%
    cols_merge(
      columns = c(fmt_mean_HD95, fmt_sd_HD95),
      hide_columns = c(fmt_sd_HD95),
      pattern = "{1} ({2})"
    ) %>%
    cols_label(
      fmt_mean_VLTPR = "VL-TPR",
      fmt_mean_VLPPV = "VL-PPV",
      fmt_mean_VLDice = "VL-Dice",
      fmt_mean_CLTPR = "CL-TPR",
      fmt_mean_CLPPV = "CL-PPV",
      fmt_mean_CLDice = "CL-Dice",
      fmt_mean_HD95 = "HD95"
    ) %>%
    gtsave(filename = glue("{fname_prefix}_comparison_metric_summary_{cohort}.html"),
           path = glue("{getwd()}/metric_comparison/"))

  return(table_dat)
}

# Now create summary tables
# 1) SHIVA vs others across All cohorts (Table 2)
main_all_table <- make_summary_table()

# 2) SHIVA vs others in each cohort separately (Table 3)
main_mrishare_table <- make_summary_table(cohort = 'MRi-Share')
main_ukb_table <- make_summary_table(cohort = 'UKB')
main_mwc_table <- make_summary_table(cohort = 'MWC')

# 3) SHIVA (FLAIR only) vs others across All cohorts (Supp Table 1)
supp_all_table <- make_summary_table(ref = 'enh270_mono')

#####

##### 4) Plot comparing CL- and VL-Dice across methods
#####
# Plot of CL- or VL-Dice scores with asterisk indicating any significant 
# difference between SHIVA and other methods

# Function to plot paired t-test results
plot_dice_comparison <- function(output_type = 'Voxel', cohort = 'MRi-Share', ref = 'enh270', pwc_step_size = 0, pwc_ypos = NULL) {
  if (ref == 'enh270') {
    fname_prefix = 'SHIVA'
    ref_color <- 'green'
    pwc_dat <- main_pwc
    
    if (cohort == 'MWC') {
      method_level = c('enh270', 'LST-LPA', 'HPM')
      method_labels = c('SHIVA', 'LST-LPA', 'HPM')
    } else {
      method_level = c('enh270', 'LST-LPA', 'WMH_pgs', 'HPM')
      method_labels = c('SHIVA', 'LST-LPA', 'PGS', 'HPM')
    }
    
  } else if (ref == 'enh270_mono') {
    fname_prefix = 'SHIVA-FLonly'
    ref_color <- 'greenyellow'
    pwc_dat <- supp_pwc
    if (cohort == 'MWC') {
      method_level = c('enh270_mono', 'LST-LPA', 'HPM')
      method_labels = c('SHIVA (FLAIR only)', 'LST-LPA', 'HPM')
    } else {
      method_level = c('enh270_mono', 'LST-LPA', 'WMH_pgs', 'HPM')
      method_labels = c('SHIVA (FLAIR only)', 'LST-LPA', 'PGS', 'HPM')
    }
  }

  # plot data selection
  plot_dat <- dat %>%
    filter(
      output_type == {{output_type}},
      cohort == {{cohort}},
      metric == 'F1',
      methods %in% method_level
    ) %>%
    convert_as_factor(SUB_ID) %>%
    mutate(methods = factor(methods, levels = method_level, labels = method_labels))
  
  # Set y position for significance indications
  if (!is.null(pwc_ypos)) {
    ypos = pwc_ypos
  } else {
    ypos <- seq(
      max(plot_dat$score, na.rm = T),
      by=0.022,
      length.out = length(method_level) - 1
      )
  }
  
  # Pairwise p values to use for significance indications
  pwc_pvals <- pwc_dat %>%
    filter(
      output_type == {{output_type}},
      cohort == {{cohort}},
      metric == 'F1'
    ) %>%
    mutate(y.position = ypos)
    
  # Palette and y labels of the graph: for MWC, there will be no PGS to compare 
  if (cohort == 'MWC') {
    pal = c(ref_color, 'magenta2', 'darkorange1')
  } else {
    pal = c(ref_color, 'magenta2', 'cyan3', 'darkorange1')
  }
  
  if (output_type == 'Voxel') {
    ylab1 = 'VL'
  } else {
    ylab1 = 'CL'
  }
  
  ylabel = glue("{ylab1}-Dice")
  
  # Plot Title
  if (cohort == 'All') {
    title = 'All test subjects (N = 31)'
  } else if (cohort == 'MRi-Share') {
    title = 'MRi-Share (n = 10)'
  } else if (cohort == 'UKB') {
    title = 'UKB (n = 11)'
  } else if (cohort == 'MWC') {
    title = 'MWC (n = 10)'
  }
  
  # Now plot!
  p <- ggpaired(
    plot_dat, x = 'methods', y = 'score', id = 'SUB_ID',
    color = 'lightsteelblue3', line.color = "gray", line.size = 0.2,
    point.size = 2, fill = 'methods',
    palette = pal,
    ylab = ylabel, xlab = '',
    ylim = c(0, 1.0), legend = "right"
  ) +
  stat_pvalue_manual(
    pwc_pvals, hide.ns = FALSE, step.increase = pwc_step_size
  ) +
  scale_x_discrete(guide = guide_axis(n.dodge=2)) +
  ggtitle(title)

  # Save
  save_fname = glue('{fname_prefix}_{ylabel}_{cohort}.png')
  ggsave(glue("metric_comparison/{save_fname}"),
         p,
         device="png",
         width=6,
         height=6,
         units="in",
         dpi=300)

  return(p)

}

# Now plot!
# 1) SHIVA vs others in each cohort separately (Figure 5)
# VL-Dice
main_vlDice_mrishare <- plot_dice_comparison(pwc_step_size = 0.02)
main_vlDice_ukb <- plot_dice_comparison(cohort = 'UKB', pwc_step_size = 0.02)
main_vlDice_mwc <- plot_dice_comparison(cohort = 'MWC', pwc_step_size = 0.02)

# CL-Dice
main_clDice_mrishare <- plot_dice_comparison(
  output_type = 'Cluster',
  pwc_step_size = 0.02,
  pwc_ypos = c(0.91, 0.93, 0.95))
main_clDice_ukb <- plot_dice_comparison(
  output_type = 'Cluster',
  cohort = 'UKB',
  pwc_step_size = 0.02,
  pwc_ypos = c(0.91, 0.93, 0.95))
main_clDice_mwc <- plot_dice_comparison(
  output_type = 'Cluster',
  cohort = 'MWC',
  pwc_step_size = 0.02,
  pwc_ypos = c(0.91, 0.93))

# 2) SHIVA (FLAIR-only) vs others in each cohort separately (Supp Figure 2)
# VL-Dice
supp_vlDice_mrishare <- plot_dice_comparison(ref = 'enh270_mono', pwc_step_size = 0.02)
supp_vlDice_ukb <- plot_dice_comparison(cohort = 'UKB', ref = 'enh270_mono', pwc_step_size = 0.02)
supp_vlDice_mwc <- plot_dice_comparison(cohort = 'MWC', ref = 'enh270_mono', pwc_step_size = 0.02)

# CL-Dice
supp_clDice_mrishare <- plot_dice_comparison(
  output_type = 'Cluster',
  ref = 'enh270_mono',
  pwc_step_size = 0.02)
supp_clDice_ukb <- plot_dice_comparison(
  output_type = 'Cluster',
  cohort = 'UKB',
  ref = 'enh270_mono',
  pwc_step_size = 0.02)
supp_clDice_mwc <- plot_dice_comparison(
  output_type = 'Cluster',
  cohort = 'MWC',
  ref = 'enh270_mono',
  pwc_step_size = 0.02,
  pwc_ypos = c(0.91, 0.93))

#####
