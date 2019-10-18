library(ggplot2)
library(tidyverse)
library(cowplot)
library(colorblindr)
theme_set(theme_cowplot())
df = read_csv('../data/estimate_min_num_wells_to_identify_viability_sensitivity_999_simulations.csv')
df$inoculum<-factor(df$inoculum, levels=c(1,1.27, 1.96,2,3,4,5))


tmp_df=data.frame(x1=460, x2=460, y1=2.9, y2=0)
p_i2 <- ggplot(df %>% filter(inoculum==2), aes(x=number_of_inoculated_wells, y=p_i*100)) +
    geom_point() +
    scale_y_log10() +
    scale_x_log10() +
    labs(x='Number of inoculated wells per experiment', y='Minimum relative abundance where no positive wells is significant (%)', title='Inoculum of 2 cells per well') +
    geom_segment(aes(x=x1, xend=x2, y=y1, yend=y2), color='blue', linetype='dashed', data=tmp_df) +
    geom_segment(aes(x=0, xend=x2, y=y1, yend=y1), color='blue', linetype='dashed', data=tmp_df) +
    annotate(geom='text', x=600, y=1.0, label='460 wells', color='blue') +
    annotate(geom='text', x=150, y=3.2, label='2.9%', color='blue')


ggsave('../figures/estimated_sensitivity_for_viability_studies.pdf', device='pdf')


inoculum_cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", '#CC79A7')


p_v123<- ggplot(df %>% filter(inoculum %in% c(1,2,3)), aes(x=number_of_inoculated_wells, y=p_i*100)) +
    geom_point(aes(color=inoculum, fill=inoculum), size=2.5, alpha=0.5, shape=21) +
    scale_y_log10() +
    scale_x_log10() +
    labs(x='Number of inoculated wells per experiment', y='Minimum relative abundance where zero positive wells is significant (%)') +
    geom_vline(xintercept=460, linetype='dashed') +
    scale_color_manual(
        values = darken(inoculum_cols, 0.3)
    ) +
    scale_fill_manual(
        values = inoculum_cols
    ) +
    theme_minimal_hgrid(16, rel_small = 1) +
    theme(
        legend.position = "top",
        legend.justification = "right",
        legend.text = element_text(size = 14),
        legend.box.spacing = unit(0, "pt")
    )
    
ggsave('../figures/estimated_sensitivity_for_viability_studies_i123.pdf', device='pdf')

p_all<- ggplot(df, aes(x=number_of_inoculated_wells, y=p_i*100)) +
    geom_point(aes(color=inoculum, fill=inoculum), size=2.5, alpha=0.5, shape=21) +
    scale_y_log10() +
    scale_x_log10() +
    labs(x='Number of inoculated wells per experiment', y='Minimum relative abundance where zero positive wells is significant (%)') +
    geom_vline(xintercept=460, linetype='dashed') +
    scale_color_manual(
        values = darken(inoculum_cols, 0.3)
    ) +
    scale_fill_manual(
        values = inoculum_cols
    ) +
    theme_minimal_hgrid(16, rel_small = 1) +
    theme(
        legend.position = "top",
        legend.justification = "right",
        legend.text = element_text(size = 14),
        legend.box.spacing = unit(0, "pt")
    )

ggsave('../figures/estimated_sensitivity_for_viability_studies_all.pdf', device='pdf')
