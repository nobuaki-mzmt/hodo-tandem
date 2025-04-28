
# -------------------------------------------------------------------------------- #
# Hod sjo
# -------------------------------------------------------------------------------- #
{
  df <- arrow::read_feather("data_fmt/2D_hist_Hod-sjo_control_FM_20_df.feather")
  df$context <- "others"
  f_leader <- df$fTip_mHead_dis <= antenna_length["Hod-sjo"] & df$m_behind
  m_leader <- df$mTip_fHead_dis <= antenna_length["Hod-sjo"] & df$f_behind
  df$context[f_leader & !m_leader] <- "tandem_f-m"
  df$context[m_leader & !f_leader] <- "tandem_m-f"
  df$context[!m_leader & !f_leader] <- "separation"
  df_plot <- df
  
  df_plot <- rbind(
    data.frame(df_plot[,c("sex", "context")],
               relative_antennatip_x = df_plot$relative_antennatipr_x,
               relative_antennatip_y = df_plot$relative_antennatipr_y,
               side = "R"),
    data.frame(df_plot[,c("sex", "context")],
               relative_antennatip_x = df_plot$relative_antennatipl_x,
               relative_antennatip_y = df_plot$relative_antennatipl_y,
               side = "L")
  )
  df_plot$context <- factor(df_plot$context, levels = c("tandem_m-f", "tandem_f-m", "separation", "others"))
  
  
  plot_hod <- function(plot_type){
    df_plot_temp <- df_plot[df_plot$context == plot_type, ]
    ggplot(df_plot_temp, aes(x=relative_antennatip_x, y=relative_antennatip_y)) +
      stat_density_2d(aes(fill = after_stat(level)), geom="polygon")+ 
      coord_fixed(xlim=c(-8, 8), ylim=c(-2.5, 6))+
      scale_y_continuous(breaks = seq(-2, 6, 2)) +
      scale_x_continuous(breaks = seq(-8, 8, 4)) +
      theme_bw() +
      scale_fill_viridis()+
      xlab("x (mm)") + 
      ylab("y (mm)") +
      facet_grid(. ~ sex,
                 switch = "y", 
                 labeller = labeller(
                   sex = c("F" = "Female",
                           "M" = "Male"))) +
      theme(axis.title = element_text(size = 8),
            axis.text = element_text(size = 6),
            legend.position = "top",
            legend.title = element_blank(),
            legend.direction = "horizontal",
            legend.key.size = unit(0.2, "cm"),
            legend.text = element_text(size = 5),
            legend.background = element_rect(fill = "transparent", color = NA),
            legend.margin = margin(0, 0, -25, 0)  ,
            legend.spacing.y = unit(0.05, "cm"),
            legend.title.align = 0,
            strip.placement = "outside",
            strip.background = element_blank(),
            strip.text = element_text(size = 7),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()   
      )  +
      annotate(geom = "point", 
               x = mean(df$relative_antennabaser_x), 
               y = mean(df$relative_antennabaser_y), 
               color = 1, size = 0.5, shape = 1) +
      annotate(geom = "point", 
               x = mean(df$relative_antennabasel_x), 
               y = mean(df$relative_antennabasel_y), 
               color = 1, size = 0.5, shape = 1) +
      annotate(geom = "point", x = 0, y = 0, color = 2, size = 0.5, shape = 1)
    #ggsave(paste0("output/2Dantenna_dist_Hod-sjo_", plot_type, ".pdf"), width = 4, height = 2.6)
  }
  plot_hod("tandem_f-m")
  plot_hod("tandem_m-f")
  plot_hod("separation")
}
# -------------------------------------------------------------------------------- #