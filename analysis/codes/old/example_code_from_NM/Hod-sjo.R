dis_change_plot <- function(df, species, yrange = 0.1, xrange = 5, reverse = F, lag = 1){
  print(paste(species, length(unique(df$video)), "pairs detected"))
  
  df <- subset(df, sex == "M")
  #df <- df[df$frame %% (30/fps) == 0, ]
  
  if(reverse){
    df$fTip_mHead_dis <- df$mTip_fHead_dis
  }
  
  tandem_pair <- unique(df$video)[tapply(df$fTip_mHead_dis < antenna_length[species],
                                         df$video, sum, na.rm = T)/length(unique(df$frame)) > 0.5]
  df <- df[df$video %in% tandem_pair,]
  print(paste(length(unique(df$video)), "pairs showed > 0.5 tandem runs"))
  floor_frame <- min(tapply(df$frame, df$video, max))
  print(paste("minimum frame length:", floor_frame))
  df <- df[df$frame < floor_frame, ]
  
  ## fTip_mHead
  df$fTip_mHead_dis_diff <- c(diff(df$fTip_mHead_dis, lag = 30), rep(NA, 30))
  df$fTip_mHead_dis_diff[df$frame > floor_frame - 30] <- NA
  
  df$fTip_mHead_dis <- round(df$fTip_mHead_dis, 1)
  df <- subset(df, fTip_mHead_dis <= xrange)
  df$approch <- (df$fTip_mHead_dis_diff < 0)*1
  
  if(F){ # this is for distance change
    p1 <- ggplot(df)+
      geom_hline(yintercept = 0, linetype = 2, linewidth = 1, alpha = .6) + 
      geom_vline(xintercept = antenna_length[species]) + 
      stat_summary(aes(x = fTip_mHead_dis, y = fTip_mHead_dis_diff),
                   fun = 'mean', geom = 'line') +
      stat_summary(aes(x = fTip_mHead_dis, y = fTip_mHead_dis_diff),
                   fun = 'mean', geom = 'point') +
      stat_summary(aes(x = fTip_mHead_dis, y = fTip_mHead_dis_diff),
                   fun.data = 'mean_cl_normal', geom = 'ribbon', alpha = 0.2) +
      theme_classic()+
      coord_cartesian(xlim=c(0, xrange), ylim=c(-yrange, yrange)) +
      scale_y_continuous(breaks = c(-yrange, 0, yrange)) +
      xlab("Leader-follower distance (mm)")+
      ylab("Distance change (mm)")+
      theme(axis.title = element_text(size = 9),
            axis.text = element_text(size = 7),
            aspect.ratio = 1/1.6180339887,
            panel.background = element_blank()
      ) 
  }
  
  p1 <- ggplot(df)+
    geom_hline(yintercept = .5, linetype = 2, linewidth = 1, alpha = .6) + 
    geom_vline(xintercept = antenna_length[species]) + 
    geom_vline(xintercept = palp_length[species]) + 
    stat_summary(aes(x = fTip_mHead_dis, y = approch),
                 fun = 'mean', geom = 'line') +
    stat_summary(aes(x = fTip_mHead_dis, y = approch),
                 fun = 'mean', geom = 'point') +
    stat_summary(aes(x = fTip_mHead_dis, y = approch),
                 fun.data = 'mean_cl_normal', geom = 'ribbon', alpha = 0.2) +
    theme_classic()+
    coord_cartesian(xlim=c(0, xrange), ylim=c(0, 1)) +
    scale_y_continuous(breaks = c(0, 0, 1)) +
    xlab("Leader-follower distance (mm)")+
    ylab("Probability to approach")+
    theme(axis.title = element_text(size = 9),
          axis.text = element_text(size = 7),
          aspect.ratio = 1/1.6180339887,
          panel.background = element_blank()
    ) 
  p2 <- ggplot(df, aes(x=fTip_mHead_dis))+
    geom_histogram(binwidth = 0.1, fill = "#333333", color="black")+
    geom_vline(xintercept = antenna_length[species]) + 
    geom_vline(xintercept = palp_length[species]) + 
    ylab("Count")+
    theme_classic()+
    scale_y_log10()+
    coord_cartesian(xlim=c(0, xrange)) +
    #scale_y_continuous(labels = scales::label_comma()) +
    theme(axis.text = element_text(size = 7),
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          aspect.ratio = 1/5,
          panel.background = element_blank()
    )
  g1 <- ggplotGrob(p1)
  g2 <- ggplotGrob(p2)
  g <- rbind(g2, g1, size = "first")
  g$widths <- unit.pmax(g1$widths, g2$widths)
  return(list(g, df))
  #return(g)
}

# Hod sjo
{
  df <- arrow::read_feather("data_fmt/antenna_angle_Hod-sjo_control_FM_20_df.feather")
  
  g <- dis_change_plot(df, "Hod-sjo", yrange = .25)
  grid.newpage()
  grid.draw(g[[1]])
  ggsave(paste0("output/dis_change_Hod-sjo_fleader.pdf"),
         g[[1]], width=3, height=4, bg = "transparent")
  
  g <- dis_change_plot(df, "Hod-sjo", yrange = .25, reverse = T)
  grid.newpage()
  grid.draw(g[[1]])
  ggsave(paste0("output/dis_change_Hod-sjo_mleader.pdf"),
         g[[1]], width=3, height=4, bg = "transparent")
}


tandem_process <- function(df, species, experiments, antenna, reverse = F, fps = 5){
  df <- subset(df, sex == "M")
  df <- df[df$frame %% (30/fps) == 0, ] # down sample into 5 FPS
  if(reverse){
    df$tandem <- df$mTip_fHead_dis < antenna_length[species]*1.1
  } else {
    df$tandem <- df$fTip_mHead_dis < antenna_length[species]*1.1
  }
  
  df_sum_temp <- data.frame(
    video = unique(df$video),
    species, experiments, antenna, reverse,
    tandem_prop = tapply(df$tandem, df$video, sum)/length(unique(df$frame))
  )
  
  tandem_pair <- unique(df$video)[tapply(df$tandem, df$video, sum, na.rm = T)/length(unique(df$frame)) > 0.5]
  df <- df[df$video %in% tandem_pair,]
  
  df_tandem = NULL; df_sep = NULL;
  if(dim(df)[1] > 0){ 
    videos <- unique(df$video)
    for(i_video in videos){
      df_temp <- subset(df, video == i_video)
      tandem = df_temp$tandem
      tandem.smoothing <- function(vec, min.sec){ 
        if(sum(vec)>0){
          timing <- which(vec)[c(T, diff(which(vec))>1)]
          end    <- which(vec)[c(diff(which(vec))>1,T)]
          for(fi in 1:length(timing)){
            if(length( vec[timing[fi]:end[fi]]) < min.sec ){
              vec[timing[fi]:end[fi]] <- F
            }
          }
        }
        return(vec)
      }
      tandem <- tandem.smoothing(tandem, fps * 2)
      if(sum(tandem)>0){
        tan.end <- which(tandem)[c(diff(which(tandem))>1,T)]
        tan.sta <- which(tandem)[c(T, diff(which(tandem))>1)]
        tan_duration <- tan.end - tan.sta + 1
        sep_end <- tan.sta[-1]
        sep_sta <- tan.end[-length(tan.end)]
        sep_duration <- sep_end - sep_sta + 1
        
        tan_cens <- rep(1, length(tan_duration))
        tan_cens[tan.sta == 1 | tan.end == max(df_temp$frame)] <- 0
        
        sep_dis <- NULL
        for(i_s in 1:length(sep_sta)){
          sep_dis <- c(sep_dis, max(df_temp[(sep_sta[i_s]):sep_end[i_s], ]$fTip_mHead_dis))
        }
        
        # data summarize  
        df_tandem_temp <- data.frame(
          name =  i_video,
          species,
          experiments,
          antenna,
          reverse,
          tan_duration = tan_duration * 6/30,
          tan_cens
        )
        df_tandem <- rbind(df_tandem, df_tandem_temp)
        
        df_sep_temp <- data.frame(
          name =  i_video,
          species,
          experiments,
          antenna,
          reverse,
          sep_duration = sep_duration * 6/30,
          sep_dis
        )
        df_sep <- rbind(df_sep, df_sep_temp)
      }
    }
  }
  
  return(list(df_sum_temp, df_tandem, df_sep))
}

df <- arrow::read_feather("data_fmt/antenna_angle_Hod-sjo_control_FM_20_df.feather")
df_list <- tandem_process(df, "Hod-sjo", "control", "20")
df_sum <- rbind(df_sum, df_list[[1]]); df_tandem <- rbind(df_tandem, df_list[[2]]); df_sep <- rbind(df_sep, df_list[[3]])
df_list <- tandem_process(df, "Hod-sjo", "control", "20", reverse = T)
df_sum <- rbind(df_sum, df_list[[1]]); df_tandem <- rbind(df_tandem, df_list[[2]]); df_sep <- rbind(df_sep, df_list[[3]])

# -------------------------------------------------------------------------------- #
# Hod sjo
# -------------------------------------------------------------------------------- #
{
  y_range <- 0.2
  
  df <- arrow::read_feather("data_fmt/antenna_angle_Hod-sjo_control_FM_20_df.feather")
  
  # focus on pairs that show tandem runs for more than 50% of the time
  #tandem_pair <- unique(df$video)[tapply(df$fTip_mHead_dis < antenna_length["Hod-sjo"]*1.1, df$video, sum)/(max(df$frame)/(30/analyze_FPS)) > 0.5]
  #df <- df[df$video %in% tandem_pair,]
  
  ## fTip_mHead
  df$antenna_r_diff <- c(abs(diff(df$angle_antenna_r)), NA)
  df$antenna_r_diff[which(df$frame==0)-1] = NA
  df$antenna_r_diff[df$antenna_r_diff>pi & !is.na(df$antenna_r_diff)] <- 
    2*pi - df$antenna_r_diff[df$antenna_r_diff>pi & !is.na(df$antenna_r_diff)]
  
  df$antenna_l_diff <- c(abs(diff(df$angle_antenna_l)), NA)
  df$antenna_l_diff[which(df$frame==0)-1] = NA
  df$antenna_l_diff[df$antenna_l_diff>pi & !is.na(df$antenna_l_diff)] <-
    2*pi - df$antenna_l_diff[df$antenna_l_diff>pi & !is.na(df$antenna_l_diff)]
  
  df$fTip_mHead_dis <- round(df$fTip_mHead_dis, 1)
  df$mTip_fHead_dis <- round(df$mTip_fHead_dis, 1)
  
  df_tandem <- subset(df, fTip_mHead_dis <= 5 & mTip_fHead_dis > 5 & m_behind)
  linetypes <- c("R" = "solid", "L" = "dashed")
  ggplot(df_tandem, aes(fill=sex))+
    geom_vline(xintercept = antenna_length[2]) + 
    stat_summary(aes(x=fTip_mHead_dis, y = antenna_r_diff, col=sex, linetype = "R"),
                 fun = 'mean', 
                 geom = 'line') +
    stat_summary(aes(x=fTip_mHead_dis, y = antenna_r_diff),
                 fun.data = 'mean_cl_normal',
                 geom = 'ribbon', alpha = 0.2) +
    stat_summary(aes(x=fTip_mHead_dis, y = antenna_l_diff, col=sex, linetype = "L"),
                 fun = 'mean', 
                 geom = 'line') +
    stat_summary(aes(x=fTip_mHead_dis, y = antenna_l_diff),
                 fun.data = 'mean_cl_normal',
                 geom = 'ribbon', alpha = 0.2) +
    coord_cartesian(xlim=c(0, 5.2), ylim=c(0, y_range))+
    scale_fill_viridis(discrete = T, labels = c("Female", "Male"), direction = -1)+
    scale_color_viridis(discrete = T, labels = c("Female", "Male"), direction = -1)+
    scale_linetype_manual(values = linetypes)+
    theme_classic() +
    theme(
      legend.margin = margin(0, 0, 0, 0)  ,
      legend.position = c(0.2, 0.85),
      legend.title = element_blank(),
      legend.text = element_text(size=5),
      legend.direction = "vertical",
      legend.key.width = unit(0.6, "cm"),
      legend.key.height = unit(0.3, "cm"),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.spacing.x = unit(0, "cm"),
      legend.spacing.y = unit(0, "cm"),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 6),
      panel.background = element_blank()
    )  +
    xlab("Leader-Follower distance (mm)") +
    ylab("Antenna movement (rad)")
  
  ggsave("output/Hod-sjo_fleader_antenna_move.pdf", width=3, height=2.5, bg = "transparent")
  
  tandem_pair <- unique(df$video)[tapply(df$fTip_mHead_dis < antenna_length["Hod-sjo"]*1.1, df$video, sum)/(max(df$frame)/(30/analyze_FPS)) > 0.5]
  df <- df[df$video %in% tandem_pair,]
  
  df$tandem <- 0
  df$tandem[df$mTip_fHead_dis > antenna_length["Hod-sjo"] & df$fTip_mHead_dis <= antenna_length["Hod-sjo"] & df$f_behind] <- "f_leader"
  df$tandem[df$mTip_fHead_dis > antenna_length["Hod-sjo"] & df$fTip_mHead_dis <= 0.7 & df$f_behind] <- "f_leader2"
  df$tandem[df$mTip_fHead_dis > antenna_length["Hod-sjo"] & df$fTip_mHead_dis > antenna_length["Hod-sjo"]] <- "sep"
  
  
  df_antenna_move_stat <- do.call(rbind, list(
    calculate_mean(df, "antenna_r_diff", c("video", "sex", "tandem"), "R", "FM"),
    calculate_mean(df, "antenna_l_diff", c("video", "sex", "tandem"), "L", "FM")
  ))
  colnames(df_antenna_move_stat)[4] <- "antenna_move"
  df_antenna_move_stat$colony <- sub("^[^_]*_[^_]*_[^_]*_([^_]*).*", "\\1", df_antenna_move_stat$video)
  
  df_antenna_move_stat <- df_antenna_move_stat[df_antenna_move_stat$tandem != "0",]
  
  df_antenna_move_stat$tandem <- factor(df_antenna_move_stat$tandem, levels = c("m_leader","m_leader2", "f_leader","f_leader2", "sep"))
  ggplot(subset(df_antenna_move_stat, tandem == "f_leader" | tandem == "f_leader2" | tandem == "sep"), 
         aes(x = tandem, y = antenna_move, color = sex, fill = sex)) +
    geom_boxplot(width = 0.4, outlier.shape = NA, position = position_dodge(0.8), alpha=.25, col= 1) +  # Boxplot with dodge
    geom_point(alpha = 0.4, size = .75, position = position_jitterdodge(0.2)) +  # Jitter with dodge
    scale_color_viridis_d(direction = -1, end = 0.95) +  # Viridis color palette
    scale_fill_viridis_d(direction = -1, end = 0.95) +  # Viridis color palette
    coord_cartesian(ylim=c(0,0.3))+
    labs(x = "", y = "Antenna Movement (rad/sec)", color = NULL) +  # Axis labels and no legend title
    theme_classic() +  # Classic theme
    theme(
      legend.position = c(0.5, 0.9),  # Position legend at the top center
      legend.direction = "horizontal",  # Make the legend horizontal
      legend.box = "horizontal"  # Ensure legend items are aligned horizontally
    )
  
  r <- lmer(antenna_move ~ tandem + (1|video), data = subset(df_antenna_move_stat, sex == "F"))
  Anova(r)
  summary(glht(r, linfct = mcp(tandem = "Tukey")))
  
  df <- arrow::read_feather("data_fmt/antenna_angle_Hod-sjo_control_FM_20_df.feather")
  # focus on pairs that show tandem runs for more than 50% of the time
  #tandem_pair <- unique(df$video)[tapply(df$mTip_fHead_dis < antenna_length["Hod-sjo"]*1.1, df$video, sum)/(max(df$frame)/(30/analyze_FPS))/2 > 0.3]
  #df <- df[df$video %in% tandem_pair,]
  
  df$antenna_r_diff <- c(abs(diff(df$angle_antenna_r)), NA)
  df$antenna_r_diff[which(df$frame==0)-1] = NA
  df$antenna_r_diff[df$antenna_r_diff>pi & !is.na(df$antenna_r_diff)] <- 
    2*pi - df$antenna_r_diff[df$antenna_r_diff>pi & !is.na(df$antenna_r_diff)]
  
  df$antenna_l_diff <- c(abs(diff(df$angle_antenna_l)), NA)
  df$antenna_l_diff[which(df$frame==0)-1] = NA
  df$antenna_l_diff[df$antenna_l_diff>pi & !is.na(df$antenna_l_diff)] <-
    2*pi - df$antenna_l_diff[df$antenna_l_diff>pi & !is.na(df$antenna_l_diff)]
  
  df$fTip_mHead_dis <- round(df$fTip_mHead_dis, 1)
  df$mTip_fHead_dis <- round(df$mTip_fHead_dis, 1)
  
  
  df_tandem <- subset(df, fTip_mHead_dis > 5 & mTip_fHead_dis <= 5 & f_behind)
  ggplot(df_tandem, aes(fill=sex))+
    geom_vline(xintercept = antenna_length[2]) + 
    stat_summary(aes(x=mTip_fHead_dis, y = antenna_r_diff, col=sex, linetype = "R"),
                 fun = 'mean', 
                 geom = 'line') +
    stat_summary(aes(x=mTip_fHead_dis, y = antenna_r_diff),
                 fun.data = 'mean_cl_normal',
                 geom = 'ribbon', alpha = 0.2) +
    stat_summary(aes(x=mTip_fHead_dis, y = antenna_l_diff, col=sex, linetype = "L"),
                 fun = 'mean', 
                 geom = 'line') +
    stat_summary(aes(x=mTip_fHead_dis, y = antenna_l_diff),
                 fun.data = 'mean_cl_normal',
                 geom = 'ribbon', alpha = 0.2) +
    coord_cartesian(xlim=c(0, 5.2), ylim=c(0, y_range))+
    scale_fill_viridis(discrete = T, labels = c("Female", "Male"), direction = -1)+
    scale_color_viridis(discrete = T, labels = c("Female", "Male"), direction = -1)+
    scale_linetype_manual(values = linetypes)+
    theme_classic() +
    theme(
      legend.margin = margin(0, 0, 0, 0)  ,
      legend.position = c(0.2, 0.85),
      legend.title = element_blank(),
      legend.text = element_text(size=5),
      legend.direction = "vertical",
      legend.key.width = unit(0.6, "cm"),
      legend.key.height = unit(0.3, "cm"),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.spacing.x = unit(0, "cm"),
      legend.spacing.y = unit(0, "cm"),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 6),
      panel.background = element_blank()
    )  +
    xlab("Leader-Follower distance (mm)") +
    ylab("Antenna movement (rad)")
  
  ggsave("output/Hod-sjo_mleader_antenna_move.pdf", width=3, height=2.5, bg = "transparent")
  
  tandem_pair <- unique(df$video)[tapply(df$mTip_fHead_dis < antenna_length["Hod-sjo"]*1.1, df$video, sum)/(max(df$frame)/(30/analyze_FPS)) > 0.5]
  df <- df[df$video %in% tandem_pair,]
  
  df$tandem <- 0
  df$tandem[df$fTip_mHead_dis > antenna_length["Hod-sjo"] & df$mTip_fHead_dis <= antenna_length["Hod-sjo"] & df$f_behind] <- "m_leader"
  df$tandem[df$fTip_mHead_dis > antenna_length["Hod-sjo"] & df$mTip_fHead_dis <= 0.7 & df$f_behind] <- "m_leader2"
  df$tandem[df$mTip_fHead_dis > antenna_length["Hod-sjo"] & df$fTip_mHead_dis > antenna_length["Hod-sjo"]] <- "sep"
  
  
  df_antenna_move_stat <- do.call(rbind, list(
    calculate_mean(df, "antenna_r_diff", c("video", "sex", "tandem"), "R", "FM"),
    calculate_mean(df, "antenna_l_diff", c("video", "sex", "tandem"), "L", "FM")
  ))
  colnames(df_antenna_move_stat)[4] <- "antenna_move"
  df_antenna_move_stat$colony <- sub("^[^_]*_[^_]*_[^_]*_([^_]*).*", "\\1", df_antenna_move_stat$video)
  
  df_antenna_move_stat <- df_antenna_move_stat[df_antenna_move_stat$tandem != "0",]
  
  df_antenna_move_stat$tandem <- factor(df_antenna_move_stat$tandem, levels = c("m_leader","m_leader2", "f_leader","f_leader2", "sep"))
  ggplot(subset(df_antenna_move_stat, tandem == "m_leader" | tandem == "m_leader2" | tandem == "sep"), 
         aes(x = tandem, y = antenna_move, color = sex, fill = sex)) +
    geom_boxplot(width = 0.4, outlier.shape = NA, position = position_dodge(0.8), alpha=.25, col= 1) +  # Boxplot with dodge
    geom_point(alpha = 0.4, size = .75, position = position_jitterdodge(0.2)) +  # Jitter with dodge
    scale_color_viridis_d(direction = -1, end = 0.95) +  # Viridis color palette
    scale_fill_viridis_d(direction = -1, end = 0.95) +  # Viridis color palette
    coord_cartesian(ylim=c(0,0.3))+
    labs(x = "", y = "Antenna Movement (rad/sec)", color = NULL) +  # Axis labels and no legend title
    theme_classic() +  # Classic theme
    theme(
      legend.position = c(0.5, 0.9),  # Position legend at the top center
      legend.direction = "horizontal",  # Make the legend horizontal
      legend.box = "horizontal"  # Ensure legend items are aligned horizontally
    )
  
  
  r <- lmer(antenna_move ~ tandem + (1|video), data = subset(df_antenna_move_stat, sex == "F"))
  Anova(r)
  summary(glht(r, linfct = mcp(tandem = "Tukey")))
}
