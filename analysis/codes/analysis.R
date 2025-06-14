# data formatting and analysis on Hodotermopsis tandem project
# created by William Chambliss
# eddited / commented by Elijah Carroll, Nobuaki Mizumoto

# --------------------------------------------------------------- #
# PACKAGES + PARAMETERS
# --------------------------------------------------------------- #
{
  rm(list = ls())
  
  library(stringr)
  library(data.table)
  library(arrow)
  library(dplyr)
  library(tidyr)
  
  library(ggplot2)
  library(viridis)
  library(patchwork)
  library(zoo)
  library(cowplot)
  
  library(knitr)
  
  library(MASS)
  library(survival)
  library(survminer)
  library(coxme)
  library(car)
  
  # parameters for defininng tandem running
  {
    # Tandem Running: Coordinated behavior in which one individual (the leader) walks ahead and is closely followed by another individual (the follower), maintaining continuous physical or near-physical contact, often via antennal or body-tip touching. 
    
    # Maximum angle difference (in radians) allowed between male and female headings during tandem
    tandemAngle = 60 * (pi / 180)  # 60 degrees converted to radians
    
    # Minimum number of consecutive frames (at 5 FPS) required to classify a period as valid tandem or separation behavior.
    # Used in the tandem.smoothing() function to remove brief, potentially spurious transitions in both tandem and close/separation states.
    # Segments shorter than tandemsmooth frames are flipped to match the surrounding state (e.g., a brief tandem blip during separation is smoothed out).
    # This ensures continuity in behavioral classification and avoids overcounting interruptions or false transitions.
    # Example: tandemsmooth = 20 → minimum duration = 4 seconds.
    tandemsmooth = 20
    
    # Same as tandemsmooth, but used for smoothing male/female leadership states within tandem runs.
    leadsmooth = 15  # 15 frames = 3 seconds minimum for a valid lead
    
    # Minimum speed (mm/s) required for both individuals to be considered "moving" in a tandem run. 
    tandemSpeed = 1.213 * 0.96
    
    # Distance threshold for defining proximity: individuals are considered “close” if their centers are within
    # 0.6 × the sum of their body lengths. This is a prerequisite for tandem classification.
    tandem_threshold_dis = 0.6
  }
  
} 
# --------------------------------------------------------------- #

# --------------------------------------------------------------- #
# Tandem Smoothing Function
# --------------------------------------------------------------- #
{
  tandem.smoothing <- function(vec, threshold) {
    r <- rle(vec)
    cs <- cumsum(r$lengths)
    starts <- c(1, head(cs, -1) + 1)
    for(i in seq_along(r$lengths)) {
      if(r$lengths[i] < threshold) {
        vec[starts[i]:cs[i]] <- !r$values[i]
      }
    }
    return(vec)
  }
}
# --------------------------------------------------------------- #

# --------------------------------------------------------------- #
# Data Formatting / Storing into rda, downsampling, and set to mm
# --------------------------------------------------------------- #
{
  ### Upload data and define columns  
  dataset <- ("data_fmt/data_raw_df.feather")
  
  # Initialize data frame and load data
  # df_all is the original data, while df is the original data along with all calculation inferences from the data
  df_all <- data.frame()
  
  print(paste("working on", dataset))
    
  df <- arrow::read_feather(dataset)
  colnames(df)[1] <- "frame"

  ### Get bodysize data 
  df_body <- fread("data_fmt/data_raw_bodysize.csv")
  for (i in 1:nrow(df_body)) {
    if (grepl("150", df_body[i, 1])) {
      df_body[i, 2:3] <- df_body[i, 2:3] / 2000 * 143
    } else {
      df_body[i, 2:3] <- df_body[i, 2:3] / 2000 * 87
    }
  }
  df_body_scaled <- df_body

  ####
  # 1. scale it into mm and second and downsample it in 5 FPS
  ## Conversion from pixels to mm##
  df$scaling_factor <- ifelse(grepl("150", df$video), 143, 87)
  df$fx <- df$fx / 2000 * df$scaling_factor
  df$fy <- df$fy / 2000 * df$scaling_factor
  df$mx <- df$mx / 2000 * df$scaling_factor
  df$my <- df$my / 2000 * df$scaling_factor
  
  ## Scaling body parts
  df$fx_abdomen <- df$fx_abdomen / 2000 * df$scaling_factor
  df$fy_abdomen <- df$fy_abdomen / 2000 * df$scaling_factor
  df$mx_abdomen <- df$mx_abdomen / 2000 * df$scaling_factor
  df$my_abdomen <- df$my_abdomen / 2000 * df$scaling_factor
  
  df$fx_headtip <- df$fx_headtip / 2000 * df$scaling_factor
  df$fy_headtip <- df$fy_headtip / 2000 * df$scaling_factor
  df$mx_headtip <- df$mx_headtip / 2000 * df$scaling_factor
  df$my_headtip <- df$my_headtip / 2000 * df$scaling_factor
  
  ## Down-sampling by specifying data at every 5 FPS##
  colnames(df)[1] <- 'time'
  df <- df[df$time %% 6 == 0, ]
  df$time <- df$time / 30
    
  df_all <- rbind(df_all, data.frame(df))

  ### Save dataframes and use to obtain speed and stability
  save(df_all, file="data_fmt/df_all.rda")
  save(df_body_scaled, file="data_fmt/df_body_scaled.rda")
}
# --------------------------------------------------------------- #

# --------------------------------------------------------------- #
# Data Set up
# --------------------------------------------------------------- #
{
  load("data_fmt/df_all.rda")
  load("data_fmt/df_body_scaled.rda")
  
  df <- df_all
  
  #### SPEED ####
  df$f_spe <- c(0, sqrt( diff(df$fx)^2 + diff(df$fy)^2 ))*5
  df$m_spe <- c(0, sqrt( diff(df$mx)^2 + diff(df$my)^2 ))*5
  df$f_spe[df$time == 0] <- 0
  df$m_spe[df$time == 0] <- 0
  
  #### TANDEM + CLOSE ####
  body_size_f <- mean(df_body_scaled$female)
  body_size_m <- mean(df_body_scaled$male)
  tandem_threshold = (body_size_f + body_size_m) * tandem_threshold_dis
  df$partner_dis = sqrt((df$fy - df$my)^2 + (df$fx - df$mx)^2)
  
  df$heading_f <- atan2(df$fy_headtip - df$fy_abdomen, df$fx_headtip - df$fx_abdomen)
  df$heading_m <- atan2(df$my_headtip - df$my_abdomen, df$mx_headtip - df$mx_abdomen)
  heading_diff <- abs((df$heading_f - df$heading_m + pi) %% (2 * pi) - pi)
  
  df$close <- (df$partner_dis < tandem_threshold)
  df$tandem <- (df$f_spe > tandemSpeed)&(df$m_spe > tandemSpeed)&(df$close) &(heading_diff < tandemAngle)
  df$tandem[is.na(df$tandem)] <- FALSE
  
  ### what is the average duration of male led tandem runs vs what is the average duration of female led tandem runs
  {
    tandem_df <- data.frame(tandem = logical(0))
    close_df <- data.frame(close = logical(0))
    lead_df <- data.frame(m_lead = logical(0),f_lead = logical(0))
    df_list <- list()
    videos <- unique(df$video)
    
    for(i_v in 1:length(videos)){ 
      
      df_temp <- subset(df, video == videos[i_v])
      
      ## Applying Tandem Smoothing to df
      repeat {
        old_rle <- rle(df_temp$tandem)
        df_temp$tandem <- tandem.smoothing(df_temp$tandem, tandemsmooth)
        new_rle <- rle(df_temp$tandem) 
        if (identical(old_rle$lengths, new_rle$lengths)) break
      }
      tandem_df <- rbind(tandem_df, data.frame(tandem = df_temp$tandem))
      
      repeat {
        old_rle <- rle(df_temp$close)
        df_temp$close <- tandem.smoothing(df_temp$close, tandemsmooth)
        new_rle <- rle(df_temp$close) 
        if (identical(old_rle$lengths, new_rle$lengths)) break
      }
      close_df <- rbind(close_df, data.frame(close = df_temp$close))
      
      ## Male and female lead for each frame/smoothing leaders
      {
        df_temp$m_lead <- FALSE
        df_temp$f_lead <- FALSE
        for (f in seq_len(nrow(df_temp))) {
          if (df_temp$tandem[f] == TRUE) {
            mTof <- sqrt((df_temp$mx_headtip[f] - df_temp$fx_abdomen[f])^2 + 
                           (df_temp$my_headtip[f] - df_temp$fy_abdomen[f])^2)
            fTom <- sqrt((df_temp$fx_headtip[f] - df_temp$mx_abdomen[f])^2 +
                           (df_temp$fy_headtip[f] - df_temp$my_abdomen[f])^2)
            if (mTof > fTom) {
              df_temp$m_lead[f] <- TRUE
              df_temp$f_lead[f] <- FALSE
            } else {
              df_temp$f_lead[f] <- TRUE
              df_temp$m_lead[f] <- FALSE
            }
          } else {
            df_temp$m_lead[f] <- FALSE
            df_temp$f_lead[f] <- FALSE
          }
        }
        
        # smoothing leaders
        repeat {
          old_rle <- rle(df_temp$m_lead)  # Store previous RLE state
          df_temp$m_lead <- tandem.smoothing(df_temp$m_lead, leadsmooth)
          new_rle <- rle(df_temp$m_lead)  # Check new state
          if (identical(old_rle$lengths, new_rle$lengths))
            break  # Stop when no changes occur
        }
        
        repeat {
          old_rle <- rle(df_temp$f_lead)  # Store previous RLE state
          df_temp$f_lead <- tandem.smoothing(df_temp$f_lead, leadsmooth)
          new_rle <- rle(df_temp$f_lead)  # Check new state
          if (identical(old_rle$lengths, new_rle$lengths))
            break  # Stop when no changes occur
        }
        
        lead_df <- rbind(lead_df,
                         data.frame(m_lead = df_temp$m_lead, f_lead = df_temp$f_lead))
      }
        
      tandem <- df_temp$tandem
      
      # Tandem calculations
      tan.end <- which(tandem)[c(diff(which(tandem)) > 1, T)]
      tan.sta <- which(tandem)[c(T, diff(which(tandem)) > 1)]
      tan_duration <- (tan.end - tan.sta) / 5
      tan_cens <- tan.end != dim(df_temp)[1] # NM: add if tandem was completed (not ended due to the end of the video)
      
      # Separation calculations
      separation <- !tandem
      sep.end <- which(separation)[c(diff(which(separation)) > 1, T)]
      sep.sta <- which(separation)[c(T, diff(which(separation)) > 1)]
      sep_duration <- (sep.end - sep.sta) / 5
      sep_cens <- sep.end != dim(df_temp)[1] # NM: add if tandem was completed (not ended due to the end of the video)
      
      # Create data frames for tandem and separation
      df_tandem <- data.frame(
        video = videos[i_v],
        tandem = TRUE,
        duration = tan_duration,
        m_lead = df_temp$m_lead[tan.sta],
        f_lead = df_temp$f_lead[tan.sta],
        cens   = tan_cens
      )
      
      df_separation <- data.frame(
        video = videos[i_v],
        tandem = FALSE,
        duration = sep_duration,
        m_lead = NA,
        f_lead = NA,
        cens    = sep_cens
      )
      
      df_list[[i_v]] <- rbind(df_tandem, df_separation)
    }
    
    df_final <- do.call(rbind, df_list)
    df$tandem <- tandem_df$tandem
    df$close  <- close_df$close
    df$close  <- ifelse(df$tandem, FALSE, df$close)
    df$m_lead <- lead_df$m_lead
    df$f_lead <- lead_df$f_lead
    df$lead <- ifelse(df$m_lead, "male",
                      ifelse(df$f_lead, "female", NA))
  }
  
  ### how frequently pairs switch leading roles during tandem running.
  {
    df$switch = FALSE
    df_switches <- data.frame(video = character(0), switches = integer(0))
    frame=0
    # NM: This is minor comment but we may stick with i_v in 1:length(videos) as other parts do so instead of v in videos. it is nice to have consistency in the same file.
    for(v in videos){
      video_rows <- which(df$video == v)
      m_lead = 0
      f_lead = 0
      switches = 0
      previous_leader = ""
      for(j in video_rows){
        frame = frame + 1
        if(df$tandem[j] == TRUE){
          current_leader = df$lead[j]
          # Increment switch count if the leader changes
          if (!is.na(previous_leader) && !is.na(current_leader)&&previous_leader != "" && current_leader != previous_leader) {
            switches = switches + 1
            df$switch[frame] = TRUE
          }
          previous_leader = current_leader
        }
      }
      df_lead_temp <- data.frame(
        video = v,
        switches = switches
      )
      df_switches <- rbind(df_switches, df_lead_temp)
    }
  }
  
  ### finding the last led for seperation events
  {
    df$lastled <- NA
    
    # NM: This is good, but I wonder why several different approaches were used for similar treatments here and there. Here we use lapply but we use simple for loop for other places. It is better to make it consistent.
    video_list <- split(df, df$video)
    processed_list <- lapply(video_list, function(sub_df) {
      sub_df <- sub_df[order(sub_df$time), ]
      last_leader <- NA
      for (i in 1:nrow(sub_df)) {
        if (sub_df$time[i] == 0) {
          sub_df$lastled[i] <- NA
        } else {
          if (is.na(sub_df$lead[i])) {
            sub_df$lastled[i] <- last_leader
          } else {
            sub_df$lastled[i] <- NA
            last_leader <- sub_df$lead[i]
          }
        }
      }
      return(sub_df)
    })
    df <- do.call(rbind, processed_list)
    df <- df[order(df$video, df$time), ]
  }
  
  ## Compare tandem run duration between different dish sizes.
  df_filtered150 <- df_final[df_final$tandem == TRUE & grepl("150", df_final$video), ]
  df_filtered90 <- df_final[df_final$tandem == TRUE & grepl("90", df_final$video), ]
  average150RunDur <- sum(df_filtered150$duration) / nrow(df_filtered150)
  average90RunDur <- sum(df_filtered90$duration) / nrow(df_filtered90)
  
  ## Compare tandem speed between different dish sizes.
  df_filteredspeed150 <- df[df$tandem == TRUE & grepl("150", df$video), ]
  df_filteredspeed90 <- df[df$tandem == TRUE & grepl("90", df$video), ]
  favg_speed90 <- mean(df_filteredspeed90$f_spe, na.rm = TRUE)
  favg_speed150 <- mean(df_filteredspeed150$f_spe, na.rm = TRUE)
  mavg_speed90 <- mean(df_filteredspeed90$m_spe, na.rm = TRUE)
  mavg_speed150 <- mean(df_filteredspeed150$m_spe, na.rm = TRUE)
  
  # Tandem runs Seperations
  df_tandem_runs <- subset(df_final, tandem == TRUE)
  df_tandem_runs$lead <- ifelse(df_tandem_runs$m_lead, "male",
                             ifelse(df_tandem_runs$f_lead, "female", NA))
  df_tandem_runs$size <- ifelse(grepl("150", df_tandem_runs$video), "150", 
                             ifelse(grepl("90", df_tandem_runs$video), "90", NA))
  
  # seperations duration by sex set up
  df_sep_runs <- df[c("video", "lastled", "time")]
  df_sep_runs$malelast   = FALSE
  df_sep_runs$femalelast = FALSE
  df_sep_runs <- df_sep_runs %>%
    mutate(lastled = replace_na(lastled, "nothing"))
  df_sep_runs <- df_sep_runs %>%
    mutate(malelast = as.character(lastled) == "male")
  df_sep_runs <- df_sep_runs %>%
    mutate(femalelast = lastled == "female")
  df_sep_runs$malelast[is.na(df_sep_runs$malelast)] <- FALSE
  df_sep_runs$malelast[is.na(df_sep_runs$femalelast)] <- FALSE
  
  df_sepfin = df_sepvid <- NULL
  for(i_v in unique(df_sep_runs$video)){
    df_temp <- subset(df_sep_runs, video == i_v)
    df_malesep = rle(df_temp$malelast)
    df_malesep <- data.frame(
      mlast   = df_malesep$values,
      lengths = df_malesep$lengths/5, # NM: convert to second
      censor  = c(rep(1, length(df_malesep$values)-1), 0), # NM: if it is the very last event. The event is censored due to the end of the observation, rather than true end of the behavior
      video   = i_v,
      size    = ifelse(grepl("150", i_v), "150mm", "90mm"),
      leader  = "Male"
    )
    df_malesep <- df_malesep %>%
      filter(mlast != FALSE)
    
    df_fmalesep = rle(df_temp$femalelast)
    df_fmalesep <- data.frame(
      mlast   = df_fmalesep$values,
      lengths = df_fmalesep$lengths/5, # NM: convert to second
      censor  = c(rep(1, length(df_fmalesep$values)-1), 0), # NM: if it is the very last event. The event is censored due to the end of the observation, rather than true end of the behavior
      video   = i_v,
      size    = ifelse(grepl("150", i_v), "150mm", "90mm"),
      leader  = "Female"
    )
    df_fmalesep <- df_fmalesep %>%
      filter(mlast != FALSE)
    
    df_sepfin <- rbind(df_sepfin, rbind(df_malesep, df_fmalesep))
    df_sepvid <- rbind(df_sepvid, rbind(df_malesep, df_fmalesep))
    
  }
  
  # Box plot of speed in tandem based on dish size Set Up
  df_tan <- df %>%
    filter(tandem == TRUE) %>%
    mutate(size = ifelse(grepl("150", video), "150 mm", "90 mm"))
  
  df_sum <- data.frame()
  for(i_v in seq_along(videos)){
    v <- videos[i_v]
    df_tem <- subset(df_tan, video == v)
    tandem_index <- df_tem$tandem
    df_sum_temp <- data.frame(
      video            = df_tem$video[i_v],
      size            = df_tem$size[i_v],
      f_speed         = mean(df_tem$f_spe, na.rm = TRUE),
      m_speed         = mean(df_tem$m_spe, na.rm = TRUE),
      leader_speed    = mean(c(df_tem$f_spe[df_tem$f_lead], df_tem$m_spe[df_tem$m_lead]))
    )
    df_sum <- rbind(df_sum, df_sum_temp)
  }
  
  # Box plot of switches based on dish size Set Up
  df_switches2 <- df_switches %>%
    mutate(
      dish_size = ifelse(grepl("150", video), "150 mm", "90 mm")
    )
  
  save(df,              # seems like we need to save df too for output
       df_sepfin,       #
       df_tandem_runs,  #
       df_sepvid,       #
       df_switches2,    # for # switch comparison
       df_sum,          # for speed comparison
       file = "data_fmt/df_processed.rda")
}
# --------------------------------------------------------------- #

# --------------------------------------------------------------- #
# Output (Data Visualization / Stats)
# --------------------------------------------------------------- #
{
  load("data_fmt/df_processed.rda")
  
  # Description of interaction pattern
  {
    df$state <- ifelse(df$m_lead, "male led tandem",
                       ifelse( df$f_lead, "female led tandem",
                               ifelse(df$close, "close", "separated")
                       ))
    df_state <- as.data.frame(table(df$video, df$state))
    colnames(df_state) <- c("video", "state", "count")
    df_state <- df_state %>%
      group_by(video) %>%
      mutate(prop = count / sum(count))
    df_state$state <- factor(df_state$state, levels = c("separated", "female led tandem",
                                                        "male led tandem", "close"))
    ggplot(df_state, aes(x = forcats::fct_rev(video), y = prop, fill = state)) +
      geom_bar(stat = "identity") +
      scale_y_continuous(breaks = c(0, .5, 1), labels = c(0, .5, 1)) +
      scale_fill_viridis_d(option = "A", begin = 0.2, direction = -1) +
      labs(x = "Video", y = "Proportion", fill = "State") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "top",
            aspect.ratio = .5)
    
    ggsave("output/state_description.pdf", width = 6, height = 6)
  }
  
  
  # Comparison of non-tandem Events Duration by sex and dish size Graph
  {
    # plot
    {
      fit_by_lead <- survfit(Surv(lengths, censor) ~ leader, data = df_sepfin)
      fit_size <- survfit(Surv(lengths, censor) ~ size, data = df_sepvid)
      
      ggsurv <- ggsurvplot(
        fit      = fit_by_lead,
        facet.by = "size",
        xlim     = c(0, 1800),
        conf.int = TRUE,
        xlab     = "Tandem interruption duration (min)",
        ylab     = "Probability",
        censor   = TRUE,
        palette = viridis(2, direction = -1, end =.5, option = "D")
      )
      ggsurv +
        scale_x_continuous(breaks = seq(0, 1800, 600), labels = c(0, 10, 20, 30)) +
        scale_y_continuous(breaks = seq(0, 1, 1)) +
        theme_classic() +
        theme(
          legend.position = c(0.7, 0.8),
          legend.title = element_blank(),
          legend.text = element_text(face = "italic", size = 12),
          text = element_text(size = 12),
          aspect.ratio = .75
        ) +
        guides(fill = "none")
      
      ggsave("output/plot_nontandem.pdf", height = 4, width = 6)
      
    }
    
    # stat
    {
      fit_sex <- coxme(Surv(lengths, censor) ~ leader * size + (1 | video), data = df_sepfin)
      Anova(fit_sex)
    }
  }
  
  # Tandem Run Duration
  {
    # plot
    {
      df_tandem_runs$lead <- ifelse(df_tandem_runs$m_lead, "male",
                                    ifelse(df_tandem_runs$f_lead, "female", NA))
      df_tandem_runs$lead <- factor(df_tandem_runs$lead, levels = c("male", "female"))
      fit_by_lead <- survfit(Surv(duration, cens) ~ lead, data = df_tandem_runs)
      
      df_tandem_runs$size <- factor(df_tandem_runs$size, levels = c("150", "90"))
      fit_by_size <- survfit(Surv(duration, cens) ~ size, data = df_tandem_runs)
      
      ggsurv <- ggsurvplot(
        fit      = survfit(Surv(duration, cens) ~ lead, data = df_tandem_runs),
        facet.by = "size",
        xlim     = c(0, 1000),
        conf.int = TRUE,
        xlab     = "Tandem run duration (min)",
        ylab     = "Probability",
        censor   = TRUE,
        palette = viridis(2, direction = -1, end =.5, option = "D")
      )
      ggsurv +
        facet_wrap( ~ size,
                    ncol = 1,
                    as.table = F,
                    strip.position = "top") +
        scale_x_continuous(breaks = c(0, 600), labels = c(0, 10)) +
        scale_y_continuous(breaks = seq(0, 1, 1)) +
        theme_classic() +
        theme(
          strip.placement = "outside",
          strip.background = element_blank(),
          panel.grid = element_blank(),
          legend.position = c(0.9, 0.9),
          legend.title = element_blank(),
          text = element_text(size = 12),
          aspect.ratio = .66
        ) +
        guides(fill = "none")
      
      ggsave("output/plot_tandem.pdf", height = 6, width = 4)
      
      #Probability Density of Tandem Run Durations
      {
        ggplot(df_tandem_runs, aes(x = duration, fill = lead)) +
          geom_density(alpha = 0.5) +
          labs(
            title = "Probability Density of Tandem Run Durations",
            x = "Duration (seconds)",
            y = "Density",
            fill = "Leading Sex"
          ) +
          theme_minimal()
      }
      
      #Tandem Runs and seperation Runs Duration Graph
      {
        breaks_seq <- seq(1.0, 2000, by = 0.1)
        truehist(df_final[!df_final$tandem, ]$duration,
                 breaks = breaks_seq,
                 xlim = c(0, 30))
        truehist(df_final[df_final$tandem, ]$duration,
                 breaks = seq(0, 1000, 0.1),
                 xlim = c(0, 30))
      }
      
    }
    
    # stat
    {
      cox_model <- coxme(Surv(duration, cens) ~ lead * size + (1 | video),
                         data = df_tandem_runs)
      Anova(cox_model)
    }
  }
  
  
  # Switches based on dish size
  {
    ggplot(df_switches2, aes(x = dish_size, y = switches)) +
      geom_dotplot(binaxis = "y", stackdir = "center") +
      labs(
        x = "Dish Size",
        y = "Number of Switches",
        title = "Number of Leading-Role Switches by Dish Size"
      ) +
      theme_minimal()
    
    r <- glm(switches ~ dish_size, family = "poisson", data = df_switches2)
    Anova(r)
    
    write.csv(df_switches2, file = "output/switch.csv", row.names = F)
    
  }
  
  # Box plot of speed in tandem based on dish size
  {
    df_sum$size <- factor(df_sum$size, levels = c("90 mm","150 mm"))
    ggplot(df_sum, aes(x = size, y = leader_speed)) +
      geom_boxplot(fill = c("darkgreen", "orange"), alpha = .25, width = .5) +
      geom_point() + 
      labs(x = "", y = "Speed (mm/s)") +
      theme_classic() +
      scale_y_continuous(limits = c(0,30))
    ggsave("output/speed_comp.pdf", width = 2.5, height = 3)
    
    
    t.test(leader_speed ~ size, data = df_sum)
    
    #Average Speeds and durations
    {
      cat("90mm male average speed:", mavg_speed90, "\n")
      cat("150mm male average speed:", mavg_speed150, "\n")
      cat("\n")
      cat("90mm female average speed:", favg_speed90, "\n")
      cat("150mm female average speed:", favg_speed150, "\n")
      cat("\n")
      cat("90mm average tandem run Durration:", average90RunDur, "\n")
      cat("150mm average tandem run Durration:",
          average150RunDur,
          "\n")
    }
  }
  
} 
# --------------------------------------------------------------- #
