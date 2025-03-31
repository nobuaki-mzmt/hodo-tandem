#### PACKAGES ########
{
  if (!requireNamespace("arrow", quietly = TRUE)) install.packages("arrow")
  library(stringr)
  library(data.table)
  library(arrow)
  library(dplyr)
  
  # parameters
  tandemAngle = 60 * (pi / 180)
} 

### Tandem Smoothing Function
{
  tandem.smoothing <- function(vec, threshold) {
    # Get run-length encoding of the vector
    r <- rle(vec)
    # Precompute start and end indices for each run
    cs <- cumsum(r$lengths)
    starts <- c(1, head(cs, -1) + 1)
    
    # Loop over each run
    for(i in seq_along(r$lengths)) {
      if(r$lengths[i] < threshold) {
        # For short runs, flip the value in that run.
        vec[starts[i]:cs[i]] <- !r$values[i]
      }
    }
    return(vec)
  }
}

### Data Formatting / Storing into rda, downsampling, and set to mm
{
  ### Upload data and define columns  
  setwd("C:/Users/Mizumoto-lab/Desktop/hodo-tandem/analysis")
  dataset<- ("data_fmt/data_raw_df.feather")
  
  # Initialize data frame and load data
  # df_all is the orgional data, while df is the origional data along with all calculation inferences from the data
  df_all <- data.frame()
  
  print(paste("working on", dataset))
    
  df <- arrow::read_feather(dataset)
  colnames(df)[1] <- "frame"

  ### Get bodysize data 
  df_body <- fread("data_fmt/data_raw_bodysize.csv")
  for (i in 1:nrow(df_body)) {  # Start from row 2 to skip the first row (it is the title)
    if (grepl("150", df_body[i, 1])) {
      # If "150" is in the video name, scale by 150
      df_body[i, 2:3] <- df_body[i, 2:3] / 2000 * 150
    } else {
      # If "150" is not in the video name, scale by 90
      df_body[i, 2:3] <- df_body[i, 2:3] / 2000 * 90
    }
  }
  df_body_scaled <- df_body

  ####
  # 1. scale it into mm and second and downsample it in 5 FPS
  ## Conversion from pixels to mm##
  df$scaling_factor <- ifelse(grepl("150", df$video), 150, 90)
  df$fx <- df$fx / 2000 * df$scaling_factor
  df$fy <- df$fy / 2000 * df$scaling_factor
  df$mx <- df$mx / 2000 * df$scaling_factor
  df$my <- df$my / 2000 * df$scaling_factor
  
  ##Scaling body parts
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
    
    df_all <- rbind(df_all,data.frame(df))

  

  ### Save dataframes and use to obtain speed and stability
  save(df_all, file="data_fmt/df_all.rda")
  save(df_body_scaled, file="data_fmt/df_body_scaled.rda")
  
}

### Data Analysis
{
  load("data_fmt/df_all.rda")
  
  #### SPEED ####
  df$f_spe <- c(0, sqrt( diff(df$fx)^2 + diff(df$fy)^2 ))*5
  df$m_spe <- c(0, sqrt( diff(df$mx)^2 + diff(df$my)^2 ))*5
  df$f_spe[df$time == 0] <- 0
  df$m_spe[df$time == 0] <- 0
  #### TANDEM ####
  
  body_size_f <- mean(df_body_scaled$female)
  body_size_m <- mean(df_body_scaled$male)
  tandem_threshold = (body_size_f + body_size_m) * 0.6
  df$partner_dis = sqrt((df$fy - df$my)^2 + (df$fx - df$mx)^2)
  
  df$heading_f <- atan2(df$fy_headtip - df$fy_abdomen, df$fx_headtip - df$fx_abdomen)
  df$heading_m <- atan2(df$my_headtip - df$my_abdomen, df$mx_headtip - df$mx_abdomen)
  heading_diff <- abs((df$heading_f - df$heading_m + pi) %% (2 * pi) - pi)
  
  summary(df$f_spe)
  
  tandemSpeed = 1.213
  df$tandem <- (df$f_spe > tandemSpeed)&(df$m_spe > tandemSpeed)&(df$partner_dis < tandem_threshold) &(heading_diff < tandemAngle)
  df$tandem[is.na(df$tandem)] <- FALSE
  
  

  
  
  
  
  
  ### I DONT THINK I ANSWERED THE CORRECT QUESTION (what is the average duration of male led tandem runs vs what is the average duration of female led tandem runs)
  df_list <- list()
  for(i_v in 1:length(videos)){ 
    df_temp <- subset(df, video == videos[i_v])
    ## Applying Tandem Smoothing to df
    {
      tandem <- df_temp$tandem
      tandem = !tandem.smoothing(!tandem, 15)
      tandem <- tandem.smoothing(tandem, 15)
      df_temp$tandem <- tandem
  }
    ##Male and female lead for each frame
    {
    df_temp$m_lead <- FALSE
    df_temp$f_lead <- FALSE
      for (f in seq_len(nrow(df_temp))) {  
      if (df_temp$tandem[f] == TRUE) {
        mTof <- sqrt((df_temp$mx_headtip[f] - df_temp$fx_abdomen[f])^2 + (df_temp$my_headtip[f] - df_temp$fy_abdomen[f])^2)
        fTom <- sqrt((df_temp$fx_headtip[f] - df_temp$mx_abdomen[f])^2 + (df_temp$fy_headtip[f] - df_temp$my_abdomen[f])^2)
        
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
    }
    tandem <- df_temp$tandem
    # Tandem calculations
    tan.end <- which(tandem)[c(diff(which(tandem)) > 1, T)]
    tan.sta <- which(tandem)[c(T, diff(which(tandem)) > 1)]
    tan_duration <- (tan.end - tan.sta) / 5
    
    # Separation calculations
    separation <- !tandem
    sep.end <- which(separation)[c(diff(which(separation)) > 1, T)]
    sep.sta <- which(separation)[c(T, diff(which(separation)) > 1)]
    sep_duration <- (sep.end - sep.sta) / 5
    sep_cens <- sep.end != length(separation)
    

    
    # Create data frames for tandem and separation
    df_tandem <- data.frame(
      video = videos[i_v],
      tandem = TRUE,
      duration = tan_duration,
      m_lead = df_temp$m_lead[tan.sta],
      f_lead = df_temp$f_lead[tan.sta]
    )
    
    df_separation <- data.frame(
      video = videos[i_v],
      tandem = FALSE,
      duration = sep_duration,
      m_lead = NA,
      f_lead = NA
    )
    
    df_list[[i_v]] <- rbind(df_tandem, df_separation)
  }
  
  df_final <- do.call(rbind, df_list)
  
  
  
  
  
  ###how frequently pairs switch leading roles during tandem running.
  
  df_switches <- data.frame(video = character(0), switches = integer(0))
  for(v in videos){
    video_rows <- which(df$video == v)
    m_lead = 0
    f_lead = 0
    switches = 0
    previous_leader = ""
    
    for(j in video_rows){
      if(df$tandem[j] == TRUE){
        male_distance = sqrt((df$mx_headtip[j] - df$fx_abdomen[j])^2 + (df$my_headtip[j] - df$fy_abdomen[j])^2)
        female_distance = sqrt((df$fx_headtip[j] - df$mx_abdomen[j])^2 + (df$fy_headtip[j] - df$my_abdomen[j])^2)
        current_leader = ifelse(male_distance > female_distance, "male", "female")
        
        # Increment switch count if the leader changes
        if (previous_leader != "" && current_leader != previous_leader) {
          switches = switches + 1
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
  
  cat("90mm male average speed:", mavg_speed90, "\n")
  cat("150mm male average speed:", mavg_speed150, "\n")
  cat("\n")
  cat("90mm female average speed:", favg_speed90, "\n")
  cat("150mm female average speed:", favg_speed150, "\n")
  cat("\n")
  cat("90mm average tandem run Durration:", average90RunDur, "\n")
  cat("150mm average tandem run Durration:", average150RunDur, "\n")
  #truehist(df_final[df_final$tandem,]$duration, breaks = seq(0,1000,0.1), xlim=c(0,30))
  
}
save(df_lead, file = "data_fmt/df_lead.rda")
dim(df)
load("data_fmt/df_lead.rda")