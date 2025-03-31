#### PACKAGES ########
{
  if (!requireNamespace("arrow", quietly = TRUE)) install.packages("arrow")
  library(stringr)
  library(data.table)
  library(arrow)
} 

### Tandem Smoothing Function
{
tandem.smoothing <- function(vec, min.frame){ 
  if (length(vec) == 0) {
    return(vec)
  }
  if(sum(vec)>0){
    timing <- which(vec)[c(T, diff(which(vec))>1)]
    end    <- which(vec)[c(diff(which(vec))>1,T)]
    for(fi in 1:length(timing)){
      if(length( vec[timing[fi]:end[fi]]) < min.frame ){
        vec[timing[fi]:end[fi]] <- F
      }
    }
  }
  return(vec)
}
}
### Data Formatting / Storing into rda, downsampling, and set to mm
{
  ### Upload data and define columns  
  setwd("C:/Users/Mizumoto-lab/Desktop/hodo-tandem/analysis")
  dataset<- ("data_fmt/data_fmt_df.feather")
  
  # Initialize data frame and load data
  # df_all is the orgional data, while df is the origional data along with all calculation inferences from the data
  df_all <- data.frame()
  
  print(paste("working on", dataset))
    
  df <- arrow::read_feather(dataset)
  colnames(df)[1] <- "frame"

  ### Get bodysize data 
  df_body <- fread("data_fmt/data_fmt_bodysize.csv")
  for (i in 2:nrow(df_body)) {  # Start from row 2 to skip the first row (it is the title)
    if (grepl("150", df_body[i, 1])) {
      # If "150" is in the video name, scale by 150
      df_body[i, 2:3] <- df_body[i, 2:3] / 2000 * 150
    } else {
      # If "150" is not in the video name, scale by 90
      df_body[i, 2:3] <- df_body[i, 2:3] / 2000 * 90
    }
  }
  df_body_scaled <- df_body
  ###
  
  ##      WHY DONT WE SCALE ALL THE ROWS HERE???
  
  
  ####
  # 1. scale it into mm and second and downsample it in 5 FPS
  ## Conversion from pixels to mm##
  scaling_factors <- ifelse(grepl("150", df[, 2]), 150, 90)
  
  # Apply scaling to columns 3 to 6 based on the scaling factor
  df[, 3:6] <- df[, 3:6] / 2000 * scaling_factors
  
  ## Down-sampling by specifying data at every 5 FPS##
  colnames(df)[1] <- 'time'
  df <- df[df$time %% 6 == 0, ]
  df[,1] = df[,1] / 30
      

    
    df_all <- rbind(df_all,data.frame(df))

  

  ### Save dataframes and use to obtain speed and stability
  save(df_all, file="data_fmt/df_all.rda")
  save(df_body_scaled, file="data_fmt/df_body_scaled.rda")
  
}

### Data Analysis
{
  load("data_fmt/df_all.rda")
  #### SPEED ####
  df$f_spe <- c(NA, sqrt( diff(df$fx)^2 + diff(df$fy)^2 ))
  df$m_spe <- c(NA, sqrt( diff(df$mx)^2 + diff(df$my)^2 ))
  
  #### TANDEM ####
  
  ###     NEED MORE STRICT DEF FOR TANDEM AND MAYBE APPLY TANDEM SMOOTHING AT THE START FOR THE ENTIRE DATASET
  
  
  body_size_f <- mean(df_body_scaled$female)
  body_size_m <- mean(df_body_scaled$male)
  tandem_threshold = (body_size_f + body_size_m) * 0.6
  df$partner_dis = sqrt((df$fy - df$my)^2 + (df$fx - df$mx)^2)
  df$tandem = df$partner_dis < tandem_threshold
  ## Applying Tandem Smoothing to df
  tandem <- df$tandem
  tandem = !tandem.smoothing(!tandem, 10)
  tandem <- tandem.smoothing(tandem, 10)
  df$tandem <- tandem
  ###Compare the duration of male-led vs. female-led tandem runs
  
  ### I DONT THINK I ANSWERED THE CORRECT QUESTION
  
  
  df_lead <- data.frame(video = character(), m_lead = numeric(), f_lead = numeric())
  videos = unique(df$video)
  for(v in videos){
    video_rows <- which(df$video == v)
    m_lead = 0
    f_lead = 0
    
    for(j in video_rows){
      if(df$tandem[j] == TRUE){
        mTof = sqrt((df$mx_headtip[j] - df$fx_abdomen[j])^2 + (df$my_headtip[j] - df$fy_abdomen[j])^2)
        fTom = sqrt((df$fx_headtip[j] - df$mx_abdomen[j])^2 + (df$fy_headtip[j] - df$my_abdomen[j])^2)
        if(mTof > fTom){
          m_lead = m_lead + 1
        } else {
          f_lead = f_lead + 1
        }
      }
    }
    df_lead_temp <- data.frame(
      video = v,
      m_lead = m_lead / 5,  
      f_lead = f_lead / 5   
    )
    df_lead <- rbind(df_lead, df_lead_temp)
  }
  
  ### how many different times termites tandem run in each video
  df_list <- list()
  for(i_v in 1:length(videos)){ 
    df_temp <- subset(df, video == videos[i_v])
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
    sep_duration[sep.sta == 1] <- NA
    
    # Create data frames for tandem and separation
    df_tandem <- data.frame(
      video = videos[i_v],
      tandem = TRUE,
      duration = tan_duration
    )
    
    df_separation <- data.frame(
      video = videos[i_v],
      tandem = FALSE,
      duration = sep_duration
    )
    df_video <- rbind(df_tandem, df_separation)
    df_list[[i_v]] <- df_video
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
    
    # Create a temporary data frame for the current video
    df_lead_temp <- data.frame(
      video = v,
      switches = switches
    )
    
    # Append the temporary result to df_switches (use a list to collect data for efficiency)
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
}
save(df_lead, file = "data_fmt/df_lead.rda")
dim(df)
load("data_fmt/df_lead.rda")