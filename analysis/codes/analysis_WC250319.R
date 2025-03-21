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
  
  # load data
  df_all <- NULL
  df_body_all <- NULL
  
  print(paste("working on", dataset))
    
  df <- arrow::read_feather(dataset)
  colnames(df)[1] <- "frame"

  ### Get bodysize data 
  d_body <- fread("data_fmt/data_fmt_bodysize.csv")
  for (i in 2:nrow(d_body)) {  # Start from row 2 to skip the first row
    if (grepl("150", d_body[i, 1])) {
      # If "150" is in the video name, scale by 150
      d_body[i, 2:3] <- d_body[i, 2:3] / 2000 * 150
    } else {
      # If "150" is not in the video name, scale by 90
      d_body[i, 2:3] <- d_body[i, 2:3] / 2000 * 90
    }
  }
  df_body_scaled <- d_body
  
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
  body_size_f <- mean(df_body_scaled$female)
  body_size_m <- mean(df_body_scaled$male)
  tandem_threshold = (body_size_f + body_size_m) * 0.6
  
  df$partner_dis = sqrt((df$fy - df$my)^2 + (df$fx - df$mx)^2)
  df$tandem = df$partner_dis < tandem_threshold
  
  ###RLE Calculations to analyze how long and how often tandem behavior occurs
  
  tempdat <- rle(df$tandem)
  tempdat$lengths
  tempdat$values
  tempdat_df <- data.frame(values = tempdat$values, lengths = tempdat$lengths)
  
  ### Sum of the entire Tandem running duration (NEED TO ADJUST TO MAKE DIFFERNT)
  tan_duration <- sum(tempdat$lengths[tempdat$values == TRUE]) / 5
  print(tan_duration)
 
  ###Compare the duration of male-led vs. female-led tandem runs
  
  df_lead = NULL
  j=0
  m_lead=0
  f_lead=0
  mTof=0
  fTom=0
  videos = unique(df$video)
  for(v in videos){
    video_rows <- which(df$video == v)
  for(j in video_rows){
    if(df$tandem[j]==TRUE){
      mTof=sqrt((df$mx_headtip[j]-df$fx_abdomen[j])^2+(df$my_headtip[j]-df$fy_abdomen[j])^2)
      fTom=sqrt((df$fx_headtip[j]-df$mx_abdomen[j])^2+(df$fy_headtip[j]-df$my_abdomen[j])^2)
      if(mTof>fTom){
        m_lead=m_lead+1
      } else{
        f_lead=f_lead+1
      }
    }
    }
    df_lead_temp <- data.frame(
      video = v,
      m_lead = sum(m_lead, na.rm = TRUE)/5,
      f_lead = sum(f_lead, na.rm = TRUE)/5
    )
    
    df_lead <- rbind(df_lead_temp, df_lead)
    m_lead=0
    f_lead=0
  }
  
  ### how many different times termites tandem run in each video
  {
    ### NOT FULLY NEEDED PER SAY
  videos = unique(df$video)
  vidtempdata_list <- list()
  for(v in videos){
    subset_data <- df$tandem[df$video == v]
    tempdat <- rle(subset_data)
    vidtempdata_df <- data.frame(video=v,values = tempdat$values, lengths = tempdat$lengths)
    vidtempdata_list[[v]] <- vidtempdata_df
  }
  video_tandem_rle_df <- do.call(rbind, vidtempdata_list)
  rownames(video_tandem_rle_df) <- 1:nrow(video_tandem_rle_df)
  
  true_counts_list <- list()
  for (v in names(vidtempdata_list)) {
    vid_df <- vidtempdata_list[[v]]
    values <- vid_df$values
    true_count <- sum(values == TRUE) 
    true_counts_list[[v]] <- data.frame(video = v, true_count = true_count)
  }
  
  video_true_counts_df <- do.call(rbind, true_counts_list)
  rownames( video_true_counts_df) <- 1:nrow( video_true_counts_df)
  }
  video_names = unique(df$video)
  df_list <- list()
  df_sum = df_tandem_sum = df_separation_sum<- NULL
  
  # for loop for each video
  for(i_v in 1:length(video_names)){
    df_temp <- subset(df, video == video_names[i_v])
    
    tandem <- df_temp$tandem
    ## obtain tandem events
    {
      # remove short separation
      tandem = !tandem.smoothing(!tandem, 10)
      #separation = !tandem.smoothing(!separation,10)
      # remove short tandem
      tandem <- tandem.smoothing(tandem, 10) 
      
      # detect tandem events
      tan.end <- which(tandem)[c(diff(which(tandem))>1,T)]
      tan.sta <- which(tandem)[c(T, diff(which(tandem))>1)]
      tan_duration <- (tan.end - tan.sta)/5
      
      ####Include separation
      separation <- !tandem
      #detect separation events
      sep.end <- which(separation)[c(diff(which(separation))>1,T)]
      sep.sta <- which(separation)[c(T, diff(which(separation))>1)]
      sep_duration <- (sep.end - sep.sta)/5
      #print(sep_duration)
      sep_cens     = sep.end != length(separation)
      sep_duration[sep.sta == 1] = NA
      
      df_tandem <- data.frame(
        video = video_names[i_v],
        tandem = TRUE,
        duration = tan_duration
      )
      
      df_separation <- data.frame(
        video = video_names[i_v],
        tandem = FALSE,
        duration = sep_duration
      )
      
      # Combine tandem and separation data for this video
      df_video <- rbind(df_tandem, df_separation)
      
      # Store in the list
      df_list[[i_v]] <- df_video
      
    }
    df_final <- do.call(rbind, df_list)
    
  
  }
  
  ###how frequently pairs switch leading roles during tandem running.

  df_switches <- NULL
  videolist = unique(df$video)
  for(v in videolist){
    video_rows <- which(df$video == v)
    m_lead = 0
    f_lead = 0
    switches = 0
    previous_leader = NULL
    for(j in video_rows){
      if(df$tandem[j] == TRUE){
        maTofe = sqrt((df$mx_headtip[j] - df$fx_abdomen[j])^2 + (df$my_headtip[j] - df$fy_abdomen[j])^2)
        feToma = sqrt((df$fx_headtip[j] - df$mx_abdomen[j])^2 + (df$fy_headtip[j] - df$my_abdomen[j])^2)
        current_leader = ifelse(maTofe > feToma, "male", "female")
        if (!is.null(previous_leader) && current_leader != previous_leader) {
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
  average150 <- sum(df_filtered150$duration) / nrow(df_filtered150)
  average90 <- sum(df_filtered90$duration) / nrow(df_filtered90)
  
  ## Compare tandem speed between different dish sizes.
  
  df_filteredspeed150_list <- list()
  df_filteredspeed90_list <- list()
  
  for(i_v in 1:length(video_names)){
    df_temp <- subset(df, video == video_names[i_v])
    
    tandem <- df_temp$tandem

    {
    
        tandem = !tandem.smoothing(!tandem, 10)
        tandem <- tandem.smoothing(tandem, 10)
    }
    
    df_temp$tandem <- tandem
    
    df_filteredspeed150 <- df_temp[df_temp$tandem == TRUE & grepl("150", df_temp$video), ]
    df_filteredspeed90 <- df_temp[df_temp$tandem == TRUE & grepl("90", df_temp$video), ]
    df_filteredspeed150_list[[i_v]] <- df_filteredspeed150
    df_filteredspeed90_list[[i_v]] <- df_filteredspeed90
  }
  
  df_filteredspeed150_combined <- do.call(rbind, df_filteredspeed150_list)
  df_filteredspeed90_combined <- do.call(rbind, df_filteredspeed90_list)
  favg_speed90 <- mean(df_filteredspeed90_combined$f_spe, na.rm = TRUE)
  favg_speed150 <- mean(df_filteredspeed150_combined$f_spe, na.rm = TRUE)
  mavg_speed90 <- mean(df_filteredspeed90_combined$m_spe, na.rm = TRUE)
  mavg_speed150 <- mean(df_filteredspeed150_combined$m_spe, na.rm = TRUE)
  
  
  
}
  
save(df_lead, file = "data_fmt/df_lead.rda")
dim(df)

load("data_fmt/df_lead.rda")

