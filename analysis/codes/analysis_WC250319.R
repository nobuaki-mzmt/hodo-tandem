#### PACKAGES ########
{
  if (!requireNamespace("arrow", quietly = TRUE)) install.packages("arrow")
  library(stringr)
  library(data.table)
  library(arrow)
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
  body_size_f <- mean(df_body$female)
  body_size_m <- mean(df_body$male)
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
}
  
save(df_lead, file = "data_fmt/df_lead.rda")
dim(df)

load("data_fmt/df_lead.rda")

