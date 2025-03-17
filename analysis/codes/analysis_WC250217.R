#### PACKAGES ########
{
library (stringr)
library(data.table)

} 
load("data_fmt/df_all_hs.rda")
### Storing into rda, downsampling, and set to mm
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
  ## NEED TO DIFFEREINATE BETWEEN 90 and 150 for formatting
  d_body[,2:3] = d_body[,2:3] / 2000 * 150
  df_body_all = d_body
  
  # 1. scale it into mm and second and downsample it in 5 FPS
  ## Conversion from pixels to mm##
  df[,3:6] = df[,3:6] / 2000 * 150
  ## Down-sampling by specifying data at every 5 FPS##
  colnames(df)[1] <- 'time'
  df <- df[df$time %% 6 == 0, ]
  df[,1] = df[,1] / 30
      

    
    df_all <- rbind(df_all,
                    data.frame(
                      df
                    ))
  
  df_body <- df_body_all

  

  ### Save dataframes and use to obtain speed and stability
  save(df_all, file="data_fmt/df_all_hs.rda")
  save(df_body, file="data_fmt/df_body_hs.rda")
  
}


{
  
  #### SPEED ####
  df$f_spe <- c(NA, sqrt( diff(df$fx)^2 + diff(df$fy)^2 ))
  df$m_spe <- c(NA, sqrt( diff(df$mx)^2 + diff(df$my)^2 ))
  
  #### TANDEM ####
  body_size_f <- mean(df_body$female)
  body_size_m <- mean(df_body$male)
  tandem_threshold = (body_size_f + body_size_m) * 0.6
  
  df$partner_dis = sqrt((df$fy - df$my)^2 + (df$fx - df$mx)^2)
  df$tandem = df$partner_dis < tandem_threshold
  
  tempdat <- rle(df$tandem)
  tempdat$lengths
  tempdat$values

  tan_duration=0
  i=0
  for(i_v in tempdat$values){
    i=i+1
    if(i_v==T){
      tan_duration=tan_duration+tempdat$lengths[i]
    } else
      next;
  }
  tan_duration <- sum(tempdat$lengths[tempdat$values == TRUE]) / 5
  tempdat_df <- data.frame(values = tempdat$values, 
                           lengths = tempdat$lengths)
  
 
  df_lead = NULL
  j=0
  m_lead=0
  f_lead=0
  mTof=0
  fTom=0
  video = unique(df$video)
  for(v in 1:video){
  for(i in v){
    j=j+1
    if(df$tandem[j]){
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
      video = unique(df$video),
      m_lead = sum(m_lead, na.rm = TRUE)/0.2,
      f_lead = sum(f_lead, na.rm = TRUE)/0.2
    )
    
    df_lead <- rbind(df_lead_temp, df_lead)
    m_lead=0
    f_lead=0
  }
}
  save(df_lead, file = "data_fmt/df_lead.rda")
dim(df)

load("data_fmt/df_lead.rda")
