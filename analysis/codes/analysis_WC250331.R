#### PACKAGES + PARAMETERS ########
{
  update.packages(ask = FALSE, checkBuilt = TRUE)
  library(stringr)
  library(data.table)
  library(arrow)
  library(dplyr)
  library(MASS)
  library(ggplot2)
  library(patchwork)
  library(knitr)
  library(survival)
  library(survminer)
  library(zoo)
  library(cowplot)
  library(coxme)
  library(tidyr)

  # parameters
  tandemAngle = 60 * (pi / 180)
  tandemsmooth =20
  leadsmooth=15
  tandemSpeed = 1.213
} 
### Tandem Smoothing Function
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
  for (i in 1:nrow(df_body)) {
    if (grepl("150", df_body[i, 1])) {
      df_body[i, 2:3] <- df_body[i, 2:3] / 2000 * 150
    } else {
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
### Data Set up
{
  load("data_fmt/df_all.rda")
  
  #### SPEED ####
  df$f_spe <- c(0, sqrt( diff(df$fx)^2 + diff(df$fy)^2 ))*5
  df$m_spe <- c(0, sqrt( diff(df$mx)^2 + diff(df$my)^2 ))*5
  df$f_spe[df$time == 0] <- 0
  df$m_spe[df$time == 0] <- 0
  
  #### TANDEM + CLOSE####
  body_size_f <- mean(df_body_scaled$female)
  body_size_m <- mean(df_body_scaled$male)
  tandem_threshold = (body_size_f + body_size_m) * 0.6
  df$partner_dis = sqrt((df$fy - df$my)^2 + (df$fx - df$mx)^2)
  
  df$heading_f <- atan2(df$fy_headtip - df$fy_abdomen, df$fx_headtip - df$fx_abdomen)
  df$heading_m <- atan2(df$my_headtip - df$my_abdomen, df$mx_headtip - df$mx_abdomen)
  heading_diff <- abs((df$heading_f - df$heading_m + pi) %% (2 * pi) - pi)
  
  df$close <- (df$partner_dis < tandem_threshold)
  df$tandem <- (df$f_spe > tandemSpeed)&(df$m_spe > tandemSpeed)&(df$partner_dis < tandem_threshold) &(heading_diff < tandemAngle)
  df$tandem[is.na(df$tandem)] <- FALSE
  
  ### what is the average duration of male led tandem runs vs what is the average duration of female led tandem runs
  tandem_df <- data.frame(tandem = logical(0))
  close_df <- data.frame(close = logical(0))
  lead_df <- data.frame(m_lead = logical(0),f_lead = logical(0))
  df_list <- list()
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
      
    ##Male and female lead for each frame/smoothing leaders
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
    #smoothing leaders
    repeat {
      old_rle <- rle(df_temp$m_lead)  # Store previous RLE state
      df_temp$m_lead <- tandem.smoothing(df_temp$m_lead, leadsmooth)
      new_rle <- rle(df_temp$m_lead)  # Check new state
      if (identical(old_rle$lengths, new_rle$lengths)) break  # Stop when no changes occur
    }
    repeat {
      old_rle <- rle(df_temp$f_lead)  # Store previous RLE state
      df_temp$f_lead <- tandem.smoothing(df_temp$f_lead, leadsmooth)
      new_rle <- rle(df_temp$f_lead)  # Check new state
      if (identical(old_rle$lengths, new_rle$lengths)) break  # Stop when no changes occur
    }
    
    lead_df <- rbind(lead_df, data.frame(m_lead = df_temp$m_lead,f_lead = df_temp$f_lead))
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
  df$tandem=tandem_df$tandem
  df$close=close_df$close
  df$close <- ifelse(df$tandem, FALSE, df$close)
  df$m_lead=lead_df$m_lead
  df$f_lead=lead_df$f_lead
  df$lead <- ifelse(df$m_lead, "male",
                    ifelse(df$f_lead, "female", NA))
  
  ###how frequently pairs switch leading roles during tandem running.
  df$switch=FALSE
  df_switches <- data.frame(video = character(0), switches = integer(0))
  frame=0
  for(v in videos){
    video_rows <- which(df$video == v)
    m_lead = 0
    f_lead = 0
    switches = 0
    previous_leader = ""
    for(j in video_rows){
      frame=frame+1
      if(df$tandem[j] == TRUE){
        current_leader = df$lead[j]
        # Increment switch count if the leader changes
        if (previous_leader != "" && current_leader != previous_leader) {
          switches = switches + 1
          df$switch[frame]=TRUE
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
  # finding the last led for seperation events
  df$lastled <- NA
  video_list <- split(df, df$video)
  processed_list <- lapply(video_list, function(sub_df) {
  sub_df <- sub_df[order(sub_df$time), ]
  last_leader <- NA
  for(i in 1:nrow(sub_df)) {
    if(sub_df$time[i] == 0) {
        sub_df$lastled[i] <- NA
      } else {
        if(is.na(sub_df$lead[i])) {
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
  
  
  #Tandem runs Seperations
  tandem_runs <- subset(df_final, tandem == TRUE)
  tandem_runs$event <- 1
  tandem_runs$lead <- ifelse(tandem_runs$m_lead, "male",
                             ifelse(tandem_runs$f_lead, "female", NA))
  tandem_runs$size <- ifelse(grepl("150", tandem_runs$video), "150", 
                             ifelse(grepl("90", tandem_runs$video), "90", NA))
  
  #seperations duration by sex set up
  sep_runs <- df[c("video", "lastled")]
  sep_runs$malelast=FALSE
  sep_runs$femalelast=FALSE
  sep_runs <- sep_runs %>%
    mutate(lastled = replace_na(lastled, "nothing"))
  sep_runs <- sep_runs %>%
    mutate(malelast = as.character(lastled) == "male")
  sep_runs <- sep_runs %>%
    mutate(femalelast = lastled == "female")
  sep_runs$malelast[is.na(sep_runs$malelast)] <- FALSE
  sep_runs$malelast[is.na(sep_runs$femalelast)] <- FALSE
  
  df_malesep=rle(sep_runs$malelast)
  df_malesep <- data.frame(
    mlast = df_malesep$values,
    lengths = df_malesep$lengths
  )
  df_malesep <- df_malesep %>%
    filter(mlast != FALSE)
  
  df_fmalesep=rle(sep_runs$femalelast)
  df_fmalesep <- data.frame(
    mlast = df_fmalesep$values,
    lengths = df_fmalesep$lengths
  )
  df_fmalesep <- df_fmalesep %>%
    filter(mlast != FALSE)
  df_fmalesep$leader="Female"
  df_malesep$leader="Male"
  df_sepfin <- rbind(df_malesep, df_fmalesep)
  df_sepfin$event=1
  
  rle_m <- rle(sep_runs$malelast)
  starts_m <- cumsum(c(1, head(rle_m$lengths, -1)))
  df_malesep$video <- df$video[starts_m[rle_m$values]]
  
  # For Female-led separations
  rle_f <- rle(sep_runs$femalelast)
  starts_f <- cumsum(c(1, head(rle_f$lengths, -1)))
  df_fmalesep$video <- df$video[starts_f[rle_f$values]]
  df_malesep$size <- ifelse(grepl("150", df_malesep$video), "150mm", "90mm")
  df_fmalesep$size <- ifelse(grepl("150", df_fmalesep$video), "150mm", "90mm")
  df_sepvid <- rbind(df_malesep, df_fmalesep)
  df_sepvid$event=1
  
  #Makes sure that the indexing is correct
  cat("Check counts:\n")
  cat("Number of male-led runs:", nrow(df_malesep), "\n")
  cat("Number of starts for male-led runs:", length(starts_m[rle_m$values]), "\n")
  
  #Box plot of speed in tandem based on dish size Set Up
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
      m_speed         = mean(df_tem$m_spe, na.rm = TRUE)
    )
    df_sum <- rbind(df_sum, df_sum_temp)
  }
  
  #Box plot of switches based on dish size Set Up
  df_switches2 <- df_switches %>%
    mutate(
      dish_size = ifelse(grepl("150", video), "150 mm", "90 mm")
    )
  
}
#### Data Visualization
{
  #Separation Events Duration by sex Graph
  {
  surv_obj <- Surv(time = df_sepfin$lengths, event = df_sepfin$event)
  fit <- survfit(surv_obj ~ leader, data = df_sepfin)
  fit_by_lead <- survfit(Surv(lengths, event) ~ leader, data = df_sepfin)
  sep_plot=plot(fit_by_lead, col = c("Red", "Blue"), lwd = 2, 
       xlab = "Separation Event Duration (seconds)", 
       ylab = "Survival Probability",
       main = "Survivorship Curve by Sex")
  legend("topright", legend = c("female last led", "male last led"), col = c("Red", "Blue"), lwd = 2)
  }
  
  #Tandem Run Duration  by Sex
  {
  tandem_runs$lead <- ifelse(tandem_runs$m_lead, "male",
                               ifelse(tandem_runs$f_lead, "female", NA))
  tandem_runs$lead <- factor(tandem_runs$lead, levels = c("male", "female"))
  fit_by_lead <- survfit(Surv(duration, event) ~ lead, data = tandem_runs)
  plot(fit_by_lead, col = c("blue", "red"), lwd = 2, 
       xlab = "Tandem Run Duration (seconds)", 
       ylab = "Survival Probability",
       main = "Survivorship Curve by Leading Sex")
  legend("topright", legend = c("Male-led", "Female-led"), col = c("blue", "red"), lwd = 2)
  }
  
  #Separation Run Duration by size of dish
  {
    fit_size <- survfit(Surv(lengths, event) ~ size, data = df_sepvid)
    
    plot(fit_size, col = c("darkgreen", "orange"), lwd = 2,
         xlab = "Separation Event Duration (seconds)",
         ylab = "Survival Probability",
         main = "Separation Duration by Dish Size")
    
    legend("topright", legend = levels(factor(df_sepvid$size)), 
           col = c("darkgreen", "orange"), lwd = 2)
  }
  
  #Tandem Run Duration by size of dish
  {
    tandem_runs$lead <- factor(tandem_runs$size, levels = c("150", "90"))
    fit_by_lead <- survfit(Surv(duration, event) ~ lead, data = tandem_runs)
    plot(fit_by_lead, col = c("blue", "red"), lwd = 2, 
         xlab = "Tandem Run Duration (seconds)", 
         ylab = "Survival Probability",
         main = "Survivorship Curve by Dish Size")
    legend("topright", legend = c("150mm", "90mm"), col = c("blue", "red"), lwd = 2)
  }
  
  #Box plot of switches based on dish size
  {
    ggplot(df_switches2, aes(x = dish_size, y = switches)) +
      geom_boxplot() +
      labs(
        x = "Dish Size",
        y = "Number of Switches",
        title = "Number of Leading-Role Switches by Dish Size"
      ) +
      theme_minimal()
  }
  
  #Box plot of speed in tandem based on dish size and Sex
  {
    
    # 1. Pivot from wide to long: gather female and male speed into one 'speed' column
    df_tandem_long <- df_sum %>%
      dplyr:: select(video, size, f_speed, m_speed) %>%
      pivot_longer(
        cols = c("f_speed", "m_speed"),
        names_to = "sex",
        values_to = "speed"
      ) %>%
      mutate(
        sex = ifelse(sex == "f_speed", "Female", "Male"),
        group = paste(size, sex)
      )
    # 2. Plot with group on x-axis and speed on y-axis
    ggplot(df_tandem_long, aes(x = group, y = speed)) +
      geom_boxplot() +
      labs(
        x = "Dish Size and Sex",
        y = "Speed (mm/s)",
        title = "Tandem Speeds: 90 Female, 90 Male, 150 Female, 150 Male"
      ) +
      theme_minimal()
  }

} 
#Check for interactions/ Other Plots
{
  #states
  {
    df$state <- ifelse(df$m_lead, 
                       "male led tandem",
                       ifelse(df$f_lead, 
                              "female led tandem",
                              ifelse(df$close, 
                                     "close", 
                                     "separated")))
    freqs <- table(df$state)
    percent <- prop.table(freqs) * 100
    cbind(Frames = freqs, Percent = round(percent, 2))
    
  }
  #P value
  {
cox_model <- coxme(Surv(duration, event) ~ lead + (1 | video), data = tandem_runs)
summary(cox_model)
  }
  # Mixed effects Cox model: survival ~ lead (i.e., sex) + random intercept for video
  {
fit_sex <- coxme(Surv(duration, event) ~ lead + (1 | video), data = tandem_runs)
summary(fit_sex)
anova(fit_sex)
  }
  # Plot tandem and separation durations
  {
  # Define lead type based on logical conditions
  tandem_runs$lead <- factor(ifelse(tandem_runs$m_lead, "Male-led", "Female-led"), 
                             levels = c("Male-led", "Female-led"))
  
  
  
  # Fit survival model for tandem runs
  fit_by_lead <- survfit(Surv(duration, event) ~ lead, data = tandem_runs)
  
  # Create survival plot for tandem runs with CIs and correct legend
  tan_plot <- ggsurvplot(fit_by_lead, 
                         data = tandem_runs,
                         xlab = "Tandem Run Duration (seconds)", 
                         ylab = "Probability of Tandem",
                         conf.int = TRUE,      # Add confidence intervals
                         palette = c("blue", "red"),
                         legend.title = "Sex of leader", 
                         legend.labs = c("Male", "Female"),  # Legend labels
                         ggtheme = theme_classic())$plot +
    theme(plot.title = element_text(hjust = -0.1, size = 16, face = "bold"))
}
  #Checking between 'size' and 'lead'
  {
  #Checking between 'size' and 'lead'
  tandem_runs$lead_size <- factor(interaction(tandem_runs$lead, tandem_runs$size), 
                                  levels = c("Male-led.90", "Male-led.150", "Female-led.90", "Female-led.150"),
                                  labels = c("Male 90", "Male 150", "Female 90", "Female 150"))
  
  # Fit survival model by both 'size' and 'lead'
  fit_by_lead_size <- survfit(Surv(duration, event) ~ lead_size, data = tandem_runs)
  
  # Create survival plot with 'size' as a grouping factor
  tan_plot <- ggsurvplot(fit_by_lead_size, 
                         data = tandem_runs,
                         xlab = "Tandem Run Duration (seconds)", 
                         ylab = "Probability of Tandem",
                         conf.int = FALSE,      # Add confidence intervals
                         palette = c("blue", "red", "green", "purple"),
                         legend.title = "Leader and Size", 
                         legend.labs = c("Male 90", "Male 150", "Female 90", "Female 150"),  # Legend labels
                         ggtheme = theme_classic())$plot +
    theme(plot.title = element_text(hjust = -0.1, size = 16, face = "bold"))
  
  # Print the plot
  print(tan_plot)
}
  #Probability Density of Tandem Run Durations
  {
    ggplot(tandem_runs, aes(x = duration, fill = lead)) +
      geom_density(alpha = 0.5) +
      labs(title = "Probability Density of Tandem Run Durations",
           x = "Duration (seconds)", y = "Density",
           fill = "Leading Sex") +
      theme_minimal()
  }
  #Tandem Runs and seperation Runs Duration Graph
  {
    breaks_seq <- seq(1.0, 2000, by = 0.1)
    truehist(df_final[!df_final$tandem, ]$duration, breaks = breaks_seq, xlim = c(0, 30))
    truehist(df_final[df_final$tandem,]$duration, breaks = seq(0,1000,0.1), xlim=c(0,30))
  }
  #Average Speeds and durations
  {
    cat("90mm male average speed:", mavg_speed90, "\n")
    cat("150mm male average speed:", mavg_speed150, "\n")
    cat("\n")
    cat("90mm female average speed:", favg_speed90, "\n")
    cat("150mm female average speed:", favg_speed150, "\n")
    cat("\n")
    cat("90mm average tandem run Durration:", average90RunDur, "\n")
    cat("150mm average tandem run Durration:", average150RunDur, "\n")
  }

}
#Extra Fun Stuff
{
  
}
#Etc
{
save(df_lead, file = "data_fmt/df_lead.rda")
dim(df)
load("data_fmt/df_lead.rda")
}