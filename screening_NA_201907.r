################ Long record stations: >= 20 yrs.
#Step 1: start from the daily snow depth time series, find the peak date for each snow year (Oct 1 - July 31)
##### This function is to find the peak date and value. It will be saved in outdf dataframe.
peakdate_finder <- function(years, snow_depth_value){
  #browser()
  ## --find the largest snow depth across the entire record
  #snow_depth_value[which(snow_depth_value$datum == max(snow_depth_value$datum)), ]
  ## --find the largest value of each year
  outdf <- as.data.frame(matrix(data=NA, ncol=ncol(snow_depth_value), nrow= years))
  colnames(outdf) <- colnames(snow_depth_value)
  year <- as.character(unique(snow_depth_value$year))
  year_extrm <- as.vector(NA)
  for (i in 1:years){
    #browser()
    mintime = strptime(paste(as.character(as.numeric(year[i])-1),"-10-01", sep = ""), format = '%Y-%m-%d')
    maxtime = strptime(paste(year[i],"-07-31", sep = ""), format = '%Y-%m-%d')
    df <- snow_depth_value[which(snow_depth_value$obs_time >= mintime & snow_depth_value$obs_time <= maxtime),]
    # handle the situation that there's no records in the snow year
    if (nrow(df) == 0) {
      print(paste('No data in the snow year', year[i], sep = ' '))
      next
    }
    # Handle the situation when all data are zero in the snow year records. Keep datum as 0; vars_id; and history_id 
    if (all(df$datum == 0)) {
      flag_zero <- unique(df$datum)
      outdf[i,5] <- flag_zero
      outdf[i,7] <- historyid
      print(paste('The records in the snow year', year[i], 'are all 0', sep = ' '))
      next
    } 
    ayear_value <- df$datum
    year_extrm[i] <- max(ayear_value)
    max_val <- df[which(df$datum == year_extrm[i]), ]
    if (nrow(max_val) > 1 | nrow(max_val == 1)) {
      max_val = max_val[nrow(max_val),] #choose the closest date as peak snow date if there's more than one peak snow values
      ################ QC: check whether max_val is valid ###################
      # if the next two days are both !NA then the SD is considered valid. If both exist, then the SD is considered not valid if both of them are 0 and SD is > 30cm
      date_max <- max_val$adj_jdate
      sdmax_n1 <- df[which(df$adj_jdate == date_max + 1),]
      sdmax_n2 <- df[which(df$adj_jdate == date_max + 2),]
      # For cases that the next two days do not have data
      #browser()
      if ((nrow(sdmax_n1) == 0 | nrow(sdmax_n2) == 0) & max_val$datum <= 30) {
        max_val <- max_val
      } else if((nrow(sdmax_n1) == 0 | nrow(sdmax_n2) == 0) & max_val$datum > 30){
        second_max <- sort(ayear_value, TRUE)[2]
        max_val <- df[which(df$datum == second_max),]
        max_val = max_val[nrow(max_val),]
      } else if (!is.na(sdmax_n1$datum & sdmax_n2$datum) & ((sdmax_n1$datum | sdmax_n2$datum) > 0 | max_val$datum <= 30)) {
        max_val <- max_val
      } else if (!is.na(sdmax_n1$datum & sdmax_n2$datum) & ((sdmax_n1$datum & sdmax_n2$datum) == 0 | max_val$datum > 30)){ #stn 3808
        second_max <- sort(ayear_value, TRUE)[2]
        max_val <- df[which(df$datum == second_max),]
        max_val = max_val[nrow(max_val),]
      } else {
        max_val <- NA
      }
    }
    outdf[i,] <- max_val
  }
  outdf$X <- year #column X is the index, not too much use in the dataframe, replace this column with the snow year
  colnames(outdf)[colnames(outdf) == "X"] <- "snow_year"
  return(outdf)
}

leap_year_check <- function(ayear){
  ndays <- rep(366,length(ayear))
  ndays[(ayear %% 4 != 0)] <- 365
  ndays[(ayear %% 100 != 0 & (ayear %% 4 == 0))] <- 366
  ndays[(ayear %% 400 != 0) & (ayear %% 100 == 0 & (ayear %% 4 == 0))] <- 365
  
  return(ndays)
  #if (ayear %% 4 != 0) {return(365) #common year
  #} else if (ayear %% 100 != 0) {return(366) #leap year
  #} else if (ayear %% 400 != 0) {return(365) #common year
  #} else {return(366) #leap year
  #}
}

QC <- function(snow_depth_value) {
  #browser()
  # Fist: remove 990 cm
  error_limit <- 900
  if (any(snow_depth_value$datum > error_limit)) {
    print(paste("station",historyid, "has error"))
    error_line <- snow_depth_value[which(snow_depth_value$datum > error_limit), ]
    snow_depth_value <- snow_depth_value[!rownames(snow_depth_value) %in% rownames(error_line),] #remove the lines with value greatert than 900
    return(snow_depth_value)
  } else {
    return(snow_depth_value)
  }
}


########### Create a dataframe that stores each station's peak snow date and values
# Read all stations in snow depth folder
#path = "/home/yaqiongw/pcic_storage/NRC/data_mscdb/snow_depth"
#path = "/storage/home/yaqiongw/NRC/data_mscdb/snow_depth"
path = "/storage/home/yaqiongw/NRC/data_mscdb/snow_depth"
metadata_20 <- read.csv("/storage/home/yaqiongw/NRC/metadata_raw_long_records_stns.csv", header = TRUE, stringsAsFactors = FALSE)
l <- nrow(metadata_20)
peaksnow_df_all <- as.data.frame(matrix(data=NA, nrow = l, ncol = 6))
names(peaksnow_df_all) <- c("station_id", "min_depth","max_depth","raw_length_yrs","length_comp_yrs_zeros","length_comp_yrs_nozeros")

for (i in 1:l) {
  #browser()
  historyid <- metadata_20$history_id[i]
  print(paste("history id is:", historyid, "& i = ", i))
  filename <- paste("/storage/home/yaqiongw/NRC/data_mscdb/snow_depth/",historyid,"_snow_depth_dly.csv", sep = "")
  snow_depth_value <- read.csv(filename, header = TRUE, stringsAsFactors = FALSE)
  # QC check
  snow_depth_value <- QC(snow_depth_value)
  ###### Data characteristics
  min_obs_time <- min(snow_depth_value$obs_time) #min obs time of the station
  max_obs_time <- max(snow_depth_value$obs_time) #max obs time of the station
  years <- length(unique(snow_depth_value$year)) #length of the record
  #### Add adjusted jdate column to snow_depth_value: Find all leap year, -366; find all non-leap year -355-------------------------------
  snow_depth_value$adj_jdate <- snow_depth_value$jdate
  snow_depth_value$adj_jdate[leap_year_check(snow_depth_value$year)==366 & snow_depth_value$jdate > 270 ] <- snow_depth_value$jdate[leap_year_check(snow_depth_value$year)==366 & snow_depth_value$jdate > 270] - 366
  snow_depth_value$adj_jdate[leap_year_check(snow_depth_value$year)==365 & snow_depth_value$jdate > 270 ] <- snow_depth_value$jdate[leap_year_check(snow_depth_value$year)==365 & snow_depth_value$jdate > 270] - 365
  ## Find the peak snow date and value
  #browser()
  peaksnow_df <- peakdate_finder(years, snow_depth_value)
  ## for station like 425, all years are either NA or 0
  if (all(is.na(peaksnow_df$year))) {
    next
  }
  #browser()
  pk_zeros_df <- peaksnow_df[which(peaksnow_df$datum == 0),] #years with annmax = 0
  ############### zeros: 100% comp #################
  if (nrow(pk_zeros_df) == 0) { #means there's no years which has annmax=0
    ############# Non-zeros: Calculate the completeness
    pk_snow_year_df <- peaksnow_df[!is.na(peaksnow_df$year),]
    pk_snow_year_df$completeness <- NA
    num_pk_years <- nrow(pk_snow_year_df)
    #browser()
    for (j in 1:num_pk_years) {
      #browser()
      pk_date_annual <- pk_snow_year_df$adj_jdate[j]
      pk_date_start <- pk_date_annual - 30
      pk_date_end <- pk_date_annual +30
      window_size <- 61
      # Got peaksnow_df (annual peak date and maxima), window of the year, find the daily data in the window
      mintime = strptime(paste(as.character(as.numeric(pk_snow_year_df$year[j])-1),"-10-01", sep = ""), format = '%Y-%m-%d')
      maxtime = strptime(paste(pk_snow_year_df$year[j],"-07-31", sep = ""), format = '%Y-%m-%d')
      asnow_year <- snow_depth_value[which(snow_depth_value$obs_time >= mintime & snow_depth_value$obs_time <= maxtime),]
      comp_records <- asnow_year[which(asnow_year$adj_jdate >= pk_date_start & asnow_year$adj_jdate <= pk_date_end),]
      num_nona <- nrow(comp_records)
      if (num_nona == 0) {
        print('No data in the window')
        next
      }
      na_rate <- (num_nona/window_size) * 100
      pk_snow_year_df$completeness[j] <- na_rate 
    }
    #browser()
  } else { #means there are years that have annmax = 0
    #years with zeros are 100% complete: line 24.
    pk_zeros_df$completeness <- 100
    ############# Non-zeros: Calculate the completeness
    pk_snow_year_df <- peaksnow_df[!is.na(peaksnow_df$year),]
    pk_snow_year_df$completeness <- NA
    num_pk_years <- nrow(pk_snow_year_df)
    #browser()
    for (j in 1:num_pk_years) {
      #browser()
      pk_date_annual <- pk_snow_year_df$adj_jdate[j]
      pk_date_start <- pk_date_annual - 30
      pk_date_end <- pk_date_annual +30
      window_size <- 61
      # Got peaksnow_df (annual peak date and maxima), window of the year, find the daily data in the window
      mintime = strptime(paste(as.character(as.numeric(pk_snow_year_df$year[j])-1),"-10-01", sep = ""), format = '%Y-%m-%d')
      maxtime = strptime(paste(pk_snow_year_df$year[j],"-07-31", sep = ""), format = '%Y-%m-%d')
      asnow_year <- snow_depth_value[which(snow_depth_value$obs_time >= mintime & snow_depth_value$obs_time <= maxtime),]
      comp_records <- asnow_year[which(asnow_year$adj_jdate >= pk_date_start & asnow_year$adj_jdate <= pk_date_end),]
      num_nona <- nrow(comp_records)
      if (num_nona == 0) {
        print('No data in the window')
        next
      }
      na_rate <- (num_nona/window_size) * 100
      pk_snow_year_df$completeness[j] <- na_rate 
    }
  }
  num_comp_years_nozeros <- nrow(pk_snow_year_df[which(pk_snow_year_df$completeness == 100),])
  pk_combine <- rbind(pk_zeros_df, pk_snow_year_df)
  pk_combine <- pk_combine[order(pk_combine$snow_year),]
  pk_combine_100 <- pk_combine[which(pk_combine$completeness == 100),]
  # Output the peak_snow_df for each station
  path_peaksnow <- "/storage/home/yaqiongw/NRC/data_mscdb/annual_peak_snow_60day_window_100comp_zerosin/"
  if (any(pk_combine_100$datum == 0)) {
    write.csv(pk_combine_100, paste(path_peaksnow, historyid,"_",i,"_annual_peak_snow_100comp_zero.csv",sep = ""), row.names = FALSE)
  } else {
    write.csv(pk_combine_100, paste(path_peaksnow, historyid,"_",i,"_annual_peak_snow_100comp.csv", sep = ""), row.names = FALSE)
  }
  # Check the number of usable stations
  peaksnow_df_all$station_id[i] <- historyid
  peaksnow_df_all$min_depth[i] <- min(peaksnow_df$datum, na.rm = TRUE)
  peaksnow_df_all$max_depth[i] <- max(peaksnow_df$datum, na.rm = TRUE)
  peaksnow_df_all$raw_length_yrs[i] <- years
  peaksnow_df_all$length_comp_yrs_zeros[i] <- nrow(pk_combine_100) #number of complete years after filtering with zeros
  peaksnow_df_all$length_comp_yrs_nozeros[i] <- num_comp_years_nozeros
  peaksnow_df_all <- peaksnow_df_all[order(as.numeric(peaksnow_df_all$station_id)),] #order the station id
}
write.csv(peaksnow_df_all, "/storage/home/yaqiongw/NRC/data_mscdb/peaksnow_df_all_60day_window_zerosin.csv", row.names = FALSE)
