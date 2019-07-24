library('crmp')
library('climdex.pcic')
library('geosphere')
library('RPostgreSQL')
library(dplyr)
library(gumbel)
library(VGAM)
#Based off of Durre et al. (2010), similar to crmp QAQC scripts but for snow
#Doesn't use crmp_objs like the other qaqc scripts since they don't exist for dbmsc data
#Now also includes outlier detection for max time series scripts as discussed in Kotte (1984)


#grabs metadata
con <- dbConnect(PostgreSQL(),user=,password=,host="dbmsc.pcic.uvic.ca",dbname="msc")
querystr <- 'SELECT meta_history.station_id, meta_history.history_id, meta_history.station_name, meta_history.lon, meta_history.lat, meta_history.elev, meta_history.freq, meta_station.network_id, meta_station.native_id FROM meta_history NATURAL JOIN meta_station;'

try.rv <- try({
  rs <- dbSendQuery(con,querystr)
  stnmeta <- fetch(rs,n=-1)
}, silent=T)
if (inherits(try.rv, 'try-error')) {
  close_crmp_connections()
  print(try.rv)
}

#grabs snow depth data for target station
historyid <- 395
snowd_ids= '1398'
con <- dbConnect(PostgreSQL(),user="cdmb",password="DNRE2019&",host="dbmsc.pcic.uvic.ca",dbname="msc")
querystr <- paste('SELECT * FROM crmp.obs_raw WHERE history_id = ',historyid, ' AND vars_id = ', snowd_ids,
                  ' order by obs_time ASC',sep='')

try.rv <- try({
  rs <- dbSendQuery(con,querystr)
  targetstndat <- fetch(rs,n=-1)
  #temp <- fetch(rs,n=-1)
}, silent=T)
if (inherits(try.rv, 'try-error')) {
  close_crmp_connections()
  print(try.rv)
}
close_crmp_connections()





#snow depth values in Durre are in mm, whereas dbmsc snow depth values are in cm
world_record_flag <- function(targetstndat){
 
  sdflags <- targetstndat$datum[2:nrow(targetstndat)] < 0 | targetstndat$datum[2:nrow(targetstndat)] > 1146
  sdincreaseflags <- diff(targetstndat$datum) >= 192.5
  
  sdincreaseflags[is.na(sdincreaseflags)] <- FALSE
  sdflags[is.na(sdflags)] <- FALSE
  
  wrflagslist <- (list(sdincreaseflags,sdflags))
  names(wrflagslist) <- c('sdincreaseflags','sdflags')
  return(wrflaglist)
  
}


repetition_flag <- function(targetstndat){
  
  sd <- targetstndat$datum[2:nrow(targetstndat)]
  #rle doesn't like lists, gotta pass it a vector
  sdflags <- rep(FALSE,nrow(targetstndat)-1)
  rsd <- rle(sd)
  reps <- which(!is.na(rsd$values) & rsd$values > 0 & rsd$lengths >= 90)
  if(any(reps)){
    rep.lengths.cumsum <- cumsum(rsd$lengths)
    ends <- rep.lengths.cumsum[reps]
    newindex = ifelse(reps>1, reps-1, 0)
    starts = rep.lengths.cumsum[newindex] + 1
    if (0 %in% newindex) starts = c(1,starts)
    for (i in 1:length(starts)){
      sdflags[starts[i]:ends[i]] <- TRUE
    }
  }
  sdflags[is.na(sdflags)] <- FALSE
  return(sdflags)
  
  
}


gap_flag <- function(targetstndat){
  
  sdflags_monthly <- list()
  for(month in 1:12){
    
    monthdat <- targetstndat$datum[as.numeric(format(targetstndat$obs_time, '%m')) == month]
    
    sortdat <- sort(monthdat)
    gaps <- diff(sortdat)
    mindex <- round(length(sortdat)/2)
    
    sdflags <- rep(FALSE,length(monthdat))
    
    tail1 <- last(which(gaps[1:mindex] >= 35))
    if(length(tail1) > 0 & !is.na(tail1)){
      sdflags[1:tail1] <- TRUE
    }
    tail2 <- (first(which(gaps[mindex:length(gaps)] >= 35)) + mindex)
    if(length(tail2) > 0 & !is.na(tail2)){
      sdflags[tail2:length(monthdat)] <- TRUE
    }
    sdflags <- data.frame(sortdat,sdflags)
    sdflags_monthly[[month]] <- sdflags
  }
  names(sdflags_monthly) <- c(1:12)
  return(sdflags_monthly)
  
}


snow_temp_consistency_flag <- function(targetstndat,historyid,con){
  
  tmin_id= '1393'
  querystr <- paste('SELECT * FROM crmp.obs_raw WHERE history_id = ',historyid, ' AND vars_id = ', tmin_id,
                    ' order by obs_time ASC',sep='')
  
  try.rv <- try({
    rs <- dbSendQuery(con,querystr)
    targetstndat_tmin <- fetch(rs,n=-1)
    #temp <- fetch(rs,n=-1)
  }, silent=T)
  if (inherits(try.rv, 'try-error')) {
    close_crmp_connections()
    print(try.rv)
  }
  close_crmp_connections()
  
  targetstndat_tmin <- targetstndat_tmin[targetstndat_tmin$obs_time %in% targetstndat$obs_time,]
 
  tmin_window <- (cbind(targetstndat_tmin$datum[2:nrow(targetstndat_tmin)-1],targetstndat_tmin$datum[2:nrow(targetstndat_tmin)],targetstndat_tmin$datum[2:nrow(targetstndat_tmin)+1]))
  tmin_window <- apply(tmin_window, 1, FUN=min)
   
  sdflags <- diff(targetstndat$datum[2:nrow(targetstndat)]) > 0 & 
             as.numeric(diff.POSIXt(targetstndat$obs_time[2:nrow(targetstndat)])) == 24 & 
             tmin_window[1:(nrow(targetstndat)-2)] >= 70
  
  return(sdflags)
  
  
}




snowfall_snowdepth_flag <- function(targetstndat,historyid,con){
  
  snf_id= '1396'
  querystr <- paste('SELECT * FROM crmp.obs_raw WHERE history_id = ',historyid, ' AND vars_id = ', snf_id,
                    ' order by obs_time ASC',sep='')
  
  try.rv <- try({
    rs <- dbSendQuery(con,querystr)
    targetstndat_snf <- fetch(rs,n=-1)
    #temp <- fetch(rs,n=-1)
  }, silent=T)
  if (inherits(try.rv, 'try-error')) {
    close_crmp_connections()
    print(try.rv)
  }
  close_crmp_connections()
  
  targetstndat_snf <- targetstndat_snf[targetstndat_snf$obs_time %in% targetstndat$obs_time,]
  
  sdflags <- diff(targetstndat$datum[2:nrow(targetstndat)]) > (targetstndat_snf$datum[2:nrow(targetstndat_snf)-1] + targetstndat_snf$datum[2:nrow(targetstndat_snf)] + 25) &
             diff(targetstndat$datum[2:nrow(targetstndat)]) > (targetstndat_snf$datum[2:nrow(targetstndat_snf)] + targetstndat_snf$datum[2:nrow(targetstndat_snf)+1] + 25)
    
  sdflags <- rep(FALSE,nrow(targetstndat))
  sdflags <- diff(targetstndat$datum[targetstndat$obs_time %in% targetstndat_snf$obs_time]) > (targetstndat_snf$datum[2:nrow(targetstndat_snf)-1] + targetstndat_snf$datum[2:nrow(targetstndat_snf)] + 25) 
  
  which(targetstndat_snf$datum[2:nrow(targetstndat_snf)-1] + targetstndat_snf$datum[2:nrow(targetstndat_snf)] > 25)
  targetstndat$datum[targetstndat$obs_time[2:nrow(targetstndat)] %in% targetstndat_snf$obs_time] - targetstndat$datum[targetstndat$obs_time[2:nrow(targetstndat)-1] %in% targetstndat_snf$obs_time]
  
  
  which(targetstndat$obs_time %in% targetstndat_snf$obs_time)
  
  
  
}



#--------------------------------------------------------------------------------------------------------------------------------
#Outlier detection for max time series

increase_check <- function(stndat){
  
  sdmax <- max(stndat$datum)
  lastsdmax <- tail((stndat$obs_time[stndat$datum == sdmax]),1)
  
  sdnextday <- stndat$datum[as.Date(stndat$obs_time) == (as.Date(lastsdmax)+1)]
  sddayplus2 <- stndat$datum[as.Date(stndat$obs_time) == (as.Date(lastsdmax)+2)]
  stndatsorted <- stndat[order(stndat$datum,decreasing=TRUE),]
  
  #teststn stuff
#  sdmax <- max(teststndat2$datum)
#  lastsdmax <- tail((teststndat2$obs_time[teststndat2$datum == sdmax]),1)
  
#  sdnextday <- teststndat2$datum[as.Date(teststndat2$obs_time) == (as.Date(lastsdmax)+1)]
#  sddayplus2 <- teststndat2$datum[as.Date(teststndat2$obs_time) == (as.Date(lastsdmax)+2)]
#  stndatsorted <- teststndat2[order(teststndat2$datum,decreasing=TRUE),]
  
  
  
  if(length(sdnextday) == 0 | length(sddayplus2) == 0){
    nacheck <- TRUE
    zerocheck <- TRUE
  } else{
    nacheck <- FALSE
    if(sdnextday == 0 & sddayplus2 == 0){
      zerocheck <- TRUE
    } else{
      zerocheck <- FALSE
    }
  }
  
  if(nacheck == TRUE){
    sdcheck <- TRUE
  } else if(sdmax >= 30 & zerocheck == TRUE){
    sdcheck <- TRUE
  } else{
    sdcheck <- FALSE
  }
  
  
  i=1
  while((nacheck | zerocheck == TRUE) & sdcheck == TRUE){
    sdmax <- stndatsorted$datum[i]
    sdnextday <- stndat$datum[as.Date(stndat$obs_time) == (as.Date(stndatsorted$obs_time[i])+1)]
    sddayplus2 <- stndat$datum[as.Date(stndat$obs_time) == (as.Date(stndatsorted$obs_time[i])+2)]
    
    if(is.na(sdmax)){
      sdmax<-NA
      break
    }
    
    if(length(sdnextday) == 0 | length(sddayplus2) == 0){
      nacheck <- TRUE
      zerocheck <- TRUE
    } else{
        nacheck <- FALSE
        if(sdnextday == 0 & sddayplus2 == 0){
        zerocheck <- TRUE
      } else{
        zerocheck <- FALSE
      }
    }
    
    if(nacheck == TRUE){
      sdcheck <- TRUE
    } else if(sdmax >= 30 & zerocheck == TRUE){
      sdcheck <- TRUE
    } else{
      sdcheck <- FALSE
    }
    
    
    
    i = i+1
    #browser()
  }
  
  return(sdmax)
}



generate_hist <- function(composite_stns){
  for (i in 1:nrow(composite_stns)){
    yeardat <- t(composite_stns[i,52:120])
    years <- c(1950:2018)
    hist(yeardat)
    browser()
  }
}

composite_stns_v3<-read.csv('composite_stns_ext_v3_withstats.csv')
composite_stns_v4<-read.csv('composite_stns_ext_v4_withstats.csv')
composite_stns_v3 <- composite_stns_v3[composite_stns_v3$record_length >= 20,]
composite_stns_v3 <- composite_stns_v3[composite_stns_v3$average_distance <= 50,]
composite_stns_v3 <- composite_stns_v3[composite_stns_v3$average_elevation_difference <= 100,]

generate_boxplot <- function(composite_stns,fid){
   yeardat <- composite_stns[composite_stns$fid == fid,]
   yeardat <- t(yeardat[52:120])
   years <- c(1950:2018)
   layout(matrix(c(2,1,2,1), 2, 2, byrow = TRUE))
   hist(yeardat)
   boxplot(yeardat)
   return(boxplot.stats(yeardat)$out)
    
}


#more plots than just the boxplot, arranged nicely
generate_outlier_plots<-function(composite_stns,fid){

  yeardat <- composite_stns[composite_stns$fid == fid,]
  yeardat <- t(yeardat[52:120])
  years <- c(1950:2018)
  layout(matrix(c(1,1,2,3,2,3), nrow=3, byrow=TRUE))
  #layout.show(n=3)
  plot(years,yeardat,type='l',main=fid)
  hist(yeardat)
  boxplot(yeardat)
  outliers <- boxplot.stats(sort(t(yeardat),decreasing=TRUE))$out
  uppertail <- outliers[outliers > mean(yeardat[!is.na(yeardat)])]
  return(uppertail)
  #return(boxplot.stats(yeardat)$out)
}



for(i in 1:nrow(composite_stns)){
  n <- generate_outlier_plots(composite_stns,composite_stns$fid[i])
  print(n)
  browser()
}

generate_outlier_plots(composite_stns,713)



generate_outlier_stats <- function(composite_stns){
  outlier_stats_list <- list()
  for (i in 1:nrow(composite_stns)){
    nlist <- list()
    nlist[[1]] <- generate_outlier_plots(composite_stns,composite_stns$fid[i])
    nlist[2] <- length(nlist[[1]])
    nlist[[3]] <- sort(t(composite_stns[i,52:120]),decreasing=TRUE)[-(1:length(nlist[[1]]))]
    names(nlist) <- c('outliers','k','N-k')
    outlier_stats_list[[i]]<-nlist
  }
  names(outlier_stats_list) <- composite_stns$fid
  return(outlier_stats_list)
}

stations_with_outliers <- outlier_stats[lapply(outlier_stats,function(x) any(x[[1]]))==TRUE]
stations_with_outliers_old <- outlier_stats[lapply(outlier_stats_old,function(x) any(x[[1]]))==TRUE]

fitdistr(stations_with_outliers[[1]][[3]],'gamma')


difflist <- vector()
y=1
compositeoverlap<- composite_stns[composite_stns$fid %in% composite_stns_old$fid,]
for(i in 1:nrow(composite_stns_old)){
  if(length(composite_stns_old[i,!is.na(composite_stns_old[i,])]) == length(compositeoverlap[i,!is.na(compositeoverlap[i,])])){
    if(!all(composite_stns_old[i,!is.na(composite_stns_old[i,])] == compositeoverlap[i,!is.na(compositeoverlap[i,])])){
      difflist[y] <- composite_stns_old$fid[i] 
      y=y+1
    }
  } else{
    difflist[y] <- composite_stns_old$fid[i] 
    y=y+1
  }
  #mega_agg_old[i,!is.na(mega_agg_old[i,])] == mega_agg[i,!is.na(mega_agg[i,])]
  
}

temp<-as.vector(composite_stns_v3[82,8:121])

temp <- temp[!is.na(temp)]
temp<-sort(temp)

dgumbelII(temp,shape=16.11938,scale=16.4939)
pgumbel(temp,scale=16.4939)
qgumbel(temp,scale=16.4939)







          