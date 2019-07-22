library('crmp')
library('climdex.pcic')
library('geosphere')
library('RPostgreSQL')
library(dplyr)
#install_github("pacificclimate/crmp")
#install.packages('C:/Users/cdmb/Documents/R Scripts/crmp_0.5.tar.gz',repos=NULL,type='source')

#source('L:/home/cdmb/crmp_fanslow/crmp_climo_utils.r')

#crmp_obj <- pull_one_crmp(con,sid,hid)
load('L:/data/projects/crmp/crmp_climdex/sid_6608_hid_8293_8110_crmp_climdex_objs.RData')

netid_dictionary <- c('5'='BCH','2'='motie','3'='motim','12'='mofr','22'='ard','16'='FRBC','19'='ecraw','9'='moe_aqn')


#for the purposes of getting the daily coverage lookup table, or something like it, going:
#SELECT history_id from obs_count_per_month_history_mv WHERE date_trunc > '1880-12-31' AND date_trunc < '1881-02-01'  ;
#2019/04/23 going to change querystr to include record start and end

con <- dbConnect(PostgreSQL(),user="fanslow",host="monsoon.pcic.uvic.ca",dbname="crmp")
querystr <- 'SELECT meta_history.station_id, meta_history.history_id, meta_history.station_name, meta_history.lon, meta_history.lat, meta_history.elev, meta_history.sdate, meta_history.edate, meta_history.freq, meta_station.network_id, meta_station.native_id FROM meta_history NATURAL JOIN meta_station;'

try.rv <- try({
    rs <- dbSendQuery(con,querystr)
      stnmeta <- fetch(rs,n=-1)
}, silent=T)
if (inherits(try.rv, 'try-error')) {
    close_crmp_connections()
  print(try.rv)
}


  #want to get the daily, windowed biweighted mean, standard deviation, and record min and max values (only max for precip)
  #use a 15 day window and for temperature a 29 day window for precip and working only on days with precip greater than zero.
  #Use 6x bistd for temp and 9x 95p for precip

twindow <- 15
pptwindow <- 29
doys <- seq(1,365)
make_doymask_matrix <- function(windowsize) {
  doys <- seq(1,365)
  doyidxs <-matrix(data=NA,nrow=365,ncol=windowsize)
  for (didx  in seq(1,windowsize)) {
    doysnow <- doys-((windowsize-1)/2)+(didx-1)
    doysnow[doysnow<=0] <- doysnow[doysnow<=0]+365
    doysnow[doysnow>365] <- doysnow[doysnow>365]-365
    doyidxs[,didx] <- doysnow
  }
  return(doyidxs)
}


#Make a matrix that yields the set of indexes for a given doy for temp and precip. Mainly, this pre-specified the wrapping over a year as needed.
tdoyidxs <- make_doymask_matrix(twindow)
pdoyidxs <- make_doymask_matrix(pptwindow)
get_doy_por_stats <- function(crmp_obj) {


  #want to get the daily, windowed biweighted mean, standard deviation, and record min and max values (only max for precip)
  #use a 15 day window and for temperature a 29 day window for precip and working only on days with precip greater than zero.
  #Use 6x bistd for temp and 9x 95p for precip
  outdf <- as.data.frame(matrix(data=NA,nrow=366,ncol=15))
  names(outdf) <- c('doy','tn_bimean','tn_bimstd','tnx','tnn','tnN','tx_bimean','tx_bimstd','txx','txn','txN','ppt_bimean','ppt_95p','pptx','pptN')

  #creat a vector of doy from the crmp object
  thetimes <- slot(crmp_obj,'time')
  datadoys <- as.numeric(format(thetimes,format='%j'))



  #loop through and pull out the pool of data for the given variable
  for (adoy in seq(1,365)) {
    obj_doysub <- crmp_obj[which(adoy == datadoys),]
    minrecs <- apply(obj_doysub,2,min,na.rm=T)
    maxrecs <- apply(obj_doysub,2,max,na.rm=T)
    bimstats_onedoy <- as.data.frame(apply(obj_doysub,2,getbimstd,minobsthresh=0.2,na.rm=T))
    #browser()
    twinddoy <- tdoyidxs[adoy,]
    pwinddoy <- pdoyidxs[adoy,]
    tdataidxs <- which(datadoys %in% twinddoy)
    pdataidxs <- which(datadoys %in% pwinddoy)

    tnloc <- which.is.tmin(crmp_obj)
    txloc <- which.is.tmax(crmp_obj)
    pptloc <- which.is.precip(crmp_obj)
    #We have the indexes that define the pool, now pull out the data for the indexes
    obj_sub <- crmp_obj[tdataidxs,]
    tndata <- obj_sub[,tnloc]
    txdata <- obj_sub[,txloc]
    tnN <- sum(!is.na(tndata))
    txN <- sum(!is.na(txdata))
    bimstats_twind <- as.data.frame(apply(obj_sub,2,getbimstd,minobsthresh=0.2,na.rm=T))
    #tnN <- sum()
    obj_sub <- crmp_obj[pdataidxs,]
    bimstats_pwind <- as.data.frame(apply(obj_sub,2,getbimstd,minobsthresh=0.2,na.rm=T))
    pptdata <- obj_sub[,pptloc]
    pptidxs <- !is.na(pptdata) & (pptdata > 0)
    pptN <- sum(pptidxs)
    ppt95p <- quantile(pptdata[pptidxs],probs=c(0.95),type=8,na.rm=T)
    outvec <- c(adoy,bimstats_twind[(2*tnloc)],bimstats_twind[(2*tnloc-1)],maxrecs[tnloc],minrecs[tnloc],tnN,bimstats_twind[(2*txloc)],bimstats_twind[(2*txloc-1)],maxrecs[txloc],minrecs[txloc],txN,bimstats_pwind[(2*pptloc)],ppt95p,maxrecs[pptloc],pptN)
    outdf[adoy,] <- outvec
  }
  outdf[366,] <- apply(outdf[c(1,365),],2,mean,na.rm=T)
  outdf[366,1] <-366
  return(outdf)
}

#World record check
#Section 3 table 1
world_record_flag <- function(crmp_obj){
  tnloc <- which.is.tmin(crmp_obj)
  tnflags <- crmp_obj[2:nrow(crmp_obj),tnloc] < -87.4 | crmp_obj[2:nrow(crmp_obj),tnloc] > 57.7
  txloc <- which.is.tmax(crmp_obj)
  txflags <- crmp_obj[2:nrow(crmp_obj),txloc] < -87.4 | crmp_obj[2:nrow(crmp_obj),txloc] > 57.7
  pptloc <- which.is.precip(crmp_obj)
  pptflags <- crmp_obj[2:nrow(crmp_obj),pptloc] < 0 | crmp_obj[2:nrow(crmp_obj),pptloc] > 1925
  
  txflags[is.na(txflags)] <- FALSE
  tnflags[is.na(tnflags)] <- FALSE
  pptflags[is.na(pptflags)] <- FALSE
  return(list(tmaxflags=txflags,tminflags=tnflags,pptflags=pptflags))
  
}

#Repetition check
#Section 3 table 1
repetition_flag_old <- function(crmp_obj){
  tnloc <- which.is.tmin(crmp_obj)
  tn <- crmp_obj[2:nrow(crmp_obj),tnloc] 
  tnflags <- rep(FALSE,nrow(crmp_obj))
  counter = 0
   for (i in 2:nrow(crmp_obj)){
    if (!is.na(tn[i]) & !is.na(tn[i-1])){
     if (tn[i] == tn[i-1]) {
      counter = counter + 1
    } else {
      counter <- 0
    }
     if (counter > 20){
       a <- i-counter
       tnflags[a:i] <- TRUE
     }
    }
   }
  
  txloc <- which.is.tmax(crmp_obj)
  tx <- crmp_obj[2:nrow(crmp_obj),txloc] 
  txflags <- rep(FALSE,nrow(crmp_obj))
  counter = 0
  for (i in 2:nrow(crmp_obj)){
    if (!is.na(tx[i]) & !is.na(tx[i-1])){
      if (tx[i] == tx[i-1]) {
        counter = counter + 1
      } else {
        counter <- 0
      }
      if (counter > 20){
        a <- i-counter
        txflags[a:i] <- TRUE
      }
    }
  }
  
  pptloc <- which.is.precip(crmp_obj)
  ppt <- crmp_obj[2:nrow(crmp_obj),pptloc] 
  pptflags <- rep(FALSE,nrow(crmp_obj))
  counter = 0
  for (i in 2:nrow(crmp_obj)){
    if ((!is.na(ppt[i]) & !is.na(ppt[i-1])) & ppt[i] != 0){
      if (ppt[i] == ppt[i-1]) {
        counter = counter + 1
      } else {
        counter <- 0
      }
      if (counter >= 20){
        a <- i-counter
        pptflags[a:i] <- TRUE
      }
    }
  }
  txflags[is.na(txflags)] <- FALSE
  tnflags[is.na(tnflags)] <- FALSE
  pptflags[is.na(pptflags)] <- FALSE
  return(list(tmaxflags=txflags,tminflags=tnflags,pptflags=pptflags))
  
  
}


#The first rep flag works but is pretty slow, gonna try and make it a bit quicker using rle
#rle will include nas no matter how many there are in a row, so make sure to ignore any reps which are nas
#identical value streak check
repetition_flag <- function(crmp_obj){
  tnloc <- which.is.tmin(crmp_obj)
  tn <- crmp_obj[2:nrow(crmp_obj),tnloc]
  #rle doesn't like lists, gotta pass it a vector
  tn <- tn@.Data[[1]]
  tnflags <- rep(FALSE,nrow(crmp_obj))
  rtn <- rle(tn)
  reps <- which(!is.na(rtn$values) & rtn$lengths >= 20)
  rep.lengths.cumsum <- cumsum(rtn$lengths)
  ends <- rep.lengths.cumsum[reps]
  newindex = ifelse(reps>1, reps-1, 0)
  starts = rep.lengths.cumsum[newindex] + 1
  if (0 %in% newindex) starts = c(1,starts)
  for (i in 1:length(starts)){
    tnflags[starts[i]:ends[i]] <- TRUE
  }
  
  txloc <- which.is.tmax(crmp_obj)
  tx <- crmp_obj[2:nrow(crmp_obj),txloc]
  tx <- tx@.Data[[1]]
  txflags <- rep(FALSE,nrow(crmp_obj))
  rtx <- rle(tx)
  reps <- which(!is.na(rtx$values) & rtx$lengths >= 20)
  rep.lengths.cumsum <- cumsum(rtx$lengths)
  ends <- rep.lengths.cumsum[reps]
  newindex = ifelse(reps>1, reps-1, 0)
  starts = rep.lengths.cumsum[newindex] + 1
  if (0 %in% newindex) starts = c(1,starts)
  for (i in 1:length(starts)){
    txflags[starts[i]:ends[i]] <- TRUE
  }
  
  pptloc <- which.is.precip(crmp_obj)
  ppt <- crmp_obj[2:nrow(crmp_obj),pptloc]
  ppt <- ppt@.Data[[1]]
  pptflags <- rep(FALSE,nrow(crmp_obj))
  rppt <- rle(ppt)
  reps <- which(!is.na(rppt$values) & rppt$values > 0 & rppt$lengths >= 20)
  rep.lengths.cumsum <- cumsum(rppt$lengths)
  ends <- rep.lengths.cumsum[reps]
  newindex = ifelse(reps>1, reps-1, 0)
  starts = rep.lengths.cumsum[newindex] + 1
  if (0 %in% newindex) starts = c(1,starts)
  for (i in 1:length(starts)){
    pptflags[starts[i]:ends[i]] <- TRUE
  }
  
  txflags[is.na(txflags)] <- FALSE
  tnflags[is.na(tnflags)] <- FALSE
  pptflags[is.na(pptflags)] <- FALSE
  return(list(tmaxflags=txflags,tminflags=tnflags,pptflags=pptflags))
  
  
}

#Section 3 table 1
#This will require the calculation of more quantile data, gonna grab some of Faron's code from climate outlier check
#Seems to work well but I kinda want to do some more testing just to make sure
#identical value frequent value check
repetition_frequency_flag <- function(crmp_obj){

#modified the climate outlier function to create a dataframe of precip data with the required percentile info
  pptoutdf <- as.data.frame(matrix(data=NA,nrow=366,ncol=9))
  names(pptoutdf) <- c('doy','ppt_bimean','ppt_95p','ppt_90p','ppt_70p','ppt_50p','ppt_30p','pptx','pptN')
  
  #creat a vector of doy from the crmp object
  thetimes <- slot(crmp_obj,'time')
  datadoys <- as.numeric(format(thetimes,format='%j'))
  
  #loop through and pull out the pool of data for the given variable
  for (adoy in seq(1,365)) {
    obj_doysub <- crmp_obj[which(adoy == datadoys),]
    minrecs <- apply(obj_doysub,2,min,na.rm=T)
    maxrecs <- apply(obj_doysub,2,max,na.rm=T)
    bimstats_onedoy <- as.data.frame(apply(obj_doysub,2,getbimstd,minobsthresh=0.2,na.rm=T))
    #browser()
    pwinddoy <- pdoyidxs[adoy,]
    pdataidxs <- which(datadoys %in% pwinddoy)
    
    pptloc <- which.is.precip(crmp_obj)
    obj_sub <- crmp_obj[pdataidxs,]
    bimstats_pwind <- as.data.frame(apply(obj_sub,2,getbimstd,minobsthresh=0.2,na.rm=T))
    pptdata <- obj_sub[,pptloc]
    pptidxs <- !is.na(pptdata) & (pptdata > 0)
    pptN <- sum(pptidxs)
    ppt95p <- quantile(pptdata[pptidxs],probs=c(0.95),type=8,na.rm=T)
    ppt90p <- quantile(pptdata[pptidxs],probs=c(0.90),type=8,na.rm=T)
    ppt70p <- quantile(pptdata[pptidxs],probs=c(0.70),type=8,na.rm=T)
    ppt50p <- quantile(pptdata[pptidxs],probs=c(0.50),type=8,na.rm=T)
    ppt30p <- quantile(pptdata[pptidxs],probs=c(0.30),type=8,na.rm=T)
    #outvec <- c(adoy,bimstats_twind[(2*tnloc)],bimstats_twind[(2*tnloc-1)],maxrecs[tnloc],minrecs[tnloc],tnN,bimstats_twind[(2*txloc)],bimstats_twind[(2*txloc-1)],maxrecs[txloc],minrecs[txloc],txN,bimstats_pwind[(2*pptloc)],ppt95p,maxrecs[pptloc],pptN)
    outvec <- c(adoy,bimstats_pwind[(2*pptloc)],ppt95p,ppt90p,ppt70p,ppt50p,ppt30p,maxrecs[pptloc],pptN)
    pptoutdf[adoy,] <- outvec
  }
  
  pptoutdf[366,] <- apply(pptoutdf[c(1,365),],2,mean,na.rm=T)
  pptoutdf[366,1] <-366
  
  pptloc <- which.is.precip(crmp_obj)
  ppt <- crmp_obj[2:nrow(crmp_obj),pptloc]
  pptflags <- rep(FALSE,nrow(crmp_obj))
  
  pptdata <- obj_sub[,pptloc]
  pptidxs <- !is.na(pptdata) & (pptdata > 0)
  testflags <- rep(FALSE,length(testarr))
  
#9/10 value test  
  for (i in 1:nrow(ppt)){
    j <- i+9
    arr <- ppt[i:j]
    datareps <- rle(arr)
    runs <- which(datareps$lengths >=9 & datareps$values > 0)
    if(any(runs) == TRUE){
      runs.lengths.cumsum = cumsum(datareps$lengths)
      ends = runs.lengths.cumsum[runs]
      newindex = ifelse(runs>1, runs-1, 0)
      starts = runs.lengths.cumsum[newindex] + 1
      if (0 %in% newindex) starts = c(1,starts)
   
    #snew and enew should be the indices of the run in the original data  
      snew <- i + starts[1] - 1
      enew <- i+ ends[1] - 1
      repdata <- ppt[snew:enew]
   
    #create an array of the 30th percentile values for the corresponding days
      percentiles <- pptoutdf$ppt_30p[datadoys[snew]:datadoys[enew]]
      a <- as.logical(which(repdata >= percentiles))
      a <- sum(a == TRUE)
      if(sum(as.logical(which(repdata >= percentiles))) >= 9){
        pptflags[snew:enew] <- TRUE
      }
      
    }
  }

#8/10 value test  
  for (i in 1:nrow(ppt)){
    j <- i+9
    arr <- ppt[i:j]
    datareps <- rle(arr)
    runs <- which(datareps$lengths >= 8 & datareps$values > 0)
    if(any(runs) == TRUE){
      runs.lengths.cumsum = cumsum(datareps$lengths)
      ends = runs.lengths.cumsum[runs]
      newindex = ifelse(runs>1, runs-1, 0)
      starts = runs.lengths.cumsum[newindex] + 1
      if (0 %in% newindex) starts = c(1,starts)
      
      #snew and enew should be the indices of the run in the original data  
      snew <- i + starts[1] - 1
      enew <- i+ ends[1] - 1
      repdata <- ppt[snew:enew]
      
      #create an array of the 30th percentile values for the corresponding days
      percentiles <- pptoutdf$ppt_50p[datadoys[snew]:datadoys[enew]]
      a <- as.logical(which(repdata >= percentiles))
      a <- sum(a == TRUE)
      if(sum(as.logical(which(repdata >= percentiles))) >= 8){
        pptflags[snew:enew] <- TRUE
      }
      
    }
  }

#7/10 value test    
  for (i in 1:nrow(ppt)){
    j <- i+9
    arr <- ppt[i:j]
    datareps <- rle(arr)
    runs <- which(datareps$lengths >= 7 & datareps$values > 0)
    if(any(runs) == TRUE){
      runs.lengths.cumsum = cumsum(datareps$lengths)
      ends = runs.lengths.cumsum[runs]
      newindex = ifelse(runs>1, runs-1, 0)
      starts = runs.lengths.cumsum[newindex] + 1
      if (0 %in% newindex) starts = c(1,starts)
      
      #snew and enew should be the indices of the run in the original data  
      snew <- i + starts[1] - 1
      enew <- i+ ends[1] - 1
      repdata <- ppt[snew:enew]
      
      #create an array of the 30th percentile values for the corresponding days
      percentiles <- pptoutdf$ppt_70p[datadoys[snew]:datadoys[enew]]
      a <- as.logical(which(repdata >= percentiles))
      a <- sum(a == TRUE)
      if(sum(as.logical(which(repdata >= percentiles))) >= 7){
        pptflags[snew:enew] <- TRUE
      }
      
    }
  }
  
#5/10 value test  
  for (i in 1:nrow(ppt)){
    j <- i+9
    arr <- ppt[i:j]
    datareps <- rle(arr)
    runs <- which(datareps$lengths >=5 & datareps$values > 0)
    if(any(runs) == TRUE){
      runs.lengths.cumsum = cumsum(datareps$lengths)
      ends = runs.lengths.cumsum[runs]
      newindex = ifelse(runs>1, runs-1, 0)
      starts = runs.lengths.cumsum[newindex] + 1
      if (0 %in% newindex) starts = c(1,starts)
      
      #snew and enew should be the indices of the run in the original data  
      snew <- i + starts[1] - 1
      enew <- i+ ends[1] - 1
      repdata <- ppt[snew:enew]
  
      #create an array of the 30th percentile values for the corresponding days
      percentiles <- pptoutdf$ppt_90p[datadoys[snew]:datadoys[enew]]
      a <- as.logical(which(repdata >= percentiles))
      a <- sum(a == TRUE)
      if(sum(as.logical(which(repdata >= percentiles))) >= 5){
        pptflags[snew:enew] <- TRUE
      
    #checks to see if there's a second run of 5 within the testing array, then redoes the test for that run  
      if(length(runs) >= 2){
        runs.lengths.cumsum = cumsum(datareps$lengths)
        ends = runs.lengths.cumsum[runs]
        newindex = ifelse(runs>1, runs-1, 0)
        starts = runs.lengths.cumsum[newindex] + 1
        if (0 %in% newindex) starts = c(1,starts)
        
        #snew and enew are the indices of the run in the original data  
        snew <- i + starts[2] - 1
        enew <- i+ ends[2] - 1
        repdata <- ppt[snew:enew]
        
        #create an array of the 30th percentile values for the corresponding days
        percentiles <- pptoutdf$ppt_90p[datadoys[snew]:datadoys[enew]]
        a <- as.logical(which(repdata >= percentiles))
        a <- sum(a == TRUE)
        if(sum(as.logical(which(repdata >= percentiles))) >= 5){
          pptflags[snew:enew] <- TRUE
        
      }  
      }
      
    }
   }
  }
  return(pptflags)
}



hid = 8293
sid = 6608
maxdist=75000
#function that ranks the closest stations to the target station 
#stnmeta file needs start and end date
get_nearby_stns <- function(crmp_obj,stnmeta,hid,sid,maxdist){
  #maxdist is in meters
  targetstation <- stnmeta[stnmeta$history_id == hid,]
  distance_ranked <- stnmeta

 for (i in 1:nrow(stnmeta)){
    
    distance_ranked$distance[i] <- distm(c(targetstation$lon,targetstation$lat),c(distance_ranked$lon[i],distance_ranked$lat[i])
                                        ,fun=distHaversine)
 
   }
 
  
  distance_ranked <- distance_ranked[order(distance_ranked$distance),]
  distance_ranked <- distance_ranked[distance_ranked$distance < maxdist,]
  distance_ranked <- distance_ranked[distance_ranked$distance > 0,]
  distance_ranked <- distance_ranked[!is.na(distance_ranked$station_id),]
  distance_ranked <- distance_ranked[]
  
  #pulling the crmp obj for nearby stations seems to have an issue with stations with no data, most of which are EV-AQN (network ID 9)
  distance_ranked <- distance_ranked[distance_ranked$network_id != 9,]
  
  for (i in 1:nrow(distance_ranked)){
    if(null_test(distance_ranked$station_id[i],distance_ranked$history_id[i]) == TRUE){
      distance_ranked$is_null[i] = TRUE
    } else{
      distance_ranked$is_null[i] = FALSE
    }
    
    
  }
  distance_ranked <- distance_ranked[distance_ranked$is_null == FALSE,]
    
  return(distance_ranked)
  
}

#since null crmp_obj issue is back (2019/04/18) for no apparent reason, this test will try and boot empty objs
#this makes the get_nearby_stns function much slower, ugh
null_test <- function(stnid,histid){
  con <- dbConnect(PostgreSQL(),user="cdmb",password="7mU@u&",host="monsoon.pcic.uvic.ca",dbname="crmp")
 # histid <- distance_ranked$history_id[i]
#  stnid <- distance_ranked$station_id[i]
  
  test_crmp_obj <- pull_one_crmp(con,stnid,histid)
  close_crmp_connections()
  
  if(is.null(test_crmp_obj)){
    return(TRUE)
  } else {
    return(FALSE)
  }
  
}  
 

#query for grabbing timestamp of the most recent observation from a station 
  con <- dbConnect(PostgreSQL(),user="cdmb",password="7mU@u&",host="monsoon.pcic.uvic.ca",dbname="crmp")
  querystr <- paste('select max(obs_time) from crmp.obs_raw where history_id = ',histid)
  try.rv <- try({
    rs <- dbSendQuery(con,querystr)
    lastobs <- fetch(rs,n=-1)
  }, silent=T)
  if (inherits(try.rv, 'try-error')) {
    close_crmp_connections()
    print(try.rv)
  }
  
  
  


#distance_ranked is list of nearest stations generated by get_nearby_stations
spatial_corroboration_temp_flag <- function(crmp_obj,distance_ranked){
  
  #calculates stats for each of the closest stations in order to calculate anomalies
  #as well as an accompanying day vector for indexing purposes
  test_stn_stats <- list()
  doyvectlist <- list()
  if (nrow(distance_ranked > 100)){
    nstns = 100
  } else {
    nstns = nrow(distance_ranked)
  }
  
   for (i in 1:nstns){
    con <- dbConnect(PostgreSQL(),user="cdmb",password="7mU@u&",host="monsoon.pcic.uvic.ca",dbname="crmp")
    histid <- distance_ranked$history_id[i]
    stnid <- distance_ranked$station_id[i]
    
    test_crmp_obj <- pull_one_crmp(con,stnid,histid)
    close_crmp_connections()
    
    #if(!is.nulltest_crmp_obj)
    
    # obj_stats <- get_doy_por_stats(test_crmp_obj)
    test_stn_stats[[i]] <- get_doy_por_stats(test_crmp_obj)
    
    #create a vector of doy from the crmp object
    thetimes <- slot(test_crmp_obj,'time')
    datadoys <- as.numeric(format(thetimes,format='%j'))
    doyvectlist[[i]] <- datadoys
  }
  close_crmp_connections()
  stnnames <- distance_ranked$history_id[1:nstns]
  names(test_stn_stats) <- c(stnnames)
  names(doyvectlist) <- c(stnnames)
  
  #calculates stats and doy list for target station
  targetstationstats <- get_doy_por_stats(crmp_obj)
  thetimes <- slot(crmp_obj,'time')
  targetdatadoys <- as.numeric(format(thetimes,format='%j'))
  
  
  
  #minimum temperature anomaly test  
  tnloc <- which.is.tmin(crmp_obj)
  tn <- crmp_obj[2:nrow(crmp_obj),tnloc] 
  tnflags <- rep(FALSE,nrow(crmp_obj))
  tn_id <- '(427,463,469,493,511,517,545,633)'
  
  
  for(i in 2:(nrow(tn)-1)){
    dayminone <- tn@time[i-1]
    doyminone <- as.numeric(format(dayminone,format='%j'))
    dayzero <- tn@time[i]
    doyzero <- as.numeric(format(dayzero,format='%j'))
    dayplusone <- tn@time[i+1]
    doyplusone <- as.numeric(format(dayplusone,format='%j'))
    
    #this grabs the days within the window at the nearby stations
    hilist <- distance_ranked$history_id[1:nstns]
    hilist <- hilist[!is.na(hilist)]
    stndat <- list()
    con <- dbConnect(PostgreSQL(),user="cdmb",password="7mU@u&",host="monsoon.pcic.uvic.ca",dbname="crmp")
    for (hid in 1:length(hilist)){
      
      querystr <- paste('SELECT * FROM crmp.obs_raw WHERE history_id = ',hilist[hid],' AND ','(','obs_time BETWEEN ',
                        '\'',dayminone,'\'',' AND ','\'',dayplusone,'\'',')', ' AND vars_id in ', tn_id,sep='')
        try.rv <- try({
        rs <- dbSendQuery(con,querystr)
        stndat[[hid]] <- fetch(rs,n=-1)
        #temp <- fetch(rs,n=-1)
      }, silent=T)
      if (inherits(try.rv, 'try-error')) {
        close_crmp_connections()
        print(try.rv)
      }
      
    }
    close_crmp_connections()
    
    #list of nearby stations with data for desired period  
    names(stndat) <- c(hilist)
    nearbystations <- Filter(length,stndat)
    
    #it looks like stations need data for all three days, so omit stations that have length < 3
    nearbystations <- nearbystations[which(lapply(nearbystations, function(x) nrow(x) >= 3) == TRUE)]
    
    #Grabs 3-7 closest stations, returns a list or skips test due to lack of appropriate stations
    if (length(nearbystations) >= 3){
      if (length(nearbystations) > 7){
        teststns <- nearbystations[1:7]
      } else {
        teststns <- nearbystations
      }
    }
    
    if(length(nearbystations) >=3){
    
      #time to calculate the anomalies for the test stations 
      dayminone_anomalies <- vector()
      dayzero_anomalies <- vector()
      dayplusone_anomalies <- vector()
      for(obj in 1:length(teststns)){
        
        hid <- names(teststns[obj])
        stnstats <- test_stn_stats[[hid]]
        
        for (dy in 1:nrow(teststns[[hid]])){
          doy <- as.numeric(format(teststns[[hid]][dy,]$obs_time,format='%j')) 
          teststns[[hid]][dy,7] <- teststns[[hid]][dy,]$datum - stnstats$tn_bimean[doy]
        }
        dayminone_anomalies[obj] <- teststns[[hid]][1,7]
        dayzero_anomalies[obj] <- teststns[[hid]][2,7]
        dayplusone_anomalies[obj] <- teststns[[hid]][3,7]
      }
      all_anomalies <- c(dayminone_anomalies,dayzero_anomalies,dayplusone_anomalies)
      
      
        #then calculate for the target station
        if(!is.na(tn[i])){
        targetdat <- c(tn[i-1],tn[i],tn[i+1])
        targetdays <- c(dayminone,dayzero,dayplusone)
        targetdat <- as.data.frame(cbind(targetdays,targetdat))
        names(targetdat) <- c('obs_time','datum')
        targetdat$obs_time <- as.Date(targetdays)
      
        doy <- as.numeric(format(targetdat$obs_time,format='%j')) 
        targetdat[,3] <- targetdat$datum - targetstationstats$tn_bimean[doy]
      
        #now we can compare anomalies
        anomalytest <- abs(targetdat$V3[2]) - abs(all_anomalies)
        anomalytest <- (abs(anomalytest)) >= 10
        
        if(all(anomalytest) == TRUE){
          tnflags[i] <- TRUE
        } 
      }
    }
  }
  
  
  
  #max temperature anomaly flag
  txloc <- which.is.tmax(crmp_obj)
  tx <- crmp_obj[2:nrow(crmp_obj),txloc] 
  txflags <- rep(FALSE,nrow(crmp_obj))
  tx_id <- c('(428,464,470,494,512,518,544,632)')
  
  
  for(i in 2:(nrow(tx)-1)){
    dayminone <- tx@time[i-1]
    doyminone <- as.numeric(format(dayminone,format='%j'))
    dayzero <- tx@time[i]
    doyzero <- as.numeric(format(dayzero,format='%j'))
    dayplusone <- tx@time[i+1]
    doyplusone <- as.numeric(format(dayplusone,format='%j'))
    
    #this grabs the days within the window at the nearby stations
    hilist <- distance_ranked$history_id[1:nstns]
    hilist <- hilist[!is.na(hilist)]
    stndat <- list()
    con <- dbConnect(PostgreSQL(),user="cdmb",password="7mU@u&",host="monsoon.pcic.uvic.ca",dbname="crmp")
    for (hid in 1:length(hilist)){
      
      querystr <- paste('SELECT * FROM crmp.obs_raw WHERE history_id = ',hilist[hid],' AND ','(','obs_time BETWEEN ',
                        '\'',dayminone,'\'',' AND ','\'',dayplusone,'\'',')', ' AND vars_id in ', tx_id,sep='')
      try.rv <- try({
        rs <- dbSendQuery(con,querystr)
        stndat[[hid]] <- fetch(rs,n=-1)
        #temp <- fetch(rs,n=-1)
      }, silent=T)
      if (inherits(try.rv, 'try-error')) {
        close_crmp_connections()
        print(try.rv)
      }
      
    }
    close_crmp_connections()
    
    #list of nearby stations with data for desired period  
    names(stndat) <- c(hilist)
    nearbystations <- Filter(length,stndat)
    
    #it looks like stations need data for all three days, so omit stations that have length < 3
    nearbystations <- nearbystations[which(lapply(nearbystations, function(x) nrow(x) >= 3) == TRUE)]
    
    #Grabs 3-7 closest stations, returns a list or skips test due to lack of appropriate stations
    if (length(nearbystations) >= 3){
      if (length(nearbystations) > 7){
        teststns <- nearbystations[1:7]
      } else {
        teststns <- nearbystations
      }
    }
    
    if(length(nearbystations) >=3){
      
      #time to calculate the anomalies for the test stations 
      dayminone_anomalies <- vector()
      dayzero_anomalies <- vector()
      dayplusone_anomalies <- vector()
      for(obj in 1:length(teststns)){
        
        hid <- names(teststns[obj])
        stnstats <- test_stn_stats[[hid]]
        
        for (dy in 1:nrow(teststns[[hid]])){
          doy <- as.numeric(format(teststns[[hid]][dy,]$obs_time,format='%j')) 
          teststns[[hid]][dy,7] <- teststns[[hid]][dy,]$datum - stnstats$tx_bimean[doy]
        }
        dayminone_anomalies[obj] <- teststns[[hid]][1,7]
        dayzero_anomalies[obj] <- teststns[[hid]][2,7]
        dayplusone_anomalies[obj] <- teststns[[hid]][3,7]
      }
      all_anomalies <- c(dayminone_anomalies,dayzero_anomalies,dayplusone_anomalies)
      
      
      #then calculate for the target station
      if(!is.na(tx[i])){
        targetdat <- c(tx[i-1],tx[i],tx[i+1])
        targetdays <- c(dayminone,dayzero,dayplusone)
        targetdat <- as.data.frame(cbind(targetdays,targetdat))
        names(targetdat) <- c('obs_time','datum')
        targetdat$obs_time <- as.Date(targetdays)
        
        doy <- as.numeric(format(targetdat$obs_time,format='%j')) 
        targetdat[,3] <- targetdat$datum - targetstationstats$tx_bimean[doy]
        
        #now we can compare anomalies
        anomalytest <- abs(targetdat$V3[2]) - abs(all_anomalies)
        anomalytest <- (abs(anomalytest)) >= 10
        
        if(all(anomalytest) == TRUE){
          txflags[i] <- TRUE
        } 
      }
    }
  }
  
  tflags <- cbind(txflags,tnflags)
  return(tflags)
  
}

#rewritting this one with the same updates as the new spatial corroboration ppt test
spatial_corroboration_temp_flag_new <- function(crmp_obj,distance_ranked){
  
  #calculates stats for each of the closest stations in order to calculate anomalies
  #as well as an accompanying day vector for indexing purposes
  test_stn_stats <- list()
  doyvectlist <- list()
  test_stn_objs <- list()
  target_crmp_obj <- crmp_obj
  if (nrow(distance_ranked) > 100){
    nstns = 100
  } else {
    nstns = nrow(distance_ranked)
  }
  
  for (i in 1:nstns){
    histid <- distance_ranked$history_id[i]
    stnid <- distance_ranked$station_id[i]
    
    #data in /crmp_climdex is hourly to daily crmp_objs, if the data isn't in there it should be daily in obs_raw
    #and we should be able to pull it from the database
    
    if(file.exists(paste('L:/data/projects/crmp/crmp_climdex/sid_',stnid,'_hid_',histid,'_8110_crmp_climdex_objs.RData',sep=''))){
      load(paste('L:/data/projects/crmp/crmp_climdex/sid_',stnid,'_hid_',histid,'_8110_crmp_climdex_objs.RData',sep=''))
      test_crmp_obj <- crmp_obj
    } else{
      con <- dbConnect(PostgreSQL(),user="cdmb",password="7mU@u&",host="monsoon.pcic.uvic.ca",dbname="crmp")
      test_crmp_obj <- pull_one_crmp(con,stnid,histid)
      close_crmp_connections()
      
    }
    
    # obj_stats <- get_doy_por_stats(test_crmp_obj)
    test_stn_objs[[i]] <- test_crmp_obj
    test_stn_stats[[i]] <- get_doy_por_stats(test_crmp_obj)
    
    #create a vector of doy from the crmp object
    thetimes <- slot(test_crmp_obj,'time')
    datadoys <- as.numeric(format(thetimes,format='%j'))
    doyvectlist[[i]] <- datadoys
    
  }
  close_crmp_connections()
  stnnames <- distance_ranked$history_id[1:nstns]
  names(test_stn_objs) <- c(stnnames)
  names(test_stn_stats) <- c(stnnames)
  names(doyvectlist) <- c(stnnames)
  crmp_obj <- target_crmp_obj
  
  #calculates stats and doy list for target station
  targetstationstats <- get_doy_por_stats(crmp_obj)
  thetimes <- slot(crmp_obj,'time')
  targetdatadoys <- as.numeric(format(thetimes,format='%j'))
  
  
  
  #minimum temperature anomaly test  
  tnloc <- which.is.tmin(crmp_obj)
  tn <- crmp_obj[2:nrow(crmp_obj),tnloc] 
  tnflags <- rep(FALSE,nrow(crmp_obj))
  tn_id <- '427'
  
  
  for(i in 2:(nrow(tn)-1)){
    
    dayzero <- tn@time[i]
    doyzero <- as.numeric(format(dayzero,format='%j'))
    
    #dayminone <- ppt@time[i-1]
    dayminone <- as.Date(dayzero) -1 
    doyminone <- as.numeric(format(dayminone,format='%j'))
    
    # dayplusone <- ppt@time[i+1]
    dayplusone <- as.Date(dayzero) +1 
    doyplusone <- as.numeric(format(dayplusone,format='%j'))
    
    
    
    
    
    #this grabs the days within the window at the nearby stations
    hilist <- distance_ranked$history_id[1:nstns]
    hilist <- hilist[!is.na(hilist)]
   # hilist <- names(test_stn_objs)
    stndat <- list()
    #con <- dbConnect(PostgreSQL(),user="cdmb",password="7mU@u&",host="monsoon.pcic.uvic.ca",dbname="crmp")
    for (hid in 1:length(hilist)){
      
      testtminloc <- which.is.tmin(test_stn_objs[[hid]])
      if(any(testtminloc)){
        testtmin <- test_stn_objs[[hid]][2:nrow(test_stn_objs[[hid]]),testtminloc]
        testtime <- test_stn_objs[[hid]]@time[2:nrow(test_stn_objs[[hid]])]
      
        temp <- c(test_stn_objs[[hid]][[testtminloc]][test_stn_objs[[hid]]@time == dayminone],
                test_stn_objs[[hid]][[testtminloc]][test_stn_objs[[hid]]@time == dayzero],
                test_stn_objs[[hid]][[testtminloc]][test_stn_objs[[hid]]@time == dayplusone])
        
        temptime <- c(dayminone,dayzero,dayplusone)
        

       if(length(temp) < 3){
         stndat[[hid]] <- list()
       } else{
         stndat[[hid]] <- as.data.frame(as.POSIXct(temptime, "%Y-%m-%d"))
         stndat[[hid]][,2] <- temp
         colnames(stndat[[hid]]) <- c("obs_time","datum")
        }
      } else{
        stndat[[hid]] <- list()
      }
     
    }
   # close_crmp_connections()
    
    #list of nearby stations with data for desired period  
    names(stndat) <- c(hilist)
    nearbystations <- Filter(length,stndat)
    
    #it looks like stations need data for all three days, so omit stations that have length < 3
    nearbystations <- nearbystations[which(lapply(nearbystations, function(x) nrow(x) >= 3) == TRUE)]
    nearbystations <- nearbystations[which(lapply(nearbystations,function(x) any(is.na(x)))==FALSE)]
    
    #Grabs 3-7 closest stations, returns a list or skips test due to lack of appropriate stations
    if (length(nearbystations) >= 3){
      if (length(nearbystations) > 7){
        teststns <- nearbystations[1:7]
      } else {
        teststns <- nearbystations
      }
    }
    
    if(length(nearbystations) >=3){
      
      #time to calculate the anomalies for the test stations 
      dayminone_anomalies <- vector()
      dayzero_anomalies <- vector()
      dayplusone_anomalies <- vector()
      for(obj in 1:length(teststns)){
        
        hid <- names(teststns[obj])
        stnstats <- test_stn_stats[[hid]]
        
        for (dy in 1:nrow(teststns[[hid]])){
          doy <- as.numeric(format(teststns[[hid]][dy,]$obs_time,format='%j')) 
          teststns[[hid]][dy,3] <- teststns[[hid]][dy,]$datum - stnstats$tn_bimean[doy]
        }
        dayminone_anomalies[obj] <- teststns[[hid]][1,3]
        dayzero_anomalies[obj] <- teststns[[hid]][2,3]
        dayplusone_anomalies[obj] <- teststns[[hid]][3,3]
      }
      all_anomalies <- c(dayminone_anomalies,dayzero_anomalies,dayplusone_anomalies)
      
      
      #then calculate for the target station
      if(!is.na(tn[i])){
        targetdat <- c(tn[i-1],tn[i],tn[i+1])
        targetdays <- c(dayminone,dayzero,dayplusone)
        targetdat <- as.data.frame(cbind(targetdays,targetdat))
        names(targetdat) <- c('obs_time','datum')
        targetdat$obs_time <- as.Date(targetdays)
        
        doy <- as.numeric(format(targetdat$obs_time,format='%j')) 
        targetdat[,3] <- targetdat$datum - targetstationstats$tn_bimean[doy]
        
        #now we can compare anomalies
        anomalytest <- abs(targetdat$V3[2]) - abs(all_anomalies)
        anomalytest <- (abs(anomalytest)) >= 10
        
        if(all(anomalytest) == TRUE){
          tnflags[i] <- TRUE
        } 
      }
    }
  }
  
  
  #max temperature anomaly flag
  txloc <- which.is.tmax(crmp_obj)
  tx <- crmp_obj[2:nrow(crmp_obj),txloc] 
  txflags <- rep(FALSE,nrow(crmp_obj))
  tx_id <- '428'
  
  
  for(i in 2:(nrow(tx)-1)){
    dayzero <- tn@time[i]
    doyzero <- as.numeric(format(dayzero,format='%j'))
    
    #dayminone <- ppt@time[i-1]
    dayminone <- as.Date(dayzero) -1 
    doyminone <- as.numeric(format(dayminone,format='%j'))
    
    # dayplusone <- ppt@time[i+1]
    dayplusone <- as.Date(dayzero) +1 
    doyplusone <- as.numeric(format(dayplusone,format='%j'))
    
    
    
    #this grabs the days within the window at the nearby stations
    hilist <- distance_ranked$history_id[1:nstns]
    hilist <- hilist[!is.na(hilist)]
    stndat <- list()
    #con <- dbConnect(PostgreSQL(),user="cdmb",password="7mU@u&",host="monsoon.pcic.uvic.ca",dbname="crmp")
    for (hid in 1:length(hilist)){
      
      testtmaxloc <- which.is.tmax(test_stn_objs[[hid]])
      if(any(testtmaxloc)){
        testtmax <- test_stn_objs[[hid]][2:nrow(test_stn_objs[[hid]]),testtmaxloc]
        testtime <- test_stn_objs[[hid]]@time[2:nrow(test_stn_objs[[hid]])]
        
        temp <- c(test_stn_objs[[hid]][[testtmaxloc]][test_stn_objs[[hid]]@time == dayminone],
                  test_stn_objs[[hid]][[testtmaxloc]][test_stn_objs[[hid]]@time == dayzero],
                  test_stn_objs[[hid]][[testtmaxloc]][test_stn_objs[[hid]]@time == dayplusone])
        
        temptime <- c(dayminone,dayzero,dayplusone)
        
        
        if(length(temp) < 3){
          stndat[[hid]] <- list()
        } else{
          stndat[[hid]] <- as.data.frame(as.POSIXct(temptime, "%Y-%m-%d"))
          stndat[[hid]][,2] <- temp
          colnames(stndat[[hid]]) <- c("obs_time","datum")
        }
      } else{
        stndat[[hid]] <- list()
      }
     
      
    }
    #close_crmp_connections()
    
    #list of nearby stations with data for desired period  
    names(stndat) <- c(hilist)
    nearbystations <- Filter(length,stndat)
    
    #it looks like stations need data for all three days, so omit stations that have length < 3
    nearbystations <- nearbystations[which(lapply(nearbystations, function(x) nrow(x) >= 3) == TRUE)]
    nearbystations <- nearbystations[which(lapply(nearbystations,function(x) any(is.na(x)))==FALSE)]
    
    #Grabs 3-7 closest stations, returns a list or skips test due to lack of appropriate stations
    if (length(nearbystations) >= 3){
      if (length(nearbystations) > 7){
        teststns <- nearbystations[1:7]
      } else {
        teststns <- nearbystations
      }
    }
    
    if(length(nearbystations) >=3){
      
      #time to calculate the anomalies for the test stations 
      dayminone_anomalies <- vector()
      dayzero_anomalies <- vector()
      dayplusone_anomalies <- vector()
      for(obj in 1:length(teststns)){
        
        hid <- names(teststns[obj])
        stnstats <- test_stn_stats[[hid]]
        
        for (dy in 1:nrow(teststns[[hid]])){
          doy <- as.numeric(format(teststns[[hid]][dy,]$obs_time,format='%j')) 
          teststns[[hid]][dy,3] <- teststns[[hid]][dy,]$datum - stnstats$tx_bimean[doy]
        }
        dayminone_anomalies[obj] <- teststns[[hid]][1,3]
        dayzero_anomalies[obj] <- teststns[[hid]][2,3]
        dayplusone_anomalies[obj] <- teststns[[hid]][3,3]
      }
      all_anomalies <- c(dayminone_anomalies,dayzero_anomalies,dayplusone_anomalies)
      
      
      #then calculate for the target station
      if(!is.na(tx[i])){
        targetdat <- c(tx[i-1],tx[i],tx[i+1])
        targetdays <- c(dayminone,dayzero,dayplusone)
        targetdat <- as.data.frame(cbind(targetdays,targetdat))
        names(targetdat) <- c('obs_time','datum')
        targetdat$obs_time <- as.Date(targetdays)
        
        doy <- as.numeric(format(targetdat$obs_time,format='%j')) 
        targetdat[,3] <- targetdat$datum - targetstationstats$tx_bimean[doy]
        
        #now we can compare anomalies
        anomalytest <- abs(targetdat$V3[2]) - abs(all_anomalies)
        anomalytest <- (abs(anomalytest)) >= 10
        
        if(all(anomalytest) == TRUE){
          txflags[i] <- TRUE
        } 
      }
    }
  }
  
  tflags <- cbind(txflags,tnflags)
  return(tflags)
  
}


spatial_corroboration_ppt_flag <- function(crmp_obj,distance_ranked){
  #calculates stats for each of the closest stations in order to calculate anomalies
  #as well as an accompanying day vector for indexing purposes
  test_stn_stats <- list()
  doyvectlist <- list()
  if (nrow(distance_ranked > 100)){
    nstns = 100
  } else {
    nstns = nrow(distance_ranked)
  }
  
  for (i in 1:nstns){
    con <- dbConnect(PostgreSQL(),user="cdmb",password="7mU@u&",host="monsoon.pcic.uvic.ca",dbname="crmp")
    histid <- distance_ranked$history_id[i]
    stnid <- distance_ranked$station_id[i]
    
    test_crmp_obj <- pull_one_crmp(con,stnid,histid)
    close_crmp_connections()
    
  #  if(!is.null(test_crmp_obj)){

    test_stn_stats[[i]] <- get_doy_por_stats(test_crmp_obj)
    
    #create a vector of doy from the crmp object
    thetimes <- slot(test_crmp_obj,'time')
    datadoys <- as.numeric(format(thetimes,format='%j'))
    doyvectlist[[i]] <- datadoys
  #  } else{
  #    test_stn_stats[[i]] <- list()
  #    doyvectlist[[i]] <- vector()
  #  }
  }
  close_crmp_connections()
  stnnames <- distance_ranked$history_id[1:nstns]
  names(test_stn_stats) <- c(stnnames)
  names(doyvectlist) <- c(stnnames)
  
  
  
  
  #calculates stats and doy list for target station
  targetstationstats <- get_doy_por_stats(crmp_obj)
  thetimes <- slot(crmp_obj,'time')
  targetdatadoys <- as.numeric(format(thetimes,format='%j'))
  
  
  
  
  pptloc <- which.is.precip(crmp_obj)
  ppt <- crmp_obj[2:nrow(crmp_obj),pptloc] 
  pptflags <- rep(FALSE,nrow(crmp_obj))
  #list of vars ids that indicate daily precip totals
  #ppt_ids <- c("519", "496", "429", "513", "501", "526", "527", "441", "495", "471", "456", "465", "546", "631", "639")
  ppt_ids <- c('(519, 429, 513, 501, 526, 527, 441, 495, 471, 456, 465, 546, 631, 639)')
  
  #vars_ppt<- read.csv('C:/Users/cdmb/Downloads/data-1555603111544.csv')             
               
  recordstart <- format(as.Date(ppt@time[2]),'%Y')
  recordend <- format(as.Date(ppt@time[length(ppt@time)]),'%Y')
  
  
  
  for(i in 2:(nrow(ppt)-1)){
   if(!is.na(ppt[i])){
    dayminone <- ppt@time[i-1]
    doyminone <- as.numeric(format(dayminone,format='%j'))
    dayzero <- ppt@time[i]
    doyzero <- as.numeric(format(dayzero,format='%j'))
    dayplusone <- ppt@time[i+1]
    doyplusone <- as.numeric(format(dayplusone,format='%j'))
   
    #dates for 29 day window
    windowstart <- format(as.Date(ppt@time[i]) - 14,'%m-%d')
    windowstartdow <- as.numeric(format(as.Date(ppt@time[i]) - 14,format='%j')) 
    windowend <- format(as.Date(ppt@time[i])+14,'%m-%d')
    windowenddow <- as.numeric(format(as.Date(ppt@time[i])+14,format='%j')) 
    
    
    #this grabs the data within the 3 day window at the nearby stations
    hilist <- distance_ranked$history_id[1:nstns]
    hilist <- hilist[!is.na(hilist)]
    stndat <- list()
    con <- dbConnect(PostgreSQL(),user="cdmb",password="7mU@u&",host="monsoon.pcic.uvic.ca",dbname="crmp")
    for (hid in 1:length(hilist)){
      querystr <- paste('SELECT * FROM crmp.obs_raw WHERE history_id = ',hilist[hid],' AND ','(','obs_time BETWEEN ',
                        '\'',dayminone,'\'',' AND ','\'',dayplusone,'\'',')', ' AND vars_id in ', ppt_ids, sep='')
      try.rv <- try({
         rs <- dbSendQuery(con,querystr)
        stndat[[hid]] <- fetch(rs,n=-1)
        #temp <- fetch(rs,n=-1)
      }, silent=T)
      if (inherits(try.rv, 'try-error')) {
         close_crmp_connections()
        print(try.rv)
      }
      
    }
    close_crmp_connections()
    names(stndat) <- c(hilist)
    
    #list of nearby stations with data for desired period  
    nearbystations <- Filter(length,stndat)
    
    #it looks like stations need data for all days, so omit stations that have length < 3 and < 29 
    nearbystations <- nearbystations[which(lapply(nearbystations, function(x) nrow(x) >= 3) == TRUE)]
   
    #extract stations with pptN >20
    test_stn_stats2 <- test_stn_stats[which(lapply(test_stn_stats,function(x) x[doyzero,'pptN'] >= 20) == TRUE)]
    
    nearbystations <- nearbystations[names(nearbystations) %in% names(test_stn_stats2)]

    #Grabs 3-7 closest stations, returns a list or skips test due to lack of appropriate stations
    if (length(nearbystations) >= 3){
      if (length(nearbystations) > 7){
        teststns <- nearbystations[1:7]
        #testwindow <- nearbystationwindow[1:7]
      } else {
        teststns <- nearbystations
        #testwindow <- nearbystationwindow
      }
    }
    
    #get the minimum absolute target neighbour difference
    if(length(nearbystations) >=3){
      #step 1 - minimum absolute target neighbour difference
      nearbyprecip <- lapply(teststns, '[[', 'datum')
      nearbyprecip <- unlist(nearbyprecip)
      diffs <- ppt[i] - nearbyprecip
      if(ppt[i] > max(nearbyprecip) | ppt[i] < min(nearbyprecip)){
        mintargetdiff <- abs(min(diffs))
      } else{
        mintargetdiff <- 0
      }
      
      #calculate climatological percent ranks for the test stations at the target day
      teststnclimpercentrank <- vector()
      for (station in 1:length(teststns)){
        dayval <- teststns[[station]][,'datum']
        meanval <- test_stn_stats2[[names(teststns[station])]][doyminone:doyplusone,'ppt_bimean']
        teststnclimpercentrank <- c(teststnclimpercentrank,(dayval * 100)/meanval)
      }
      
      #calculate climatological percent ranks for the target station at the target day
      dayval <- ppt@.Data[[1]][i]
      meanval <- targetstationstats[doyzero,'ppt_bimean']
      targetclimpercentrank <- (dayval * 100)/meanval
      
      mintargetdiff_prank <- min(abs(targetclimpercentrank - teststnclimpercentrank))
      test_threshold <- -45.72*log10(mintargetdiff_prank) + 269.24
      
      if(mintargetdiff > test_threshold){
        pptflags[i] <- TRUE
      }
      
    } else{
      #part 4 - test threshold
      nearbyprecip <- lapply(teststns, '[[', 'datum')
      nearbyprecip <- unlist(nearbyprecip)
      diffs <- ppt[i] - nearbyprecip
      if(ppt[i] > max(nearbyprecip) | ppt[i] < min(nearbyprecip)){
        mintargetdiff <- abs(min(diffs))
      } else{
        mintargetdiff <- 0
      }
      test_threshold <- 269.24
      if(mintargetdiff > test_threshold){
        pptflags[i] <- TRUE
      }
    }
   }
  }
  return(pptflags)
  
}



#going to try and redo spatial corrob using test station crmp objs instead of obs_raw
#this is faster and grabs more stations than old spatial corroboration ppt flag
spatial_corroboration_ppt_flag_new <- function(crmp_obj,distance_ranked){
  #calculates stats for each of the closest stations in order to calculate anomalies
  #as well as an accompanying day vector for indexing purposes
  targetstnobj <- crmp_obj
  test_stn_stats <- list()
  test_stn_objs <- list()
  doyvectlist <- list()
  if (nrow(distance_ranked) > 100){
    nstns = 100
  } else {
    nstns = nrow(distance_ranked)
  }
  
  for (i in 1:nstns){
    histid <- distance_ranked$history_id[i]
    stnid <- distance_ranked$station_id[i]
    
    #data in /crmp_climdex is hourly to daily crmp_objs, if the data isn't in there it should be daily in obs_raw
    #and we should be able to pull it from the database
    
    if(file.exists(paste('L:/data/projects/crmp/crmp_climdex/sid_',stnid,'_hid_',histid,'_8110_crmp_climdex_objs.RData',sep=''))){
      load(paste('L:/data/projects/crmp/crmp_climdex/sid_',stnid,'_hid_',histid,'_8110_crmp_climdex_objs.RData',sep=''))
      test_crmp_obj <- crmp_obj
    } else{
      con <- dbConnect(PostgreSQL(),user="cdmb",password="7mU@u&",host="monsoon.pcic.uvic.ca",dbname="crmp")
      test_crmp_obj <- pull_one_crmp(con,stnid,histid)
      close_crmp_connections()
      
    }
    
    # obj_stats <- get_doy_por_stats(test_crmp_obj)
    test_stn_objs[[i]] <- test_crmp_obj
    test_stn_stats[[i]] <- get_doy_por_stats(test_crmp_obj)
      
    #create a vector of doy from the crmp object
    thetimes <- slot(test_crmp_obj,'time')
    datadoys <- as.numeric(format(thetimes,format='%j'))
    doyvectlist[[i]] <- datadoys

  }
  
  
  close_crmp_connections()
  stnnames <- distance_ranked$history_id[1:nstns]
  names(test_stn_objs) <- c(stnnames)
  names(test_stn_stats) <- c(stnnames)
  names(doyvectlist) <- c(stnnames)
  
  
  crmp_obj <- targetstnobj
  #calculates stats and doy list for target station
  targetstationstats <- get_doy_por_stats(crmp_obj)
  thetimes <- slot(crmp_obj,'time')
  targetdatadoys <- as.numeric(format(thetimes,format='%j'))
  
  
  
  
  pptloc <- which.is.precip(crmp_obj)
  ppt <- crmp_obj[2:nrow(crmp_obj),pptloc] 
  pptflags <- rep(FALSE,nrow(crmp_obj))
  #list of vars ids that indicate daily precip totals
  ppt_ids <- c("519", "496", "429", "513", "501", "514", "526", "527", "547", "441", "430", "495", "471", "456", "465", "466", "546", "631", "639")
  
  #vars_ppt<- read.csv('C:/Users/cdmb/Downloads/data-1555603111544.csv')             
  
  recordstart <- format(as.Date(ppt@time[2]),'%Y')
  recordend <- format(as.Date(ppt@time[length(ppt@time)]),'%Y')
  
  
  
  for(i in 2:(nrow(ppt)-1)){
    if(!is.na(ppt[i])){
      dayzero <- ppt@time[i]
      doyzero <- as.numeric(format(dayzero,format='%j'))
      
      #dayminone <- ppt@time[i-1]
      dayminone <- as.Date(dayzero) -1 
      doyminone <- as.numeric(format(dayminone,format='%j'))
     
     # dayplusone <- ppt@time[i+1]
      dayplusone <- as.Date(dayzero) +1 
      doyplusone <- as.numeric(format(dayplusone,format='%j'))
      
      
      
      #dates for 29 day window
      windowstart <- format(as.Date(ppt@time[i]) - 14,'%m-%d')
      windowstartdow <- as.numeric(format(as.Date(ppt@time[i]) - 14,format='%j')) 
      windowend <- format(as.Date(ppt@time[i])+14,'%m-%d')
      windowenddow <- as.numeric(format(as.Date(ppt@time[i])+14,format='%j')) 
      
      
      #this grabs the data within the 3 day window at the nearby stations
      hilist <- distance_ranked$history_id[1:nstns]
      hilist <- hilist[!is.na(hilist)]
      stndat <- list()
      for (hid in 1:length(hilist)){
        
        testpptloc <- which.is.precip(test_stn_objs[[hid]])
        if(any(testpptloc)){
          testppt <- test_stn_objs[[hid]][2:nrow(test_stn_objs[[hid]]),testpptloc]
          testtime <- test_stn_objs[[hid]]@time[2:nrow(test_stn_objs[[hid]])]
          
          temp <- c(test_stn_objs[[hid]][[testpptloc]][test_stn_objs[[hid]]@time == dayminone],
                         test_stn_objs[[hid]][[testpptloc]][test_stn_objs[[hid]]@time == dayzero],
                         test_stn_objs[[hid]][[testpptloc]][test_stn_objs[[hid]]@time == dayplusone])
          
          if(length(temp) == 0){
            stndat[[hid]] <- list()
          } else{
            temp <- as.list(temp)
            stndat[[hid]] <- temp
          }
        }else{
          temp <- as.list(temp)
          stndat[[hid]] <- temp
        }
        
      }

      names(stndat) <- c(hilist)
      
      #list of nearby stations with data for desired period  
      nearbystations <- Filter(length,stndat)
      
      #it looks like stations need data for all days, so omit stations that have length < 3
      nearbystations <- nearbystations[which(lapply(nearbystations, function(x) length(x) >= 3) == TRUE)]
      nearbystations <- nearbystations[which(lapply(nearbystations,function(x) any(is.na(x)))==FALSE)]
      
      #extract stations with pptN >20
      test_stn_stats2 <- test_stn_stats[which(lapply(test_stn_stats,function(x) x[doyzero,'pptN'] >= 20) == TRUE)]
      
      nearbystations <- nearbystations[names(nearbystations) %in% names(test_stn_stats2)]
      
      #Grabs 3-7 closest stations, returns a list or skips test due to lack of appropriate stations
      if (length(nearbystations) >= 3){
        if (length(nearbystations) > 7){
          teststns <- nearbystations[1:7]
          #testwindow <- nearbystationwindow[1:7]
        } else {
          teststns <- nearbystations
          #testwindow <- nearbystationwindow
        }
      }
      
      #get the minimum absolute target neighbour difference
      if(length(nearbystations) >=3){
        #step 1 - minimum absolute target neighbour difference
        #nearbyprecip <- lapply(teststns, '[[', 'datum')
        nearbyprecip <- unlist(teststns)
        diffs <- ppt[i] - nearbyprecip
        if(ppt[i] > max(nearbyprecip) | ppt[i] < min(nearbyprecip)){
          mintargetdiff <- abs(min(diffs))
        } else{
          mintargetdiff <- 0
        }
        
        #calculate climatological percent ranks for the test stations at the target day
        teststnclimpercentrank <- vector()
        for (station in 1:length(teststns)){
          dayval <- unlist(teststns[[station]])
          meanval <- c(test_stn_stats2[[names(teststns[station])]][doyminone,'ppt_bimean'],
                       test_stn_stats2[[names(teststns[station])]][doyzero,'ppt_bimean'],
                       test_stn_stats2[[names(teststns[station])]][doyplusone,'ppt_bimean'])
          teststnclimpercentrank <- c(teststnclimpercentrank,(dayval * 100)/meanval)
        }
        
        #calculate climatological percent ranks for the target station at the target day
        dayval <- ppt@.Data[[1]][i]
        meanval <- targetstationstats[doyzero,'ppt_bimean']
        targetclimpercentrank <- (dayval * 100)/meanval
        
        mintargetdiff_prank <- min(abs(targetclimpercentrank - teststnclimpercentrank))
        test_threshold <- -45.72*log10(mintargetdiff_prank) + 269.24
        
        if(mintargetdiff > test_threshold){
          pptflags[i] <- TRUE
        }
        
      } else{
        #part 4 - test threshold
        nearbyprecip <- lapply(teststns, '[[', 'datum')
        nearbyprecip <- unlist(nearbyprecip)
        diffs <- ppt[i] - nearbyprecip
        if(ppt[i] > max(nearbyprecip) | ppt[i] < min(nearbyprecip)){
          mintargetdiff <- abs(min(diffs))
        } else{
          mintargetdiff <- 0
        }
        testthreshold <- 269.24
        if(mintargetdiff > test_threshold){
          pptflags[i] <- TRUE
        }
      }
    }
  }
  return(pptflags)
  
}


#So, now have the stats needed to look at climatological quality control. Need to next implement several procedures.
#section 4 table 2
climo_outlier_flag <- function(crmp_obj,statsdf) {
  #looking for tmax greater than
  tnloc <- which.is.tmin(crmp_obj)
  txloc <- which.is.tmax(crmp_obj)
  pptloc <- which.is.precip(crmp_obj)
  thetimes <- slot(crmp_obj,'time')
  datadoys <- as.numeric(format(thetimes,format='%j'))
  tnflags <- rep(FALSE,nrow(crmp_obj))
  txflags <- rep(FALSE,nrow(crmp_obj))
  pptflags <- rep(FALSE,nrow(crmp_obj))
  for (adoy in seq(1,366)) {
    doyidxs <- which(datadoys == adoy)
    crmp_obj_sub <- crmp_obj[doyidxs,]
    tndata <- crmp_obj_sub[,tnloc]
    tnzscore <- (tndata-statsdf$tn_bimean[adoy])/statsdf$tn_bimstd[adoy]
    tndoyflags <- abs(tnzscore) >= 6.0 & rep(statsdf$tnN[adoy],length(tnzscore)) >= 100
    tnflags[doyidxs] <- tndoyflags
    txdata <- crmp_obj_sub[,txloc]
    txzscore <- (txdata-statsdf$tx_bimean[adoy])/statsdf$tx_bimstd[adoy]
    txdoyflags <- abs(txzscore) >= 6.0 & rep(statsdf$txN[adoy],length(txzscore)) >= 100
    txflags[doyidxs] <- txdoyflags
    pptdata <- crmp_obj_sub[,pptloc]
    pptdoyflags <- rep(FALSE,length(pptdata))
    pptthresh <- (txdata+tndata)*.05 > 0
    pptthresh[is.na(pptthresh)] <- FALSE
    pptdoyflags[pptthresh] <- pptdata[pptthresh] >= 9*statsdf$ppt_95p[adoy] & rep(statsdf$pptN[adoy],length(pptdata[pptthresh])) >= 20
    pptdoyflags[!pptthresh] <- pptdata[!pptthresh] >= 5*statsdf$ppt_95p[adoy] & rep(statsdf$pptN[adoy],length(pptdata[!pptthresh])) >= 20
    pptflags[doyidxs] <- pptdoyflags
  }
  txflags[is.na(txflags)] <- FALSE
  tnflags[is.na(tnflags)] <- FALSE
  pptflags[is.na(pptflags)] <- FALSE
  return(list(tmaxflags=txflags,tminflags=tnflags,pptflags=pptflags))
}

#section 5 table 3
temperature_spike_flag <- function(crmp_obj) {
  tnloc <- which.is.tmin(crmp_obj)
  tnflags <- abs(diff(crmp_obj[2:nrow(crmp_obj),tnloc])) >= 25 & as.numeric(diff.POSIXt(crmp_obj@time[2:nrow(crmp_obj)])) == 24 & abs(diff(crmp_obj[1:(nrow(crmp_obj)-1),tnloc])) >= 25 & as.numeric(diff.POSIXt(crmp_obj@time[1:(nrow(crmp_obj)-1)])) == 24
  txloc <- which.is.tmax(crmp_obj)
  txflags <- abs(diff(crmp_obj[2:nrow(crmp_obj),txloc])) >= 25 & as.numeric(diff.POSIXt(crmp_obj@time[2:nrow(crmp_obj)])) == 24 & abs(diff(crmp_obj[1:(nrow(crmp_obj)-1),txloc])) >= 25 & as.numeric(diff.POSIXt(crmp_obj@time[1:(nrow(crmp_obj)-1)])) == 24
  txflags[is.na(txflags)] <- FALSE
  tnflags[is.na(tnflags)] <- FALSE
  return(list(tmaxflags=txflags,tminflags=tnflags))
}

#section 5 table 3 
internal_consistency_flag <- function(crmp_obj) {
  #tmax(0) < tmin(0) - 1
  #tmax(0) < tmin(1) - 1
  #tmin(0) > tmax(1) + 1
  txloc <- which.is.tmax(crmp_obj)
  tnloc <- which.is.tmin(crmp_obj)
  tndata <- crmp_obj[,tnloc]
  txdata <- crmp_obj[,txloc]
  flagstat <- TRUE
  while (flagstat) {
    tnflags1 <- txdata < tndata - 1 #have to flag the current tmin and tmax
    tnflags1[is.na(tnflags1)] <- FALSE
    txflags1 <- tnflags1
    txflags2 <- (txdata[1:(length(txdata)-1)] < tndata[2:length(tndata)] - 1) & as.numeric(diff.POSIXt(crmp_obj@time)) == 24
    txflags2[is.na(txflags2)] <- FALSE
    tnflags2 <- rep(FALSE,length(txflags2))
    tnflags2[2:length(txflags2)] <- txflags2[1:(length(txflags2)-1)]
    tnflags3 <- (tndata[1:(length(tndata)-1)] > txdata[2:length(txdata)] + 1) & as.numeric(diff.POSIXt(crmp_obj@time)) == 24
    tnflags3[is.na(tnflags3)] <- FALSE
    txflags3 <- rep(FALSE,length(tnflags3))
    txflags3[2:length(tnflags3)] <- tnflags3[1:(length(tnflags3)-1)]
    tn1count <- as.numeric(tnflags1)
    tn2count <- as.numeric(c(tnflags2,FALSE))
    tn3count <- as.numeric(c(tnflags3,FALSE))
    tx1count <- as.numeric(txflags1)
    tx2count <- as.numeric(c(txflags2,FALSE))
    tx3count <- as.numeric(c(txflags3,FALSE))
    tncount <- tn1count + tn2count + tn3count
    txcount <- tx1count + tx2count + tx3count
    if (any(tncount > 1 | txcount > 1)) {
      browser()
    } else {
      flagstat=FALSE
      break
    }
  }
  return(list(tmaxflags=as.logical(txcount),tminflags=as.logical(tncount)))
}


gianttime <- make_giant_daily_timescale()
giantdoytime <- format(gianttime,'%Y-%j')

#section 5 table 3
lagged_temp_range_flag <- function(crmp_obj) {
  #tmax(0) >= max(tmin(-1:1))+40
  #tmin(0) <= min(tmax(-1:1))-40
  txloc <- which.is.tmax(crmp_obj)
  tnloc <- which.is.tmin(crmp_obj)
  txdata <- crmp_obj[,txloc]
  tndata <- crmp_obj[,tnloc]
  txflags <- rep(FALSE,nrow(crmp_obj))
  tnflags <- rep(FALSE,nrow(crmp_obj))
  #Well, if we throw the data onto a common timeline with NA infilling the actual gaps, then taking max and min over windows of appropriate lags becomes arbitrary as the NAs can either be excluded
  obsdoytime <- format(crmp_obj@time,'%Y-%j')
  gianttx <- rep(NA,length(giantdoytime))
  gianttn <- rep(NA,length(giantdoytime))
  giantdf <- data.frame(time=giantdoytime,tx=gianttx,tn=gianttn)
  #Drop the replicated time rows, then append the rows with data onto the bottom and then resort the time column.
  giantdf <- giantdf[which(!(giantdf$time %in% obsdoytime)),]
  datadf <- data.frame(time=obsdoytime,tx=txdata,tn=tndata)
  giantdf <- rbind(giantdf,datadf)
  giantdf <- giantdf[order(giantdf$time),]
  #trim down the giantdf to reduce later computation over nothingness
  minidx <- which(giantdf$time == obsdoytime[1])
  maxidx <- which(giantdf$time == obsdoytime[length(obsdoytime)])
  giantdf <- giantdf[minidx:maxidx,]
  tdf <- giantdf[,c('time','tn','tx')]
  nrows <- nrow(tdf)
  tdf$tnp <- rep(NA,nrows)
  tdf$tnp[1:(nrows-1)] <- tdf$tn[2:nrows]
  tdf$tnm <- rep(NA,nrows)
  tdf$tnm[2:nrows] <- tdf$tn[1:(nrows-1)]
  tdf$txp <- rep(NA,nrows)
  tdf$txp[1:(nrows-1)] <- tdf$tx[2:nrows]
  tdf$txm <- rep(NA,nrows)
  tdf$txm[2:nrows] <- tdf$tx[1:(nrows-1)]
  #all set up, now just calculate the max and mins and compare
  tdf$tx3n <- apply(tdf[,c('tx','txp','txm')],1,min,na.rm=T)-40
  tdf$tn3x <- apply(tdf[,c('tn','tnp','tnm')],1,max,na.rm=T)+40
  #Alright, have everything, make some flags!
  tdf$tnflags <- tdf$tn <= tdf$tx3n
  tdf$tnflags[which(is.infinite(tdf$tnflags)| is.na(tdf$tnflags))] <- FALSE
  tdf$tnflags[which(is.infinite(tdf$tx3n))] <- FALSE
  tdf$txflags <- tdf$tx >= tdf$tn3x
  tdf$txflags[which(is.infinite(tdf$txflags)| is.na(tdf$txflags))] <- FALSE
  tdf$txflags[which(is.infinite(tdf$tn3x))] <- FALSE
  #remap the flags to the original length of the crmp object and etc.
  flagsdf <- data.frame(time=obsdoytime,txflags=txflags,tnflags=tnflags)
  tdf <- tdf[which(tdf$time %in% obsdoytime),]
  flagsdf <- flagsdf[which(!(flagsdf$time %in% tdf$time)),]
  flagsdf <- rbind(flagsdf,tdf[,c('time','txflags','tnflags')])
  return(list(tmaxflags=flagsdf$txflags,tminflags=flagsdf$tnflags))
}

#table 2
get_onevar_dist_flag <- function(datavec_raw,twotail=TRUE,gapthresh=10) {
  flags <- rep(FALSE,length(datavec_raw))
  datavec <- datavec_raw[which(!is.na(datavec_raw))]
  datavec <- datavec[order(datavec)]
  datavec_high <- datavec[datavec >= median(datavec,na.rm=T)]
  diffdata <- diff(datavec_high)
  maxdif_high <- max(diffdata,na.rm=T)
  maxdifval <- datavec_high[(which(diffdata == maxdif_high)+1)]
  datavec_low <- datavec[which(datavec < median(datavec,na.rm=TRUE))]
  diffdata <- diff(datavec_low)
  maxdif_low <- max(diffdata,na.rm=T)
  mindifval <- datavec_low[(which(diffdata == maxdif_low))]
  if (maxdif_high >= gapthresh) {
    flags[which(datavec_raw >= maxdifval)] <- TRUE
  }
  if (twotail & maxdif_low >= gapthresh) {
    flags[which(datavec_raw <= mindifval)] <- TRUE
  }
  return(flags)
}



distribution_flag <- function(crmp_obj) {
  tdata <- crmp_obj[,which.is.tmax(crmp_obj)]
  tmaxflags <- get_onevar_dist_flag(tdata,TRUE,10)
  tdata <- crmp_obj[,which.is.tmin(crmp_obj)]
  tminflags <- get_onevar_dist_flag(tdata,TRUE,10)
  tdata <- crmp_obj[,which.is.precip(crmp_obj)]
  pptflags <- get_onevar_dist_flag(tdata,FALSE,300)
  return(list(tmaxflags=tmaxflags,tminflags=tminflags,pptflags=pptflags))
}

get_all_month_dist_mega_flags <- function(crmp_obj) {
  for (amonth in seq(1,12)) {
    print(sprintf('%02d',amonth))
    month_crmp_obj <- get_month_subset(crmp_obj,amonth)
    monthDistFlagsList <- distribution_flag(month_crmp_obj)
    monthMegaconsistFlagsList <- megaconsist_flag(month_crmp_obj)
    if (any(unlist(lapply(monthDistFlagsList,any)))) {
      #have some flags that we need to worry about
      print('month dist flags found. need to process.')
      browser()
    }
    if (any(unlist(lapply(monthMegaconsistFlagsList,any)))) {
      #have some flags that we need to worry about
      print('month megaconsistency flags found. need to process.')
      browser()
    }
  }
}

#table 5
megaconsist_flag <- function(crmp_obj) {
  tndata <- crmp_obj[,which.is.tmin(crmp_obj)]
  tnflags <- rep(FALSE,length(tndata))
  txdata <- crmp_obj[,which.is.tmax(crmp_obj)]
  txflags <- rep(FALSE,length(txdata))
  txflags <- txdata < min(tndata,na.rm=T)
  txflags[is.na(txflags)] <- FALSE
  tnflags <- tndata > max(txdata,na.rm=T)
  tnflags[is.na(tnflags)] <- FALSE
  if (any(txflags)) {print('tx megaconsistency fail')}
  if (any(tnflags)) {print('tn megaconsistency fail')}
  return(list(tmaxflags=txflags,tminflags=tnflags))
}

get_histids_withdata_by_month_year <- function() {
  yearmonthcovlist <- list()
  for (year in seq(1880,2017)) {
    for (month in seq(1,12)) {
      con <- dbConnect(PostgreSQL(),user="fanslow",host="monsoon.pcic.uvic.ca",dbname="crmp")
      yearmonthstr <- paste(as.character(year),'-',sprintf('%02d',month),sep='')
      minstr <- paste(yearmonthstr,'-01',sep='')
      maxstr <- paste(yearmonthstr,'-',sprintf('%02d',get_mdays(month,year)),sep='')
      if (year == 1900 & month == 2) {
        maxstr <- paste(yearmonthstr,'-',sprintf('%02d',(get_mdays(month,year)-1)),sep='')
      }
      querystr <- paste("SELECT history_id from obs_count_per_month_history_mv WHERE date_trunc >= '",minstr,"' AND date_trunc <= '",maxstr,"'	;",sep='')
      try.rv <- try({
        rs <- dbSendQuery(con,querystr)
        monthyearcov <- fetch(rs,n=-1)
      }, silent=T)
      if (inherits(try.rv, 'try-error')) {
        close_crmp_connections()
        print(try.rv)
      }
      close_crmp_connections()
      yearmonthcovlist[[yearmonthstr]] <- monthyearcov
    }
    print(paste('done with year: ',as.character(year),sep=''))
  }
  return(yearmonthcovlist)
}

load('/home/fanslow/Work/code/crmp_fanslow/distmat.RData')
load('/home/fanslow/Work/code/crmp_fanslow/monthcovlist.RData')
library('lubridate')

kill_allna <- function(df) {
  allisna <- function(x) {
    all(is.na(x))
  }
  keeprows <-  !apply(df[,seq(2,ncol(df))],1,allisna)
  keepcols <- !apply(df[,seq(2,ncol(df))],2,allisna)
  df <- df[which(keeprows),(c(1,(which(keepcols)+1)))]
  return(df)
}

get_cov_pct <- function(df) {
  count_datarows <- function(x) {
    return(sum(!is.na(x)))
  }
  return(100*(apply(df,2,count_datarows))/nrow(df))
}

drop_insufficient_overlap_count <- function(refdf,targtime,targdata) {
  refdf_orig <- refdf
  overlap_one <- function(refdata,targdata) {
    return(sum(!is.na(refdata) & !is.na(targdata)))
  }
  targdf <- data.frame(date=as.POSIXlt(targtime),datum=targdata)
  refdf <- refdf[which(refdf$date %in% as.POSIXlt(targdf$date)),]
  targdf <- targdf[which(as.POSIXlt(targdf$date) %in% refdf$date),]
  #With the data subsetted. should be able to compare the data themselves
  counts <- apply(refdf[,2:(ncol(refdf))],2,overlap_one,targdata=targdf$datum)
  return(list(df=refdf_orig[,c(1,(which(counts >= 40)+1))],df_short=refdf[,c(1,(which(counts >= 40)+1))]))
}

get_ioa_lmstats <- function(df,targdata) {
  get_lm_slopeint <- function(x,y) {
    lmfit <- lm(y~x)
    return(c(lmfit$coefficients[2],lmfit$coefficients[1]))
  }
  ioas <- apply(df[,2:ncol(df)],2,idxOfAgreement,y=targdata)
  slopeints <- apply(df[,2:ncol(df)],2,get_lm_slopeint,y=targdata)
  return(rbind(ioas,slopeints))
}

jitter_the_dater <- function(targ,refdf) {
  jitter_one <- function(ref,targ) {
    ref <- cbind(ref,c(ref[2:(length(ref))],NA),c(NA,ref[1:(length(ref)-1)]))
    refdif <- abs(targ-ref)
    mindifs <- matrix(apply(refdif,1,min,na.rm=T),ncol=1)
    keeps <- cbind(mindifs == refdif[,1],mindifs == refdif[,2],mindifs == refdif[,3])
    keeps[!keeps] <- NA
    return(apply(matrix(ref[keeps],ncol=3),1,mean,na.rm=T))
  }
  return(apply(refdf,2,jitter_one,targ=targ))
}

crmp_obj_path <- '/home/fanslow/scratch/crmp_climdex/'
buddy_check_one_stn_txtn <- function(history_id,radius=75000,minNeighbours=3,maxNeighbours=7) {
  station_id <- stnmeta$station_id[which(stnmeta$history_id == history_id)]
  try.rv <- try(load(paste(crmp_obj_path,'sid_',as.character(station_id),'_hid_',as.character(history_id),'_8110_crmp_climdex_objs.RData',sep='')),silent=T); if(inherits(try.rv,'try-error')) {return(NULL)}
  tnloc_targ <- which.is.tmin(crmp_obj)
  txloc_targ <- which.is.tmax(crmp_obj)
  if (length(txloc_targ) == 0 | length(tnloc_targ)==0) {browser(); return(NA)}
  if (all(is.na(crmp_obj[,tnloc_targ])) | all(is.na(crmp_obj[,txloc_targ]))) {return(NULL)}
  crmp_obj <- crmp_obj[which(!duplicated(format(crmp_obj@time,'%Y-%m-%d'))),]
  targ_obj <- crmp_obj
  dist_histids <- row.names(distmat)
  distances <-data.frame(history_id=dist_histids,distance=distmat[which(dist_histids==as.character(history_id)),],stringsAsFactors=F)
  distances <- distances[order(distances$distance,decreasing=F),]
  distances <- distances[which(distances$distance <= radius),]
  #get the stations that have temporal coverage during the POR of the target station
  #We want to refine this further by seeking the number of stations in the first third and last third of the record.
  #Make sure that in the first third of the record, there are lots of stations with data available.
  #yearmo_strs <- unique(format(targ_obj@time[c((1:as.integer(length(targ_obj@time)/3)),(max((length(targ_obj@time)-100),1):length(targ_obj@time)))],'%Y-%m'))
  yearmo_strs <- unique(format(targ_obj@time,'%Y-%m'))
  ref_histids <- yearmonthcovlist[yearmo_strs]
  ref_histids <- unique(unlist(ref_histids))
  distances <- distances[which(distances$history_id %in% ref_histids),]
  ndist <- nrow(distances)-1
  #now the distances dataframe has everything that has SOME data within the range of the target series
  #load in all the data in the distances dataframe and assemble into a giant data frame...?
  nrefstns <- nrow(distances)
  if (nrefstns < minNeighbours) {print(paste('Too few stations for QC on entire station: ',targ_obj@station.name,sep='')); return(NULL)}
  timescale <- make_giant_daily_timescale()
  tndf <- as.data.frame(matrix(data=NA,nrow=length(timescale),ncol=(nrefstns+1)))
  names(tndf) <- c('date',as.character(distances$history_id))
  tndf$date <- format(timescale,'%Y-%m-%d')
  txdf <- tndf
  print(paste('loading ',as.character(nrefstns),' reference stations for buddy checking station: ',targ_obj@station.name,sep=''))
  for (stnidx in seq(1,nrefstns)) {
    #print(100*((stnidx-1)/nrefstns))
    #load in the data from the flat file
    hist_id <- distances$history_id[stnidx]
    stn_id <- stnmeta$station_id[which(stnmeta$history_id == hist_id)]
    colidx <- which(names(tndf) == as.character(hist_id))

    try.rv <- try(load(paste(crmp_obj_path,'sid_',as.character(stn_id),'_hid_',as.character(hist_id),'_8110_crmp_climdex_objs.RData',sep='')),silent=T); if(inherits(try.rv,'try-error')) {next}

    crmp_obj <- crmp_obj[which(!duplicated(format(crmp_obj@time,'%Y-%m-%d'))),]
    tnloc <- which.is.tmin(crmp_obj)
    txloc <- which.is.tmax(crmp_obj)
    indf <- data.frame(date=format(crmp_obj@time,'%Y-%m-%d'),tn=crmp_obj[,tnloc],tx=crmp_obj[,txloc])
    #drop dates that are not in the master data frame. Some crmp objects have intermediate date data such as those from the ARDA network.
    rowidxs <- which(tndf$date %in% indf$date)
    if (length(rowidxs) == nrow(indf)) {
      try.rv <- try({
        if (length(tnloc) == 0) {
          indf$tn <- rep(NA,length(rowidxs))
        }
        if (length(txloc) == 0) {
          indf$tx <- rep(NA,length(rowidxs))
        }
        tndf[rowidxs,colidx] <- indf$tn
        txdf[rowidxs,colidx] <- indf$tx
      })
      if (inherits(try.rv,'try-error')) {browser()}
    } else {
      print('date alignment mismatch'); browser()
    }
  }
  #trim off the excess by looking for all row-wise NA and drop those rows and cols with all NA as well.
  #This has dangerous consequences as may delete rows in the midst of data coverage which should otherwise be there. Better approach may be to keep track of the
  #max and min dates as the data frame is being assembled and only trim off the ends of the data frame that are older or younger. Still want to pull out the all na columns.
  #This is further coomplicated by continued removal of all NA rows in subsequent code...
  browser()
  tndf <- kill_allna(tndf)
  tndf$date <- strptime(tndf$date,format='%Y-%m-%d')
  txdf <- kill_allna(txdf)
  txdf$date <- strptime(txdf$date,format='%Y-%m-%d')
  #Data are all loaded and prepped. Time to run through the target station record on a window by window basis and do the index of agreement calculations and etc.
  for (yearmo in yearmo_strs) {
    yearmoidxs <- which(yearmo_strs == yearmo)
    windowlist <- get_regression_datewindow(yearmo)

    #pull out the data in the window and winnow by completeness
    targ_sub <- targ_obj[which(targ_obj@time >= windowlist$windowstart & targ_obj@time <= windowlist$windowend),]
    if (all(is.na(targ_sub[,tnloc_targ])) & all(is.na(targ_sub[,txloc_targ]))) {print(paste('No data in target station to QC for ',yearmo,sep='')); next}
    tndf_sub <- tndf[which(tndf$date >= windowlist$windowstart & tndf$date <= windowlist$windowend),]
    tndf_sub <- kill_allna(tndf_sub)
    if (length(tndf_sub) == 0 | is.null(ncol(tndf_sub))) {print(paste('Not enough overlap in data for QC ',yearmo,sep='')); next}
    try.rv <- try(if (ncol(tndf_sub) <= 2) {print(paste('Not enough overlap in data for QC ',yearmo,sep='')); next}); if (inherits(try.rv,'try-error')) {browser()}
    tndf_sub <- tndf_sub[,which(get_cov_pct(tndf_sub) >= 69)]
    txdf_sub <- txdf[which(txdf$date >= windowlist$windowstart & txdf$date <= windowlist$windowend),]
    txdf_sub <- kill_allna(txdf_sub)
    if (length(txdf_sub) == 0) {print(paste('Not enough overlap in data for QC ',yearmo,sep='')); next}
    if (ncol(txdf_sub) <= 2) {print(paste('Not enough overlap in data for QC ',yearmo,sep='')); next}
    txdf_sub <- txdf_sub[,which(get_cov_pct(txdf_sub) >= 69)]
    #Finally, have to make sure that there are at least 40 common days between the
    if (length(tndf_sub) == 0 | is.null(ncol(tndf_sub))) {print(paste('Not enough overlap in data for QC ',yearmo,sep='')); next}
    try.rv <- try(if (ncol(tndf_sub) <= 2) {print(paste('Not enough overlap in data for QC ',yearmo,sep='')); next}); if (inherits(try.rv,'try-error')) {browser()}
    if (ncol(txdf_sub) <= 2) {print(paste('Not enough overlap in data for QC ',yearmo,sep='')); next}

    try.rv <- try(dflist <- drop_insufficient_overlap_count(tndf_sub,targ_sub@time,targ_sub[,which.is.tmin(targ_sub)])); if(inherits(try.rv,'try-error')) {browser()}
    if (is.null(ncol(dflist$df))) {print(paste('Not enough overlap in data for QC ',yearmo,sep='')); next}
    tndf_sub <- dflist['df'][[1]]
    tndf_sub_short <- dflist['df_short'][[1]]
    nstns <- ncol(tndf_sub)
    nstnsuse <- min((nstns-2),maxNeighbours)
    if (nstnsuse < minNeighbours) {print(paste('Not enough reference data for QC ',yearmo,sep='')); next}
    nobs <- nrow(tndf_sub)
    #nobs_short <- nrow(tndf_sub_short)
    if (yearmo == '1998-10') {browser()}
    try.rv <- try(dflist <- drop_insufficient_overlap_count(txdf_sub,targ_sub@time,targ_sub[,which.is.tmax(targ_sub)])); if(inherits(try.rv,'try-error')) {browser()}
    txdf_sub <- dflist['df'][[1]]
    txdf_sub_short <- dflist['df_short'][[1]]
    tnstats <- t(get_ioa_lmstats(tndf_sub_short,targ_sub[,which.is.tmin(targ_sub)]))
    thecolnames <- colnames(tnstats); therownames <- rownames(tnstats); tnstats <- data.frame(tnstats)
    names(tnstats) <- thecolnames; rownames(tnstats) <- therownames; tnstats$history_id <- therownames
    tndf_sub <- tndf_sub[,c(1,(1+order(tnstats$ioas,decreasing=T)))]
    #REORDER THE DATAFRAME!!!!!!!!!!!!!!!!!!!
    tnstats <- tnstats[order(tnstats$ioas,decreasing=T),]
    #Somewhere in here have to introduce the jitter to enable the comparison between the target data and reference data among the centred 3-day window. Looking for the minimum absolute difference
    jittered <- jitter_the_dater(tndf_sub[,2],as.matrix(tndf_sub[,3:nstns]))
    pred_raw <- t(as.matrix(tnstats[['(Intercept)']][2:nstns]))[rep(1,nobs),seq(1,(nstns-2))]+t(as.matrix(tnstats$x[2:nstns]))[rep(1,nobs),seq(1,(nstns-2))] * as.matrix(jittered)

    ioasmat <- t(as.matrix(tnstats$ioas[2:(nstns-1)]))[rep(1,nobs),seq(1,(nstns-2))]
    namask <- is.na(pred_raw)
    ioasmat[namask] <- NA
    keepidxs <- which(format(tndf_sub$date,'%Y-%m') == yearmo)
    predictions <- apply((pred_raw[keepidxs,1:nstnsuse] * ioasmat[keepidxs,1:nstnsuse]),1,sum,na.rm=T) / apply((ioasmat[keepidxs,1:nstnsuse]),1,sum,na.rm=T)
    try.rv <- try(thecor <- cor(predictions,tndf_sub[keepidxs,2],use='complete.obs')); if (inherits(try.rv,'try-error')) {browser()}
    if (thecor >= 0.8) {
      #have a good enough fit to try to QC the data!
      resids <- predictions-tndf_sub[keepidxs,2]
      stdresids <- (resids-mean(resids,na.rm=T))/sd(resids,na.rm=T)
      if (max(abs(resids),na.rm=T) >= 8 & max(abs(stdresids),na.rm=T) >= 4) {print('found one!')}
      print(paste(as.character(history_id),' regression checked with ',as.character(nstnsuse),' stns for: ',yearmo,'. Predicted series correlation: ',sprintf('%4.2f',thecor),' with max res. ',sprintf('%3.1f',max(abs(resids),na.rm=T)),' and stdize resid. ',sprintf('%3.1f',max(abs(stdresids),na.rm=T)),sep=''))

    }
    #print(paste('done with ',yearmo,sep=''))
  }
  #browser()
}
testrun_all <- function(stnmeta,radius=75000,minNeighbours=3,maxNeighbours=7) {
  for (histid in stnmeta$history_id) {
    print(histid)
    buddy_check_one_stn_txtn(histid,radius=radius,minNeighbours=minNeighbours,maxNeighbours=maxNeighbours)
  }
}

get_regression_datewindow <- function(yearmo) {
  year <- as.numeric(strsplit(yearmo,'-')[[1]][1])
  month <- as.numeric(strsplit(yearmo,'-')[[1]][2])
  yearmonth <- paste(sprintf('%04d',year),'-',sprintf('%02d',month),sep='')
  minyearmonth <- paste(sprintf('%04d',year),'-',sprintf('%02d',(month-1)),sep='')
  maxyearmonth <- paste(sprintf('%04d',year),'-',sprintf('%02d',(month+1)),sep='')
  if (month == 1) {
    minyearmonth <- paste(sprintf('%04d',(year-1)),'-',sprintf('%02d',12),sep='')
  } else if (month == 12) {
    maxyearmonth <- paste(sprintf('%04d',year+1),'-',sprintf('%02d',1),sep='')
  }
  yearmonthstart <- strptime(paste(yearmonth,'-01',sep=''),format='%Y-%m-%d')
  yearmonthend <- strptime(paste(yearmonth,'-',sprintf('%02d',get_mdays(month,year)),sep=''),format='%Y-%m-%d')
  windowstart <- yearmonthstart - ddays(15)
  windowstart <- strptime(windowstart,format='%Y-%m-%d')
  windowend <- yearmonthend + ddays(15)
  windowend <- strptime(windowend,format='%Y-%m-%d')
  ndayswindow <- as.numeric(windowend-windowstart)+1
  return(list(yearmonth=yearmonth,minyearmonth=minyearmonth,maxyearmonth=maxyearmonth,windowstart=windowstart,windowend=windowend,ndayswindow=ndayswindow))
}


order.crmp <- function(obj) {
  stopifnot(inherits(obj, 'crmp'))
  t <- slot(obj, 'time')
  i <- order(t)
  obj[i,,drop=F]
}

get_month_subset <- function(crmp_obj,month) {
  monthstr <- sprintf('%02d',month)
  months <- format(crmp_obj@time,'%m')
  return(crmp_obj[which(months == monthstr),])
}


process_auto_baserange <- function(crmp_obj) {
  lows <- c(1961,1971,1981)
  highs <- lows+29
  minyear <- min(as.numeric(format(slot(crmp_obj,'time'),'%Y')))
}
