#devtools::install_github('buzacott/bomWater')
#devtools::install_github('WillemMaetens/standaRdized')

library(bomWater)
library(data.table)
library(dplyr)
library(forecast)
library(fst)
library(geodist)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(lubridate)
library(magrittr)
library(raster)
library(readxl)
library(rnaturalearth)
library(rprojroot)
library(standaRdized)
library(sf)
library(SPEI)
library(stringr)

#Define directory structure
rootdir <- find_root(has_dir("src"))
resdir <- file.path(rootdir, "results")
figdir <- file.path(resdir, 'figures')
datdir <- file.path(rootdir, "data")

#To get rainfall data (get station number first)
#bomrang #https://docs.ropensci.org/bomrang/articles/bomrang.html
#bomrang::get_historical

source('SamZipper.R')

#------------------- Functions ------------------------------------------------
#------ zero_lomf -----------------
#' Last \[non-zero\] Observation Moved Forward (lomf)
#'
#' Finds the index, for each row, of the previous row with a non-zero value
#'
#' @param x Numeric vector.
#' @param first (logical) Whether to consider first value as a non-zero value
#'   whose index is moved forward even if it is zero. This prevents having NAs
#'   in the results and somewhat assumes that, for a time series, the day prior
#'   to the first value is non-zero.
#'
#' @return Numeric vector of the indices of the previous non-zero for each
#'   element of the input vector.
#'
#' @examples
#' test1 <- c(1,1,1,0,0,0,0,1,1)
#' zero_lomf(test1)
#' test2 <- c(0,0,0,0,0,1,1,0,1)
#' zero_lomf(test2, first=FALSE)
#' zero_lomf(test2, first=TRUE)
#'
#' @export
zero_lomf <- function(x, first=TRUE) {
  if (length(x) > 0) {
    non.zero.idx <- which(x != 0)
    if(first==T & x[1]==0) {
      non.zero.idx=c(1,non.zero.idx)
      #Repeat index of previous row with non-zero as many times gap until next non-zero values
      out_ind <- rep.int(non.zero.idx, diff(c(non.zero.idx, length(x) + 1)))
    } else if (first==F & x[1]==0) {
      out_ind <- c(rep(NA, min(non.zero.idx)-1), 
                   rep.int(non.zero.idx, diff(c(non.zero.idx, length(x) + 1))))
    }
    
    return(out_ind)
  }
}



#------------------- Import and format data ------------------------------------------------
#Download daily data for Ashburton River at Nanutarra
ash_raw <- get_daily(parameter_type = 'Water Course Discharge',
                     station_number = '706003', #'706209'
                     start_date = '1900-01-01',
                     end_date = '2021-06-27') %>%
  as.data.table
min(ash_raw$Timestamp) #Check beginning of time series

#Keep only period included in Van Huijegoort
ash_sub <- ash_raw[
  Timestamp >= as.POSIXct('1973-01-01 00:00:00') & 
    Timestamp < as.POSIXct('2006-01-01 00:00:00'),]

#Plot time series
#ash_raw[is.na(Value), Value:= 9999]
ggplot(ash_sub[Timestamp >= as.POSIXct('1988-01-01 00:00:00') & 
                 Timestamp < as.POSIXct('1993-01-01 00:00:00'),], 
       aes(x=Timestamp, y=Value+0.1)) + 
  geom_line() +
  geom_point(data=ash_sub[Value==9999,], color='red') +
  scale_y_log10() + 
  theme_classic()



ash_preddates <- ash_sub[
  Timestamp >= as.POSIXct('1988-01-01 00:00:00') & 
    Timestamp < as.POSIXct('1993-01-01 00:00:00'),Timestamp]

#-------------- Compute Standardized Streamflow Index (SSI) with standaRdized package--------------------
#Convert streamflow data to xts
ashxts <- as.xts.data.table(ash_sub[, .(Timestamp, Value)])

#As implemented by https://github.com/WillemMaetens/standaRdized
#Compute 6-month SSI
SSI1 <- standardized.index(data= (ashxts/100)+0.01, #Numerical stability problem with algorithm due to large values — https://stackoverflow.com/questions/53557022/error-code-100-fitting-exp-distribution-using-fitdist-in-r
                           index.out = ash_preddates,
                           agg.length = 365,
                           distr = 'gev',
                           agg.fun = "sum",
                           method = 'lmom',
                           agg.na.thres = 10,
                           agg.interpolation = 'linear',
                           output.attrs = 'all'
)

ggplot(as.data.table(SSI1), aes(x=index, y=value)) + 
  geom_line(color='red') + 
  geom_line(data=ash_sub[as.POSIXct(Timestamp) %in% ash_preddates,], 
            aes(x=as.Date(Timestamp), y=Value)) + 
  scale_y_continuous(
    trans=scales::pseudo_log_trans(base = 10))


#-------------- Compute Standardized Streamflow Index (SSI) with SPEI package --------------------




#-------------- Compute Variable threshold level method (Van Huijgevoort) ------------------------
#As implemented by van Huijgevoort
compute_tlm <- function(in_dt, datecol, qcol, prob, na_mothresh, 
                        IDdrought = TRUE, tdrought = c('daily', 'motnhly')) {
  in_dt <- copy(in_dt)
  
  #Create a month and a year-month column for processing
  in_dt[, `:=`(
    month=format(get(datecol), '%m'),
    ym=format(get(datecol), '%Y-%m')
  )]
  
  #Aggregate to monthly time steps
  dtmo <- in_dt[, list(
    Qmo = mean(get(qcol), na.rm=T), #Compute mean monthly discharge
    Nmiss = .SD[is.na(get(qcol)), .N], #Compute monthly number of missing days of data
    date = .SD[.N, as.Date(get(datecol))] #Keep date column, converted to date
  ), by=.(ym, month)]

  #Compute percentiles
  dtperct <- dtmo[, quantmo := .SD[Nmiss < na_mothresh, #Only consider months with number of NA discharge below threshold
                                 quantile(Qmo, prob, na.rm=T)], #Compute quantile
                  by= month]

  if (tdrought == 'daily') {
    dtperct <- in_dt[dtperct[, -'month', with=F], on='ym']
  }
  
  #Identify droughts
  if (IDdrought) {
    if (tdrought == 'monthly') {
      dtperct <- dtperct[Qmo <= quantmo, drought := Qmo] 
    } else if (tdrought == 'daily') {
      dtperct[get(qcol) <= quantmo, drought := get(qcol)] 
    }
    
    dtperct[, drought_id := rleid(is.na(drought))]
  }

  return(dtperct)
}

ashperct <- compute_tlm(
  in_dt = ash_sub,
  datecol = 'Timestamp',
  qcol = 'Value',
  IDdrought = TRUE,
  probs = 0.05,
  na_mothresh = 4,
  tdrought = 'daily')

ggplot(ashperct, aes(x=as.Date(Timestamp), y=quantmo, group=1)) + 
  geom_line(color='blue', size=1.2) + 
  #geom_point(data=ashperct[Value==0,], aes(x=as.Date(Timestamp), y=Value), color='green', size=1.5) +
  geom_line(data=ash_sub, aes(y=Value)) + 
  geom_line(aes(y=drought, group=drought_id), color='red', size=1.2) +
  scale_y_continuous(
    #trans=scales::pseudo_log_trans(base = 10), 
    breaks=c(0,0.1,0.5,1,10,100,1000, 2500)) + 
  scale_x_date(limits=c(as.Date('1988-01-01'), 
                        as.Date('1993-01-01'))) + 
  coord_cartesian(ylim=c(0,0.2)) + 
  theme_classic()

#-------------- Compute mixed threshold level and consecutive dry period method (Van Huijgevoort) ------------------------
compute_tlm_cdpm <- function(in_dt, datecol, qcol, prob, na_mothresh, qperc_window = 1,
                             IDdrought = TRUE, tdrought = c('daily', 'monthly')
) {
  #Compute drought percentile based
  out_dt <- copy(in_dt)
  
  #Compute percentile statistics for each DOY (step 1 in Van Huijgevoort)
  out_dt[, doy := as.integer(format(Timestamp, '%j'))]
  
  ecdf_list <- lapply(1:366, function(j) {
    out_dt[doy %in% seq(j - floor(qperc_window),
                        j + floor(qperc_window)), 
           ecdf(Value)]
  })
  
  out_dt[, q_percentile := mapply(function(d, v) {ecdf_list[[d]](v)}, 
                                  d=doy,
                                  v=Value)]
  
  #Get proportion of zero-flow values
  Fwet <- out_dt[, ecdf(get(qcol))(0)]
  
  # #Get overall empirical cumulative distribution with moving window mean 
  # qecdf <- out_dt[,ecdf(frollmean(get(qcol), n=qperc_window, fill=NA,
  #                                 align="center", na.rm=T))]
  
  if (out_dt[, quantile(get(qcol), 0.05, na.rm=T) > 0]) { #If fifth percentile > 0, return TLM method.
    out_dt[, drought_final := as.numeric(qpercentile <= prob)] #(step 2 in Van Huijgevoort)
    return(out_dt)
  } 
  else { #Otherwise compute Van Huijgoort method
    
    #For each record, compute date of last non-zero flow day (step 3 in Van Huijgevoort)
    out_dt[, prevflowdate := out_dt[zero_lomf(get(qcol), first=T),
                                  datecol, 
                                  with=F]
    ] %>% 
      .[(get(qcol) != 0) | is.na(get(qcol)), prevflowdate:=NA] %>% #If non-zero flow, set prevflowdate to NA
      .[!is.na(prevflowdate), `:=`(
        zero_consec = as.integer(difftime(get(datecol),  #Number of consecutive zero-flow days
                                          prevflowdate,
                                          units='days')
        ),
        dryordrought = 1)] 
    
    #Compute percentile statistics for each daily discharge value (step 4 in Van Huijgevoort)
    out_dt[, qpercentile := qecdf(frollmean(get(qcol), n=qperc_window, fill=NA,
                                            align="center", na.rm=T))]
    #Compute drought for Q > 0 (step 4 in Van Huijgevoort)
    out_dt[(get(qcol) > 0 & qpercentile <= prob & !is.na(get(qcol))), 
           dryordrought := 1] 
    
    #Compute consecutive days of drought with Q >0 | dry days (i.e., Ndry,drought; step 5 in Van Huijgevoort)
    out_dt[, drought_consec := cumsum(dryordrought), by=rleid(dryordrought)] %>%
      .[is.na(drought_consec), drought_consec := 0] %>% 
      .[is.na(zero_consec), zero_consec := 0]
    
    #Compute percentile value for each drought day (in terms of consecutive # of drought days)
    #based on empirical cumulative distribution function of zero-flow days only (step 6 in Van Huijgevoort)
    zero_ecdf <- out_dt[zero_consec > 0, ecdf(zero_consec)]
    out_dt[, drought_percentile := (1-zero_ecdf(drought_consec))*Fwet]
    
    #Combine percentile values for wet and dry periods
    out_dt[, final_percentile := fifelse(get(qcol) > 0, 
                                         qpercentile,
                                         drought_percentile)]
    
    #Determine whether a drought is occurring (step 7 in Van Huijgevoort)
    out_dt[, drought_final := as.numeric(final_percentile <= prob)]
    
  }
  
  return(out_dt)
}

ash_tlmcdpm <- compute_tlm_cdpm(in_dt = ash_sub[Timestamp %in% ash_preddates,],
                                datecol = 'Timestamp',
                                qcol = 'Value',
                                qperc_window = 31,
                                IDdrought = FALSE,
                                prob = 0.2,
                                na_mothresh = 4,
                                tdrought = 'daily')

ggplot(ash_tlmcdpm, aes(x=Timestamp)) + 
  geom_line(aes(y=Value)) +
  geom_point(data=ash_tlmcdpm[drought_final==1,], aes(y=Value), color='red') +
  #geom_line(aes(y=zero_consec), color='green', size=1.1) +
  #geom_line(aes(y=quantzero), color='green') +
  geom_line(aes(y=100*final_percentile), color='blue') + 
  coord_cartesian(ylim=c(0,100))
  

#Compute percentile statistics threshold across all consecutive zero-flow day 
#out_dt[, zerothresh := quantile(zero_consec, 1-prob, na.rm=T)] #To reproduce Fig. 6

#------------Compare with Sam Zipper's function --------------------------------
dv_in <- ash_sub[,
                 list(Date = mdy(as.character(format(Timestamp, "%m/%d/%Y"))),
                      discharge = forecast::na.interp(Value))] %>%
  as.data.frame

ash_droughti <- HydrologicalDrought_Daily_VanHuijgevoortEtAl2012HESS(
  dv=dv_in, M=31, thres.Q.prc=0.2)
setDT(ash_droughti)

ggplot(ash_droughti[Date %in% as.Date(ash_preddates),], aes(x=Date)) + 
  geom_line(aes(y=discharge)) +
  geom_point(data=ash_droughti[drought==1 & Date %in% as.Date(ash_preddates),], 
             aes(y=discharge), color='red') +
  #geom_line(aes(y=zero_consec), color='green', size=1.1) +
  #geom_line(aes(y=quantzero), color='green') +
  geom_step(aes(y=100*prc), color='blue') + 
  coord_cartesian(ylim=c(0,100))
