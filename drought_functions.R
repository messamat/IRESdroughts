#------ diny -----------------
#' Number of days in the year
#'
#' Computes the number of days in a given year, taking in account leap years
#'
#' @param year Numeric or integer vector of the year (format: 4-digit year, %Y).
#'
#' @return Number of days in that year.
#'
#' @examples
#' diny(1999)
#' diny(2000)
#' diny(2004)
#' diny(2100)
#' diny(1600)
#'
#' @export
diny <- function(year) {
  365 + (year %% 4 == 0) - (year %% 100 == 0) + (year %% 400 == 0)
}

#------ readformatGRDC -----------------
#' Read and pre-format GRDC data
#'
#' Reads text file of daily discharge data for a single GRDC station.
#' Creates columns for year, month, and date of last non-zero flow day +
#' computes yearly number of days of missing data
#'
#' @param path (character) path to the text file of daily discharge data in
#'   standard GRDC format.
#'
#' @return \link[data.table]{data.table} of daily discharge data with additional columns
#'
#' @export
readformatGRDC<- function(path) {
  #extract GRDC unique ID by formatting path
  gaugeno <- strsplit(basename(path), '[.]')[[1]][1]
  
  #Read GRDC text data
  gaugetab <- cbind(fread(path, header = T, skip = 40, sep=";",
                          colClasses = c('character', 'character', 'numeric',
                                         'numeric', 'integer')),
                    GRDC_NO = gaugeno)%>%
    setnames('YYYY-MM-DD', 'dates') %>%
    setorder(GRDC_NO, dates)
  
  #Format data
  gaugetab[, `:=`(year = as.numeric(substr(dates, 1, 4)), #Create year column
                  month = as.numeric(substr(dates, 6, 7)), #create month column
                  integervalue = fifelse(Original == round(Original), 1, 0) #Flag integer discharge values]
  )]
  
  #For each record, compute date of last non-zero flow day
  gaugetab[, prevflowdate := gaugetab[zero_lomf(x=Original),'dates', with=F]] %>% #Get previous date with non-zero flow
    .[Original != 0, prevflowdate:=NA] #If non-zero flow, set prevflowdate to NA
  
  #Compute number of missing days per year, excluding NoData values
  gaugetab[!(Original %in% c(-999, -99, -9999, 99, 999, 9999) | is.na(Original)),
           `:=`(missingdays = diny(year)-.N,
                datadays = .N),
           by= 'year']
  
  return(gaugetab)
}

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
#' 
zero_lomf <- function(x, first=TRUE) {
  if (length(x) > 0) {
    non.zero.idx <- which(x != 0)
    if(first==T & x[1]==0)
      non.zero.idx=c(1,non.zero.idx)
    #Repeat index of previous row with non-zero as many times gap until next non-zero values
    rep.int(non.zero.idx, diff(c(non.zero.idx, length(x) + 1)))
  }
}

# zero_lomf <- function(x, first=TRUE) {
#   if (length(x) > 0) {
#     non.zero.idx <- which(x != 0)
#     if(first==T & x[1]==0) {
#       non.zero.idx=c(1,non.zero.idx)
#       #Repeat index of previous row with non-zero as many times gap until next non-zero values
#       out_ind <- rep.int(non.zero.idx, diff(c(non.zero.idx, length(x) + 1)))
#     } else if (first==F & x[1]==0) {
#       out_ind <- c(rep(NA, min(non.zero.idx)-1), 
#                    rep.int(non.zero.idx, diff(c(non.zero.idx, length(x) + 1))))
#     } else if (first==F & x[1]!=0) {
#       out_ind <- rep.int(non.zero.idx, diff(c(non.zero.idx, length(x) + 1)))
#     }
#     return(out_ind)
#   }
# }


#------ compute_zero_consec -----------------
compute_zero_consec <- function(in_dt, datecol, qcol) {
  in_dt[, prevflowdate := in_dt[zero_lomf(get(qcol), first=T),
                                datecol, 
                                with=F]
  ] %>% 
    .[(get(qcol) != 0) | is.na(get(qcol)), prevflowdate:=NA] %>% #If non-zero flow, set prevflowdate to NA
    .[!is.na(prevflowdate), `:=`(
      zero_consec = as.integer(difftime(get(datecol),  #Number of consecutive zero-flow days
                                        prevflowdate,
                                        units='days')
      )
    )]
}

#---------------- compute_tlm: variable threshold level method (Van Huijgevoort) -------------
compute_tlm <- function(in_dt, datecol, qcol, probs, qperc_window=1, na_mothresh, 
                        IDdrought = TRUE, tdrought = c('daily', 'monthly')) {
  process_dt <- copy(in_dt)
  
  #Create a month and a year-month column for processing
  process_dt[, `:=`(
    month=as.numeric(format(get(datecol), '%m')),
    yrmo=format(get(datecol), '%Y-%m'),
    doy=as.integer(format(get(datecol), '%j'))
  )]
  
  if (tdrought == 'daily') {
    #Compute percentile statistics for each DOY
    ecdf_list <- lapply(1:366, function(j) {
      process_dt[doy %in% seq(j - floor(qperc_window/2),
                              j + floor(qperc_window/2)), 
                 ecdf(get(qcol))]
    })
    
    doy_qthresh_dt <- data.table(doy=1:366,
                                 doy_qthresh = unlist(
                                   lapply(1:366, function(j) {
                                     process_dt[doy %in% seq(j - floor(qperc_window/2),
                                                             j + floor(qperc_window/2)), 
                                                quantile(get(qcol), probs, na.rm=T)]
                                   }))
    )
    
    #Compute percentile statistics for each daily record (Step 4)
    process_dt[, quantdoyr := mapply(function(d, v) {ecdf_list[[d]](v)}, 
                                       d=doy,
                                       v=get(qcol))]
    
    out_dt <- process_dt[doy_qthresh_dt, on='doy']
    
    out_dt[get(qcol) <= doy_qthresh, drought := 1]
    
    out_dt[, drought_id := rleid(is.na(drought))]
  } 
  else if (tdrought == 'monthly') {
    #Aggregate to monthly time steps
    dtmo <- process_dt[, list(
      Qmo = mean(get(qcol), na.rm=T), #Compute mean monthly discharge
      Nmiss = .SD[is.na(get(qcol)), .N], #Compute monthly number of missing days of data
      date = .SD[.N, as.Date(get(datecol))] #Keep date column, converted to date
    ), by=.(yrmo, month)]
    
    
    ecdf_list <- lapply(1:12, function(m) {
      dtmo[month == m, 
           ecdf(Qmo)]
    })
    
    dtmo[, qpercentile := mapply(function(m, v) {ecdf_list[[m]](v)}, 
                                 m=month,
                                 v=Qmo)]
    
    #Compute percentiles
    out_dt <- dtmo[, quantmo := .SD[Nmiss < na_mothresh, #Only consider months with number of NA discharge below threshold
                                    quantile(Qmo, probs, na.rm=T)], #Compute quantile
                   by= month]
    
    out_dt[Qmo <= quantmo, drought := 1] 
  }
  
  return(out_dt)
}

#-------------- compute_tlm_cdpm: mixed threshold level and consecutive dry period method (Van Huijgevoort) ------------------------
compute_tlm_cdpm <- function(in_dt, datecol, qcol, prob, na_mothresh, qperc_window = 1,
                             IDdrought = TRUE, tdrought = c('daily', 'monthly')
) {
  #Compute drought percentile based
  out_dt <- copy(in_dt)
  
  #Compute percentile statistics for each DOY (step 1 in Van Huijgevoort)
  out_dt[, doy := as.integer(format(get(datecol), '%j'))]
  
  ecdf_list <- lapply(1:366, function(j) {
    out_dt[doy %in% seq(j - floor(qperc_window/2),
                        j + floor(qperc_window/2)), 
           ecdf(get(qcol))]
  })
  
  #Compute percentile statistics for each daily record (Step 4)
  out_dt[, qpercentile := mapply(function(d, v) {ecdf_list[[d]](v)}, 
                                 d=doy,
                                 v=get(qcol))]
  
  # #Get overall empirical cumulative distribution with moving window mean 
  # qecdf <- out_dt[,ecdf(frollmean(get(qcol), n=qperc_window, fill=NA,
  #                                 align="center", na.rm=T))]
  
  if (out_dt[, quantile(get(qcol), 0.05, na.rm=T) > 0]) { #If fifth percentile > 0, return TLM method.
    out_dt[, drought_final := as.numeric(qpercentile <= prob)] #(step 2 in Van Huijgevoort)
    return(out_dt)
  } 
  else { #Otherwise compute Van Huijgoort method
    
    #For each record, compute date of last non-zero flow day (step 3 in Van Huijgevoort)
    compute_zero_consec(in_dt=out_dt, datecol=datecol, qcol=qcol) 
    out_dt[!is.na(prevflowdate), dryordrought := 1] 
    
    #Compute drought for Q > 0 (step 4 in Van Huijgevoort)
    out_dt[(get(qcol) > 0 & qpercentile <= prob & !is.na(get(qcol))), 
           dryordrought := 1] 
    
    #Compute consecutive days of drought with Q >0 | dry days (i.e., Ndry,drought; step 5 in Van Huijgevoort)
    out_dt[, drought_consec := cumsum(dryordrought), by=rleid(dryordrought)] %>%
      .[is.na(drought_consec), drought_consec := 0] %>% 
      .[is.na(zero_consec), zero_consec := 0]
    
    
    #Get proportion of zero-flow values for each date (with moving window)
    Fdry_dt <- data.table(doy = 1:366,
                          Fdry = unlist(lapply(1:366, function(j) ecdf_list[[j]](0))))
    out_dt <- out_dt[Fdry_dt, on='doy']
    
    #Compute percentile value for each drought day (in terms of consecutive # of drought days)
    #based on empirical cumulative distribution function of zero-flow days only 
    #(i.e., normalized by proportion of zero-flow days) — step 6 in Van Huijgevoort
    zero_ecdf <- out_dt[zero_consec > 0, ecdf(zero_consec)]
    out_dt[, drought_percentile := (1-zero_ecdf(drought_consec))*Fdry]
    
    #Combine percentile values for wet and dry periods
    out_dt[, final_percentile := fifelse(get(qcol) > 0, 
                                         qpercentile,
                                         drought_percentile)]
    
    #Determine whether a drought is occurring (step 7 in Van Huijgevoort)
    out_dt[, drought_final := as.numeric(final_percentile <= prob)]
    
  }
  
  return(out_dt)
}
