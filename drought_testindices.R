#devtools::install_github('buzacott/bomWater')
#devtools::install_github('WillemMaetens/standaRdized')

source('drought_packages.R')
source('SamZipper.R')

#Define directory structure
rootdir <- find_root(has_dir("src"))
resdir <- file.path(rootdir, "results")
figdir <- file.path(resdir, 'figures')
datdir <- file.path(rootdir, "data")

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

#-------------- Compute Standardized Streamflow Index (SSI) with SPEI package --------------------
ssi_df <- spei(data=ts(na.interp(ash_sub$Value), frequency=12), scale=365, na.rm=T)
plot(ssi_df)

plot(ecdf((na.interp(ash_sub$Value))))

#-------------- Compute Standardized Streamflow Index for IR (SSIIR) with SPEI package --------------------
ssiir_indt <- copy(ash_sub)
compute_zero_consec(ssiir_indt, datecol = 'Timestamp', qcol = 'Value')
ssiir_indt[is.na(zero_consec), zero_consec := 0]
ssiir_df <- spei(data=ts(na.interp(ssiir_indt$zero_consec), frequency=12), scale=365, na.rm=T)

plot(ssiir_df)
ssiir_indt[, ssiir := ssiir_df$fitted]

ggplot(ssiir_indt, aes(x=as.Date(Timestamp), y=-ssiir)) + 
  geom_line(color='red') + 
  geom_line(aes(y=Value)) +
  geom_hline(yintercept = 0) +
  coord_cartesian(ylim=c(-3,3)) +
  scale_x_date(limits=c(as.Date('1988-01-01'), 
                        as.Date('1993-01-01')))

#-------------- Compute Standardized Streamflow Index  with SCI package --------------------
in_scidt <- ash_sub[, mean(Value, na.rm=T), by=format(Timestamp, '%Y-%m')]
sci_para <- fitSCI(na.interp(in_scidt$V1),
                   distr='gamma', time.scale = 12, first.mon=1, p0 = TRUE, 
                   p0.center.mass=TRUE)
in_scidt[, sci := transformSCI(na.interp(in_scidt$V1),first.mon=1,obj=sci_para)]

ggplot(in_scidt, aes(x=as.Date(paste0(format, '-01')), y=sci)) + 
  geom_line(color='blue') +
  geom_line(data=ssiir_indt, aes(x=as.Date(Timestamp), y=-ssiir), color='red') +
  geom_hline(yintercept=0) +
  coord_cartesian(ylim=c(-3,3)) +
  scale_x_date(limits=c(as.Date('1988-01-01'), 
                        as.Date('1993-01-01')))

#-------------- Compute Variable threshold level method (Van Huijgevoort) ------------------------
ashperct_mwindow <- compute_tlm(
  in_dt = ash_sub,
  datecol = 'Timestamp',
  qcol = 'Value',
  IDdrought = TRUE,
  probs = 0.2,
  na_mothresh = 4,
  qperc_window = 31,
  tdrought = 'daily')


ggplot(ashperct_mwindow,
       aes(x=as.Date(Timestamp), y=qpercentile, group=1)) + 
  geom_line() +
  geom_line(aes(y=doy_qthresh)) +
  coord_cartesian(ylim=c(0,0.4)) +
  scale_x_date(limits=c(as.Date('1988-01-01'), 
                        as.Date('1993-01-01')))

ashperct <- compute_tlm(
  in_dt = ash_sub,
  datecol = 'Timestamp',
  qcol = 'Value',
  IDdrought = TRUE,
  probs = 0.05,
  na_mothresh = 4,
  tdrought = 'monthly')

ggplot(ashperct, aes(x=as.Date(paste0(ym, '-01')), y=quantmo, group=1)) + 
  geom_line(color='blue', size=1.2) + 
  #geom_point(data=ashperct[Value==0,], aes(x=as.Date(Timestamp), y=Value), color='green', size=1.5) +
  geom_line(data=ash_sub, aes(x=as.Date(Timestamp), y=Value)) + 
  #geom_line(aes(y=drought, group=drought_id), color='red', size=1.2) +
  scale_y_continuous(
    #trans=scales::pseudo_log_trans(base = 10), 
    breaks=c(0,0.1,0.5,1,10,100,1000, 2500)) + 
  scale_x_date(limits=c(as.Date('1988-01-01'), 
                        as.Date('1993-01-01'))) + 
  coord_cartesian(ylim=c(0,0.2)) + 
  theme_classic()

#-------------- Compute mixed threshold level and consecutive dry period method (Van Huijgevoort) ------------------------
ash_tlmcdpm <- compute_tlm_cdpm(in_dt = ash_sub,
                                datecol = 'Timestamp',
                                qcol = 'Value',
                                qperc_window = 30,
                                IDdrought = FALSE,
                                prob = 0.2,
                                tdrought = 'daily')

ggplot(ash_tlmcdpm[Timestamp %in% ash_preddates,], aes(x=Timestamp)) + 
  geom_line(aes(y=Value)) +
  geom_point(data=ash_tlmcdpm[Timestamp %in% ash_preddates & drought_final==1,], aes(y=Value), color='red') +
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
