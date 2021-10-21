source('drought_packages.R')
source('drought_functions.R')


#Define directory structure
rootdir <- find_root(has_dir("src"))
resdir <- file.path(rootdir, "results")
figdir <- file.path(resdir, 'figures')
datdir <- file.path(rootdir, "data")
grdcdir <- file.path(datdir, 'GRDCdat')

#------------- Get streamflow gauging stations ---------------------------------------------------
#Toise River at Forkroad
perennial_dt <- readformatGRDC(file.path(grdcdir, "1160775.txt")) %>%
  .[Original == -999, Original := NA]

#GQUNUBE RIVER @ OUTSPAN (27672101) -32.801667	27.855278
iresdrought_dt <- readformatGRDC(file.path(grdcdir, "1160660.txt")) %>%
  .[Original == -999, Original := NA] 

#BUFFELSRIVIER @ SCHURVEPOORT (27661402)
shift_dt <- readformatGRDC(file.path(grdcdir, "1160851.txt")) %>%
  .[Original == -999, Original := NA]  %>%
  .[2:.N, Original := zoo::na.locf(Original)]


#------------ Search a weather station around 1160635 and 1160660 --------------
#For perennial station 1160635
out <- isd_stations_search(lat=-32.5, lon=27, rad=50) 

dohne <- lapply(1973:2013, function(y) isd(usaf=out[out$station_name=='DOHNE', 'usaf'],
                                           wban="99999", year=y, cores=4, parallel=T)) %>% 
  rbindlist(fill=T, use.names=T) %>%
  setDT %>%
  .[as.numeric(AA1_depth) >= 1000, AA1_depth := NA] %>%
  .[, yrmo := as.Date(paste0(format(as.Date(date, '%Y%m%d'), '%Y-%m'), '-01'))] %>%
  .[yrmo >= as.Date('1976-09-01') & yrmo <= as.Date('1996-09-01'), list(precip_mo = sum(as.numeric(AA1_depth), na.rm=T)), by=yrmo]


ggplot(dohne, aes(x=yrmo, y=precip_mo, group=1)) + 
  geom_line() + 
  geom_smooth()

#------------ Compute SPI ------------------------------------------------------
spi_para <- fitSCI(na.interp(dohne$precip_mo),
                   distr='gamma', time.scale = 24, first.mon=1, p0 = FALSE)
dohne[, SPI := transformSCI(na.interp(dohne$precip_mo),first.mon=1,obj=spi_para)]

ggplot(dohne, aes(x=yrmo, y=SPI)) + 
  geom_line(color='blue') +
  #geom_line(aes(x=yrmo, y=-ssiir), color='red') +
  geom_hline(yintercept=0) +
  coord_cartesian(ylim=c(-3,3))

#------------ Format and join data.table for plotting  -------------------------
disprecip_dt <- rbindlist(list(perennial_dt, iresdrought_dt, shift_dt), fill=T) %>%
  .[as.Date(dates) %in% seq(as.Date('1978-09-01'), as.Date('1990-01-01'), by='day'),
    `:=`(plotgrp =  rleid(as.numeric(is.na(Original))),
           recmax = 120), by=GRDC_NO] %>%
  .[, yrmo := as.Date(paste0(format(as.Date(dates, '%Y-%m-%d'), '%Y-%m'), '-01'))] %>%
  .[dohne, on='yrmo']

#Create labels for each plot
anno_dt <- data.table(GRDC_NO = c(1160660, 1160775, 1160851),
                      GRDC_name = c("a. Gqunube River at Outspan",
                                    "b. Toise River at Forkroad",
                                    "c. Buffelsrivier at Schurvepoort"
                      ),
                      dates = rep('1978-12-01',3),
                      Original = 100
)


#------------ Plot hydrographs and meteorological drought for three gauges  -------------------------
ZA1983_droughtplot <- ggplot(disprecip_dt, aes(x=as.Date(dates), y=na.interp(Original))) +
  geom_rect(aes(xmin=as.Date(dates), xmax=as.Date(dates) + 30, 
                ymin=0, ymax=recmax, 
                fill = SPI), alpha=0.1) +
  geom_line(color='#737373') +
  geom_point(data=disprecip_dt[Original==0,], color='black', size=1) +
  geom_text(data=anno_dt, aes(label=GRDC_name), hjust = 0) +
  scale_x_date(name='Date', 
               date_breaks='2 years', date_labels = '%Y',
               limits=as.Date(c('1978-09-01', '1990-01-01')), expand=c(0,0)) +
  scale_y_sqrt(name=expression(Discharge~m^3~s^-1), limits=c(0, 120), expand=c(0,0)) + 
  scale_fill_distiller(name='SPI', palette='RdBu', direction=1) +
  facet_wrap(~GRDC_NO, ncol=1, scales='free_y') + 
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(0.6, "lines"),
        plot.margin = unit(c(0.2,0.1,0.1,0.1), unit='cm'),
        legend.margin = margin(l = 0, unit='cm'))


# pdf(paste0(file.path(figdir, 'example_hydrographs_drought'),
#            format(Sys.Date(), '%Y%m%d'), '.pdf'),
#     width=6.654, height=4)
# ZA1983_droughtplot
# dev.off()

png(paste0(file.path(figdir, 'example_hydrographs_drought'),
           format(Sys.Date(), '%Y%m%d'), '.png'),
    width=6.654, height=4, unit='in', res=600)
ZA1983_droughtplot
dev.off()



#------------ Plot drought indices  --------------------------------------------
#Compute SSI
gqunube_dis <- iresdrought_dt[, dates := as.Date(dates, "%Y-%m-%d")] %>%
  .[, yrmo := as.Date(paste0(format(dates, '%Y-%m'), '-01'))] %>%
  .[min(which(!is.na(Original))):max(which(!is.na(Original))),]

gqunube_dismo <- gqunube_dis[, list(qmo = mean(Original, na.rm=T)), 
                             by=yrmo]

ssi_para <- fitSCI(na.interp(gqunube_dismo$qmo),
                   distr='gamma', time.scale = 24, first.mon=1, p0 = TRUE, 
                   p0.center.mass=TRUE)
gqunube_dismo[, SSI := transformSCI(gqunube_dismo$qmo,
                                    first.mon=1,obj=ssi_para)]


#Compute SSIIR
# compute_zero_consec(gqunube_dis, datecol = 'dates', qcol = 'Original')
# gqunube_dis[is.na(zero_consec), zero_consec := 0]

gqunube_irmo <- gqunube_dis[, list(irmo = sum(as.numeric(Original == 0), 
                                              na.rm=T)), 
                            by=yrmo]

ssiir_para <- fitSCI(gqunube_irmo$irmo,
                     distr='gamma', time.scale = 24, first.mon=1, p0 = TRUE, 
                     p0.center.mass=TRUE)
gqunube_irmo[, SFI := -transformSCI(gqunube_irmo$irmo,
                                   first.mon=1,obj=ssiir_para)]

sci_dtformat <- rbindlist(list(gqunube_irmo, gqunube_dismo, dohne), 
                          use.names = T, fill = T) %>%
  melt(id.vars = 'yrmo', measure.vars=c('SPI','SSI','SFI')) %>%
  .[!is.na(value),]


#Compute Variable threshold level method (Van Huijgevoort)
gqunube_tlm_window <- compute_tlm(
  in_dt = gqunube_dis,
  datecol = 'dates',
  qcol = 'Original',
  IDdrought = TRUE,
  probs = 0.1,
  qperc_window = 31,
  tdrought = 'daily') %>%
  .[, yrmo := as.Date(paste0(yrmo, '-01'))]%>%
  .[order(dates), plotgrp := rleid(is.na(drought))]

gqunube_tlm_monthly <- compute_tlm(
  in_dt = gqunube_dis,
  datecol = 'dates',
  qcol = 'Original',
  IDdrought = TRUE,
  probs = 0.1,
  na_mothresh = 4,
  tdrought = 'monthly') %>%
  .[, yrmo := as.Date(paste0(yrmo, '-01'))] %>%
  .[, plotgrp := rleid(is.na(drought))]


gqunube_tlmcdpm <- compute_tlm_cdpm(in_dt = gqunube_dis,
                                    datecol = 'dates',
                                    qcol = 'Original',
                                    qperc_window = 60,
                                    IDdrought = FALSE,
                                    prob = 0.10,
                                    tdrought = 'daily') %>%
  .[order(dates), plotgrp := rleid(drought_final)]


#Plot
sci_tsplot <- ggplot(sci_dtformat, aes(x=yrmo, y=value, color=variable)) + 
  geom_hline(yintercept=c(-2, 2), alpha=0.1) +
  geom_hline(yintercept=0, alpha=0.5) +
  geom_vline(xintercept=seq(as.Date('1980-01-01'), 
                            as.Date('1996-01-01'), '2 years'), alpha=0.1) +
  geom_line(size=1) +
  scale_color_manual(name='', values=c('gray', '#1f78b4', '#ff7f00')) +
  scale_y_continuous(name='Standardized index') +
  scale_x_date(name='Date', 
               date_breaks='2 years', date_labels = '%Y',
               limits=as.Date(c('1978-09-01', '1996-09-01')), expand=c(0,0)) +
  theme_classic() +
  theme(plot.margin = unit(c(t=-2,r=0,b=0,l=0), unit='cm'))

ybks <- c(0, 1, 10, 25, 50, 100)
threshold_tsplot <- ggplot(gqunube_tlmcdpm, aes(x=dates, y=sqrt(Original))) + 
  geom_vline(xintercept=seq(as.Date('1980-01-01'), 
                            as.Date('1996-01-01'), '2 years'), alpha=0.1) +
  geom_line(color='#737373', aes(group=1)) + 
  scale_y_continuous(name=expression(Discharge~m^3~s^-1),
                     breaks=sqrt(ybks), labels=ybks, expand=c(0.01,0)) +
  #geom_point(data=gqunube_tlmcdpm, aes(y=final_percentile), color='red') +
  geom_line(data=gqunube_tlm_window[drought==1,], 
            aes(group=plotgrp), 
            color='#d6604d', size=1.2) +
  geom_line(data=gqunube_tlmcdpm[drought_final==1,], 
            aes(group=plotgrp, y=Original+0.2), 
            color='#92c5de', size=0.75) +
  scale_x_date(date_breaks='2 years', date_labels = '%Y',
               limits=as.Date(c('1978-09-01', '1996-09-01')), expand=c(0,0)) +
  theme_classic() +
  theme(plot.margin = unit(c(t=0,r=0,b=0,l=-1), unit='cm'),
        axis.title.x = element_blank())
  


pdf(paste0(file.path(figdir, 'example_drought_indices'),
           format(Sys.Date(), '%Y%m%d'), '.pdf'),
    width=6.654, height=4)
threshold_tsplot / sci_tsplot 
dev.off()

