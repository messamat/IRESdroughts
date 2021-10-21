
#Define directory structure
rootdir <- find_root(has_dir("src"))
resdir <- file.path(rootdir, "results")
figdir <- file.path(resdir, 'figures')
datdir <- file.path(rootdir, "data")

GSIMy <- file.path(datdir, 'GSIM', 'GSIM_indices', 'TIMESERIES', 'yearly')

GSIMy_list <-  file.path(GSIMy, list.files(GSIMy))


read_GSIMy <- function(path) {
  print(path)
  #Read GSIM text data
  gaugetab <- fread(path, header=T, skip = 21, sep=",",
                    colClasses=c('Date', rep('numeric', 23),
                                 'integer', 'integer')) %>%
    setnames(new=gsub('[\\\t]|["]', '', names(.))) %>% #Remove tab in field names
    setorder(date) %>%
    .[, year := as.numeric(substr(date, 1, 4))] #Create month column
  
  #Compute minimum number of zero-flow days for each month
  gaugetab[, mDur_miny := fcase(
    P90==0, as.integer(floor(0.9*n.available)),
    P80==0, as.integer(floor(0.8*n.available)),
    P70==0, as.integer(floor(0.7*n.available)),
    P60==0, as.integer(floor(0.6*n.available)),
    P50==0, as.integer(floor(0.5*n.available)),
    P40==0, as.integer(floor(0.4*n.available)),
    P30==0, as.integer(floor(0.3*n.available)),
    P20==0, as.integer(floor(0.2*n.available)),
    P10==0, as.integer(floor(0.1*n.available)),
    MIN7==0L, 7L,
    MIN==0L, 1L,
    default = 0L
    )]
  
  outtab <- gaugetab[n.missing < 20, list(mDur_min = mean(mDur_miny),
                                          nyears_keep = .N)
  ]
  
  return(outtab)
} 

GSIMstats <- lapply(GSIMy_list, read_GSIMy) %>% rbindlist
GSIMstats[, IRES := fifelse(mDur_min >= 7, '1', '0')]

#Get stats for paper
GSIMstats[, .N, by=IRES]
GSIMstats[, .SD[IRES=='1',.N]/.N]
GSIMstats[, mean(nyears_keep), by=IRES]

ggplot(GSIMstats, aes(x=mDur_min, y=nyears_keep)) + 
  geom_point() +
  geom_quantile(quantiles=c(0.25, 0.5, 0.75, 0.9))
