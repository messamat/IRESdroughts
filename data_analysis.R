library(dataRetrieval)
library(data.table)
library(ggplot2)
library(magrittr)

######################### SCOPPETTONE ET AL 2015 ############################
#Get hydrologic data for Scoppettone et al. 2015
scoppettone_Q <- readNWISdv(siteNumber = "10351700", 
                        parameterCd = '00060',
                        startDate = '1800-01-01',
                        endDate = Sys.Date())
min(scoppettone_Q$X_00060_00003)

######################### GOLLADAY ET AL. 2004 ##################################
#Get hydrologic data for Golladay et al. 2004
golladay_Q1 <- readNWISdv(siteNumber = "02357000", 
                            parameterCd = '00060',
                            startDate = '1800-01-01',
                            endDate = Sys.Date()) %>%
  as.data.table
min(golladay_Q1$X_00060_00003)
ggplot(golladay_Q1, aes(x=Date, y=X_00060_00003)) + 
  geom_point() +
  geom_point(data=golladay_Q1[X_00060_00003==0,], color='red') +
  scale_y_sqrt()

golladay_Q2 <- readNWISdv(siteNumber = "02354500", 
                          parameterCd = '00060',
                          startDate = '1800-01-01',
                          endDate = Sys.Date()) %>%
  as.data.table
readNWISsite('02354500')

ggplot(golladay_Q2, aes(x=Date, y=X_00060_00003)) + 
  geom_point() +
  geom_point(data=golladay_Q2[X_00060_00003==0,], color='red') +
  scale_y_sqrt()

########################## STEFERRUD ET AL. 2011 #################################
stefferud_Q <- readNWISdv(siteNumber = "09430500",  #Gila near Gila
                          parameterCd = '00060',
                          startDate = '1800-01-01',
                          endDate = Sys.Date()) %>%
  as.data.table
min(stefferud_Q$X_00060_00003)

stefferud_Q2 <- readNWISdv(siteNumber = "09430030",  #Gila River nr Gila Hot Springs, NM
                           parameterCd = '00060',
                           startDate = '1800-01-01',
                           endDate = Sys.Date()) %>%
  as.data.table
min(stefferud_Q2$X_00060_00003)
  

stefferud_Q3 <- readNWISdv(siteNumber = "09430010", #WEST FORK GILA RIVER AT GILA CLIFF DWELLINGS, NM
                           parameterCd = '00060',
                           startDate = '1800-01-01',
                           endDate = Sys.Date()) %>%
  as.data.table
min(stefferud_Q3$X_00060_00003, na.rm=T)

stefferud_Q4 <- readNWISdv(siteNumber = "09430020", #W FORK GILA R BLW MDL FORK NR GILA HOT SPRINGS, NM
                           parameterCd = '00060',
                           startDate = '1800-01-01',
                           endDate = Sys.Date()) %>%
  as.data.table
min(stefferud_Q4$X_00060_00003, na.rm=T)
  



vaughn_Q <- readNWISdv(siteNumber = "07336200", #W FORK GILA R BLW MDL FORK NR GILA HOT SPRINGS, NM
                           parameterCd = '00060',
                           startDate = '1800-01-01',
                           endDate = Sys.Date()) %>%
  as.data.table
min(vaughn_Q$X_00060_00003, na.rm=T)

ggplot(vaughn_Q, aes(x=Date, y=X_00060_00003)) + 
  geom_point() +
  geom_point(data=vaughn_Q[X_00060_00003==0,], color='red') +
  scale_y_sqrt()

