source('drought_packages.R')

#Define directory structure
rootdir <- find_root(has_dir("src"))
resdir <- file.path(rootdir, "results")
figdir <- file.path(resdir, 'figures')
datdir <- file.path(rootdir, "data")

#Spreadsheet
data_spreadsheet <- file.path(datdir, 'Review Ecohydro_20210530.xlsx')
dataraw <- read_excel(data_spreadsheet, sheet = "Final_data") %>%
  setDT
misuse_ini <- read_excel(data_spreadsheet, sheet = "misuse drought_drying") %>%
  setDT

IRESstats <- fread(file.path(datdir, 'IRES_globalstats_Messageretal2021.csv')) %>%
  .[, Climate := factor(Climate, levels=.SD[, Climate])]

#Download climate zone data
gensdl <- file.path(datdir, 'GEnSv3')
if (!(dir.exists(gensdl))) {
  download.file("https://datashare.ed.ac.uk/bitstream/handle/10283/3089/GEnSv3.zip?sequence=7&isAllowed=y",
                destfile = paste0(gensdl,'.zip'))
  unzip(gensdl)
}

gens_shp <- file.path(gensdl, 'GEnS_v3.shp')
gens_ras <- file.path(gensdl, 'gens_v3.tif')


#---------------- Overall statistics -------------------------------------------
#Eligible (without excluding for misusing drought)
dataraw[Suitable %in% c('Y', 'NA') | `Mix-up drought and drying`=='Y', .N] + nrow(misuse_ini)

#Eligible that misuse drought
dataraw[`Mix-up drought and drying`=='Y', .N] + nrow(misuse_ini)
dataraw[Suitable == 'Y' & `Mix-up drought and drying`=='Y', .N] #For the sake of simplicity, keep those 2 out of the count of misusing as it's unsure

#Remaining studies
datasub <- dataraw[Suitable %in% c('Y', 'NA'),] %>%
  .[, `:=`(Latitude = as.numeric(Latitude),
           Longitude = as.numeric(Longitude))]
nrow(datasub)

#---------------- Plot studies -------------------------------------------
#Intersect studies with climate zones
datap <- st_as_sf(x = datasub[!is.na(Longitude) & Suitable == 'Y',], #Exclude Resh et al. 2013 because re-analysis of data from multiple studies
                  coords = c('Longitude', 'Latitude'),
                  crs = 4326)

#Intersect with climate
tempfst <- file.path(resdir, 'datap_gens.fst')
if (!file.exists(tempfst)) {
  gens <- read_sf(gens_shp) %>%
    st_transform(crs=4326)
  datap_gens <- st_join(datap, gens) %>%
    as.data.table %>%
    cbind(st_coordinates(datap))
  write.fst(datap_gens[,-"geometry", with=F], tempfst)
} else {
  datap_gens <- read.fst(tempfst) %>%
    as.data.table
}

# cluster all points using a hierarchical clustering approach
studyclus <- geodist(st_coordinates(datap),
                     measure='geodesic') %>%
  as.dist %>%
  hclust(method="complete")

datap_gens[, `:=`(cluster = cutree(studyclus, h=100000),
                  Climate = gsub('^[A-Z][.]\\s', '', GEnZname)
)] %>%
  .[, Climate := factor(Climate, levels=.SD[order(GEnZname), unique(Climate)])]
datap_gens[, .N, by=.(cluster, Climate)]

clustcentro <- datap_gens[, list(lon = mean(X), lat = mean(Y), N = .N), 
                          by=.(cluster, Climate)] %>%
  st_as_sf(coords = c('lon', 'lat'), crs = 4326)


#Map it
crs_wintri = "+proj=wintri +datum=WGS84 +no_defs +over"

# gens_ras <- raster(gens_ras) %>%
#   projectRaster(crs=crs_wintri)
cluscentro_gens_wintri <- st_transform(clustcentro, crs_wintri)

wcountries <- rnaturalearth::ne_countries(
  scale = "medium", returnclass = "sf") %>%
  st_as_sf %>%
  st_transform(crs_wintri, use_gdal = FALSE)


names(gens_ras)
studiesmap <- ggplot() +
  geom_sf(data=wcountries, color='#bdbdbd', alpha=1/2) +
  geom_sf(data = cluscentro_gens_wintri, 
          aes(size=N, color=Climate), alpha=0.75) + #color = alpha('#3182bd', 1/2)
  scale_size_continuous(name=str_wrap('Number of studies within 100 km', 20),
                        range=c(1, 5), breaks=c(1, 2, 5, 9)) +
  scale_color_manual(name='', values=c('#08519c', '#2171b5', '#67a9cf', '#006d2c',
                                       '#41ab5d', '#fd8d3c', '#e31a1c',
                                       '#800026')) +
  coord_sf(datum = NA, expand=F) +
  theme_map() +
  theme(legend.position=c(0, 0.25),
        #legend.direction="horizontal",
        legend.text = element_text(size=9),
        legend.title = element_text(size=9),
        legend.key.height = unit(0.1,"in"),
        legend.background = element_blank())


#Plot histogram of study subject (fish, plant, amphibian)
datasub[, c("organism1", "organism2", "organism3") := tstrsplit(`Orgnism/process`, ";", fixed=TRUE)]

orgadat <- melt(datasub[Suitable == 'Y',], id.vars = c('Author Full Names', 'Article Title'), 
                measure.vars=c('organism1', 'organism2', 'organism3'))[!is.na(value),] %>%
  .[, organism:= gsub('^ ', '', str_to_title(value))] %>%
  .[, list(nstudies = .N), by=organism] %>%
  .[, organism := factor(organism, levels=.SD[order(-nstudies,),organism])]

orgap <- ggplot(orgadat, aes(x=organism, y=nstudies)) + 
  geom_bar(stat='identity', fill='#006d2c', alpha=1/3) + 
  scale_y_reverse(name='Number of studies', expand=c(0,0),
                  breaks=c(0,2, 5, 10, 20, 30)) +
  scale_x_discrete(name='Ecosystem component of study', position='top') +
  coord_flip() +
  theme_classic()
#Map of location with climate zones background
#Mention gaps in coverage based on climate zone
#Point to proportion of studies that are either about the Millennial Drought or the California Drought of the 2010s

#Plot histogram by climate
IRESstats[, totperc := (Intermittence_mapped*Length_mapped)/(0.41*sum(Length_mapped))] 

climcomp <- merge(IRESstats,
                  as.data.table(datap_gens)[,list(studyperc = 100*.N/nrow(datap_gens)), 
                                            by=Climate], 
                  by='Climate', all.x=T) %>%
  melt(id.vars='Climate', measure.vars=c('studyperc', 'totperc'), 
       value.name='perc')

climp <- ggplot(climcomp, aes(x=Climate, y=perc, fill=variable)) +
  geom_bar(stat='identity', position='dodge', alpha=0.6)+ 
  scale_fill_manual(labels = c('Studies on drought in IRES', 'Global length of IRES'),
                    values = c('#fa9fb5', '#c51b8a')) +
  scale_y_continuous(name='Percentage of total', expand=c(0,0)) +
  scale_x_discrete(name='Climate zone') +
  coord_flip() +
  theme_classic() + 
  theme(legend.position=c(0.65,0.15),
        legend.title = element_blank(),
        legend.background = element_blank())
#Map of location with climate zones 


pdf(file.path(figdir, paste0('summary_fig', format(Sys.Date(), '%Y%m%d'), '.pdf')), 
    height = 8, width = 8)
grid.arrange(studiesmap, climp, orgap, layout_matrix = rbind(c(1,1,1,1,1,1,1),
                                                             c(1,1,1,1,1,1,1),
                                                             c(2,2,2,2,3,3,3)))
dev.off()

#Check type of studies for written statistics
datasub[, c("type1", "type2") := tstrsplit(Theme, ";", fixed=TRUE)]
typedat <- melt(datasub[!(Ecosystem %in% c('Review', 'Experiment', 'River, experiment')),],
                id.vars='DOI', measure.vars = c('type1', 'type2'))
datasub[!(Ecosystem %in% c('Review', 'Experiment', 'River, experiment')), .N]
typedat[, .N/datasub[!(Ecosystem %in% c('Review', 'Experiment', 'River, experiment')), .N], by=str_to_lower(value)]

#Check organism and climate of study for written statistics (excludes reviews, includes experiments)
orgadat[, nstudies/datasub[Suitable == 'Y',.N], by=organism]

#Studies in mediterranean, semi-arid and arid areas
targetclims <- c('Cool temperate and dry', 'Cool temperate and xeric',
                 'Warm temperate and mesic', 'Warm temperate and xeric',
                 'Hot and dry', 'Hot and arid',
                 'Extremely hot and arid', 'Extremely hot and xeric')
climcomp[Climate %in% targetclims & variable == 'studyperc', sum(perc, na.rm=T)]

climcomp[!(Climate %in% targetclims) & variable == 'totperc', sum(perc, na.rm=T)]

#---------------- Analyze drought reporting of studies -------------------------
drought_measures <- datasub[!(`Type of drought measure` %in% c(NA, 'NA')),
                            c('Drought measure', 'Type of drought measure', 'Transferability of drought measure')]

drought_measures[, .N, by='Type of drought measure']
drought_measures[, .N, by=c('Type of drought measure','Transferability of drought measure')]







