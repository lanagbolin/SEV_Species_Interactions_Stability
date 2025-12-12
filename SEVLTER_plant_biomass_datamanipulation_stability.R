####### To reshape SEV estimated plant biomass data to wide format for analysis #############
#######  J Rudgers July 2023 ##############################################################
rm(list=ls(all=TRUE)) #give R a blank slate
library(reshape2)
library(vegan)
library(tidyverse)
library(dplyr)

# Set working directory
setwd("C:/Users/jrcassvc/Desktop/SEV Data/NPP/")

################ import data from EDI
# Or use this function to import
mass<-read.csv("C:/Users/jrcassvc/Desktop/SEV Data/NPP/Sevilleta_allbiomass_05Feb2023.csv",stringsAsFactors = T) #
str(mass)

# create unique identifier for quadrat - need to adjust this line based on the focal site or experiment
# for core sites need site, web, plot, quad to ID the quadrat - this includes everything, clunky, but comprehensive
mass$quad_ID<-as.factor(paste(mass$site,mass$web,mass$transect,mass$treatment,mass$block,mass$plot,mass$subplot,mass$quad,sep="_"))

# cast the data to create columns for each plant species so zeros are added for dates where species was absent  ---------------------
# for all quads
mass_species<-dcast(mass,year+site+season+season.precip+GDD+SPEI.comp+SiteCluster+MetStation+web+transect+treatment+block+plot+subplot+quad+quad_ID ~ kartez,sum,value.var="biomass.BM",fill=0)
# can do the same for photosynthetic pathway or life history: annual, perennial, grass, forb
#mass_path<-dcast(mass,year+site+treat+season+web+plot+subplot+quad_ID+quad~path,sum,value.var="weight",fill=0)
#mass_apgf<-dcast(mass_cntrl,year+site+treat+season+web+plot+subplot+quad_ID+quad~apgf,sum,value.var="weight",fill=0)
summary(mass_species)
#get rid of column called EMPTY - this was a placeholder for quadrats that had bare ground, not a real species
mass_species <- mass_species %>% select(-EMPTY)
mass_species$LATR2 <-mass_species$LATR2+mass_species$STEM
mass_species <- mass_species %>% select(-STEM)
# calculate total biomass
mass_species$totmass<-rowSums(mass_species[,17:338])
summary(mass_species$totmass)

# add richness, diversity, evenness
mass_species$richness<-specnumber(mass_species[,17:338])
mass_species$shannonH<-diversity(mass_species[,17:338])
mass_species$evenness<-mass_species$shannonH/log(mass_species$richness+1)
write.csv(mass_species,"SEV_all_quadrat_biomass_cast_2022.csv")

# reshape again so each species in each quadrat is a row, cols are timepoints
# create a continuous time variable
mass_species$seas<-as.numeric(as.character(recode_factor(mass_species$season,fall="0.5",spring="0")))
summary(mass_species$seas)
mass_species$year_seas<-as.factor(mass_species$year+mass_species$seas )                            
summary(mass_species$year_seas)
mass_species$time<-as.numeric(as.factor(mass_species$year_seas))
summary(mass_species$time)
mass_species$year.f<-as.factor(mass_species$year)

##### create web scale data for mammals and hoppers-----------
# core sites only
summary(mass_species$site)
mass_core<-subset(mass_species, site=="core_black"|
                                site=="core_blue"|
                                site=="core_creosote")
# add ecosystem
mass_core$ecosystem<- recode_factor(mass_core$site,
                                           core_black="Desert grassland", 
                                           core_blue="Plains grassland",
                                           core_creosote="Desert shrubland")
#remove extraneous columns so we can average the rest by web
mass_core_select <- mass_core %>% select(-MetStation,-transect,-site,
                                         -SiteCluster,-transect,-treatment,
                                         -block, -plot, -subplot,-quad, -quad_ID)
summary(mass_core)
# take average per web so data are in g/m2
plant_mass_web<- mass_core_select %>%  group_by(ecosystem,year.f,year,season,year_seas,time,web,season.precip,GDD,SPEI.comp) %>% summarise_at(vars(ACNE:evenness),mean, na.rm = TRUE)
summary(plant_mass_web)
# calculate diversity at the scale of the web rather than avg per quadrat
plant_mass_web$web_shannonH<-diversity(plant_mass_web[,11:332])
plant_mass_web$web_richness<-specnumber(plant_mass_web[,11:332])
plant_mass_web$web_evenness<-plant_mass_web$shannonH/log(plant_mass_web$richness+1)

#export data
SEVPLANTMASS_WEB<-plant_mass_web
summary(SEVPLANTMASS_WEB)
write.csv(SEVPLANTMASS_WEB,"C:/Users/jrcassvc/Desktop/SEV Data/NPP/SEVPLANTMASS_WEB.csv")

  
###### total biomass stability ----------------------
#tot mass
totmass<-mass_species %>% select(year.f,site,season,seas,year_seas,season.precip,GDD,SPEI.comp,
                                 SiteCluster,MetStation,web,transect,treatment,block,plot,subplot,quad,quad_ID,totmass)

totmass_cast<-dcast(totmass,site+season+seas+SiteCluster+MetStation+web+transect+treatment+block+plot+subplot+quad+quad_ID ~ year.f , sum, value.var="totmass",fill=-9999)
# fix problem with fill in dcast
totmass_cast[totmass_cast== "-9999"] <- "NA"
summary(totmass_cast)
totmass_cast <- totmass_cast %>% mutate_if(is.character, as.numeric)

# calculate totmass sd
totmass_cast$sd<-transform(totmass_cast, totmass_stdev=apply(totmass_cast[,14:37], 1, sd, na.rm=TRUE))


# subset for particular experiments
summary(mass_species$site)
mass_MVE_Plains<-subset(mass_species, site=="meanvar_blue")
write.csv(mass_MVE_Plains, "plants_mass__MVE_Plains.csv")




#cast to make a row for each species in each quad, columns are times of observation
#fill missing data with NA because that quad was not observed on that date and zeros have been filled in
mass_species_melt<-melt(mass_species,id=c("year","site","season","seas","time","year_seas","season.precip","GDD","SPEI.comp",
                             "SiteCluster","MetStation","web","transect","treatment","block","plot","subplot","quad","quad_ID"),variable.name="kartez",value.name="biomass.BM")
summary(mass_species_melt$biomass.BM)
summary(mass_species_melt$kartez)

mass_species_cast<-dcast(mass_species_melt,site+SiteCluster+MetStation+web+transect+treatment+block+plot+subplot+quad+quad_ID+kartez ~ year_seas,sum,value.var = "biomass.BM",fill="NA")

#remove all species that were never observed in a quad
mass_species_cast$sum<-rowSums(mass_species_cast[,17:340])

mass_species_stability <-mass_species_cast %>%  filter(sum>0)

#count up number of observations for each species/quad
mass_species_stability$count<-rowSums(!is.na(df[,2:4]))

#calculate temporal SD

#calculate stability as CV


write.csv(mass_species, "plants_mass_stability.csv")

