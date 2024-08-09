# ## PREAMBULE
# SDISph = 0.375*((C1)/(C0*C1))*((a)/(0.5*MD))
# SDIExp = 0.317*((C1)/(C0*C1))*((a)/(0.5*MD))
# SDIGau = 0.504*((C1)/(C0*C1))*((a)/(0.5*MD))

rm(list=ls());graphics.off();cat("\014")  # Initialization 
setwd("C:/Users/AGON/Documents/Tâches/Dr. AGBANGBA/Article/Manucrit Dr/Good data")
R = 10000
Rexp = 5000


################################################################################
##################### SPHERIQUE DATA SOMULATION ################################
################################################################################

######################################################################
### Libraries 
library(readxl)
library(gstat)# for kriging and simulation
library(sp)   # spplot
library(e1071)# for skew_excess and kurt_excess
#library(rgdal)# for shapefiles reading 
library(writexl) # for saving and reading xls file. 
library(raster)
library(sf)
library(dplyr)
library(tidyverse)
library(units)
library(distances)
library(tcltk)
cat("\014") 


## Import coordinate
# CoordCRD <- as.data.frame(read_excel("coord_5000_CRD.xlsx",sheet = 1,col_names = T))

                            # Seed
# rm(list=ls());graphics.off();cat("\014")  # Initialization 

CoordCRD1 <- read.table("Final_coordinates_CRD_50000.txt", header = T,sep = ";")
CoordCRD2 <- CoordCRD1[,c(2:5)]
str(CoordCRD2)

library(dplyr)
# Filtrer les sites ayant exactement 96 observations
CoordCRD3 <- CoordCRD2 %>%
  group_by(Site) %>%
  filter(n() == 96) %>%
  ungroup()

# Renuméroter les sites de 1 à n
CoordCRD <- CoordCRD3 %>%
  mutate(Site = as.integer(factor(Site, levels = unique(Site))))

CoordCRD <- as.data.frame(CoordCRD)
remove(CoordCRD1,CoordCRD2,CoordCRD3)

###############################################################################
## Skewness and Kurtosis function
skew_excess=function (x, na.rm = FALSE, type = 3) 
{
  if (any(ina <- is.na(x))) {
    if (na.rm) 
      x <- x[!ina]
    else return(NA)
  }
  if (!(type %in% (1:3))) 
    stop("Invalid 'type' argument.")
  n <- length(x)
  x <- x - mean(x)
  m2 <- (1/n)*(sum(x^2))
  m3 <- (1/n)*(sum(x^3))
  y <- sqrt(n*(n-1))/(n-2) * m3/m2^(3/2)
  y
}

kurt_excess=function (x, na.rm = FALSE, type = 3) 
{
  if (any(ina <- is.na(x))) {
    if (na.rm) 
      x <- x[!ina]
    else return(NA)
  }
  if (!(type %in% (1:3))) 
    stop("Invalid 'type' argument.")
  n <- length(x)
  x <- x - mean(x)
  m2 <- (1/n)*(sum(x^2))
  m4 <- (1/n)*(sum(x^4))
  y <- (n-1)/((n-2)*(n-3))*(((n+1)*((m4/m2^2)-3))+6)
  y
}



## Simulation argument
R = 10000            # Repetition
Rexp = 5000
nbAtt=100           # Number of attribute
nbAtt1=nbAtt+1      # Number of attribute + 1
nbAtt2=nbAtt+2      # Number of attribute + 2
n = 96              # Number of coordinate per Site of 1ha (100m x 100m)

## Variogram parameters and Spatial Dependence level
nugget11=0.9;sill11=1 # Weak
nugget12=0.5;sill12=1 # Medium
nugget13=0.1;sill13=1 # Strong

## Spatial variogram structure
model1="Sph"

## Import coordinate
# CoordCRD <- as.data.frame(read_excel("coord_5000_CRD.xlsx",sheet = 1,col_names = T))
# CoordCBD <- as.data.frame(read_excel("coordonnees_points_cbd.xlsx",sheet = 1,col_names = T))
cat("\014") 


########################## SPHERIQUE SIMULATION ###########################################
pb <- tkProgressBar(title = " WEAK SPHERIQUE DATA SIMULATION",      # Window title
                    label = "Percentage completed", # Window label
                    min = 0,      # Minimum value of the bar
                    max = Rexp, # Maximum value of the bar
                    initial = 0,  # Initial value of the bar
                    width = 500)  # Width of the window


############# LEVEL 1
Hdata = list()
mtrx <- NULL


for(s in 1:R){
  Site <- filter(CoordCRD, Site == paste(s))
  ShpPtsDat <- Site[,c(4:5)]
  colnames(ShpPtsDat) <- c("x","y")
  
  data_sf <- st_as_sf(ShpPtsDat, coords = c("x", "y")) %>%
    st_set_crs(32631)
  MD <- drop_units(max(st_distance(data_sf)))
  
  ## Varigram parameters
  psill = sill11
  range = (MD/10)+5
  nugget = nugget11
  model1="Sph"      ## Model
  spa.d <- "Weak"   ## Spatial dependence
  
  SDISph = 0.375*((psill)/(nugget+psill))*((range)/(0.5*MD))*100
  
  if (SDISph < 7) {
    g.dummy1 <- gstat(formula=z~1, locations=~x+y, dummy=T,beta=0,model=vgm(psill,model1,range,nugget), nmax=20)
    
      yy <- predict(g.dummy1, newdata=ShpPtsDat, nsim=nbAtt);cat("\014")
        if(nrow(yy)==n){
          Zy=as.data.frame(yy)
          Zyy = Zy[,-c(1:2)]+range
          skew.yy= apply(Zyy,2,skew_excess)  # not returning all values 6/10
          kurt.yy= apply(Zyy,2,kurt_excess)    # not returning all values
          skew.yy=round(skew.yy,0);kurt.yy=round(kurt.yy,0)
          
          #nSimSK=which(skew.yy == 0 & kurt.yy == 0)   ;nSimSK     ## Find
          nSimSK=which(round(skew.yy,0) == 0 & abs(round(skew.yy,1))< 0.1 & 
                         round(kurt.yy,0) == 0 & abs(round(kurt.yy,1))< 0.1)
          
          AtSimSK= Zyy[nSimSK]
          Attribut <- AtSimSK[,1]
          
          data_Site <- cbind(Site,Attribut)
          colnames(data_Site) <- c("Site","Bloc","Unite","Longitude","Latitude","Attribut")

          Hdata[[s]]=data_Site
          Hdata <- compact(Hdata)
          
          pctg <- paste(round(length(Hdata)/Rexp *100, 0), "% completed with",length(Hdata),"Number of site")
          setTkProgressBar(pb, length(Hdata), label = pctg)

          if(length(Hdata)==Rexp){
            mtrx <- data.frame(NA_col = rep(NA, n)) # Crée un dataframe avec une colonne remplie de NA

            # Ajoute 7 colonnes supplémentaires remplies de NA
            for (k in 2:6) {
              mtrx[[paste0("col_", k)]] <- rep(NA, n)
            }

            for(j in 1:length(Hdata)){
              colnames(mtrx) <- c("Site","Bloc","Unite","Longitude","Latitude","Attribut")
              mtrx <- rbind(mtrx, Hdata[[j]])
              
              if(nrow(mtrx)==(((Rexp*n/n+1))*n)){
                mtrx_vf <- mtrx[c((n+1):(n*(Rexp+1))),]
                write_xlsx(mtrx_vf,paste(spa.d,model1,"CRD_data_Sites",Rexp,".xlsx",sep = "_"))
                }
              }
         break
            }
      }
  } 

}





pb <- tkProgressBar(title = " MODERATE SPHERIQUE DATA SIMULATION",      # Window title
                    label = "Percentage completed", # Window label
                    min = 0,      # Minimum value of the bar
                    max = Rexp, # Maximum value of the bar
                    initial = 0,  # Initial value of the bar
                    width = 500)  # Width of the window


############# LEVEL 1
Hdata = list()
mtrx <- NULL


for(s in 1:R){
  Site <- filter(CoordCRD, Site == paste(s))
  ShpPtsDat <- Site[,c(4:5)]
  colnames(ShpPtsDat) <- c("x","y")
  
  data_sf <- st_as_sf(ShpPtsDat, coords = c("x", "y")) %>%
    st_set_crs(32631)
  MD <- drop_units(max(st_distance(data_sf)))
  
  ## Varigram parameters
  psill = sill12
  range = (MD/8)+5
  nugget = nugget12
  model1="Sph"      ## Model
  spa.d <- "Moderate"   ## Spatial dependence
  
  SDISph = 0.375*((psill)/(nugget+psill))*((range)/(0.5*MD))*100
  
  if (SDISph > 7 & SDISph <=15) {
    
    
    g.dummy1 <- gstat(formula=z~1, locations=~x+y, dummy=T,beta=0,model=vgm(psill,model1,range,nugget), nmax=20)
    
    yy <- predict(g.dummy1, newdata=ShpPtsDat, nsim=nbAtt);cat("\014")
    if(nrow(yy)==n){
      Zy=as.data.frame(yy)
      Zyy = Zy[,-c(1:2)]+range
      skew.yy= apply(Zyy,2,skew_excess)  # not returning all values 6/10
      kurt.yy= apply(Zyy,2,kurt_excess)    # not returning all values
      skew.yy=round(skew.yy,0);kurt.yy=round(kurt.yy,0)
      
      #nSimSK=which(skew.yy == 0 & kurt.yy == 0)   ;nSimSK     ## Find
      nSimSK=which(round(skew.yy,0) == 0 & abs(round(skew.yy,1))< 0.1 & 
                     round(kurt.yy,0) == 0 & abs(round(kurt.yy,1))< 0.1)
      
      AtSimSK= Zyy[nSimSK]
      Attribut <- AtSimSK[,1]
      
      data_Site <- cbind(Site,Attribut)
      colnames(data_Site) <- c("Site","Bloc","Unite","Longitude","Latitude","Attribut")
      
      Hdata[[s]]=data_Site
      Hdata <- compact(Hdata)
      
      pctg <- paste(round(length(Hdata)/Rexp *100, 0), "% completed with",length(Hdata),"Number of site")
      setTkProgressBar(pb, length(Hdata), label = pctg)
      
        if(length(Hdata)==Rexp){
          mtrx <- data.frame(NA_col = rep(NA, n)) # Crée un dataframe avec une colonne remplie de NA
          
          # Ajoute 7 colonnes supplémentaires remplies de NA
          for (k in 2:6) {
            mtrx[[paste0("col_", k)]] <- rep(NA, n)
          }
          
          for(j in 1:length(Hdata)){
            colnames(mtrx) <- c("Site","Bloc","Unite","Longitude","Latitude","Attribut")
            mtrx <- rbind(mtrx, Hdata[[j]])
            
              if(nrow(mtrx)==(((Rexp*n/n+1))*n)){
                mtrx_vf <- mtrx[c((n+1):(n*(Rexp+1))),]
                write_xlsx(mtrx_vf,paste(spa.d,model1,"CRD_data_Sites",Rexp,".xlsx",sep = "_"))
              }
            }
          break
          }
    }
  } 
  
}




pb <- tkProgressBar(title = " STRONG SPHERIQUE DATA SIMULATION",      # Window title
                    label = "Percentage completed", # Window label
                    min = 0,      # Minimum value of the bar
                    max = Rexp, # Maximum value of the bar
                    initial = 0,  # Initial value of the bar
                    width = 500)  # Width of the window


############# LEVEL 1
Hdata = list()
mtrx <- NULL


for(s in 1:R){
  Site <- filter(CoordCRD, Site == paste(s))
  ShpPtsDat <- Site[,c(4:5)]
  colnames(ShpPtsDat) <- c("x","y")
  
  data_sf <- st_as_sf(ShpPtsDat, coords = c("x", "y")) %>%
    st_set_crs(32631)
  MD <- drop_units(max(st_distance(data_sf)))
  
  ## Varigram parameters
  psill = sill13
  range = (MD/4)+5
  nugget = nugget13
  model1="Sph"      ## Model
  spa.d <- "Strong"   ## Spatial dependence
  
  SDISph = 0.375*((psill)/(nugget+psill))*((range)/(0.5*MD))*100
  
  if (SDISph > 15) {
    
    
    g.dummy1 <- gstat(formula=z~1, locations=~x+y, dummy=T,beta=0,model=vgm(psill,model1,range,nugget), nmax=20)
    
    yy <- predict(g.dummy1, newdata=ShpPtsDat, nsim=nbAtt);cat("\014")
    if(nrow(yy)==n){
      Zy=as.data.frame(yy)
      Zyy = Zy[,-c(1:2)]+range
      skew.yy= apply(Zyy,2,skew_excess)  # not returning all values 6/10
      kurt.yy= apply(Zyy,2,kurt_excess)    # not returning all values
      skew.yy=round(skew.yy,0);kurt.yy=round(kurt.yy,0)
      
      #nSimSK=which(skew.yy == 0 & kurt.yy == 0)   ;nSimSK     ## Find
      nSimSK=which(round(skew.yy,0) == 0 & abs(round(skew.yy,1))< 0.1 & 
                     round(kurt.yy,0) == 0 & abs(round(kurt.yy,1))< 0.1)
      
      AtSimSK= Zyy[nSimSK]
      Attribut <- AtSimSK[,1]
      
      data_Site <- cbind(Site,Attribut)
      colnames(data_Site) <- c("Site","Bloc","Unite","Longitude","Latitude","Attribut")
      
      Hdata[[s]]=data_Site
      Hdata <- compact(Hdata)
      
      pctg <- paste(round(length(Hdata)/Rexp *100, 0), "% completed with",length(Hdata),"Number of site")
      setTkProgressBar(pb, length(Hdata), label = pctg)
      
      if(length(Hdata)==Rexp){
        mtrx <- data.frame(NA_col = rep(NA, n)) # Crée un dataframe avec une colonne remplie de NA
        
        # Ajoute 7 colonnes supplémentaires remplies de NA
        for (k in 2:6) {
          mtrx[[paste0("col_", k)]] <- rep(NA, n)
        }
        
        for(j in 1:length(Hdata)){
          colnames(mtrx) <- c("Site","Bloc","Unite","Longitude","Latitude","Attribut")
          mtrx <- rbind(mtrx, Hdata[[j]])
          
          if(nrow(mtrx)==(((Rexp*n/n+1))*n)){
            mtrx_vf <- mtrx[c((n+1):(n*(Rexp+1))),]
            write_xlsx(mtrx_vf,paste(spa.d,model1,"CRD_data_Sites",Rexp,".xlsx",sep = "_"))

          }
        }
        break
      }
    }
  } 
  
}












################################################################################
##################### EXPONENTIAL DATA SOMULATION ##############################
################################################################################

set.seed(2024)                            # Seed
tokeep <- c("CoordCRD")
rm(list=setdiff(ls(), tokeep))

###############################################################################
## Skewness and Kurtosis function
skew_excess=function (x, na.rm = FALSE, type = 3) 
{
  if (any(ina <- is.na(x))) {
    if (na.rm) 
      x <- x[!ina]
    else return(NA)
  }
  if (!(type %in% (1:3))) 
    stop("Invalid 'type' argument.")
  n <- length(x)
  x <- x - mean(x)
  m2 <- (1/n)*(sum(x^2))
  m3 <- (1/n)*(sum(x^3))
  y <- sqrt(n*(n-1))/(n-2) * m3/m2^(3/2)
  y
}

kurt_excess=function (x, na.rm = FALSE, type = 3) 
{
  if (any(ina <- is.na(x))) {
    if (na.rm) 
      x <- x[!ina]
    else return(NA)
  }
  if (!(type %in% (1:3))) 
    stop("Invalid 'type' argument.")
  n <- length(x)
  x <- x - mean(x)
  m2 <- (1/n)*(sum(x^2))
  m4 <- (1/n)*(sum(x^4))
  y <- (n-1)/((n-2)*(n-3))*(((n+1)*((m4/m2^2)-3))+6)
  y
}



## Simulation argument
R = 10000            # Repetition
Rexp = 5000
nbAtt=100           # Number of attribute
nbAtt1=nbAtt+1      # Number of attribute + 1
nbAtt2=nbAtt+2      # Number of attribute + 2
n = 96              # Number of coordinate per Site of 1ha (100m x 100m)

## Variogram parameters and Spatial Dependence level
nugget11=0.9;sill11=1 # Weak
nugget12=0.5;sill12=1 # Medium
nugget13=0.1;sill13=1 # Strong

## Spatial variogram structure
model1="Exp"




########################## EXPONENTIAL SIMULATION ###########################################
pb <- tkProgressBar(title = " WEAK EXPONENTIAL DATA SIMULATION",      # Window title
                    label = "Percentage completed", # Window label
                    min = 0,      # Minimum value of the bar
                    max = Rexp, # Maximum value of the bar
                    initial = 0,  # Initial value of the bar
                    width = 500)  # Width of the window


############# LEVEL 1
Hdata = list()
mtrx <- NULL


for(s in 1:R){
  Site <- filter(CoordCRD, Site == paste(s))
  ShpPtsDat <- Site[,c(4:5)]
  colnames(ShpPtsDat) <- c("x","y")
  
  data_sf <- st_as_sf(ShpPtsDat, coords = c("x", "y")) %>%
    st_set_crs(32631)
  MD <- drop_units(max(st_distance(data_sf)))
  
  ## Varigram parameters
  psill = sill11
  range = (MD/10)+5
  nugget = nugget11
  model1="Exp"      ## Model
  spa.d <- "Weak"   ## Spatial dependence
  
  
  SDIExp = 0.317*((psill)/(nugget+psill))*((range)/(0.5*MD))*100
  
  if (SDIExp < 7) {
    
    
    g.dummy1 <- gstat(formula=z~1, locations=~x+y, dummy=T,beta=0,model=vgm(psill,model1,range,nugget), nmax=20)
    
    yy <- predict(g.dummy1, newdata=ShpPtsDat, nsim=nbAtt);cat("\014")
    if(nrow(yy)==n){
      Zy=as.data.frame(yy)
      Zyy = Zy[,-c(1:2)]+range
      skew.yy= apply(Zyy,2,skew_excess)  # not returning all values 6/10
      kurt.yy= apply(Zyy,2,kurt_excess)    # not returning all values
      skew.yy=round(skew.yy,0);kurt.yy=round(kurt.yy,0)
      
      #nSimSK=which(skew.yy == 0 & kurt.yy == 0)   ;nSimSK     ## Find
      nSimSK=which(round(skew.yy,0) == 0 & abs(round(skew.yy,1))< 0.1 & 
                     round(kurt.yy,0) == 0 & abs(round(kurt.yy,1))< 0.1)
      
      AtSimSK= Zyy[nSimSK]
      Attribut <- AtSimSK[,1]
      
      data_Site <- cbind(Site,Attribut)
      colnames(data_Site) <- c("Site","Bloc","Unite","Longitude","Latitude","Attribut")
      
      Hdata[[s]]=data_Site
      Hdata <- compact(Hdata)
      
      pctg <- paste(round(length(Hdata)/Rexp *100, 0), "% completed with",length(Hdata),"Number of site")
      setTkProgressBar(pb, length(Hdata), label = pctg)
      
      if(length(Hdata)==Rexp){
        mtrx <- data.frame(NA_col = rep(NA, n)) # Crée un dataframe avec une colonne remplie de NA
        
        # Ajoute 7 colonnes supplémentaires remplies de NA
        for (k in 2:6) {
          mtrx[[paste0("col_", k)]] <- rep(NA, n)
        }
        
        for(j in 1:length(Hdata)){
          colnames(mtrx) <- c("Site","Bloc","Unite","Longitude","Latitude","Attribut")
          mtrx <- rbind(mtrx, Hdata[[j]])
          
          if(nrow(mtrx)==(((Rexp*n/n+1))*n)){
            mtrx_vf <- mtrx[c((n+1):(n*(Rexp+1))),]
            write_xlsx(mtrx_vf,paste(spa.d,model1,"CRD_data_Sites",Rexp,".xlsx",sep = "_"))

          }
        }
        break
      }
    }
  } 
  
}







pb <- tkProgressBar(title = " MODERATE EXPONENTIAL DATA SIMULATION",      # Window title
                    label = "Percentage completed", # Window label
                    min = 0,      # Minimum value of the bar
                    max = Rexp, # Maximum value of the bar
                    initial = 0,  # Initial value of the bar
                    width = 500)  # Width of the window


############# LEVEL 1
Hdata = list()
mtrx <- NULL


for(s in 1:R){
  Site <- filter(CoordCRD, Site == paste(s))
  ShpPtsDat <- Site[,c(4:5)]
  colnames(ShpPtsDat) <- c("x","y")
  
  data_sf <- st_as_sf(ShpPtsDat, coords = c("x", "y")) %>%
    st_set_crs(32631)
  MD <- drop_units(max(st_distance(data_sf)))
  
  ## Varigram parameters
  psill = sill12
  range = (MD/8)+5
  nugget = nugget12
  model1="Exp"      ## Model
  spa.d <- "Moderate"   ## Spatial dependence
  
  SDIExp = 0.317*((psill)/(nugget+psill))*((range)/(0.5*MD))*100
  
  if (SDIExp > 7 & SDIExp <=15) {
    
    
    g.dummy1 <- gstat(formula=z~1, locations=~x+y, dummy=T,beta=0,model=vgm(psill,model1,range,nugget), nmax=20)
    
    yy <- predict(g.dummy1, newdata=ShpPtsDat, nsim=nbAtt);cat("\014")
    if(nrow(yy)==n){
      Zy=as.data.frame(yy)
      Zyy = Zy[,-c(1:2)]+range
      skew.yy= apply(Zyy,2,skew_excess)  # not returning all values 6/10
      kurt.yy= apply(Zyy,2,kurt_excess)    # not returning all values
      skew.yy=round(skew.yy,0);kurt.yy=round(kurt.yy,0)
      
      #nSimSK=which(skew.yy == 0 & kurt.yy == 0)   ;nSimSK     ## Find
      nSimSK=which(round(skew.yy,0) == 0 & abs(round(skew.yy,1))< 0.1 & 
                     round(kurt.yy,0) == 0 & abs(round(kurt.yy,1))< 0.1)
      
      AtSimSK= Zyy[nSimSK]
      Attribut <- AtSimSK[,1]
      
      data_Site <- cbind(Site,Attribut)
      colnames(data_Site) <- c("Site","Bloc","Unite","Longitude","Latitude","Attribut")
      
      Hdata[[s]]=data_Site
      Hdata <- compact(Hdata)
      
      pctg <- paste(round(length(Hdata)/Rexp *100, 0), "% completed with",length(Hdata),"Number of site")
      setTkProgressBar(pb, length(Hdata), label = pctg)
      
      if(length(Hdata)==Rexp){
        mtrx <- data.frame(NA_col = rep(NA, n)) # Crée un dataframe avec une colonne remplie de NA
        
        # Ajoute 7 colonnes supplémentaires remplies de NA
        for (k in 2:6) {
          mtrx[[paste0("col_", k)]] <- rep(NA, n)
        }
        
        for(j in 1:length(Hdata)){
          colnames(mtrx) <- c("Site","Bloc","Unite","Longitude","Latitude","Attribut")
          mtrx <- rbind(mtrx, Hdata[[j]])
          
          if(nrow(mtrx)==(((Rexp*n/n+1))*n)){
            mtrx_vf <- mtrx[c((n+1):(n*(Rexp+1))),]
            write_xlsx(mtrx_vf,paste(spa.d,model1,"CRD_data_Sites",Rexp,".xlsx",sep = "_"))

          }
        }
        break
      }
    }
  } 
  
}





pb <- tkProgressBar(title = " STRONG EXPONENTIAL DATA SIMULATION",      # Window title
                    label = "Percentage completed", # Window label
                    min = 0,      # Minimum value of the bar
                    max = Rexp, # Maximum value of the bar
                    initial = 0,  # Initial value of the bar
                    width = 500)  # Width of the window


############# LEVEL 1
Hdata = list()
mtrx <- NULL


for(s in 1:R){
  Site <- filter(CoordCRD, Site == paste(s))
  ShpPtsDat <- Site[,c(4:5)]
  colnames(ShpPtsDat) <- c("x","y")
  
  data_sf <- st_as_sf(ShpPtsDat, coords = c("x", "y")) %>%
    st_set_crs(32631)
  MD <- drop_units(max(st_distance(data_sf)))
  
  ## Varigram parameters
  psill = sill13
  range = (MD/4)+5
  nugget = nugget13
  model1="Exp"      ## Model
  spa.d <- "Strong"   ## Spatial dependence
  
  SDIExp = 0.317*((psill)/(nugget+psill))*((range)/(0.5*MD))*100
  
  if (SDIExp > 15) {
    
    
    g.dummy1 <- gstat(formula=z~1, locations=~x+y, dummy=T,beta=0,model=vgm(psill,model1,range,nugget), nmax=20)
    
    yy <- predict(g.dummy1, newdata=ShpPtsDat, nsim=nbAtt);cat("\014")
    if(nrow(yy)==n){
      Zy=as.data.frame(yy)
      Zyy = Zy[,-c(1:2)]+range
      skew.yy= apply(Zyy,2,skew_excess)  # not returning all values 6/10
      kurt.yy= apply(Zyy,2,kurt_excess)    # not returning all values
      skew.yy=round(skew.yy,0);kurt.yy=round(kurt.yy,0)
      
      #nSimSK=which(skew.yy == 0 & kurt.yy == 0)   ;nSimSK     ## Find
      nSimSK=which(round(skew.yy,0) == 0 & abs(round(skew.yy,1))< 0.1 & 
                     round(kurt.yy,0) == 0 & abs(round(kurt.yy,1))< 0.1)
      
      AtSimSK= Zyy[nSimSK]
      Attribut <- AtSimSK[,1]
      
      data_Site <- cbind(Site,Attribut)
      colnames(data_Site) <- c("Site","Bloc","Unite","Longitude","Latitude","Attribut")
      
      Hdata[[s]]=data_Site
      Hdata <- compact(Hdata)
      
      pctg <- paste(round(length(Hdata)/Rexp *100, 0), "% completed with",length(Hdata),"Number of site")
      setTkProgressBar(pb, length(Hdata), label = pctg)
      
      if(length(Hdata)==Rexp){
        mtrx <- data.frame(NA_col = rep(NA, n)) # Crée un dataframe avec une colonne remplie de NA
        
        # Ajoute 7 colonnes supplémentaires remplies de NA
        for (k in 2:6) {
          mtrx[[paste0("col_", k)]] <- rep(NA, n)
        }
        
        for(j in 1:length(Hdata)){
          colnames(mtrx) <- c("Site","Bloc","Unite","Longitude","Latitude","Attribut")
          mtrx <- rbind(mtrx, Hdata[[j]])
          
          if(nrow(mtrx)==(((Rexp*n/n+1))*n)){
            mtrx_vf <- mtrx[c((n+1):(n*(Rexp+1))),]
            write_xlsx(mtrx_vf,paste(spa.d,model1,"CRD_data_Sites",Rexp,".xlsx",sep = "_"))

          }
        }
        break
      }
    }
  } 
  
}









################################################################################
################### GAUSSIAN DATA SOMULATION ###################################
################################################################################

set.seed(2024)                            # Seed
tokeep <- c("CoordCRD")
rm(list=setdiff(ls(), tokeep))


###############################################################################
###############################################################################
## Skewness and Kurtosis function
skew_excess=function (x, na.rm = FALSE, type = 3) 
{
  if (any(ina <- is.na(x))) {
    if (na.rm) 
      x <- x[!ina]
    else return(NA)
  }
  if (!(type %in% (1:3))) 
    stop("Invalid 'type' argument.")
  n <- length(x)
  x <- x - mean(x)
  m2 <- (1/n)*(sum(x^2))
  m3 <- (1/n)*(sum(x^3))
  y <- sqrt(n*(n-1))/(n-2) * m3/m2^(3/2)
  y
}

kurt_excess=function (x, na.rm = FALSE, type = 3) 
{
  if (any(ina <- is.na(x))) {
    if (na.rm) 
      x <- x[!ina]
    else return(NA)
  }
  if (!(type %in% (1:3))) 
    stop("Invalid 'type' argument.")
  n <- length(x)
  x <- x - mean(x)
  m2 <- (1/n)*(sum(x^2))
  m4 <- (1/n)*(sum(x^4))
  y <- (n-1)/((n-2)*(n-3))*(((n+1)*((m4/m2^2)-3))+6)
  y
}



## Simulation argument
R = 10000            # Repetition
Rexp = 5000
nbAtt=100           # Number of attribute
nbAtt1=nbAtt+1      # Number of attribute + 1
nbAtt2=nbAtt+2      # Number of attribute + 2
n = 96              # Number of coordinate per Site of 1ha (100m x 100m)

## Variogram parameters and Spatial Dependence level
nugget11=0.9;sill11=1 # Weak
nugget12=0.5;sill12=1 # Medium
nugget13=0.1;sill13=1 # Strong

## Spatial variogram structure
model1="Gau"




########################## GAUSSIAN SIMULATION ###########################################
pb <- tkProgressBar(title = " WEAK GAUSSIAN DATA SIMULATION",      # Window title
                    label = "Percentage completed", # Window label
                    min = 0,      # Minimum value of the bar
                    max = Rexp, # Maximum value of the bar
                    initial = 0,  # Initial value of the bar
                    width = 500)  # Width of the window


############# LEVEL 1
Hdata = list()
mtrx <- NULL


for(s in 1:R){
  Site <- filter(CoordCRD, Site == paste(s))
  ShpPtsDat <- Site[,c(4:5)]
  colnames(ShpPtsDat) <- c("x","y")
  
  data_sf <- st_as_sf(ShpPtsDat, coords = c("x", "y")) %>%
    st_set_crs(32631)
  MD <- drop_units(max(st_distance(data_sf)))
  
  ## Varigram parameters
  psill = sill11
  range = (MD/10)+2
  nugget = nugget11
  model1="Gau"      ## Model
  spa.d <- "Weak"   ## Spatial dependence
  
  
  SDIGau = 0.504*((psill)/(nugget+psill))*((range)/(0.5*MD))*100
  
  if (SDIGau < 7) {
    
    
    g.dummy1 <- gstat(formula=z~1, locations=~x+y, dummy=T,beta=0,model=vgm(psill,model1,range,nugget), nmax=20)
    
    yy <- predict(g.dummy1, newdata=ShpPtsDat, nsim=nbAtt);cat("\014")
    if(nrow(yy)==n){
      Zy=as.data.frame(yy)
      Zyy = Zy[,-c(1:2)]+range
      skew.yy= apply(Zyy,2,skew_excess)  # not returning all values 6/10
      kurt.yy= apply(Zyy,2,kurt_excess)    # not returning all values
      skew.yy=round(skew.yy,0);kurt.yy=round(kurt.yy,0)
      
      #nSimSK=which(skew.yy == 0 & kurt.yy == 0)   ;nSimSK     ## Find
      nSimSK=which(round(skew.yy,0) == 0 & abs(round(skew.yy,1))< 0.1 & 
                     round(kurt.yy,0) == 0 & abs(round(kurt.yy,1))< 0.1)
      
      AtSimSK= Zyy[nSimSK]
      Attribut <- AtSimSK[,1]
      
      data_Site <- cbind(Site,Attribut)
      colnames(data_Site) <- c("Site","Bloc","Unite","Longitude","Latitude","Attribut")
      
      Hdata[[s]]=data_Site
      Hdata <- compact(Hdata)
      
      pctg <- paste(round(length(Hdata)/Rexp *100, 0), "% completed with",length(Hdata),"Number of site")
      setTkProgressBar(pb, length(Hdata), label = pctg)
      
      if(length(Hdata)==Rexp){
        mtrx <- data.frame(NA_col = rep(NA, n)) # Crée un dataframe avec une colonne remplie de NA
        
        # Ajoute 7 colonnes supplémentaires remplies de NA
        for (k in 2:6) {
          mtrx[[paste0("col_", k)]] <- rep(NA, n)
        }
        
        for(j in 1:length(Hdata)){
          colnames(mtrx) <- c("Site","Bloc","Unite","Longitude","Latitude","Attribut")
          mtrx <- rbind(mtrx, Hdata[[j]])
          
          if(nrow(mtrx)==(((Rexp*n/n+1))*n)){
            mtrx_vf <- mtrx[c((n+1):(n*(Rexp+1))),]
            write_xlsx(mtrx_vf,paste(spa.d,model1,"CRD_data_Sites",Rexp,".xlsx",sep = "_"))

          }
        }
        break
      }
    }
  } 
  
}






pb <- tkProgressBar(title = " MODERATE GAUSSIAN DATA SIMULATION",      # Window title
                    label = "Percentage completed", # Window label
                    min = 0,      # Minimum value of the bar
                    max = Rexp, # Maximum value of the bar
                    initial = 0,  # Initial value of the bar
                    width = 500)  # Width of the window


############# LEVEL 1
Hdata = list()
mtrx <- NULL


for(s in 1:R){
  Site <- filter(CoordCRD, Site == paste(s))
  ShpPtsDat <- Site[,c(4:5)]
  colnames(ShpPtsDat) <- c("x","y")
  
  data_sf <- st_as_sf(ShpPtsDat, coords = c("x", "y")) %>%
    st_set_crs(32631)
  MD <- drop_units(max(st_distance(data_sf)))
  
  ## Varigram parameters
  psill = sill12
  range = (MD/8)+2
  nugget = nugget12
  model1="Gau"      ## Model
  spa.d <- "Moderate"   ## Spatial dependence
  
  SDIGau = 0.504*((psill)/(nugget+psill))*((range)/(0.5*MD))*100
  
  if (SDIGau > 7 & SDIGau <=15) {
    
    
    g.dummy1 <- gstat(formula=z~1, locations=~x+y, dummy=T,beta=0,model=vgm(psill,model1,range,nugget), nmax=20)
    
    yy <- predict(g.dummy1, newdata=ShpPtsDat, nsim=nbAtt);cat("\014")
    if(nrow(yy)==n){
      Zy=as.data.frame(yy)
      Zyy = Zy[,-c(1:2)]+range
      skew.yy= apply(Zyy,2,skew_excess)  # not returning all values 6/10
      kurt.yy= apply(Zyy,2,kurt_excess)    # not returning all values
      skew.yy=round(skew.yy,0);kurt.yy=round(kurt.yy,0)
      
      #nSimSK=which(skew.yy == 0 & kurt.yy == 0)   ;nSimSK     ## Find
      nSimSK=which(round(skew.yy,0) == 0 & abs(round(skew.yy,1))< 0.1 & 
                     round(kurt.yy,0) == 0 & abs(round(kurt.yy,1))< 0.1)
      
      AtSimSK= Zyy[nSimSK]
      Attribut <- AtSimSK[,1]
      
      data_Site <- cbind(Site,Attribut)
      colnames(data_Site) <- c("Site","Bloc","Unite","Longitude","Latitude","Attribut")
      
      Hdata[[s]]=data_Site
      Hdata <- compact(Hdata)
      
      pctg <- paste(round(length(Hdata)/Rexp *100, 0), "% completed with",length(Hdata),"Number of site")
      setTkProgressBar(pb, length(Hdata), label = pctg)
      
      if(length(Hdata)==Rexp){
        mtrx <- data.frame(NA_col = rep(NA, n)) # Crée un dataframe avec une colonne remplie de NA
        
        # Ajoute 7 colonnes supplémentaires remplies de NA
        for (k in 2:6) {
          mtrx[[paste0("col_", k)]] <- rep(NA, n)
        }
        
        for(j in 1:length(Hdata)){
          colnames(mtrx) <- c("Site","Bloc","Unite","Longitude","Latitude","Attribut")
          mtrx <- rbind(mtrx, Hdata[[j]])
          
          if(nrow(mtrx)==(((Rexp*n/n+1))*n)){
            mtrx_vf <- mtrx[c((n+1):(n*(Rexp+1))),]
            write_xlsx(mtrx_vf,paste(spa.d,model1,"CRD_data_Sites",Rexp,".xlsx",sep = "_"))

          }
        }
        break
      }
    }
  } 
  
}





pb <- tkProgressBar(title = " STRONG GAUSSIAN DATA SIMULATION",      # Window title
                    label = "Percentage completed", # Window label
                    min = 0,      # Minimum value of the bar
                    max = Rexp, # Maximum value of the bar
                    initial = 0,  # Initial value of the bar
                    width = 500)  # Width of the window


############# LEVEL 1
Hdata = list()
mtrx <- NULL


for(s in 1:R){
  Site <- filter(CoordCRD, Site == paste(s))
  ShpPtsDat <- Site[,c(4:5)]
  colnames(ShpPtsDat) <- c("x","y")
  
  data_sf <- st_as_sf(ShpPtsDat, coords = c("x", "y")) %>%
    st_set_crs(32631)
  MD <- drop_units(max(st_distance(data_sf)))
  
  ## Varigram parameters
  psill = sill13
  range = (MD/4)+5
  nugget = nugget13
  model1="Gau"      ## Model
  spa.d <- "Strong"   ## Spatial dependence
  
  SDIGau = 0.504*((psill)/(nugget+psill))*((range)/(0.5*MD))*100
  
  if (SDIGau > 15) {
    
    
    g.dummy1 <- gstat(formula=z~1, locations=~x+y, dummy=T,beta=0,model=vgm(psill,model1,range,nugget), nmax=20)
    
    yy <- predict(g.dummy1, newdata=ShpPtsDat, nsim=nbAtt);cat("\014")
    if(nrow(yy)==n){
      Zy=as.data.frame(yy)
      Zyy = Zy[,-c(1:2)]+range
      skew.yy= apply(Zyy,2,skew_excess)  # not returning all values 6/10
      kurt.yy= apply(Zyy,2,kurt_excess)    # not returning all values
      skew.yy=round(skew.yy,0);kurt.yy=round(kurt.yy,0)
      
      #nSimSK=which(skew.yy == 0 & kurt.yy == 0)   ;nSimSK     ## Find
      nSimSK=which(round(skew.yy,0) == 0 & abs(round(skew.yy,1))< 0.1 & 
                     round(kurt.yy,0) == 0 & abs(round(kurt.yy,1))< 0.1)
      
      AtSimSK= Zyy[nSimSK]
      Attribut <- AtSimSK[,1]
      
      data_Site <- cbind(Site,Attribut)
      colnames(data_Site) <- c("Site","Bloc","Unite","Longitude","Latitude","Attribut")
      
      Hdata[[s]]=data_Site
      Hdata <- compact(Hdata)
      
      pctg <- paste(round(length(Hdata)/Rexp *100, 0), "% completed with",length(Hdata),"Number of site")
      setTkProgressBar(pb, length(Hdata), label = pctg)
      
      if(length(Hdata)==Rexp){
        mtrx <- data.frame(NA_col = rep(NA, n)) # Crée un dataframe avec une colonne remplie de NA
        
        # Ajoute 7 colonnes supplémentaires remplies de NA
        for (k in 2:6) {
          mtrx[[paste0("col_", k)]] <- rep(NA, n)
        }
        
        for(j in 1:length(Hdata)){
          colnames(mtrx) <- c("Site","Bloc","Unite","Longitude","Latitude","Attribut")
          mtrx <- rbind(mtrx, Hdata[[j]])
          
          if(nrow(mtrx)==(((Rexp*n/n+1))*n)){
            mtrx_vf <- mtrx[c((n+1):(n*(Rexp+1))),]
            write_xlsx(mtrx_vf,paste(spa.d,model1,"CRD_data_Sites",Rexp,".xlsx",sep = "_"))

          }
        }
        break
      }
    }
  } 
  
}


tokeep <- c("CoordCRD")
rm(list=setdiff(ls(), tokeep))
