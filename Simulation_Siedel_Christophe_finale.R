################################################################################
##################### SPHERIQUE DATA SiMULATION ################################
################################################################################
rm(list=ls());graphics.off();cat("\014")  # Initialization 

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
rm(list=ls());graphics.off();cat("\014")





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


# Fonction qui définie les coordonnées 
Coordonnée <- function(x_int, y_int)
  {
  
  # Définir les intervalles et le nombre de points
  x_interval <- seq(x_int[1], x_int[2], length.out = x_int[3])  # Intervalle pour x
  y_interval <- seq(y_int[1], y_int[2], length.out = y_int[3]) # Intervalle pour y
  
  # Créer le meshgrid
  CoordCRD <- expand.grid(x = x_interval, y = y_interval)
  #str(CoordCRD)
  # Afficher le meshgrid
  # print(CoordCRD)
  
  CoordCRD <- as.data.frame(CoordCRD)
  #remove(CoordCRD1,CoordCRD2,CoordCRD3)
  n=nrow(CoordCRD)
  
  # Valeurs pour filtrer
  filter_col1 = c(max(CoordCRD$x),  min(CoordCRD$x)) 
  filter_col2 = c(max(CoordCRD$y),  min(CoordCRD$y))
  
  # Filtrage
  Coord = CoordCRD[CoordCRD$x %in% filter_col1 | CoordCRD$y %in% filter_col2,]
  nrow(Coord)
  df <- CoordCRD[CoordCRD$x != filter_col1[1] & CoordCRD$x != filter_col1[2] & CoordCRD$y != filter_col2[1] & CoordCRD$y != filter_col2[2],]
  
  ShpPtsDat <- CoordCRD[,c(1:2)]
  colnames(ShpPtsDat) <- c("x","y")
  
  data_sf <- st_as_sf(ShpPtsDat, coords = c("x", "y")) %>%
    st_set_crs(32631)
  MD <- drop_units(max(st_distance(data_sf)))
  
  
  return(list(
    df1 = Coord ,
    df2 = df,
    df3 = CoordCRD ,
    df3 = MD
  ))
}



cat("\014") 



intermediare <- function(bord,psill,nugget,model1,spa.d,number,ratio,skwenes)
  {

  ShpPtsDat <- bord[,c(1:2)]
  colnames(ShpPtsDat) <- c("x","y")
  
  data_sf <- st_as_sf(ShpPtsDat, coords = c("x", "y")) %>%
    st_set_crs(32631)
  MD <- drop_units(max(st_distance(data_sf)))
  range = (MD/ratio)
  
  
# Simulation des données fixes pour les conditions au bords
SDISph = 0.375*((psill)/(nugget+psill))*((range)/(0.5*MD))*100;SDISph
g.dummy <- gstat(formula=z~1, locations=~x+y, dummy=T,beta=0,model=vgm(psill,model1,range,nugget), nmax=20)
    
yy <- predict(g.dummy, newdata=ShpPtsDat, nsim=number);cat("\014")
      Zy=as.data.frame(yy)
      Zyy = Zy[,-c(1:2)]+range
      skew.yy= apply(Zyy,2,skew_excess)  # not returning all values 6/10
      #kurt.yy= apply(Zyy,2,kurt_excess)    # not returning all values
      skew.yy=round(skew.yy,0);#kurt.yy=round(kurt.yy,0)
      
      #nSimSK=which(skew.yy == 0 & kurt.yy == 0)   ;nSimSK     ## Find
      #nSimSK=which(round(skew.yy,0) == skwenes & abs(round(skew.yy,1))> abs(skwenes)-0.1 & abs(round(skew.yy,1))< abs(skwenes)+0.1 )
      nSimSK=which(  abs(skew.yy - skwenes ) < 0.1 )
      nSimSK
      #length(nSimSK)
      AtSimSK= Zyy[nSimSK]
      
      return(list(
        df1 = AtSimSK ,
        df2 = SDISph ,
        df3 = MD ,
        df4 = range 
      ))
      
}

concatener_dataframe <- function(liste,row,numb)
{
    mtrx <- data.frame(NA_col = rep(NA, row)) # Crée un dataframe avec une colonne remplie de NA
    
    for(j in 1:length(liste)){
      mtrx <- cbind(mtrx, liste[[j]])
    }
    mtrx <- cbind(CoordCRD,mtrx[,c(2:(numb+1))])
    
    return(mtrx)
}


intern_donne <- function(CoordCRD,psill,nugget,model1,spa.d,Rexp,R)
  
{
    ############# LEVEL 1
    Hdata = list()
    donnee <- NULL
    
    
    for(s in 1:R){
      
      sortie=intermediare(CoordCRD,psill,nugget,model1,spa.d,Rexp,ratio,skwenes)
      Attribut=sortie[[1]] 
      SDISph12=sortie[[2]]
      MD=sortie[[3]]
      range=sortie[[4]]
      row=nrow(CoordCRD)
      data_Site <- Attribut
      
      Hdata[[s]]=data_Site
      Hdata <- compact(Hdata)
      
      # Nombre d'observations Hdata
      nHdata <- sapply(Hdata, ncol);nHdata
      sHdata <- sum(nHdata)
      
      pctg <- paste(round(sHdata/Rexp*100, 0), "% completed with")
      setTkProgressBar(pb, sHdata, label = pctg)
      
      
      if(sHdata>Rexp){
        donnee=concatener_dataframe(Hdata, row ,Rexp)
        #write_xlsx(donnee,paste(spa.d,model1,ske,"data",n,".xlsx",sep = "_"))
        break
      }
      
    }
    
    return(list(
      df1 = donnee ,
      df2 = SDISph12 ,
      df3 = MD ,
      df4 = range 
    ))

}


donneebord=function(bordcoord, dbord,choix)
{
  dbord1=dbord[choix]
  for(j in 1:Rexp-1)
  {
    dbord1 <- cbind(dbord1,dbord[choix])
  }
  donnefinbord <- cbind(bordcoord, dbord1)
  donnefinbord=donnefinbord[,c(1:(Rexp+2))]
  return(donnefinbord)
}



general=function(bordcoord,CoordCRD,psill,nugget,model1,spa.d,ratio,skwenes,Rexp,R)
    
  {
    
    
    bord=intermediare(bordcoord,psill,nugget,model1,spa.d,1000,ratio,skwenes)
    dbord=bord[[1]]
    SDISph11=bord[[2]]
    MD1=bord[[3]]
    range1=bord[[4]]
    
    
    interieur=intern_donne(CoordCRD,psill,nugget,model1,spa.d,Rexp,R)
    donnefin_inter=interieur[[1]]
    SDISph12=interieur[[2]]
    MD2=interieur[[3]]
    range2=interieur[[4]]
    
    # Séléctionner une condition au bord
    choix=1
    donnefinbord=donneebord(bordcoord, dbord , choix)
    
    names(donnefinbord)[3:(Rexp+2)] <- paste0("col", 3:(Rexp+2))
    names(donnefin_inter)[3:(Rexp+2)] <- paste0("col", 3:(Rexp+2))
    donnee_finale=rbind(donnefinbord,donnefin_inter)
    
    write_xlsx(donnee_finale,paste("C:/Users/Christophe/Documents/",spa.d,model1,ske,"data",n,".xlsx",sep = "_"))
    
    return(list(
      df1 =donnee_finale  ,
      df2 = SDISph11 ,
      df3 = SDISph12 ,
      df4 = MD1 ,
      df5 = range1 ,
      df6 = MD2 ,
      df7 = range2 
    ))
    
  }




## Simulation argument
R = 10000           # Repetition
Rexp = 1000
nbAtt=100  


x_int= c(1,10,8) 
y_int= c(1,10,8)
result <- Coordonnée(x_int,y_int)

bordcoord=result[[1]]
CoordCRD=result[[2]]
toutes_coord=result[[3]]
Distance_Max=result[[4]]
n = nrow(CoordCRD)+nrow(bordcoord)


## Variogram parameters and Spatial Dependence level
nugget11=0.9;sill11=0.1 # Weak
nugget12=0.6;sill12=0.4 # Medium
nugget13=0.1;sill13=0.9 # Strong
## Spatial variogram structure
model1="Sph"

asymetrie=c("negative","symetric","positive")
degré_asymetrie=c(-1,0,1)

numero=2
ske=asymetrie[numero]
skwenes=degré_asymetrie[numero]



cat("\014") 


########################## SPHERIQUE SIMULATION ###########################################
pb <- tkProgressBar(title = " WEAK SPHERIQUE SYMETRIC DATA SIMULATION",      # Window title
                    label = "Percentage completed", # Window label
                    min = 0,      # Minimum value of the bar
                    max = Rexp, # Maximum value of the bar
                    initial = 0,  # Initial value of the bar
                    width = 500)  # Width of the window

## Variogram parameters
spa.d <- "weak"   ## Spatial dependence
ratio=8
range_weak=Distance_Max/ratio

resulta_final_1=general(bordcoord,CoordCRD,sill11,nugget11,model1,spa.d,ratio,skwenes,Rexp,R)



########################## SPHERIQUE SIMULATION ###########################################
pb <- tkProgressBar(title = " MODERATE SPHERIQUE SYMETRIC DATA SIMULATION",      # Window title
                    label = "Percentage completed", # Window label
                    min = 0,      # Minimum value of the bar
                    max = Rexp, # Maximum value of the bar
                    initial = 0,  # Initial value of the bar
                    width = 500)  # Width of the window

## Vadiogram parameters
spa.d <- "moderate"   ## Spatial dependence
ratio=4
range_moderate=Distance_Max/ratio

resulta_final_2=general(bordcoord,CoordCRD,sill12,nugget12,model1,spa.d,ratio,skwenes,Rexp,R)



cat("\014") 


########################## SPHERIQUE SIMULATION ###########################################
pb <- tkProgressBar(title = " STRONG SPHERIQUE SYMETRIC DATA SIMULATION",      # Window title
                    label = "Percentage completed", # Window label
                    min = 0,      # Minimum value of the bar
                    max = Rexp, # Maximum value of the bar
                    initial = 0,  # Initial value of the bar
                    width = 500)  # Width of the window


## Variogram parameters
spa.d <- "strong"   ## Spatial dependence
ratio=2
range_strong=Distance_Max/ratio

resulta_final_3=general(bordcoord,CoordCRD,sill13,nugget13,model1,spa.d,ratio,skwenes,Rexp,R)

