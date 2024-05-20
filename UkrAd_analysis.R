

# Author:    Ruairidh Macleod
# Copyright: Copyright 2024, University of Cambridge
# Email:     r@palaeome.org
# License:   MIT


#####this way up#####

require(dplyr)
require(sp)
require(stringr)
require(ggplot2)

##Set the following to your machine’s file directory details for the two files ‘ukrad_data.csv’ (see paper SI) and data.evec (eigenvectors output from smartpca, will add to github) 

meta1240K <- read.csv('https://reichdata.hms.harvard.edu/pub/datasets/amh_repo/curated_releases/V54/V54.1.p1/SHARE/public.dir/v54.1.p1_HO_public.anno', sep="\t")
data_evec <- read_table("Downloads/data.evec", 
                         +     col_names = FALSE)
data_evec$Master_ID <- str_match(data_evec$X1, ".*:(.*)")[,2]
ukrad_data <- read.csv(“Downloads/ukkrad_data.csv”)

names(meta1240K) [8] <- "meandate"
names(meta1240K) [7] <- "datemethod"
names(meta1240K) [9] <- "dateSD"
names(meta1240K) [5] <- "yearpublished"
names(meta1240K) [10] <- "FullDate"
meta1240K <- meta1240K[!duplicated(meta1240K$`Master ID`), ] ##remove duplicates under genetic ID

####to run analysis using W. Eusrasion AADR subset in paper
data_evec$Master_ID <- str_match(data_evec$X1, ".*:(.*)")[,2]
meta1240K <- meta1240K[match(data_evec$Master_ID, meta1240K$`Master ID`), nomatch=0]

PSEI <- subset(meta1240K, meta1240K$meandate <5000 & meta1240K$meandate >2700)
SEIA <- subset(meta1240K, meta1240K$meandate <2699 & meta1240K$meandate >2300)
IAEM <- subset(meta1240K, meta1240K$meandate <2299 & meta1240K$meandate >1100)
LMEM <- subset(meta1240K, meta1240K$meandate <1099 & meta1240K$meandate >200)
rm(meta1240K) ##quite a big table to leave lying around

PSEI$Lat. <- as.numeric(PSEI$Lat.)
PSEI$Long. <- as.numeric(PSEI$Long.)
PSEI <- PSEI[!is.na(PSEI$Lat.),]
SEIA$Lat. <- as.numeric(SEIA$Lat.)
SEIA$Long. <- as.numeric(SEIA$Long.)
SEIA <- SEIA[!is.na(SEIA$Lat.),]
IAEM$Lat. <- as.numeric(IAEM$Lat.)
IAEM$Long. <- as.numeric(IAEM$Long.)
IAEM <- IAEM[!is.na(IAEM$Lat.),]
LMEM$Lat. <- as.numeric(LMEM$Lat.)
LMEM$Long. <- as.numeric(LMEM$Long.)
LMEM <- LMEM[!is.na(LMEM$Lat.),]


#####get reference area data for dategroups

dategroups <- c("PSEI", "SEIA", "IAEM", "LMEM") 

dist_haversine <- function(lat1, lon1, lat2, lon2) {
  dlat <- (lat2 - lat1) * pi / 180
  dlon <- (lon2 - lon1) * pi / 180
  a <- sin(dlat/2)^2 + cos(lat1 * pi / 180) * cos(lat2 * pi / 180) * sin(dlon/2)^2
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  R <- 6371 # Radius of Earth in km
  return(R * c)
}

get_smallest_enclosing_circle <- function(latitudes, longitudes) {
  min_lat <- min(latitudes)
  max_lat <- max(latitudes)
  min_long <- min(longitudes)
  max_long <- max(longitudes)
  midpoint_lat <- (min_lat + max_lat) / 2
  midpoint_long <- (min_long + max_long) / 2
  distances <- sqrt((latitudes - midpoint_lat)^2 + (longitudes - midpoint_long)^2)
  distances_km <- dist_haversine(midpoint_lat, midpoint_long, latitudes, longitudes)
  radius <- max(distances_km)
  return(list(midpoint_lat, midpoint_long, radius))
}

refdata = data.frame()
for(dategroup in dategroups) {
  df <- unique(cbind(as.numeric(ukrad_data$Latitude[ukrad_data$dategroup==dategroup]), as.numeric(ukrad_data$Longitude[ukrad_data$dategroup==dategroup])))
  circ <- get_smallest_enclosing_circle(df[,1], df[,2])
  output <- as.vector(c(dategroup, circ[[1]], circ[[2]], circ[[3]]))
  refdata=rbind(refdata, output)
}
colnames(refdata) <- c("dataset", "midpoint_lat", "midpoint_long", "radius_km")

#####get random coordinate sets and corresponding sampleIDs from AADR

##  Sample Ns below; test sample data range randomly downsampled to same N.
## IAEM LMEM PSEI SEIA 
##  26   19   15   27  

dist_haversine <- function(lat1, lon1, lat2, lon2) {
  dlat <- (lat2 - lat1) * pi / 180
  dlon <- (lon2 - lon1) * pi / 180
  a <- sin(dlat/2)^2 + cos(lat1 * pi / 180) * cos(lat2 * pi / 180) * sin(dlon/2)^2
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  R <- 6371 # Radius of Earth in km
  return(R * c)
}

filter_coordinates_by_radius <- function(dataset, midpoint_lat, midpoint_long, radius_km) {
  distances <- dist_haversine(as.vector(rep(midpoint_lat, times=length(dataset$Lat.))), as.vector(rep(midpoint_long, times=length(dataset$Lat.))), dataset$Lat., dataset$Long.)
  filtered_data <- dataset[distances <= radius_km, ]
  return(filtered_data)
}


PSEI_testdata = data.frame() ##for loop uses this to make a big dataframe with all sample IDs within the same radius of the 100 random PSEI samples
random100 <- PSEI[sample(nrow(PSEI), 100, replace = FALSE), ]
for(i in 1:nrow(random100)) {
  midpoint_lat <- random100$Lat.[i]
  midpoint_long <- random100$Long.[i]
  dg_radius <- as.numeric(refdata$radius_km[refdata$dataset=="PSEI"])
  filtered_coords <- filter_coordinates_by_radius(PSEI, midpoint_lat, midpoint_long, dg_radius)
  filtered_data <- filtered_coords[sample(nrow(filtered_coords), 15, replace = TRUE), ] ## downsample to UkrAd DateGroup N 
  filtered_data <- filtered_data[!duplicated(filtered_data$`Genetic ID`), ] ## remove duplicates from previous replacement being true
  output <- as.matrix(cbind(as.vector(rep(paste("PSEI_",i,sep=""), times=length(filtered_data$`Genetic ID`))), filtered_data$`Genetic ID`, filtered_data$`Master ID`, filtered_data$meandate, filtered_data$`Group ID`, filtered_data$Lat., filtered_data$Long.))
  print(filtered_data)
  PSEI_testdata=rbind(PSEI_testdata, output)
}
colnames(PSEI_testdata) <- c("iteration", "Genetic_ID", "Master_ID", "meandate", "group_ID", "Lat.", "Long.")
PSEI_testdata <- PSEI_testdata %>% group_by(iteration) %>% filter(n() >= 10) ##remove any with N < 10

SEIA_testdata = data.frame() ##for loop uses this to make a big dataframe with all sample IDs within the same radius of the 100 random SEIA samples
random100 <- SEIA[sample(nrow(SEIA), 100, replace = FALSE), ]
for(i in 1:nrow(random100)) {
  midpoint_lat <- random100$Lat.[i]
  midpoint_long <- random100$Long.[i]
  dg_radius <- as.numeric(refdata$radius_km[refdata$dataset=="SEIA"])
  filtered_coords <- filter_coordinates_by_radius(SEIA, midpoint_lat, midpoint_long, dg_radius)
  filtered_data <- filtered_coords[sample(nrow(filtered_coords), 27, replace = TRUE), ] ## downsample to UkrAd DateGroup N 
  filtered_data <- filtered_data[!duplicated(filtered_data$`Genetic ID`), ] ## remove duplicates from previous replacement being true
  output <- as.matrix(cbind(as.vector(rep(paste("SEIA_",i,sep=""), times=length(filtered_data$`Genetic ID`))), filtered_data$`Genetic ID`, filtered_data$`Master ID`, filtered_data$meandate, filtered_data$`Group ID`, filtered_data$Lat., filtered_data$Long.))
  print(filtered_data)
  SEIA_testdata=rbind(SEIA_testdata, output)
}
colnames(SEIA_testdata) <- c("iteration", "Genetic_ID", "Master_ID", "meandate", "group_ID", "Lat.", "Long.")
SEIA_testdata <- SEIA_testdata %>% group_by(iteration) %>% filter(n() >= 10) ##remove any with N < 10

IAEM_testdata = data.frame() ##for loop uses this to make a big dataframe with all sample IDs within the same radius of the 100 random IAEM samples
random100 <- IAEM[sample(nrow(IAEM), 100, replace = FALSE), ]
for(i in 1:nrow(random100)) {
  midpoint_lat <- random100$Lat.[i]
  midpoint_long <- random100$Long.[i]
  dg_radius <- as.numeric(refdata$radius_km[refdata$dataset=="IAEM"])
  filtered_coords <- filter_coordinates_by_radius(IAEM, midpoint_lat, midpoint_long, dg_radius)
  filtered_data <- filtered_coords[sample(nrow(filtered_coords), 26, replace = TRUE), ] ## downsample to UkrAd DateGroup N 
  filtered_data <- filtered_data[!duplicated(filtered_data$`Genetic ID`), ] ## remove duplicates from previous replacement being true
  output <- as.matrix(cbind(as.vector(rep(paste("IAEM_",i,sep=""), times=length(filtered_data$`Genetic ID`))), filtered_data$`Genetic ID`, filtered_data$`Master ID`, filtered_data$meandate, filtered_data$`Group ID`, filtered_data$Lat., filtered_data$Long.))
  print(filtered_data)
  IAEM_testdata=rbind(IAEM_testdata, output)
}
colnames(IAEM_testdata) <- c("iteration", "Genetic_ID", "Master_ID", "meandate", "group_ID", "Lat.", "Long.")
IAEM_testdata <- IAEM_testdata %>% group_by(iteration) %>% filter(n() >= 10) ##remove any with N < 10

LMEM_testdata = data.frame() ##for loop uses this to make a big dataframe with all sample IDs within the same radius of the 100 random LMEM samples
random100 <- LMEM[sample(nrow(LMEM), 100, replace = FALSE), ]
for(i in 1:nrow(random100)) {
  midpoint_lat <- random100$Lat.[i]
  midpoint_long <- random100$Long.[i]
  dg_radius <- as.numeric(refdata$radius_km[refdata$dataset=="LMEM"])
  filtered_coords <- filter_coordinates_by_radius(LMEM, midpoint_lat, midpoint_long, dg_radius)
  filtered_data <- filtered_coords[sample(nrow(filtered_coords), 19, replace = TRUE), ] ## downsample to UkrAd DateGroup N 
  filtered_data <- filtered_data[!duplicated(filtered_data$`Genetic ID`), ] ## remove duplicates from previous replacement being true
  output <- as.matrix(cbind(as.vector(rep(paste("LMEM_",i,sep=""), times=length(filtered_data$`Genetic ID`))), filtered_data$`Genetic ID`, filtered_data$`Master ID`, filtered_data$meandate, filtered_data$`Group ID`, filtered_data$Lat., filtered_data$Long.))
  print(filtered_data)
  LMEM_testdata=rbind(LMEM_testdata, output)
}
colnames(LMEM_testdata) <- c("iteration", "Genetic_ID", "Master_ID", "meandate", "group_ID", "Lat.", "Long.")
LMEM_testdata <- LMEM_testdata %>% group_by(iteration) %>% filter(n() >= 10) ##remove any with N < 10


#####use to analyse multidimensional PCA distances

plot(data_evec$X2[2:nrow(data_evec)],data_evec$X3[2:nrow(data_evec)], col="grey", pch=16) ##eyeballing check of PCA
points(data_evec$X2[16465:20578],data_evec$X3[16465:20578], col="lightgrey", pch=16) ##HO individuals
points(data_evec$X2[20579:20674],data_evec$X3[20579:20674], col="blue", pch=16) ##UkrAd individuals

scaled_evec <- data_evec
names(scaled_evec) [1] <- "Master_ID"
scaled_evec$Master_ID <- str_match(data_evec$X1, ".*:(.*)")[,2]
for(i in 2:ncol(data_evec)) {
  scaled_evec[, i] <- as.matrix(data_evec[, i])*as.numeric(data_evec[1,i])
  print(scaled_evec)
}
scaled_evec <- scaled_evec[2:nrow(scaled_evec),]

dategroups <- c("PSEI", "SEIA", "IAEM", "LMEM")
for(dg in dategroups) {
  mat <- scaled_evec[match(ukrad_data$Individual[ukrad_data$dategroup == dg], scaled_evec$Master_ID), nomatch=0]
  assign(paste(dg, "_distances", sep=""), c(dist(cbind(mat[,c(2:26)])))) ##no of colummns selected corresponds to no. of PCs used
}

ggplot() + ##plot samples by dategroup
  geom_boxplot(aes(x="1.PSEI", y=PSEI_distances)) +
  stat_summary(aes(x="1.PSEI", y=PSEI_distances), geom = "errorbar", fun.min = mean, fun.max = mean, fun = mean, width = .75, linetype = "dashed") +
  geom_boxplot(aes(x="2.SEIA", y=SEIA_distances)) +
  stat_summary(aes(x="2.SEIA", y=SEIA_distances), geom = "errorbar", fun.min = mean, fun.max = mean, fun = mean, width = .75, linetype = "dashed") +
  geom_boxplot(aes(x="3.IAEM", y=IAEM_distances)) +
  stat_summary(aes(x="3.IAEM", y=IAEM_distances), geom = "errorbar", fun.min = mean, fun.max = mean, fun = mean, width = .75, linetype = "dashed") +
  geom_boxplot(aes(x="4.LMEM", y=LMEM_distances)) +
  stat_summary(aes(x="4.LMEM", y=LMEM_distances), geom = "errorbar", fun.min = mean, fun.max = mean, fun = mean, width = .75, linetype = "dashed") +
  ylab("Euclidean Distance") +
  xlab("DateGroup") +
  theme_minimal()

##run PSEI comparison

PSEI_iterations <- as.vector(unique(PSEI_testdata$iteration))
PSEI_testdistres = data.frame()
for(it in PSEI_iterations) {
  mat <- data_evec[match(PSEI_testdata$Master_ID[PSEI_testdata$iteration == it], data_evec$Master_ID), nomatch=0]
  outdist <- c(dist(cbind(mat[, c(2:21)])))
  output <- as.vector(c(it, summary(outdist)))
  PSEI_testdistres=rbind(PSEI_testdistres, output)
}
colnames(PSEI_testdistres) <- c("iteration", "Min.", "1st_Qu.",  "Median", "Mean", "3rd_Qu.", "Max.")

PSEI_testdists = data.frame(matrix(ncol = length(PSEI_iterations), nrow = 105)) ###15 choose 2 = 105
colnames(PSEI_testdists) = PSEI_iterations
for(it in PSEI_iterations) {
  mat <- data_evec[match(PSEI_testdata$Master_ID[PSEI_testdata$iteration == it], data_evec$Master_ID), nomatch=0]
  outdist <- c(dist(cbind(mat[, c(2:21)])))
  #print(length(outdist))
  output <- c(outdist, rep(NA, 105 - length(outdist)))
  PSEI_testdists[it] <- output
}
par(mfrow = c(2, 2))
plot(density(PSEI_distances), col="red", ylim=c(0,80), xlim=c(0,0.6), xaxs="i", yaxs="i", main = "PSEI Density Comparison", xlab = "Euclidean distance")
for(i in 1:length(PSEI_iterations)) {
  lines(density(na.omit(PSEI_testdists[,i])), col="lightgrey")
}
lines(density(PSEI_distances), col="red")
     

##run SEIA comparison
SEIA_iterations <- as.vector(unique(SEIA_testdata$iteration))
SEIA_testdists = data.frame(matrix(ncol = length(SEIA_iterations), nrow = 351)) ###27 choose 2
colnames(SEIA_testdists) = SEIA_iterations
for(it in SEIA_iterations) {
  mat <- data_evec[match(SEIA_testdata$Master_ID[SEIA_testdata$iteration == it], data_evec$Master_ID), nomatch=0]
  outdist <- c(dist(cbind(mat[, c(2:21)])))
  #print(length(outdist))
  output <- c(outdist, rep(NA, 351 - length(outdist)))
  SEIA_testdists[it] <- output
}

plot(density(SEIA_distances), col="blue", ylim=c(0,80), xlim=c(0,0.6), xaxs="i", yaxs="i", main = "SEIA Density Comparison", xlab = "Euclidean distance")
for(i in 1:length(SEIA_iterations)) {
  lines(density(na.omit(SEIA_testdists[,i])), col="lightgrey")
}
lines(density(SEIA_distances), col="blue")


##run IAEM comparison
IAEM_iterations <- as.vector(unique(IAEM_testdata$iteration))
IAEM_testdists = data.frame(matrix(ncol = length(IAEM_iterations), nrow = 325)) ###26 choose 2
colnames(IAEM_testdists) = IAEM_iterations
for(it in IAEM_iterations) {
  mat <- data_evec[match(IAEM_testdata$Master_ID[IAEM_testdata$iteration == it], data_evec$Master_ID), nomatch=0]
  outdist <- c(dist(cbind(mat[, c(2:21)])))
  #print(length(outdist))
  output <- c(outdist, rep(NA, 325 - length(outdist)))
  IAEM_testdists[it] <- output
}

plot(density(IAEM_distances), col="purple", ylim=c(0,80), xlim=c(0,0.6), xaxs="i", yaxs="i", main = "IAEM Density Comparison", xlab = "Euclidean distance")
for(i in 1:length(IAEM_iterations)) {
  lines(density(na.omit(IAEM_testdists[,i])), col="lightgrey")
}
lines(density(IAEM_distances), col="purple")


##run LMEM comparison
LMEM_iterations <- as.vector(unique(LMEM_testdata$iteration))
LMEM_testdists = data.frame(matrix(ncol = length(LMEM_iterations), nrow = 171)) ###19 choose 2
colnames(LMEM_testdists) = LMEM_iterations
for(it in LMEM_iterations) {
  mat <- data_evec[match(LMEM_testdata$Master_ID[LMEM_testdata$iteration == it], data_evec$Master_ID), nomatch=0]
  outdist <- c(dist(cbind(mat[, c(2:21)])))
  #print(length(outdist))
  output <- c(outdist, rep(NA, 171 - length(outdist)))
  LMEM_testdists[it] <- output
}

plot(density(LMEM_distances), col="green", ylim=c(0,80), xlim=c(0,0.6), xaxs="i", yaxs="i", main = "LMEM Density Comparison", xlab = "Euclidean distance")
for(i in 1:length(LMEM_iterations)) {
  lines(density(na.omit(LMEM_testdists[,i])), col="lightgrey")
}
lines(density(LMEM_distances), col="green")







