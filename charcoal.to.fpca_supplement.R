#  R Script for fPCA analysis in Walsh et al.
#  The script requires several files for plotting that are not neccessary to
#  actually run the analysis.  These are noted below, at the point they are
#  first loaded.  Otherwise csv files used should be placed in the same
#  directory as this R script.
#  Coding by: S. Goring, contact - goring@wisc.edu or simon.j.goring@gmail.com

library(fda)
library(rgdal)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(mgcv)

# Either replace this and set the working directory manually, or add the home
#  directory within the quotes.

SET CURRENT WORKING DIRECTORY

################################################################################
#  These datasets require downloading from an external source.
#  The assumption in this code is that you will download and install into a folder
#  called 'Maps' in your parent directory.

#  This file can be acceded through the interface here:
#  http://wwwn.cdc.gov/epiinfo/script/shapefiles.aspx
usa <- readOGR('Maps/us.shp', 'us')

#  This file can be accesed through the GeoBase website here:
#  http://www.geobase.ca/geobase/en/search.do?produit=cgb&language=en
canada <- readOGR('Maps/PROVINCE.SHP', 'PROVINCE')

#  This file is available using this link: 
#   http://www.cec.org/atlas/files/Terrestrial_Ecoregions_L3/TerrestrialEcoregions_L3_Shapefile.zip
#  Change file name and location accordingly.
ecozones <- readOGR('Maps/TerrestrialEcoregions_L3_Shapefile/NA_Terrestrial_Ecoregions_Level_III/data/terrestrial_Ecoregions_updated//terrestrial_Ecoregions_updated.shp', 'terrestrial_Ecoregions_updated')

################################################################################

#  All csv files are included in the suplemental ZIP file.
cores <- list.files('data', full.names=TRUE)

#  The files aren't quite formed correctly, so the column names are off by
#  one.  Here all the data goes into a 
core.data <- lapply(cores, read.csv, row.names=NULL, stringsAsFactors=FALSE)

#  At issue here is how to smooth the data into reasonable bins:
interval <- sapply(core.data, function(x)mean(diff(x$est_age[x$est_age < 3000])))

apr.core <- function(x, interval, time.range){
  
  xout <-  seq(time.range[1], time.range[2], by = interval)
  
  z <- x
  
  x <- x[x$est_age < max(xout) + 400 & x$est_age > min(xout) - 200, ]
  
  x$influx.scale <- scale(x$influx)
  
  if(nrow(x) > 30 & any(z$est_age > max(xout))){
  
    model <- gam(influx.scale ~ s(est_age, k = nrow(x)), data = x)
    
    output <- data.frame(core = rep(as.numeric(substr(colnames(x)[1], 5, 8)), length(xout)),
                         age = xout,
                         char = predict(model, 
                                        newdata = data.frame(est_age = xout), 
                                        type='response'),
                         class = 'smooth', stringsAsFactors = FALSE)
    origin <- data.frame(core = rep(as.numeric(substr(colnames(x)[1], 5, 8)), nrow(x)),
                         age = x$est_age,
                         char = x$influx.scale,
                         class = 'raw', stringsAsFactors = FALSE)
    
    output <- rbind(output, origin)
  }
  else{
    origin <- data.frame(core = rep(as.numeric(substr(colnames(x)[1], 5, 8)), nrow(x)),
                         age = x$est_age,
                         char = x$influx.scale,
                         class = 'raw', stringsAsFactors = FALSE)
    
    output <- data.frame(core = rep(as.numeric(substr(colnames(x)[1], 5, 8)), length(xout)),
                         age = NA,
                         char = NA,
                         class = 'smooth', stringsAsFactors = FALSE)
  }
  
  output
}

time_range <- c(200, 3000)

core.approx <- lapply(core.data, apr.core, interval = 100, time.range = time_range)
#  Each column is one record.
core.mat <- sapply(core.approx, function(x) x$char[x$class == 'smooth'])

base.sites <- read.csv('PNW master sitelist_Jan 30.csv', stringsAsFactors=FALSE)
base.sites <- base.sites[match(as.numeric(substr(cores, 6, 9)), base.sites[,1]), ]

pts <- SpatialPoints(cbind(base.sites$Long, base.sites$Lat), 
                     proj4string=CRS('+proj=longlat +ellps=WGS84'))
pt.ecozones <- over(spTransform(pts, CRSobj=CRS(proj4string(ecozones))), ecozones)

#  These points are slightly mis-specified and show up in the water so they don't
#  get assigned values when the overlay command is used.  Here we assign them
#  their proper ecozone.
pt.ecozones$NAME_L1[is.na(pt.ecozones$NAME_L1)] <- 'Marine West Coast Forests'
pt.ecozones$NAME_L2[is.na(pt.ecozones$NAME_L2)] <- 'Marine West Coast Forests'

pt.ecozones$NAME_L1 <- factor(pt.ecozones$NAME_L1)
pt.ecozones$NAME_L2 <- factor(pt.ecozones$NAME_L2)

base.sites$ecozone <- pt.ecozones$NAME_L2

#  Find out the longest interval between samples in the data.  
#  We want less than 100 years.
good.sites <- colSums(is.na(core.mat)) == 0

sites <- base.sites[good.sites, ]

#  Turn the raw data into a functional data object and run fPCA.
char.mat <- core.mat[ , good.sites]

#  We want to be able to compare each of the cores just for our own edification,
#  and to make sure that this smoothing is okay:
big.comp.table <- ldply(core.approx, function(x)x)
big.comp.table$site <- base.sites$Site.Name[match(big.comp.table$core, base.sites[,1])]

big.comp.table <- big.comp.table[!is.na(big.comp.table$age),]

#  This is part of the fPCA process, first we need to define the basis functions
#  for the b-splines.
char.basis <- create.bspline.basis(rangeval=time_range, nbasis = 30)

char.model <- Data2fd(argvals = seq(min(time_range), max(time_range), by=100),
                        y = (char.mat),
                        basisobj= char.basis, 
                        lambda = 1e6)

char.fpca <- pca.fd(char.model, nharm=6, centerfns = FALSE); plot(char.fpca)


#  fPCA is done.  The rest is figures
################################################################################

col.pal <- colorRampPalette(c('#FF0000', '#FFFFFF', '#0000FF'), space='Lab')

fax.cols <- findInterval(char.fpca$scores[,1], 
                         seq(from=min(char.fpca$scores[,1]),
                             to = max(char.fpca$scores[,1]), 
                                      length.out=100))

sax.cols <- findInterval(char.fpca$scores[,2], 
                         seq(from=min(char.fpca$scores[,2]),
                             to = max(char.fpca$scores[,2]), 
                             length.out=100))

###############################
#  Some analysis:
#  char.frame contains all the factor levels, geographic data and 
fx.na <- function(x) x[!is.na(x)]

site.test <- sites

char.frame <- data.frame(scores = char.fpca$scores[,1:2],
                         lat = site.test$Lat,
                         long = site.test$Long,
                         elev = site.test$Elev,
                         inland = site.test$geog_15km,
                         high.low = site.test$elevation_500m,
                         forest = site.test$vegtype,
                         type = site.test$ecozone,
                         site = site.test$Site.Name)

summary(lm(scores.1 ~ type, data=char.frame))
summary(lm(scores.1 ~ high.low, data=char.frame))

summary(lm(scores.1 ~ inland*high.low, data=char.frame)) # neither are significant
anova(lm(scores.2 ~ inland*high.low, data=char.frame))

anova(lm(scores.1 ~ high.low, data=char.frame))
anova(lm(scores.2 ~ high.low, data=char.frame)) 

anova(lm(scores.1 ~ forest, data=char.frame))
anova(lm(scores.2 ~ forest, data=char.frame))

################################################################################
#  This plots the actual axis scores for the sites over time.  It uses code
#  modified from the original fPCA code.
source('get.fpca.R')

#  This returns an object with two elements, of which 'sets' contains the high, low
#  and mean curves for each of the first four fPCA harmonics.  'yrange' contains the
#  y ranges for each of the harmonics to assist in plotting.
plotting.fpca <- get.pos.min(char.fpca)
axis.one <- with(plotting.fpca, 
                 data.frame(age = rep(sets[[1]]$x, 3),
                            score = unlist(sets[[1]][,c(3,2,4)]),
                            line  = rep(c('high', 'mean', 'low'), each=nrow(sets[[1]]))))

axis.two <- with(plotting.fpca, 
                 data.frame(age = rep(sets[[2]]$x, 3),
                            score = unlist(sets[[2]][,c(3,2,4)]),
                            line  = rep(c('high', 'mean', 'low'), each=nrow(sets[[2]]))))

output <- rbind(axis.one, axis.two)
output$axis <- rep(c('fPCA Axis One', 'fPCA Axis Two'), each=nrow(axis.one))

big.comp.table$fpca1 <- factor(c('fPCA One < 0', 'fPCA One > 0')[as.numeric(char.fpca$scores[match(big.comp.table$site, sites$Site.Name),1]>1) + 1], levels = c('fPCA One > 0', 'fPCA One < 0'))
big.comp.table$fpca2 <- factor(c('fPCA Two < 0', 'fPCA Two > 0')[as.numeric(char.fpca$scores[match(big.comp.table$site, sites$Site.Name),2]>1) + 1], levels = c('fPCA Two > 0', 'fPCA Two < 0'))

big.comp.table <- big.comp.table[findInterval(big.comp.table$age, time_range) == 1, ]

axis.plot <- ggplot(data=output) + 
            #geom_line(data = big.comp.table, aes(x = age, y = char), alpha = 0.2) +
            geom_line(aes(x = age, y = score, group=line, color=line, size=line)) +
            scale_size_manual(values=c(0.5, 2, 2)) +
            scale_color_manual(values=c('black', 'blue', 'red')) +
            ylab('fPCA Axis Score') +
            xlab('Calendar Years Before Present') +
            theme_bw() +
            theme(axis.text = element_text(size = 14, family='serif'),
                  axis.title = element_text(size = 24, family='serif', face='bold'),
                  strip.text = element_text(size = 24, family='serif', face='bold'),
                  rect = element_rect(color='white'),
                  legend.text = element_text(size=14, family='serif'),
                  legend.position = 'none') +
            facet_wrap(facets= ~axis) +
            scale_x_reverse(expand=c(0,0), limits=rev(time_range), breaks = seq(0, 2500, by = 500)) +
      scale_y_continuous(expand=c(0,0), limits = c(-.8, 1))

bb <- ggplot(big.comp.table) +
  geom_line(data = big.comp.table, aes(x = age, y = char)) +
  facet_wrap(~fpca1, ncol = 1) +
  scale_x_reverse(expand=c(0,0), limits=rev(time_range), breaks = c(500, 1500, 2500)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() + xlab('') + ylab('Influx') +
  theme(axis.text = element_text(size = 14, family='serif'),
        axis.title = element_text(size = 18, family='serif', face='bold'),
        strip.text = element_text(size = 18, family='serif', face='bold'),
        rect = element_rect(color='white'),
        legend.text = element_text(size=14, family='serif'),
        legend.position = 'none')

cc <- ggplot(big.comp.table) +
  geom_line(data = big.comp.table, aes(x = age, y = char)) +
  facet_wrap(~fpca2, ncol = 1) +
  scale_x_reverse(expand=c(0,0), limits=rev(time_range), breaks = c(500, 1500, 2500)) +
  theme_bw() + xlab('') + ylab('Influx') +
  theme(axis.text = element_text(size = 14, family='serif'),
        axis.title = element_text(size = 18, family='serif', face='bold'),
        strip.text = element_text(size = 18, family='serif', face='bold'),
        rect = element_rect(color='white'),
        legend.text = element_text(size=14, family='serif'),
        legend.position = 'none')

#
################################################################################
#
#  Plotting countries with colors:
#  
crop.df <- function(x, field, ext){
  #  crops a spatial polygon and turns it into a data.frame for ggplot2:
  #  This code is modified from code posted to R-Sig-Geo, beginning here:
  #  https://stat.ethz.ch/pipermail/r-sig-geo/2012-June/015337.html
  
  require(PBSmapping); require(maptools)
  if(is.na(proj4string(x)))   proj4string(x) <- CRS('+proj=longlat +ellps=WGS84')
  x.ps <- SpatialPolygons2PolySet(x)
  x.c <- clipPolys(x.ps,xlim=ext[1:2],ylim=ext[3:4])
  x.sp <- PolySet2SpatialPolygons(x.c, close_polys=TRUE)
  x.df <- fortify(x.sp, field)

  x.df
}

us.df <- crop.df(usa, 'STATE_NAME', c(-125, -115, 43, 52))
ca.df <- crop.df(canada, 'NAME', c(-125, -115, 43, 52))

char.folded <- with(char.frame, data.frame(score = c(scores.1, scores.2),
                                           lat = rep(lat, 2),
                                           long = rep(long, 2),
                                           axis = rep(c('Axis One Score', 'Axis Two Score'), each=nrow(char.frame))))

countries <- ggplot(char.folded, aes(x = long, y = lat)) +
        geom_polygon(data= us.df, aes(x = long, y = lat, group=factor(group)), fill='gray50') +
        geom_path(data= us.df, aes(x = long, y = lat, group=factor(group)), color='white') +
        geom_polygon(data= ca.df, aes(x = long, y = lat, group=factor(group)), fill='gray50') +
        geom_path(data= ca.df, aes(x = long, y = lat, group=factor(group)), color='white') +
        geom_point(color='black', size=6) +
        geom_point(aes(color=score), size=4) +
        facet_wrap(~axis) +
        scale_color_gradient2(low='blue', mid = 'white', high='red', limits=c(-40,40)) +
        scale_x_continuous(breaks=c(-116, -120, -124), limits = c(-125, -115), expand=c(0,0)) +
        scale_y_continuous(limits = c(43, 52), expand=c(0,0)) +
        xlab('Longitude') + ylab('Latitude') +
        theme_bw() +
        theme(axis.text = element_text(size = 14, family='serif'),
              axis.title = element_text(size = 18, family='serif', face='bold'),
              strip.text = element_text(size = 18, family='serif', face='bold'),
              rect = element_rect(color='white'),
              legend.text = element_text(size=14, family='serif'))


grid.arrange(aa, arrangeGrob(bb, cc, nrow = 1), countries, ncol=1)
ggsave(filename='fPCAplot.png', width = 6, height = 8, dpi = 300)
#  This is the plot of all sites, it's ordered by quadrant
big.comp.table$site <- factor(big.comp.table$site, 
                              levels= c('Fish Lake', 'Lost Lake', 'Mt Constitution C11',
                                 'Frozen Lake', 'Taylor Lake OR',
                                 'Mt Constitution C32', 'Shadow Lake', 'Sunrise Lake',
                                 'Beaver Lake', 'Cooley Lake', 
                                 'Battle Ground Lake', 'Lake Oswego', 'Little Sunrise Lake',
                                 'Little Lake','Mount Barr Cirque',
                                 'Mt Constitution C38', 'Rockslide Lake', 'Tipsoo Pond',
                                 'Panther', 'Yahoo Lake'))

ggplot(big.comp.table) + 
  geom_line(aes(x = age, y = char, color = class)) + 
  facet_wrap(~site) +
  xlab('Calibrated Years Before Present') +
  ylab('Scaled Charcoal Influx') +
  scale_x_continuous(limits = c(100, 3000), expand=c(0,0), breaks = c(500, 2500)) +
  scale_y_continuous(breaks = c(0,5,10)) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.title = element_text(size = 12))

ggsave(filename='multipSitePanel.png', width = 6, height = 6)
