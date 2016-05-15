library(qtl)
library(shiny)
library(markdown)
library(ggplot2)
library(qtl)      
library(ASMap)
library(shinythemes)
library(ggthemes)
library(ggiraph)
library(gtools)
library(dplyr)
library(DT)
library(reshape2)
mstresult <- NULL
# mydata <- read.cross(format = "csv",
#            file = "data/cross.csv" ,
#            genotypes = c("A","H","B"),
#            alleles = c("A","B"),
#            estimate.map = FALSE,
#            BC.gen = 0,
#            F.gen = 7
# )
# mydata <- mydata %>% convert2riself()

sliderInput.custom <- function(inputId="placeholder", label="placeholder", ticks=TRUE, value=c(0,0), min=0, max=0, custom.ticks=c("placeholder")){
  args <- list(inputId=inputId, label=label, ticks=ticks, value=value, min=min, max=max)
  html <- do.call('sliderInput', args)
  ##<MAD HACKS>
  html$children[[2]]$attribs[['data-values']] <- paste(custom.ticks,collapse=',')
  ##</MAD HACKS>
  return(html)
}

LGChrom.facetplot <- function(posmap, min, max, cross){
  #expand ranges
  range <- seq(min+1,max+1,1)
  #map to chromosome
  # range <- c("1a","1b","2","3","4","5","6","7","8","9","10","11","12")[range]
  range <- names(cross$geno)[range]
  
  plot <- ggplot(filter(posmap,Chrom %in% range),aes(x=(GenPos/10^6) %>% round(digits = 0),y=CM, tooltip=LG, data_id=LG))
  plot <- plot + theme_dark() + facet_grid(LG.map ~ Chrom)
  plot <- plot + theme(axis.text.x=element_text(angle = 45, vjust=0.5)) + labs(x="Genome Position (Mbp)",y="Genetic position (cM)")
  plot <- plot + geom_point_interactive(size=1, col="orange", alpha=0.4)
  return(plot)
}

#Genome position mapped to centimorgan file
posmap <- read.delim("genpos_cm_map.csv", sep=",", stringsAsFactors=F)
posmap <- transform(posmap , LG.map = factor(LG.map, levels = mixedsort(posmap$LG.map, decreasing = T) %>% unique))
posmap <- transform(posmap , Chrom = factor(Chrom, levels = mixedsort(posmap$Chrom) %>% unique))

