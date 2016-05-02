library(ASMap)
library(dplyr)
library(stringr)
library(ggplot2)

mydata <- read.cross(format = "csv",
           file = "data/cross.csv" ,
           genotypes = c("a","h","b"),
           alleles = c("a","b"),
           estimate.map = FALSE,
           BC.gen = 0,
           F.gen = 6
)

test1 <- mstmap.cross(mydata, bychr = FALSE, dist.fun = "kosambi", 
                      trace = FALSE, id = "Index",
                      p.value = 4e-7)


### Rename linkage groups
themap <- pull.map(test1, as.table = T)
themap <- cbind(marker = row.names(themap), themap)
genome <- read.table("genpos_cm_map.csv", sep = ",", header = T)


themap_compare <- inner_join(genome, themap, by = c("LG"="marker") ) %>% arrange(chr,pos)
filter(themap_compare, chr == "LG8")

test2 <- breakCross(test1, split = list("LG8" = "62534652-12"))
# filter(pull.map(test2, as.table = T), chr %in% c("LG8.1","LG8.2"))
# 
# filter(themap_compare, chr %in% c("LG8.1","LG8.2"))

direction <- themap_compare %>% group_by(chr) %>% summarise(mycor = cor(GenPos,pos))
to_invert <- NULL
to_invert <- direction[direction[,"mycor"] < 0,]
themap_compare %>% head
for(i in 1:nrow(to_invert)){
    max_position <- (themap_compare[which(themap_compare$chr == to_invert[i,]$chr),])$pos %>% max
    themap_compare[which(themap_compare$chr == to_invert[i,]$chr),]$pos <- ((themap_compare[which(themap_compare$chr == to_invert[i,]$chr),])$pos - max_position)*-1
  }
p <- ggplot(aes(x = CM, y = pos), data = themap_compare)
p + geom_point() + facet_grid(Chrom~chr) + theme_bw()








### Segregation distortion
gt <- geno.table(test1, scanone.output=TRUE)
par(mfrow=c(2,1))
plot(gt, ylab=expression(paste(-log[10], " P-value")))
plot(gt, lod=3:5, ylab="Genotype frequency")
abline(h=c(0.25, 0.5), lty=2, col="gray")

## Cross over and Double XO
pg1 <- profileGen(test1, bychr = FALSE, stat.type = c("xo", "dxo",
                                                       "miss"), id = "Index", xo.lambda = 14, layout = c(1, 3), lty = 2, cex =
                    0.7)



