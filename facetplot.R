### This function takes a dataframe with columns po1,po2,chr1,chr2, and marker names. 
### po1 = position 1 (cM) po2 the bp position, chr1 name of linkage group, chr2, name of chromosomes 

### Example:
# ggdata <- data.frame(marker = c("m1","m2","m3"),
# po1 = c(0,10,25),
# po2 = c(0,100000,200000),
# chr1 = c(2,2,2),
# chr2 = c(2,2,2))

output$scatterplot <- renderPlot({
  data <- cdata()
  plotdata <- data$links %>% as.data.frame
  plotdata  <- transform(plotdata , chr2 = factor(chr2, levels = mixedsort(plotdata$chr2, decreasing = TRUE) %>% unique))
  plotdata$po1 <- plotdata$po1 %>% as.character %>% as.numeric
  plotdata$po2 <- plotdata$po2 %>% as.character %>% as.numeric
  p <- ggplot(aes(x = po2, y = po1), data = plotdata)
  p <- p + geom_point() + facet_grid(chr2~chr1) + theme_bw()
  print(p)
})
