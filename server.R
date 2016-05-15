library(shiny)
shinyServer(function(input, output, session) {
  
  ##################################
  ###### Upload data
  ##################################
  geno <- reactive({
    inFile <- input$file1
    validate(
      need(input$file1 != "", "Upload a cross file to begin")
    )
    cross <- read.cross(
      format = "csv",
      file = inFile$datapath,
      genotypes = c("A","H","B"),
      alleles = c("A","B"),
      estimate.map = FALSE,
      BC.gen = 0,
      F.gen = 6
    )
    # cross <- cross %>% convert2riself()
    cross <- est.rf(cross)
    cross
  })
  
  output$inputSummary <- renderPrint({
    validate(
      need(input$file1 != "", "Upload a cross file to begin")
    )
    summary(geno())
  })
  
###############################
###### Data exploration
###############################
	
## Plot of raw data
output$raw_plot <- renderggiraph({
  validate(
    need(input$file1 != "", "Upload a cross file to begin")
  )
  if(input$mapactivator == 0){
  mymap <- pull.map(geno(), as.table = T)
  mymap <- data.frame(marker = row.names(mymap), mymap)
  mapplot <- ggplot(aes(x = factor(chr), y = pos, tooltip=marker, data_id=marker), data = mymap)  +
    geom_path(data = mymap, aes(x = factor(chr), y = pos, group = chr),size=4.5, lineend="round", col = colors()[284]) +
    theme_dark(base_size = 12)  + ylab("cM") +
    xlab("Linkage Group") + scale_y_reverse() +
    geom_point_interactive(size=4.5, col="orange", shape = 95) 
  ggiraph(code = {print(mapplot)}, zoom_max = 2, tooltip_offx = 20, tooltip_offy = -10, hover_css = "fill:black;stroke-width:1px;stroke:wheat;cursor:pointer;alpha:1;")
  } else {
    mymap <- pull.map(mstresult(), as.table = T)
    mymap <- data.frame(marker = row.names(mymap), mymap)
    mapplot <- ggplot(aes(x = factor(chr), y = pos, tooltip=marker, data_id=marker), data = mymap)  +
      geom_path(data = mymap, aes(x = factor(chr), y = pos, group = chr),size=4.5, lineend="round", col = colors()[284]) +
      theme_dark(base_size = 12)  + ylab("cM") +
      xlab("Linkage Group") + scale_y_reverse() +
      geom_point_interactive(size=4.5, col="orange", shape = 95)
    ggiraph(code = {print(mapplot)}, zoom_max = 2, tooltip_offx = 20, tooltip_offy = -10, hover_css = "fill:black;stroke-width:1px;stroke:wheat;cursor:pointer;alpha:1;")
  }
  })

output$selectorgenoImage <- renderUI({
  if(input$mapactivator == 0){
    dataset <- geno() 
  } else {
    dataset <- mstresult()
  }
  selectizeInput("genoImagesub", "Subset graphical genotype", multiple = TRUE, choices = names(dataset$geno) %>% as.list,
                 selected = names(dataset$geno) %>% as.list)
})

output$selectorHeatmap <- renderUI({
  if(input$mapactivator == 0){
    dataset <- geno() 
  } else {
    dataset <- mstresult()
  }
  selectizeInput("selectorHeatmap", "Subset heatmap", multiple = TRUE, choices = names(dataset$geno) %>% as.list,
                 selected = names(dataset$geno) %>% as.list)
})

Heatmap <- eventReactive(input$actionheatmap,{
  validate(
    need(input$file1 != "", "Upload a cross file to begin"),
    need(input$selectorHeatmap != "", "")
  )
  if(input$mapactivator == 0){
    heatMap(geno(), main = "", chr = input$selectorHeatmap)
  } else {
    heatMap(mstresult(), main = "", chr = input$selectorHeatmap)
  }
})

output$mst <- reactive({
  if (is.null(mstresult))
    return(NULL)
  print(mstresult())
})

output$Heatmap2 <-  renderPlot({
  Heatmap()
})

tablePlot <- eventReactive(input$startPlotTable,{
  validate(
    need(input$genoImagesub != "","Please select a set of chromosomes"))
  
  if (is.null(geno))
    return(NULL)
  
  if(input$mapactivator == 0){
    table <- geno()$geno
  } else {
    table <- mstresult()$geno
  }
  
  if (!is.null(input$genoImagesub)){
    table <- lapply(input$genoImagesub, function(chr){
      table <- melt(t(table[[chr]]$data))
      table <- cbind(chr,table)
    })
    table <- as.data.frame(do.call(rbind, table))
    table$value <- factor(table$value)
  }
  
  ggplot(table, aes(x=Var1,y=Var2, fill=value)) + 
    geom_tile() + facet_grid(.~chr,scales="free") + 
    theme_dark() + theme(axis.ticks.x=element_blank(),axis.text.x=element_blank()) + 
    scale_y_discrete(expand = c(0, 0)) + 
    labs(x="Markers",y="Individuals") + 
    scale_fill_manual(values=c("#004CC7","#008A05","#C70D00"),na.value="white",labels=c("A","H","B"), guide_legend(title="Genotype"))  
})

output$genoImage <- renderPlot({
  print(tablePlot())
})

##############
### Marker QC
##############

statmark <- reactive({
  statMark(mstresult(), stat.type = c("marker","interval"), map.function = 'kosambi')
})

statgen <- reactive({
  statGen(mstresult(), stat.type = c("xo","dxo"), id = "RILs", bychr = FALSE)
})

output$qcType1 <- renderUI({
  if(input$mapactivator == 0)
    return(NULL)
  selectInput("qcType1", "QC at marker or interval level", multiple = FALSE, choices = c("marker","interval"))
})

output$qc_plot1 <- renderPlot({
  if(input$mapactivator == 0)
    return(NULL)
  if(input$qcType1 == "marker"){
    p <- ggplot(aes_string(x = "pos", y = input$qcType2), data = statmark()$marker)
    p <- p + geom_line(lwd = 1.2) + theme_bw(base_size = 14) + facet_wrap(~chr) + xlab("Position (cM)")
    print(p)
  }
  if(input$qcType1 == "interval"){
    p <- ggplot(aes_string(x = "pos", y = input$qcType3), data = statmark()$interval)
    p <- p + geom_line(lwd = 1.2) + theme_bw(base_size = 14) + facet_wrap(~chr) + xlab("Position (cM)")
    print(p)
  }
})

output$qc_plot2 <- renderPlot({
  if(input$mapactivator == 0)
    return(NULL)
  if(input$qcType4 == "number of crossovers"){
    p <- ggplot(aes(x = index, y = xo), data = data.frame(index = attributes(statgen()$xo)$names, xo = statgen()$xo))
    p <- p + geom_point(size = 1.2) + theme_bw(base_size = 10) + xlab("RIL number") + ylab("Nr.Crossovers") +
      geom_text(aes(label=index),hjust=0, vjust=0) + scale_x_discrete(breaks=NULL)
    print(p)

  }
  if(input$qcType4 == "number of double crossovers"){
    p <- ggplot(aes(x = index, y = dxo), data = data.frame(index = attributes(statgen()$dxo)$names, dxo = statgen()$dxo))
    p <- p + geom_point(size = 1.2) + theme_bw(base_size = 10)  + xlab("RIL number") + ylab("Nr.DoubleCrossovers") +
      geom_text(aes(label=index),hjust=0.1, vjust=0.1) + scale_x_discrete(breaks=NULL)
    print(p)
  }
})


##############
### Facet plot
##############

output$chromFacetPlot <- renderggiraph({
  if(is.null(input$chromInput[1]))
    return(NULL)
  validate(
    need(input$file1 != "", "Upload a cross file to begin")
  )
  min <- as.numeric(input$chromInput[1])
  max <- as.numeric(input$chromInput[2])
  if(input$mapactivator == 0){
    plot <- LGChrom.facetplot(posmap,min,max, cross = geno())
  } else {
    plot <- LGChrom.facetplot(posmap,min,max, cross = mstresult())
  }
  ggiraph(code = {print(plot)}, zoom_max = 2, tooltip_offx = 20, tooltip_offy = -10, hover_css = "fill:black;stroke-width:1px;stroke:wheat;cursor:pointer;alpha:1;")
})

output$chromInput <- renderPrint({input$chromInput})

output$chromSlider <- renderUI({
  if(input$mapactivator == 0){
    sliderInput.custom(inputId="chromInput", label="Chromosome :", value=c(2,4), min=0, max=13, custom.ticks=names(geno()$geno))
  } else {
    sliderInput.custom(inputId="chromInput", label="Chromosome :", value=c(2,4), min=0, max=13, custom.ticks=names(mstresult()$geno))
  }
})

output$chromSelect <- renderUI({
  if(input$mapactivator == 0){
    selectizeInput(inputId="chromSelect", 
                   label = "Select Chromosome", 
                   choices = names(geno()$geno)
    )
  } else {
    selectizeInput(inputId="chromSelect",
                   label = "Select Chromosome",
                   choices = names(mstresult()$geno)
    )
  }
})


output$markerSelect <- renderUI({
  if(input$mapactivator == 0){
    selectizeInput(inputId="markerSelect", 
                   label = 'Select a marker to compare', 
                   choices = colnames(geno()$geno[[input$chromSelect]]$data), 
                   multiple = TRUE
    )
  } else {
    selectizeInput(inputId="markerSelect", 
                   label = 'Select a marker to compare', 
                   choices = colnames(mstresult()$geno[[input$chromSelect]]$data), 
                   multiple = TRUE
    )
  }
})

alleleTable <- reactive({
  validate(
    need(input$file1 != "", "Upload a cross file to begin"),
    need(input$markerSelect != "NULL" , "Please select a set of markers")
  )
  
  if(input$mapactivator == 0){
    table <- t(geno()$geno[[input$chromSelect]]$data[,input$markerSelect])
  } else {
    table <- t(mstresult()$geno[[input$chromSelect]]$data[,input$markerSelect])
  }
  colnames(table) <- seq(ncol(table))
  table[is.na(table)] <- "NA"
  table[which(table==1)] <- "A"
  table[which(table==2)] <- "H"
  table[which(table==3)] <- "B"
  return(table)
})

# Filter data based on selections
output$genotable <- DT::renderDataTable(DT::datatable(alleleTable(),
                                                      extensions = "FixedColumns",
                                                      options = list(
                                                        dom = 't',
                                                        scrollX = TRUE,
                                                        scrollY = "600px",
                                                        fixedColumns = TRUE,
                                                        scrollCollapse = TRUE,
                                                        paging = FALSE
                                                      )
) 
%>% 
  formatStyle(colnames(alleleTable()),
              fontWeight = 'bold',
              backgroundColor = styleEqual(c("A","H","B","NA"),c("#004CC7","#008A05","#C70D00","#737373")),#BLUE,GREEN,RED,GRAY
              color = styleEqual(c("A","H","B","NA"),c("black","black","black","white"))
  )
)


##################################
###### Genetic map construction
##################################

mstresult <- eventReactive(input$mapactivator,{
  if (is.null(geno()))
    return(NULL)
  mapobject <- mstmap.cross(geno(), bychr = FALSE, dist.fun = input$mapping, 
                            trace = FALSE, id = "RILs",
                            p.value = 10^-input$split)
  mymap <- pull.map(mapobject, as.table = T)
  mymap <- data.frame(marker = row.names(mymap), mymap)
  mymap$refgroup <- lapply(mymap$marker %>% as.character %>% strsplit(split = "-", fixed = T),"[",2) %>% as.numeric
  mymap$bp <-lapply(mymap$marker %>% as.character %>% strsplit(split = "-", fixed = T),"[",1) %>% as.numeric
  
  ### Invert linkage groups if necessary
  cor_genfys <- mymap %>% group_by(chr) %>% summarise(cor = cor(pos, bp))
  invert <- filter(cor_genfys, cor < 0)
  #Invert groups
  if(nrow(invert) != 0){
    for(i in 1:nrow(invert)){
      invert_index <- which(mymap$chr == invert[i,]$chr)
      max_pos <- mymap[invert_index ,]$pos %>% max
      mymap[invert_index,]$pos <- (mymap[invert_index,]$pos - max_pos)*-1
    }
  }
  the_map <- split(mymap[,3], mymap$chr %>% as.factor)
  lgindex <- mymap$chr %>% as.factor %>% unique
  class(the_map) <- "map"
  for(i in seq_along(the_map)){
    names(the_map[[i]]) <- filter(mymap, chr == lgindex[i])$marker
    class(the_map[[i]]) <- "A"
  }
  mapobject <- replace.map(mapobject, the_map)
  
  ## Rename linkage groups
  ref_dictionary <- mymap %>% group_by(refgroup, chr) %>% summarise(a = n()) %>% group_by(chr) %>% top_n(wt = a,n = 1)
  subgroup_index <- which(duplicated(ref_dictionary$refgroup))
  ref_dictionary$refgroup[subgroup_index - 1] <- ref_dictionary$refgroup[subgroup_index - 1] + 0.1
  ref_dictionary$refgroup[subgroup_index] <- ref_dictionary$refgroup[subgroup_index]  + 0.2
  for(i in seq_along(ref_dictionary$chr)){
    LG_name <- which(names(mapobject$geno) == ref_dictionary$chr[i])
    names(mapobject$geno)[LG_name] <- ref_dictionary$refgroup[i]
  }
  
  ### Rearrange linkage groups name order in cross object
  chr_names_old <- names(mapobject$geno) 
  chr_names_new <- names(mapobject$geno) %>% mixedsort()
  mapobject2 <- mapobject
  mapobject2$geno <- list()
  ### Best piece of coding ever..
  for(i in seq_along(mapobject$geno)){
    k <- which(chr_names_new == chr_names_old[i])
    mapobject2$geno[[k]] <- mapobject$geno[[i]]
  }
  names(mapobject2$geno) <- chr_names_new
  mapobject2 <- est.rf(mapobject2)
  mapobject2
})
## Plot of mst ordered data
output$map_plot <- renderggiraph({
  validate(
    need(input$file1 != "", "Upload a cross file to begin")
  )
  mymap <- pull.map(mstresult(), as.table = T)
  mymap <- data.frame(marker = row.names(mymap), mymap)
  mymap$chr <- mymap$chr %>% as.character %>% as.numeric()
  mymap <- mymap %>% arrange(chr)
  mapplot <- ggplot(aes(x = factor(chr), y = pos, tooltip=marker, data_id=marker), data = mymap)  +
    geom_path(data = mymap, aes(x = factor(chr), y = pos, group = chr),size=4.5, lineend="round", col = colors()[284]) +
    theme_dark(base_size = 12)  + ylab("cM") +
    xlab("Linkage Group") + scale_y_reverse() +
    geom_point_interactive(size=4.5, col="orange", shape = 95) 
  ggiraph(code = {print(mapplot)}, zoom_max = 2, tooltip_offx = 20, tooltip_offy = -10, hover_css = "fill:black;stroke-width:1px;stroke:wheat;cursor:pointer;alpha:1;")
})

########################
#####Export results ####
########################

# Ui to save files
output$save <- renderUI({
  list(
    selectInput("datatype", label = "Select data to export", 
                choices = c("Genotype file","Genetic map","Segregation distortion","Recombination fractions"), 
                selected = 1),
    downloadButton('downloadData', 'Save')
  )
  
})

output$DownloadData <- downloadHandler(
  filename = function() {
    paste(input$Download, format(Sys.time(), "%a %b %d %X"), '.csv', sep='', col.names = NA)
  },
  content = function(con) {
    write.csv(data, con)
  }
)

# Download file handler
observeEvent(input$datatype, {
  if (input$datatype == 'Genotype file'){
    output$downloadData <- {
      downloadHandler(
        filename = function() {paste0('Genotype_', format(Sys.time(), "%a %b %d %X"), '.csv') },
        content = function(file) {
          genodata <- mstresult() %>% pull.geno %>% as.data.frame
          row.names(genodata) <- mstresult()$pheno$RILs
          genodata[genodata == 1] <- "A"
          genodata[genodata == 2] <- "H"
          genodata[genodata == 3] <- "B"
          genodata <- cbind(Rils = rownames(genodata), genodata)
          row.names(genodata) <- NULL
          write.table(genodata, file, sep = ",", row.names = FALSE)
        }
      )}
  }
  if (input$datatype == 'Genetic map'){
    output$downloadData <- {
      downloadHandler(
        filename = function() {paste0('Genetic map_', format(Sys.time(), "%a %b %d %X"), '.csv') },
        content = function(file) {
          mymap <- pull.map(mstresult(), as.table = T)
          mymap <- data.frame(marker = row.names(mymap), mymap)
          mymap$bp <-lapply(mymap$marker %>% as.character %>% strsplit(split = "-", fixed = T),"[",1) %>% as.numeric
          mymap <- mymap %>% arrange(chr,pos)
          write.table(mymap, file, sep = ",", row.names = FALSE)
        }
      )}
  }
  if (input$datatype == 'Segregation distortion'){
    output$downloadData <- {
      downloadHandler(
        filename = function() {paste0('SegDist_', format(Sys.time(), "%a %b %d %X"), '.csv') },
        content = function(file) {
          segdist_data <- mydata %>% geno.table
          colnames(segdist_data)[3:4] <- c("A","B")
          segdist_data <- cbind(marker = rownames(segdist_data), segdist_data )
          map <- pull.map(mydata, as.table = T)
          map <- cbind(marker = row.names(map),map)
          result <- left_join(map[,c(1,3)], segdist_data, by = "marker") %>% arrange(chr,pos)
          write.table(result[,1:6], file, sep = ",", row.names = FALSE)
        }
      )}
  }
  if (input$datatype == 'Recombination fractions'){
    output$downloadData <- {
      downloadHandler(
        filename = function() {paste0('RecFraction_', format(Sys.time(), "%a %b %d %X"), '.csv') },
        content = function(file) {
          rf_data <- mstresult()$rf %>% as.data.frame()
          rf_data <- cbind(marker = row.names(rf_data), rf_data )
          write.table(rf_data, file, sep = ",", row.names = FALSE)
        }
      )}
  }
})
})


