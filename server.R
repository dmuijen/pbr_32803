shinyServer(function(input, output, session) {
  
##################################
###### Upload data
##################################
  	geno <- reactive({
		inFile <- input$file1
		if (is.null(inFile))
		return(NULL)
		cross <- read.cross(
			format = "csv",
			file = inFile$datapath,
			genotypes = c("A","H","B"),
			alleles = c("A","B"),
			estimate.map = FALSE,
			BC.gen = 0,
			F.gen = 7
		)
		cross <- cross %>% convert2riself()
	  cross
	})

  	output$inputSummary <- renderPrint({
  	  if (is.null(geno())) {
  	    cat('No input data.')
  	  } else {
  	    summary(geno())
  	  }
  	})
  	
##################################
###### Data exploration
##################################
	
## Plot of raw data
output$raw_plot <- renderggiraph({
  validate(
    need(input$file1 != "", "Upload a cross file to begin")
  )
  mymap <- pull.map(geno(), as.table = T)
  mymap <- data.frame(marker = row.names(mymap), mymap)
  mapplot <- ggplot(aes(x = factor(chr), y = pos, tooltip=marker, data_id=marker), data = mymap)  +
    geom_path(data = mymap, aes(x = factor(chr), y = pos, group = chr),size=4.5, lineend="round", col = colors()[284]) +
    theme_dark(base_size = 12)  + ylab("cM") +
    xlab("Linkage Group") + scale_y_reverse() +
    geom_point_interactive(size=4.5, col="orange", shape = 95) 
  ggiraph(code = {print(mapplot)}, zoom_max = 2, tooltip_offx = 20, tooltip_offy = -10, hover_css = "fill:black;stroke-width:1px;stroke:wheat;cursor:pointer;alpha:1;")
})

output$selectorgenoImage <- renderUI({
  selectizeInput("genoImagesub", "Subset graphical genotype", multiple = TRUE, names(geno()$geno) %>% as.list)
})


output$rf_plot <- renderPlot({
	if (is.null(mstresult()))
		return(NULL)
	## A very crappy 'progressbar' for 5 secs.. 
	progress <- shiny::Progress$new(session, min=1, max=5)
	on.exit(progress$close())
	progress$set(message = 'Calculation in progress')
	for (i in 1:5) {
		progress$set(value = i)
		Sys.sleep(0.5)
	}
	rf.cross <- est.rf(mstresult())
	heatMap(rf.cross)
})

output$genoImage <- renderPlot({
  if (is.null(geno))
    return(NULL)
  if (!is.null(input$genoImagesub)){
    geno.image(geno(), chr = input$genoImagesub)
  } else {
    geno.image(geno())
    }
})

output$qc_plot <- renderPlot({
	suppressWarnings(profileMark(mstresult(), stat.type = input$qcType, id =
							"RILs", type = "a"))
})


output$chromFacetPlot <- renderggiraph({
  if (is.null(input$chromInput[1]))
    return(NULL)
  validate(
    need(input$file1 != "", "Upload a cross file to begin")
  )
  min <- as.numeric(input$chromInput[1])
  max <- as.numeric(input$chromInput[2])
  plot <- LGChrom.facetplot(posmap,min,max, cross = geno())
  ggiraph(code = {print(plot)}, zoom_max = 2, tooltip_offx = 20, tooltip_offy = -10, hover_css = "fill:black;stroke-width:1px;stroke:wheat;cursor:pointer;alpha:1;")
})

output$chromInput <- renderPrint({input$chromInput})

output$chromSlider <- renderUI({
  sliderInput.custom(inputId="chromInput", label="Chromosome :", value=c(2,4), min=0, max=13, custom.ticks=c("1a","1b","2","3","4","5","6","7","8","9","10","11","12"))
})

output$chromSelect <- renderUI({
  selectizeInput(inputId="chromSelect", 
                 label = "Select Chromosome", 
                 choices = names(geno()$geno)
  )
})

output$markerSelect <- renderUI({
  selectizeInput(inputId="markerSelect", 
                 label = 'Select a marker to compare', 
                 choices = colnames(geno()$geno[[input$chromSelect]]$data), 
                 multiple = TRUE
  )
})

alleleTable <- reactive({
  validate(
    need(input$file1 != "", "Upload a cross file to begin"),
    need(input$markerSelect != "NULL" , "Please select a set of markers")
  )
  table <- t(geno()$geno[[input$chromSelect]]$data[,input$markerSelect])
  colnames(table) <- seq(ncol(table))
  table[is.na(table)] <- "NA"
  table[which(table==1)] <- "A"
  table[which(table==2)] <- "B"
  # table[which(table==3)] <- "B"
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

mstresult <- reactive({
  if (is.null(geno()))
    return(NULL)
  mapobject <- mstmap.cross(geno(), bychr = FALSE, dist.fun = input$mapping, 
                            trace = FALSE, id = "RILs",
                            p.value = 10^-input$split)
  # names(mapobject$geno) <- paste0("LG",seq_along(mapobject$geno))
  mymap <- pull.map(mapobject, as.table = T)
  mymap <- data.frame(marker = row.names(mymap), mymap)
  mymap$refgroup <- lapply(mymap$marker %>% as.character %>% strsplit(split = "-", fixed = T),"[",2) %>% as.numeric
  ref_dictionary <- mymap %>% group_by(refgroup, chr) %>% summarise(a = n()) %>% group_by(chr) %>% top_n(wt = a,n = 1)
  subgroup_index <- which(duplicated(ref_dictionary$refgroup))
  ref_dictionary$refgroup[subgroup_index - 1] <- ref_dictionary$refgroup[subgroup_index - 1] + 0.1
  ref_dictionary$refgroup[subgroup_index] <- ref_dictionary$refgroup[subgroup_index]  + 0.2
  for(i in seq_along(ref_dictionary$chr)){
    LG_name <- which(names(mapobject$geno) == ref_dictionary$chr[i])
    names(mapobject$geno)[LG_name] <- ref_dictionary$refgroup[i]
  }
  mapobject
})

## Selectize selector to combine linkage groups
output$MapCombine <- renderUI({
  if (is.null(mstresult()))
    return(NULL)
  selectizeInput("MapCombine", "Manual combine", multiple = TRUE, names(mstresult()$geno) %>% as.list)
})

## Plot of mst ordered data
output$map_plot <- renderggiraph({
  validate(
    need(input$file1 != "", "Upload a cross file to begin")
  )
  # plotMap(mstresult())
  mymap <- pull.map(mstresult(), as.table = T)
  mymap <- data.frame(marker = row.names(mymap), mymap)
  mymap$chr <- mymap$chr %>% as.character %>% as.numeric()
  mymaps <- mymap %>% arrange(chr)
  mapplot <- ggplot(aes(x = factor(chr), y = pos, tooltip=marker, data_id=marker), data = mymap)  +
    geom_path(data = mymap, aes(x = factor(chr), y = pos, group = chr),size=4.5, lineend="round", col = colors()[284]) +
    theme_dark(base_size = 12)  + ylab("cM") +
    xlab("Linkage Group") + scale_y_reverse() +
    geom_point_interactive(size=4.5, col="orange", shape = 95) 
  ggiraph(code = {print(mapplot)}, zoom_max = 2, tooltip_offx = 20, tooltip_offy = -10, hover_css = "fill:black;stroke-width:1px;stroke:wheat;cursor:pointer;alpha:1;")
})

  ####################
	#####QTL Mapping####
	####################
	output$selectUI1 <- renderUI({ 
	  if (is.null(geno()))
	    return(NULL)
	  selectInput("trait", "Select trait", names(geno()$pheno))
	})
	
	s1 <- reactive({
	  if (is.null(geno()))
	    return(NULL)
	  scanone(geno(), pheno.col = which(names(geno()$pheno)==input$trait))
	})
	
	output$scanone <- renderPlot({
	  validate(
	    need(input$file1 != "", "Upload a cross file to begin")
	  )
	  plot(s1())
	})
	
	########################
	#####Export results ####
	########################

	# Ui to save files
	output$save <- renderUI({
	  validate(
	    need(mstresult() != "", "Run genetic map construction first")
	  )
	  if (!is.null(mstresult())) {
	    list(
	      selectInput("datatype", label = "Select data to export", 
	                  choices = c("Genotype file","Genetic map","Segregation distortion"), 
	                  selected = 1),
	      downloadButton('downloadData', 'Save')
	    )
	  }
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
	          write.table(genodata, file, sep = ",", row.names = TRUE, col.names=NA)
	          }
	      )}  
	  }
	  if (input$datatype == 'Genetic map'){
	    output$downloadData <- {
	      downloadHandler(
	        filename = function() {paste0('Genetic map_', format(Sys.time(), "%a %b %d %X"), '.csv') },
	        content = function(file) {
	          write.table(data.frame(marker = row.names(mstresult() %>% pull.map(as.table = T)), mstresult() %>% pull.map(as.table = T)) %>%
	                        as.data.frame(), file, sep = ",", row.names = FALSE)
	        }
	      )}
	  }
	  if (input$datatype == 'Segregation distortion'){
	    output$downloadData <- {
	      downloadHandler(
	        filename = function() {paste0('SegDist_', format(Sys.time(), "%a %b %d %X"), '.csv') },
	        content = function(file) {
	          write.table(mstresult() %>% geno.table, file, sep = ",", col.names = NA)
	        }
	      )}
	  }
	})
	
})


