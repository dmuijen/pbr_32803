library(shiny)
shinyUI(
  navbarPage("", inverse = TRUE, theme = shinytheme("cerulean"),
    tabPanel(div(h4("Upload data")),
          sidebarLayout(
           sidebarPanel(
                 fileInput('file1', 'Choose file to upload',
                    accept = c(
                      'text/csv',
                      'text/comma-separated-values',
                      'text/tab-separated-values',
                      'text/plain',
                      '.csv',
                      '.tsv'
                    )
                    )),
            mainPanel(
              verbatimTextOutput(outputId = "inputSummary")
                    )
          )
      ),
    tabPanel(div(h4("Data exploration")),
             tabsetPanel(id = "mapvis",
              tabPanel("Genetic Map",
               wellPanel(
               ggiraphOutput('raw_plot', width = "100%", height = "600px")
               )
               ),
              tabPanel("Graphical genotype",
                       sidebarLayout(
                         sidebarPanel(
                           uiOutput("selectorgenoImage"),
                           actionButton("startPlotTable","Plot graphical genotype")                           
                         ),
                         mainPanel(
                           plotOutput('genoImage'),
                           helpText("Graphical genotype representation - red corresponds to marker genotypes of parent A, blue to parent B, green to heterozygous markers, and white to missing genotypes")
                         )
                       )
              ), 
              tabPanel("Recombination counting",
                        sidebarLayout(
                          sidebarPanel(
                            uiOutput("chromSelect"),
                            uiOutput("markerSelect")
                          ),
                          mainPanel(
                            DT::dataTableOutput("genotable")
                          )
                        )
               ),
               tabPanel("Genomic vs. genetic distance",
                        sidebarLayout(
                          sidebarPanel(
                            uiOutput("chromSlider")
                          ),
                          mainPanel(
                            ggiraphOutput("chromFacetPlot", width="100%", height="800px")
                          )
                        )
               )
             )
    ),
    tabPanel(div(h4("Genetic map")),
             tabsetPanel(id = "mapcreate",
               tabPanel("Genetic map",
                        sidebarLayout(
                          sidebarPanel(
                            selectInput('mapping',
                                        label = 'Mapping function',
                                        choices = c("kosambi","haldane")
                            ),
                            br(),
                            sliderInput('split',
                                        label = 'Clustering -log10(P)',
                                        min = 2, max = 6,
                                        value = 4, step = 0.1
                            ),
                            actionButton("mapactivator", "Create map")
                          ),
                          mainPanel(
                            ggiraphOutput("map_plot", width="100%", height="800px")
                          )
                        )
               ),
               tabPanel("Recombination Frequency",
                        sidebarLayout(
                          sidebarPanel(
                            uiOutput("selectorHeatmap"),
                            actionButton("actionheatmap", "Plot heatmap")
                          ),
                          mainPanel(
                            helpText(h4("Pairwise recombination fractions and LOD scores for linkage", align = "center")),
                            plotOutput('Heatmap2', height = "800px")
                          )
                        )
               ), 
               tabPanel("Marker level QC",
                        sidebarLayout(
                          sidebarPanel(
                            conditionalPanel(condition = "input.qcType1 == 'marker'",
                                             h4("Marker statistics"),
                                             helpText("neglog10P = -log10 p-value from a test of segregation distortion"),
                                             helpText("missing = proportion of missing values"),
                                             helpText("AA = allele proportion of homozygous A allele"),
                                             helpText("BB = allele proportion of homozygous B allele"),
                                             selectInput("qcType1", "QC at marker or interval level", 
                                                         multiple = FALSE, choices = c("marker","interval")),
                                             selectInput("qcType2", "Select QC plot", 
                                                         multiple = FALSE, choices = c("neglog10P","missing","AA","BB"))
                            ),
                            conditionalPanel(condition = "input.qcType1 == 'interval'",
                                             h4("Interval statistics"),
                                             helpText("erf = estimated recombination fraction"),
                                             helpText("lod = LOD score for linkage"),
                                             helpText("dist = interval map distance"),
                                             helpText("mrf = map recombination fraction"),
                                             helpText("recomb = number of recombinations."),
                                             selectInput("qcType3", "Select QC plot", 
                                                         multiple = FALSE, choices = c("erf","lod","dist","mrf","recomb"))
                            )
                          ),
                          mainPanel(
                            plotOutput('qc_plot1')
                          )
                        )
               ),
               tabPanel("Genotype level QC",
                        sidebarLayout(
                          sidebarPanel(
                            h4("Genotype statistics"),
                            selectInput("qcType4", "QC genotype level", 
                                        multiple = FALSE, choices = c("number of crossovers","number of double crossovers"))
                          ),
                          mainPanel(
                            plotOutput('qc_plot2')
                          )
                        )
               )
             )
    ),
    tabPanel(div(h4("Export Results")),
             sidebarLayout(
               sidebarPanel(
                 helpText("If there issues opening in Excel, select the semicolon as delimiter"),
                                  selectInput("sep_type", "Select delimiter",
                             multiple = FALSE, choices = c(",",";")),
                 htmlOutput("save")
               ),
               mainPanel(
                 helpText("Select a dataset to export and click download")
               )
             )
    ),
    tabPanel(div(h4("About")),
             mainPanel(
               includeMarkdown("about.Rmd")
             )
    )
  )
)
