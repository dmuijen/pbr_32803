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
             tabsetPanel(
              tabPanel("Genetic Map",
               wellPanel(
               ggiraphOutput('raw_plot', width = "100%", height = "600px")
               )
               ),
              tabPanel("Graphical genotype",
                       sidebarLayout(
                         sidebarPanel(
                           uiOutput("selectorgenoImage")
                         ),
                         mainPanel(
                           plotOutput('genoImage'),
                           helpText("Graphical genotype representation - red corresponds to marker genotypes of parent A, blue to parent B, and white to missing or heterozygous genotype calls")
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
             tabsetPanel(
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
                    )
                   ),
                  mainPanel(
                    ggiraphOutput("map_plot", width="100%", height="800px")
                  )
                )
               ),
               tabPanel("Recombination Frequency",
                        helpText(h4("Pairwise recombination fractions and LOD scores", align = "center")),
                         plotOutput('rf_plot', height = "800px")
               ),
               tabPanel("Break/merge Linkage Groups",
                        sidebarLayout(
                          sidebarPanel(
                            uiOutput("breakcombine"),
                            uiOutput("LGCombine"),
                            uiOutput("MarkCombine"),
                            actionButton("button", "Go!")
                            ),
                          mainPanel(
                            # ggiraphOutput("map_plot", width="100%", height="800px")
                            helpText("blabla")
                          )
                        )
               ),
               tabPanel("Map QC",
                        sidebarLayout(
                          sidebarPanel(
                            selectizeInput('qcType',
                                    label = 'Mapping function',
                                    choices = c("dxo","seg.dist")
                            )
                          ),
                          mainPanel(
                           plotOutput('qc_plot')
                          )
                        )
                )
             )
             ),
    tabPanel(div(h4("Export Results")),
      sidebarLayout(
        sidebarPanel(
             htmlOutput("save")
             ),
        mainPanel(
          helpText("Select a dataset to export and click download")
        )
        )
    ),
    tabPanel(div(h4("QTL mapping")),
      sidebarLayout(
        sidebarPanel(
          selectInput('e2','Interval mapping',
            choices = c("Interval mapping","Stepwise")
          ),
          htmlOutput("selectUI1")
        ),
        mainPanel(
          # iplotScanone_output('scanone', width = "100%", height = "580")
          plotOutput('scanone')
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