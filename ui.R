
library(shiny)
ui <- fluidPage(
  tags$head(
    tags$title("IGEDA")
  ),
  tabsetPanel(
    tabPanel("Home", 
             tags$img(src = "https://wallpapercave.com/uwp/uwp4290623.png",
                      width = 600,
                      height = 600,
                      style = "display: block; margin: 50px auto;")
    ),
    tabPanel("RNA-Seq Analysis",
             titlePanel("RNA-Seq Analysis "),
             sidebarLayout(
               sidebarPanel(
                 fileInput("countFile", "Upload the Count Data File", accept = c(".txt")),
                 fileInput("sampleInfoFile", "Upload the Sample Information File", accept = c(".txt")),
                 selectInput("organismInput", "Select the Organism",
                             choices = c("Human" = "Homo sapiens",
                                         "Mouse" = "Mus musculus",
                                         "Rat" = "Rattus norvegicus",
                                         "Fly" = "Drosophila melanogaster",
                                         "Zebrafish" = "Danio rerio",
                                         "Arabidopsis" = "Arabidopsis thaliana",
                                         "Yeast" = "Saccharomyces cerevisiae",
                                         "Worm" = "Caenorhabditis elegans",
                                         "Bovine" = "Bos taurus",
                                         "Chicken" = "Gallus gallus",
                                         "Pig" = "Sus scrofa",
                                         "Canine" = "Canis lupus familiaris",
                                         "Rhesus" = "Macaca mulatta",
                                         "Xenopus" = "Xenopus tropicalis",
                                         "E. coli strain K12" = "Escherichia coli K-12",
                                         "Chimp" = "Pan troglodytes",
                                         "Anopheles" = "Anopheles gambiae",
                                         "E. coli strain Sakai" = "Escherichia coli Sakai",
                                         "Myxococcus xanthus DK 1622" = "Myxococcus xanthus DK 1622"
                             ),
                             selected = "Homo sapiens"),
                 actionButton("runAnalysisBtn", "Perform Analysis"),
                 uiOutput("contrastDropdown1"),
                 uiOutput("contrastDropdown2"),
                 actionButton("runContrastBtn", "Run Contrast")
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Design Matrix", tableOutput("orig_data"), downloadButton("downloaddesign", "Download",style="position: absolute; right:0; top:0; ")),
                   tabPanel("Top-Table", dataTableOutput("top_data"), downloadButton("downloadtop1", "Download",style="position: absolute; right:0; top:0; ")),
                   navbarMenu("Visualization",
                              tabPanel("MDS Plot", plotOutput("mds_map")),
                              tabPanel("HeatMap", plotOutput("heat_map"),downloadButton("downloadHeatmap", label = "Download ",style="position: absolute; right:0; top:0; ")),
                              tabPanel("Voom Plot", plotOutput("voom"),downloadButton("downloadVoom", label = "Download ",style="position: absolute; right:0; top:0; ")),
                              tabPanel("Box Plot", plotOutput("box"),downloadButton("downloadBoxplot", label = "Download ",style="position: absolute; right:0; top:0; ")),
                              tabPanel("Volcano Plot", plotOutput("vol"))
                   )
                 )
               )
             )
    ),
    navbarMenu("Microarray Analysis",   
               tabPanel("Affymetrix Microarray Analysis",
                        titlePanel("Affymetrix Microarray Analysis"),
                        sidebarLayout(
                          sidebarPanel(
                            # Input: Choose directory containing CEL files
                            fileInput("celDir", "Upload CEL Files", multiple = TRUE),
                            textInput("groups", "Enter Group Labels (comma-separated)"),
                            #Tumor,Tumor,Tumor,Tumor,Tumor,Normal,Normal,Normal,Normal,Normal
                            # Action button to trigger analysis
                            actionButton("analyzeBtn", "Perform Analysis")
                          ),
                          mainPanel(
                            # Display analysis results
                            tabsetPanel(
                              tabPanel("MA Plot", plotOutput("maPlot_ma"),downloadButton("downloadmaPlot_ma", label = "Download",style="position: absolute; right:0; top:0; ")),
                              tabPanel("BOX Plot", plotOutput("boxPlot_affy"),downloadButton("downloadboxPlot_affy", label = "Download",style="position: absolute; right:0; top:0; ")),
                              tabPanel("PCA Plot", plotOutput("pcaPlot_affy"),downloadButton("downloadpcaPlot_affy", label = "Download",style="position: absolute; right:0; top:0; ")),
                              tabPanel("Volcano Plot", plotOutput("vol_affy"),downloadButton("downloadvol_affy", label = "Download",style="position: absolute; right:0; top:0; ")),
                              tabPanel("Top-Table",dataTableOutput("top_res"),downloadButton("downloadtop_res", label = "Download",style="position: absolute; right:0; top:0; "))
                            )
                          )
                        )),
               tabPanel("Spotted Microarray Analysis",
                        titlePanel("Spotted Microarray Analysis"),
                        
                        sidebarLayout(
                          sidebarPanel(
                            fileInput("file", "Upload the Sample Information File"),
                            fileInput("gal_file", "Upload the GAL File"),
                            textInput("ref_level", "Enter the Reference Level", placeholder = "wild type"),
                            actionButton("analyze_button", "Perform Analysis")
                          ),
                          mainPanel(
                            tabsetPanel(
                              navbarMenu("Image Plot",
                                         tabPanel("Image Plot [RED] ",plotOutput("image1b"),downloadButton("downloadimage1b", label = "Download",style="position: absolute; right:0; bottom:0; ")),
                                         tabPanel("Image Plot [RED] Background",plotOutput("image1"),downloadButton("downloadimage1", label = "Download",style="position: absolute; right:0; bottom:0; ")),
                                         tabPanel("Image Plot [GREEN] ",plotOutput("image2b"),downloadButton("downloadimage2b", label = "Download",style="position: absolute; right:0; bottom:0; ")),
                                         tabPanel("Image Plot [GREEN] Background",plotOutput("image2"),downloadButton("downloadimage2", label = "Download",style="position: absolute; right:0; bottom:0; "))),
                              navbarMenu("Scatter Plot",
                                         tabPanel("Raw Scatter Plot",plotOutput("rawsp"),downloadButton("downloadrawsp", label = "Download",style="position: absolute; right:0; bottom:0; ")),
                                         tabPanel("Log2 Trasformed Scatter Plot",plotOutput("logsp"),downloadButton("downloadlogsp", label = "Download",style="position: absolute; right:0; bottom:0; "))),
                              
                              tabPanel("MA Plot",plotOutput("ma_plot"),downloadButton("downloadma_plot", label = "Download",style="position: absolute; right:0; bottom:0; ")),
                              
                              navbarMenu("Density plot",
                                         tabPanel("Density plot (BEFORE Within array normalization)",plotOutput("wdbefore"),downloadButton("downloadwdbefore", label = "Download",style="position: absolute; right:0; bottom:0; ")),
                                         tabPanel("Density plot (AFTER Within array normalization)",plotOutput("wdafter"),downloadButton("downloadwdafter", label = "Download",style="position: absolute; right:0; bottom:0; "))),
                              navbarMenu("Tiploess plot",
                                         tabPanel("Tiploess plot (BEFORE Within array normalization)",plotOutput("tpbefore"),downloadButton("downloadtpbefore", label = "Download",style="position: absolute; right:0; bottom:0; ")),
                                         tabPanel("Tiploess plot (AFTER Within array normalization)",plotOutput("tpafter"),downloadButton("downloadtpafter", label = "Download",style="position: absolute; right:0; bottom:0; "))),
                              navbarMenu("Box plot",
                                         tabPanel("Box plot (BEFORE Between array normalization)",plotOutput("wbox"),downloadButton("downloadwbox", label = "Download",style="position: absolute; right:0; bottom:0; ")),
                                         tabPanel("Box plot (AFTER Between array normalization)",plotOutput("bbox"),downloadButton("downloadbbox", label = "Download",style="position: absolute; right:0; bottom:0; "))),
                              
                              tabPanel("Top-Table",dataTableOutput("top_data_spot"),
                                       downloadButton("downloadtopspot", "Download")),
                              
                              tabPanel("Volcano Plot",plotOutput("volspot"),downloadButton("downloadvolspot", label = "Download",style="position: absolute; right:0; bottom:0; "))
                            ))
                        )),
    ),
    navbarMenu("Enrichment Analysis",
               tabPanel("Gene Ontology (GO)",
                        titlePanel("GO Pathway Enrichment Analysis"),
                        
                        sidebarLayout(
                          sidebarPanel(
                            fileInput("gofile", "Upload the Data File",accept = c(".txt",".csv")),
                            selectInput("ontology", "Select the Ontology:",
                                        choices = c("Biological Process" = "BP", "Molecular Function" = "MF", "Cellular Component" = "CC")),
                            selectInput("organismInputgo", "Select the Organism",
                                        choices = c("Human" = "Homo sapiens",
                                                    "Mouse" = "Mus musculus",
                                                    "Rat" = "Rattus norvegicus",
                                                    "Fly" = "Drosophila melanogaster",
                                                    "Zebrafish" = "Danio rerio",
                                                    "Arabidopsis" = "Arabidopsis thaliana",
                                                    "Yeast" = "Saccharomyces cerevisiae",
                                                    "Worm" = "Caenorhabditis elegans",
                                                    "Bovine" = "Bos taurus",
                                                    "Chicken" = "Gallus gallus",
                                                    "Pig" = "Sus scrofa",
                                                    "Canine" = "Canis lupus familiaris",
                                                    "Rhesus" = "Macaca mulatta",
                                                    "Xenopus" = "Xenopus tropicalis",
                                                    "E. coli strain K12" = "Escherichia coli K-12",
                                                    "Chimp" = "Pan troglodytes",
                                                    "Anopheles" = "Anopheles gambiae",
                                                    "E. coli strain Sakai" = "Escherichia coli Sakai",
                                                    "Myxococcus xanthus DK 1622" = "Myxococcus xanthus DK 1622"
                                        ),
                                        selected = "Homo sapiens"),
                            actionButton("goButtongo", "Perform Enrichment Analysis")
                          ),
                          
                          mainPanel(
                            tabPanel("Go Table",dataTableOutput("goTable"),downloadButton("downloadgoTable", label = "Download"))
                          )
                        )),
               tabPanel("Kyoto Encyclopedia of Genes and Genomes (KEGG)",
                        titlePanel("KEGG Pathway Enrichment Analysis"),
                        
                        sidebarLayout(
                          sidebarPanel(
                            fileInput("keggfile", "Upload the Data File :",accept = c(".txt",".csv")),
                            textInput("organism_code", "Enter the Organism Code (e.g., hsa, mmu, ath):", "",placeholder = "ath"),
                            actionButton("goButton", "Perform Enrichment Analysis")
                          ),
                          
                          mainPanel(
                            tabPanel("kEGG Table",dataTableOutput("enrichmentTable"),downloadButton("downloadenrichmentTable", label = "Download"))
                          )
                        )),
    ),
    navbarMenu("Clustering",
               tabPanel("Hierarchical Clustering",
                        titlePanel("Hierarchical Clustering"),
                        sidebarLayout(
                          sidebarPanel(
                            fileInput("hfile", "Upload the Data File",
                                      accept = c(
                                        "text/csv",
                                        "text/comma-separated-values,text/plain",
                                        ".csv")),
                            actionButton("run_analysis", "Run Clustering")
                          ),
                          mainPanel(
                            tabPanel("HC",plotOutput("hclust_plot"),downloadButton("downloadhclust_plot", label = "Download",style="position: absolute; right:0; top:0; "))
                          )
                        )),
               tabPanel("K-means Clustering", 
                        titlePanel("K-Means Clustering"),
                        
                        sidebarLayout(
                          sidebarPanel(
                            fileInput("kfile", "Upload the Data File"),
                            numericInput("clusters", "Enter the Number of Clusters", value = 3, min = 1, max = 10, step = 1),
                            selectInput("gene_col1", "Select the Gene Names Column", choices = NULL, selected = NULL),
                            actionButton("analyzekmeans", "Run Clustering")
                          ),
                          mainPanel(
                            tabsetPanel(
                              tabPanel("Cluster Plot", plotOutput("cluster_plot"),downloadButton("downloadcluster_plot", label = "Download",style="position: absolute; right:0; top:0; ")),
                              tabPanel("Kmeans cluster Table", tableOutput("ktable"))
                            )
                          )
                        )),
               tabPanel("PCA Clustering", 
                        
                        titlePanel("Principal Component Analysis (PCA) "),
                        
                        sidebarLayout(
                          sidebarPanel(
                            fileInput("pcafile1", "Upload the Data File"),
                            selectInput("gene_col", "Select the Gene Names Column", choices = NULL, selected = NULL),
                            actionButton('pcaanalyze','Run Clustering')
                          ),
                          
                          mainPanel(
                            tabsetPanel(
                              tabPanel("Bi Plot", plotOutput("biplotpca"),downloadButton("downloadbiplotpca", label = "Download",style="position: absolute; right:0; top:0; ")),
                              tabPanel("Scree plot", plotOutput("scree_plotpca"),downloadButton("downloadscree_plotpca", label = "Download",style="position: absolute; right:0; top:0; ")),
                              tabPanel("PCA Plot", plotOutput("variable_importance_plotpca"),downloadButton("downloadvariable_importance_plotpca", label = "Download",style="position: absolute; right:0; top:0; ")),
                              tabPanel("Individuals plot", plotOutput("individuals_plotpca"),downloadButton("downloadindividuals_plotpca", label = "Download",style="position: absolute; right:0; top:0; "))
                            )
                          )   ) )
    ),
    navbarMenu("Machine Learning",
               
               tabPanel("Artificial Neural Network(ANN)",
                        titlePanel("Artificial Neural Network(ANN)"),
                                               sidebarLayout(
                          sidebarPanel(
                            fileInput("annfile", "Upload the Data File",
                                      accept = c(".csv")),
                            numericInput("p", "Enter the Percentage for data partitioning", value = 60),
                            textInput("hidden_layers", "Enter the Number of neurons in hidden layers (comma-separated)"),
                            actionButton("train_model", "Run Analysis")
                          ),
                          mainPanel(
                            tabsetPanel(
                              tabPanel("ANN Plot", plotOutput("neuralnet_plot"),downloadButton("downloadneuralnet_plot", label = "Download",style="position: absolute; right:0; top:0; ")),
                              tabPanel("Confusion Matrix", verbatimTextOutput("confusion_matrix22"),downloadButton("downloadconfusion_matrix22", label = "Download",style="position: absolute; right:0; top:0; "))
                            ) )  )),
               tabPanel("Naive Bayes",
                        titlePanel("Naive Bayes Classifier"),
                        
                        sidebarLayout(
                          sidebarPanel(
                            fileInput("nvfile", "Upload the Data File",
                                      accept = c( "text/csv", "text/comma-separated-values,text/plain",".csv")),
                            numericInput("pp", "Enter the Percentage for data partitioning", value = 60),
                            numericInput("features", "Enter the Feature to display", value = 1),
                            actionButton("run_naive", "Run Analysis")
                          ),
                          mainPanel(
                            tabsetPanel(
                              tabPanel("Missingness Plot",plotOutput("miss_plot"),downloadButton("downloadmiss_plot", label = "Download",style="position: absolute; right:0; top:0; ")),
                              tabPanel("Pairs Plot",plotOutput("ggp_plot"),downloadButton("downloadggp_plot", label = "Download",style="position: absolute; right:0; top:0; ")),
                              tabPanel("Feature Plot",plotOutput("feature_plot"),downloadButton("downloadfeature_plot", label = "Download",style="position: absolute; right:0; top:0; ")),
                              tabPanel("Roc Plot", plotOutput("roc_plot"),downloadButton("downloadroc_plot", label = "Download",style="position: absolute; right:0; top:0; ")),
                              tabPanel("Confusion Matrix", verbatimTextOutput("confusion_matrixx3"),downloadButton("downloadconfusion_matrixx3", label = "Download",style="position: absolute; right:0; top:0; "))
                            )
                          )
                        )
               ),
               tabPanel("Regression",
                        titlePanel("LASSO and Ridge Regression"),
                        
                        sidebarLayout(
                          sidebarPanel(
                            fileInput("regfile", "Upload the Data File",
                                      accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
                            actionButton("regButton", "Run Analysis")
                          ),
                          mainPanel(
                            tabsetPanel(
                              tabPanel("LASSO Coefficients", plotOutput("lasso_plot"),downloadButton("downloadlasso_plot", label = "Download",style="position: absolute; right:0; top:0; ")),
                              tabPanel("Ridge Coefficients", plotOutput("ridge_plot"),downloadButton("downloadridge_plot", label = "Download",style="position: absolute; right:0; top:0; "))
                            ) ) )),
    )
    
  )
)