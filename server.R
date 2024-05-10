# server.R

library(shiny)
options(shiny.maxRequestSize = 30*1024^2) # Set max request size to 30MB

options(repos = "https://cran.rstudio.com/")

#if (!requireNamespace("BiocManager"))
#install.packages("BiocManager")
#BiocManager::install(c("limma", "edgeR","org.Mm.eg.db","gplots", "Glimma" ,"RColorBrewer", "NMF", "BiasedUrn", "GO.db"), ask = FALSE)
#install.packages("DT",ask=FALSE)
#if (!requireNamespace("BiocManager", quietly = TRUE)) {
#  install.packages("BiocManager")}
#BiocManager::install("AnnotationDbi")


# Function to install the package if not already installed
install_annotation_package <- function(org_name) {
  package_name <- switch(org_name,
                         "Homo sapiens" = "org.Hs.eg.db",
                         "Mus musculus" = "org.Mm.eg.db",
                         "Rattus norvegicus" = "org.Rn.eg.db",
                         "Drosophila melanogaster" = "org.Dm.eg.db",
                         "Danio rerio" = "org.Dr.eg.db",
                         "Arabidopsis thaliana" = "org.At.tair.db",
                         "Saccharomyces cerevisiae" = "org.Sc.sgd.db",
                         "Caenorhabditis elegans" = "org.Ce.eg.db",
                         "Bos taurus" = "org.Bt.eg.db",
                         "Gallus gallus" = "org.Gg.eg.db",
                         "Sus scrofa" = "org.Ss.eg.db",
                         "Canis lupus familiaris" = "org.Cf.eg.db",
                         "Macaca mulatta" = "org.Mmu.eg.db",
                         "Xenopus tropicalis" = "org.Xl.eg.db",
                         "Escherichia coli K-12" = "org.EcK12.eg.db",
                         "Pan troglodytes" = "org.Pt.eg.db",
                         "Anopheles gambiae" = "org.Ag.eg.db",
                         "Escherichia coli Sakai" = "org.EcSakai.eg.db",
                         "Myxococcus xanthus DK 1622" = "org.Mxanthus.db",
                         stop("Annotation database not available for the specified organism.")
  )
  
  if (!package_name %in% rownames(installed.packages())) {
    install.packages(package_name, repos = "https://bioconductor.org/packages/3.14/bioc")
  }
  
  # Load the library
  library(package_name, character.only = TRUE)
}

library(shiny)
library(AnnotationDbi)
library(DT)
library(affy)
library(oligo)
library(edgeR)
library(limma)
library(Glimma)
#library(org.Mm.eg.db)
library(gplots)
library(RColorBrewer)
library(NMF)
library(GO.db)
library(caTools) 
library(randomForest) 
library(cluster)
library(ggplot2)
library(dplyr)
library(tidyr)
library(neuralnet)
library(caret)
library(tidyverse)
library(factoextra) 
library(glmnet) 
library(e1071)
library(caret)
library(pROC)
library(GGally)
library(Amelia)
library(mice)
library(KEGGREST)
library(png)
library(topGO) 

server <- function(input, output, session) { 
  
  #rna seq
  observeEvent(input$runAnalysisBtn, {
    # Check if both files are uploaded
    if (is.null(input$countFile) || is.null(input$sampleInfoFile)) {
      return()
    }
    organism <- input$organismInput
    install_annotation_package(organism)
    
    org_db <- switch(organism,
                     "Homo sapiens" = org.Hs.eg.db,
                     "Mus musculus" = org.Mm.eg.db,
                     "Rattus norvegicus" = org.Rn.eg.db,
                     "Drosophila melanogaster" = org.Dm.eg.db,
                     "Danio rerio" = org.Dr.eg.db,
                     "Arabidopsis thaliana" = org.At.tair.db,
                     "Saccharomyces cerevisiae" = org.Sc.sgd.db,
                     "Caenorhabditis elegans" = org.Ce.eg.db,
                     "Bos taurus" = org.Bt.eg.db,
                     "Gallus gallus" = org.Gg.eg.db,
                     "Sus scrofa" = org.Ss.eg.db,
                     "Canis lupus familiaris" = org.Cf.eg.db,
                     "Macaca mulatta" = org.Mmu.eg.db,
                     "Xenopus tropicalis" = org.Xl.eg.db,
                     "Escherichia coli K-12" = org.EcK12.eg.db,
                     "Pan troglodytes" = org.Pt.eg.db,
                     "Anopheles gambiae" = org.Ag.eg.db,
                     "Escherichia coli Sakai" = org.EcSakai.eg.db,
                     "Myxococcus xanthus DK 1622" = org.Mxanthus.db,
                     stop("Annotation database not available for the specified organism.")
    )
    
    # Read count data
    seqdata <- read.delim(input$countFile$datapath, stringsAsFactors = FALSE)
    processed_sampleinfo <- read.delim(input$sampleInfoFile$datapath, stringsAsFactors = TRUE)
    
     # Process the sampleinfo data
    processed_sampleinfo <- processed_sampleinfo[!duplicated(processed_sampleinfo$FileName), ]
    
    processed_sampleinfo <- processed_sampleinfo %>%
      setNames(c("NewNameOfFile", "NewNameOfSample", "NewGroup1", "NewGroup2"))
    
    # Process the sampleinfo data
    #processed_sampleinfo <- sampleinfo[-1, ] %>%
    #  setNames(c("NewNameOfFile", "NewNameOfSample", "NewGroup1", "NewGroup2"))
    
    countdata <- seqdata[,-1]
    rownames(countdata) <- seqdata[,1]
    colnames(countdata) <- substr(colnames(countdata), start = 1, stop = 7)
    y <- DGEList(countdata) 
    
    #group <- paste(sampleinfo$CellType, sampleinfo$Status, sep = ".")
    group <- paste(processed_sampleinfo$NewGroup1, processed_sampleinfo$NewGroup2, sep = ".")
    group <- factor(group)
    y$samples$group <- group
    
    # Annotation
    keys <- rownames(y$counts)
    columns <- c("ENTREZID", "GENENAME")
    ann <- AnnotationDbi::select(org_db, keys = keys, columns = columns, keytype = "ENTREZID")
    y$genes <- ann
    
    myCPM <- cpm(countdata)
    thresh <- myCPM > 0.5
    keep <- rowSums(thresh) >= 2
    y <- y[keep, keep.lib.sizes=FALSE]
    logcounts <- cpm(y,log=TRUE)
    var_genes <- apply(logcounts, 1, var)
    select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
    highly_variable_lcpm <- logcounts[select_var,]
    y <- calcNormFactors(y)
    design <- model.matrix(~ 0 + group)
    colnames(design) <- levels(group)
    
    v <- voom(y,design)
    fit <- lmFit(v)
    
    output$contrastDropdown1 <- renderUI({
      selectInput("contrastColumn1", "Select First Column for Contrast:",
                  choices = colnames(design))
    })
    
    output$contrastDropdown2 <- renderUI({
      selectInput("contrastColumn2", "Select Second Column for Contrast:",
                  choices = colnames(design))
    }) 
   
    # Define a reactive expression to handle contrast calculation
    observeEvent(input$runContrastBtn, { 
      req(input$contrastColumn1)
      req(input$contrastColumn2)
      
      astr <- paste(input$contrastColumn1 ,"-", input$contrastColumn2)
      
      prestr = "makeContrasts("
      poststr = ",levels=design)" 
      commandstr = paste(prestr,astr,poststr,sep="") 
      
      cont.matrix=eval(parse(text=commandstr))
      
      print(cont.matrix) 
      fit.cont <- contrasts.fit(fit, contrast = cont.matrix)
      #print(fit.cont) 
      fit.cont <- eBayes(fit.cont)
      summa.fit <- decideTests(fit.cont)
      top <- topTable(fit.cont, coef = 1, sort.by = "p", n = Inf)   
      #print(top) 
      
      
      #HEATMAP 
      plotheat = function(){
        var_genes <- apply(logcounts, 1, var)
        select_var <- names(sort(var_genes, decreasing = TRUE))[1:500]
        highly_variable_lcpm <- logcounts[select_var,]
        mypalette <- brewer.pal(11, "RdYlBu")
        morecols <- colorRampPalette(mypalette)
        #col.cell <- c("purple", "orange")[sampleinfo$CellType]
        #group_labels <- paste(sampleinfo$CellType, sampleinfo$Status, sep = ".")
        col.cell <- c("purple", "orange")[processed_sampleinfo$NewGroup1]
        group_labels <- paste(processed_sampleinfo$NewGroup1, processed_sampleinfo$NewGroup2, sep = ".")
        heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most variable genes across samples",ColSideColors=col.cell,scale="row")
        #heatmap.2(highly_variable_lcpm, col = rev(morecols(50)), trace = "none", main = "Top 500 most variable genes across samples", ColSideColors = col.cell, labCol = group_labels, scale = "row")
      }
      
      #BOX PLOT 
      
      plotboxrna= function(){ 
        par(mfrow=c(1,2))
        boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2,main="Unnormalised logCPM")
        abline(h=median(logcounts),col="blue")
        boxplot(v$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
        abline(h=median(v$E),col="blue")
      }
      
      output$box <- renderPlot({ 
        plotboxrna()
      }, height = 600)
      
      
      #VOLCANO PLOT
      #status=summa.fit[,1],
      output$vol <- renderPlot({
        group2 <- group
        #levels(group2) <- c("basal.lactate","basal.preg","basal.virgin","lum.lactate", "lum.preg", "lum.virgin")
        # Concatenate column names with commas
        levels_string <- paste(colnames(design), collapse = ",") 
        # Split the concatenated string into separate elements
        levels_vector <- unlist(strsplit(levels_string, ","))
        # Assign the vector of separate elements as levels to group2
        levels(group2) <- levels_vector
        glXYPlot(x=fit.cont$coefficients[,1], y=fit.cont$lods[,1],
                 xlab="logFC", ylab="-log10(p-value)", main="Volcano Plot",
                 counts=v$E, groups=group2, status=summa.fit[,1],
                 anno=fit.cont$genes, side.main="ENTREZID", folder="volcano")
      })
      
      
      output$cont_data <- renderTable({
        cont.matrix
      })
      output$top_data <- DT::renderDataTable({
        DT::datatable(top, options = list(orderClasses = TRUE))
      })
      
      output$voom <- renderPlot({
        v <- voom(y,design, plot = T)
      }, height = 600)
      
      output$heat_map <- renderPlot({
        plotheat()
      }, height = 600)
      
      output$mds_map <- renderPlot({
        #labels <- paste(sampleinfo$SampleName, sampleinfo$CellType, sampleinfo$Status)
        labels <- paste(processed_sampleinfo$NewNameOfSample, processed_sampleinfo$NewGroup1, processed_sampleinfo$NewGroup2)
        glMDSPlot(y, labels=labels, groups=group, folder="mds")
        
      })
      
      
      output$downloadHeatmap <- downloadHandler(
        filename = "heatmap.png",
        content = function(file) {
          png(file)
          plotheat()
          dev.off()
        })   
      
      
      output$downloadVoom <- downloadHandler(
        filename = "voomplot.png",
        content = function(file) {
          png(file)
          v <- voom(y,design, plot = T)
          dev.off()
        }
      )
      
      output$downloadBoxplot <- downloadHandler(
        filename = "boxplot.png",
        content = function(file) {
          png(file)
          plotboxrna()
          dev.off()
        }
      )
      
      
      output$downloadtop1 <- downloadHandler(
        filename = function() {
          "top_table.csv"
        },
        content = function(file) {
          write.csv(top, file)
        }
      ) 
      
      
      
      output$top_data <- DT::renderDataTable({
        DT::datatable(top, options = list(orderClasses = TRUE)) 
      })
      
    })
    
    output$orig_data <- renderTable({
      design
    })
    
    output$downloaddesign <- downloadHandler(
      filename = function() {
        paste(design, Sys.Date(), ".csv", sep='')
      },
      content = function(file) {
        write.csv(design, file)
      }
    )
    
    
    
  })
  
  #spotted array
  
  observeEvent(input$analyze_button, {
    
    if (is.null(input$file) || is.null(input$gal_file)) {
      return()
    }
    
    targets <- readTargets(input$file$datapath)
    rg <- read.maimages(targets, source = "spot")
    rg$genes <- readGAL(input$gal_file$datapath)
    rg$printer <- getLayout(rg$genes)
    red <- (rg$R[, 1] - rg$Rb[, 1])
    green <- (rg$G[, 1] - rg$Gb[, 1])
    MA <- normalizeWithinArrays(rg, method = "none")
    MA1 <- normalizeWithinArrays(rg)
    MA2 <- normalizeBetweenArrays(MA1, method = "scale")
    #design <- modelMatrix(targets, ref = "wild type")
    # Get reference level from user input
    ref_level <- input$ref_level
    
    # Use reference level in modelMatrix function
    design <- modelMatrix(targets, ref = ref_level)
    fit <- lmFit(MA2, design)
    fit2 <- eBayes(fit)
    summa.fit <- decideTests(fit2)
    topspot <- topTable(fit2,adjust='BH',sort.by="p",n = Inf) 
    
    plotimage11 = function(){
      imageplot(rg$R[, 1], rg$printer, low = "white", high = "red",main = "Image Plot [RED] with Background")
    } 
    
    
    # Plot image
    output$image1b <- renderPlot({
      plotimage11()
    })
    
    output$downloadimage1b <- downloadHandler(
      filename = "imageplot1b.png",
      content = function(file) {
        png(file)
        plotimage11()
        dev.off()
      }
    )
    
    plotimage22 = function(){
      imageplot(rg$Rb[, 1], rg$printer, low = "white", high = "red",main = "Image Plot [RED]")
    } 
    
    output$image1 <- renderPlot({
      plotimage22()
    })
    
    output$downloadimage1 <- downloadHandler(
      filename = "imageplot1.png",
      content = function(file) {
        png(file)
        plotimage22()
        dev.off()
      }
    )
    
    
    plotimage33 = function(){
      imageplot(rg$G[, 1], rg$printer, low = "white", high = "green", main = "Image Plot [GREEN] with Background")
    } 
    
    output$image2b <- renderPlot({
      plotimage33()
    })
    
    output$downloadimage2b <- downloadHandler(
      filename = "imageplot2.png",
      content = function(file) {
        png(file)
        plotimage33()
        dev.off()
      }
    )
    
    plotimage44 = function(){
      imageplot(rg$Gb[, 1], rg$printer, low = "white", high = "green" , main = "Image Plot [GREEN]")
    } 
    
    output$image2 <- renderPlot({
      plotimage44()
    })
    
    output$downloadimage2 <- downloadHandler(
      filename = "imageplot2b.png",
      content = function(file) {
        png(file)
        plotimage44()
        dev.off()
      }
    )
    #scatter plot 
    plotrawsp = function(){
      plot(green,red,main="Raw scatter plot")
    } 
    
    output$rawsp <- renderPlot({
      plotrawsp()
    }) 
    
    output$downloadrawsp <- downloadHandler(
      filename = "rawscatter.png",
      content = function(file) {
        png(file)
        plotrawsp()
        dev.off()
      }
    )
    
    
    plotlogsp = function(){
      plot(log2(rg$R),log2(rg$G),main="Log Transform scatter plot")
    } 
    
    output$logsp <- renderPlot({
      plotlogsp()
    })
    
    output$downloadlogsp <- downloadHandler(
      filename = "LOGScatter.png",
      content = function(file) {
        png(file)
        plotlogsp()
        dev.off()
      }
    ) 
    
    # Plot MA plot
    
    plotma_plot = function(){
      plot((green + red) / 2, red - green, main = "MA Plot")
    } 
    output$ma_plot <- renderPlot({
      plotma_plot()
    })
    
    output$downloadma_plot <- downloadHandler(
      filename = "MA.png",
      content = function(file) {
        png(file)
        plotma_plot()
        dev.off()
      }
    ) 
    
    # Plot density plot
    # before
    plotwdbefore= function(){
      plotDensities(MA , main = "Within Array Density plot [BEFORE]")
    }
    
    output$wdbefore <- renderPlot({
      plotwdbefore()
    })
    
    output$downloadwdbefore <- downloadHandler(
      filename = "density plot.png",
      content = function(file) {
        png(file)
        plotwdbefore()
        dev.off()
      }
    ) 
    
    output$tpbefore <- renderPlot({
      plotPrintTipLoess(MA)
    }) 
    
    
    output$downloadtpbefore <- downloadHandler(
      filename = "tiploess1.png",
      content = function(file) {
        png(file)
        plotPrintTipLoess(MA)
        dev.off()
      }
    ) 
    
    # after
    
    output$wdafter <- renderPlot({
      plotDensities(MA1 , main = "Within Array Density plot [AFTER]")
    }) 
    
    output$downloadwdafter <- downloadHandler(
      filename = "Within Array Density plot.png",
      content = function(file) {
        png(file)
        plotDensities(MA1 , main = "Within Array Density plot [AFTER]")
        dev.off()
      }
    ) 
    
    output$tpafter <- renderPlot({
      plotPrintTipLoess(MA1)
    }) 
    
    output$downloadtpafter <- downloadHandler(
      filename = "after plot2.png",
      content = function(file) {
        png(file)
        plotPrintTipLoess(MA1)
        dev.off()
      }
    ) 
    
    # Plot boxplot
    output$wbox <- renderPlot({
      boxplot(MA1$M , main = "Within Array Box plot")
    })
    
    output$downloadwbox <- downloadHandler(
      filename = "boxplot.png",
      content = function(file) {
        png(file)
        boxplot(MA1$M , main = "Within Array Box plot")
        dev.off()
      }
    ) 
    
    
    output$bbox <- renderPlot({
      boxplot(MA2$M~col(MA2$M),names=colnames(MA2$M) , main = "Between Array  Box plot")
    })
    
    output$downloadbbox <- downloadHandler(
      filename = "boxplot2.png",
      content = function(file) {
        png(file)
        boxplot(MA2$M~col(MA2$M),names=colnames(MA2$M) , main = "Between Array  Box plot")
        dev.off()
      }
    ) 
    
    output$top_data_spot <- DT::renderDataTable({
      DT::datatable(topspot, options = list(orderClasses = TRUE))
    })
    
    output$downloadtopspot <- downloadHandler(
      filename = function() {
        "top_table_spot.csv"
      },
      content = function(file) {
        write.csv(topspot, file)
      }
    ) 
    
    
    
    
    
    # Plot volcano plot
    output$volspot <- renderPlot({
      volcanoplot(fit2, highlight = 20, names = fit2$genes$Name , main = "Volcano Plot")
      # })
      
      #glXYPlot(x=fit2$coefficients[,1], y=fit2$lods[,1],
      # xlab="logFC", ylab="-log10(p-value)", main="B.PregVsLac",
      #status=summa.fit[,1],folder="volcano")
    }) 
    
    output$downloadvolspot <- downloadHandler(
      filename = "volcano.png",
      content = function(file) {
        png(file)
        volcanoplot(fit2, highlight = 20, names = fit2$genes$Name , main = "Volcano Plot")
        dev.off()
      }
    ) 
    
  }) 
  
  #affymetrix array
  
  observeEvent(input$analyzeBtn, {
    # Check if CEL files are provided
    req(input$celDir)
    req(input$groups)
    #req(input$groups)
    # Read CEL files
    affyRaw <- read.celfiles(input$celDir$datapath)
    eset <- oligo::rma(affyRaw)
    # Fetch phenotype data
    phenotype_data <- pData(eset)
    # Perform differential expression analysis
    # Split user input into individual group labels
    groups <- unlist(strsplit(input$groups, ","))
    #Groups <- c("Tumor","Tumor","Tumor","Tumor","Tumor","Normal","Normal","Normal","Normal","Normal")
    # Create a design matrix based on user-provided groups
    design <- model.matrix(~factor(groups))
    #colnames(design) <- c("Tumor","TumorvsNormal")
    colnames(design) <- unique(groups)
    # Fit linear model and perform differential expression analysis
    fit <- lmFit(eset, design)
    fit <- eBayes(fit)
    summa.fit <- decideTests(fit)
    res <- topTable(fit, number = Inf, adjust.method = "none", coef = 1)
    
    # Output analysis results
    output$top_res <- DT::renderDataTable({
      DT::datatable(res, options = list(orderClasses = TRUE))
    }) 
    
    output$downloadtop_res <- downloadHandler(
      filename = function() {
        "top_table.csv"
      },
      content = function(file) {
        write.csv(res, file)
      }
    ) 
    
    plotmaPlot_ma= function(){ 
      plot(res$AveExpr, res$logFC, 
           col = ifelse(abs(res$logFC) > 1 & res$P.Value < 0.05, "red", "black"), 
           pch = 20, 
           cex = 0.5,
           main = "MA Plot",
           xlab = "Average Expression",
           ylab = "Log2 Fold Change")
      abline(h = 0, col = "red", lty = 2)
      abline(v = 0, col = "blue", lty = 2)
    }
    
    
    output$maPlot_ma<- renderPlot({
      plotmaPlot_ma()
    })
    
    output$downloadmaPlot_ma <- downloadHandler(
      filename = "ma.png",
      content = function(file) {
        png(file)
        plotmaPlot_ma()
        dev.off()
      }
    )
    
    
    plotboxPlot_affy <- function(groups, eset) {
      # Get unique groups
      unique_groups <- unique(groups)
      
      # Generate a sequence of colors based on the number of unique groups
      colors <- rainbow(length(unique_groups))
      
      # Create a color vector corresponding to each sample's group
      sample_colors <- colors[match(groups, unique_groups)]
      
      # Plot boxplot
      boxplot(exprs(eset),
              col = sample_colors,
              main = "Boxplot of Expression Data",
              xlab = "Groups",
              ylab = "Expression")
    }
    
    
    
    output$boxPlot_affy<- renderPlot({ 
      plotboxPlot_affy(unlist(strsplit(input$groups, ",")), eset)
    }) 
    
    output$downloadboxPlot_affy <- downloadHandler(
      filename = "boxplot.png",
      content = function(file) {
        png(file)
        plotboxPlot_affy()
        dev.off()
      }
    )
    
    
    plotpcaPlot_affy= function(){ 
      pca_data <- prcomp(exprs(eset))
      # Plot PCA
      plot(pca_data$x[,1], pca_data$x[,2], 
           xlab = "Principal Component 1", 
           ylab = "Principal Component 2", 
           main = "PCA Plot")
      
    }
    
    output$pcaPlot_affy<- renderPlot({ 
      plotpcaPlot_affy()
    }) 
    
    output$downloadpcaPlot_affy <- downloadHandler(
      filename = "boxplot.png",
      content = function(file) {
        png(file)
        plotpcaPlot_affy()
        dev.off()
      }
    )
    
    
    plotvol_affy= function(){ 
      group2 <- groups
      #levels(group2) <- c("Tumor","TumorvsNormal")
      levels(group2) <- unique(groups)
      main_title <- paste(unique(groups), collapse = " vs ")
      glXYPlot(x=fit$coefficients[,1], y=fit$lods[,1],
               xlab="logFC", ylab="-log10(p-value)", main=main_title,
               groups=group2, status=summa.fit[,1],
               anno=fit$genes, side.main="ID", folder="Volcano")
      
    }
    
    
    # Volcano Plot using glXYPlot from Glimma
    output$vol_affy <- renderPlot({
      plotvol_affy()
    })
    
    
    output$downloadvol_affy <- downloadHandler(
      filename = "boxplot.png",
      content = function(file) {
        png(file)
        plotvol_affy()
        dev.off()
      }
    )
    
  }) 
  
  #enrichment
  #go 
  
  org_packages <- list(
    "Homo sapiens" = "org.Hs.eg.db",
    "Mus musculus" = "org.Mm.eg.db",
    "Rattus norvegicus" = "org.Rn.eg.db",
    "Drosophila melanogaster" = "org.Dm.eg.db",
    "Danio rerio" = "org.Dr.eg.db",
    "Arabidopsis thaliana" = "org.At.tair.db",
    "Saccharomyces cerevisiae" = "org.Sc.sgd.db",
    "Caenorhabditis elegans" = "org.Ce.eg.db",
    "Bos taurus" = "org.Bt.eg.db",
    "Gallus gallus" = "org.Gg.eg.db",
    "Sus scrofa" = "org.Ss.eg.db",
    "Canis lupus familiaris" = "org.Cf.eg.db",
    "Macaca mulatta" = "org.Mmu.eg.db",
    "Xenopus tropicalis" = "org.Xl.eg.db",
    "Escherichia coli K-12" = "org.EcK12.eg.db",
    "Pan troglodytes" = "org.Pt.eg.db",
    "Anopheles gambiae" = "org.Ag.eg.db",
    "Escherichia coli Sakai" = "org.EcSakai.eg.db",
    "Myxococcus xanthus DK 1622" = "org.Mxanthus.db"
  )
  
  observeEvent(input$goButtongo, {
    req(input$gofile) 
    req(input$organismInputgo)  
    
    organismgo <- input$organismInputgo
    
    # Retrieve annotation package name based on selected organism
    org_dbgo <- org_packages[[organismgo]]
    print(org_dbgo)  # Print the selected annotation package name to the console for debugging
    
    # Read the data file
    data <- read.delim(input$gofile$datapath)
    
    # Extract P.Values and Gene names
    geneList <- data$P.Value
    names(geneList) <- data$Gene
    
    # Create topGOdata object
    GOdata <- new("topGOdata",
                  ontology = input$ontology,
                  allGenes = geneList,
                  geneSelectionFun = function(x) x,
                  annot = annFUN.org,
                  mapping = org_dbgo)  # Use get() to retrieve the package by name
    
    # Run Kolmogorov-Smirnov testing
    resultKS <- runTest(GOdata, algorithm = "weight01", statistic = "ks")
    
    # Generate a table of significant GO terms
    tab <- GenTable(GOdata, raw.p.value = resultKS, topNodes = length(resultKS@score), numChar = 120)
    
    # Display the table
    output$goTable <- DT::renderDataTable({
      DT::datatable(tab, options = list(orderClasses = TRUE))
    }) 
    
    
    output$downloadgoTable <- downloadHandler(
      filename = function() {
        "go_table.csv"
      },
      content = function(file) {
        write.csv(tab, file)
      }
    ) 
    
    
  }) 
  
  #kegg
  
  observeEvent(input$goButton, {
    req(input$keggfile)
    req(input$organism_code)
    # Read the DE data file
    print("Reading DE data file...")
    DEtable <- read.delim(input$keggfile$datapath)
    
    
    geneList <- DEtable$P.Value
    names(geneList) <- DEtable$Gene
    # Pull all pathways for Arabidopsis thaliana
    print("Pulling pathways...")
    pathways.list <- keggList("pathway", input$organism_code)
    
    pathway.codes <- sub("path:", "", names(pathways.list)) 
    genes.by.pathway <- sapply(pathway.codes,
                               function(pwid){
                                 pw <- keggGet(pwid)
                                 if (is.null(pw[[1]]$GENE)) return(NA)
                                 pw2 <- pw[[1]]$GENE[c(TRUE,FALSE)] # may need to modify this to c(FALSE, TRUE) for other organisms
                                 pw2 <- unlist(lapply(strsplit(pw2, split = ";", fixed = T), function(x)x[1]))
                                 return(pw2)
                               }
    )
    pVals.by.pathway <- t(sapply(names(genes.by.pathway),
                                 function(pathway) {
                                   pathway.genes <- genes.by.pathway[[pathway]]
                                   list.genes.in.pathway <- intersect(names(geneList), pathway.genes)
                                   list.genes.not.in.pathway <- setdiff(names(geneList), list.genes.in.pathway)
                                   scores.in.pathway <- geneList[list.genes.in.pathway]
                                   scores.not.in.pathway <- geneList[list.genes.not.in.pathway]
                                   if (length(scores.in.pathway) > 0){
                                     p.value <- wilcox.test(scores.in.pathway, scores.not.in.pathway, alternative = "less")$p.value
                                   } else{
                                     p.value <- NA
                                   }
                                   return(c(p.value = p.value, Annotated = length(list.genes.in.pathway)))
                                 }))
    
    # Assemble output table
    print("Assembling output table...")
    outdat <- data.frame(pathway.code = rownames(pVals.by.pathway))
    
    outdat$pathway.name <- pathways.list[outdat$pathway.code]
    outdat$p.value <- pVals.by.pathway[,"p.value"]
    outdat$Annotated <- pVals.by.pathway[,"Annotated"]
    outdat <- outdat[order(outdat$p.value),]
    #print(outdat)
    # Display the enrichment table
    output$enrichmentTable <-DT::renderDataTable({ 
      
      DT::datatable(outdat, options = list(orderClasses = TRUE))
    }) 
    
    output$downloadenrichmentTable <- downloadHandler(
      filename = function() {
        "KEGG_table.csv"
      },
      content = function(file) {
        write.csv(outdat, file)
      }
    ) 
    
    
  })
  
  
  
  #clustering
  
  data_reactive <- reactive({
    req(input$hfile)
    inFile <- input$hfile
    if (is.null(inFile)) return(NULL)
    df <- read.csv(inFile$datapath, row.names = 1, header = FALSE, sep = ",")
    df <- t(df)
    rownames(df) <- df[, 1]
    df <- df[, -1]
    mode(df) <- "numeric"
    return(df)
  }) 
  
  observeEvent(input$run_analysis, {
    
    plothclust_plot= function(){
      data_df <- data_reactive()
      if (!is.null(data_df)) {
        d <- dist(as.matrix(data_df))
        hc <- hclust(d)
        plot(hc)
      }
    }
    
    
    output$hclust_plot <- renderPlot({
      plothclust_plot()
    }) 
    
    output$downloadhclust_plot <- downloadHandler(
      filename = "hc_plot.png",
      content = function(file) {
        png(file)
        plothclust_plot()
        dev.off()
      }
    )
    
    
    
  })
  
  #kmeans
  observeEvent(input$kfile, {
    req(input$kfile)
    
    # Read the uploaded CSV file
    data <- read.csv(input$kfile$datapath)
    
    # Update choices for gene_col selectInput
    updateSelectInput(session, "gene_col1", choices = colnames(data))
  })
  observeEvent(input$analyzekmeans, {
    req(input$kfile)
    req(input$clusters)
    data = read.csv(input$kfile$datapath, header = TRUE, stringsAsFactors = FALSE)
    
    
    gene_col1 <- input$gene_col1
    # Extract gene names from the selected column
    gene_names <- data[[gene_col1]]
    
    # Remove the selected column from the data
    data <- data[, -which(names(data) %in% c(gene_col1))]
    
    # Transpose the data
    data <- t(data)
    
    # Perform k-means clustering
    set.seed(123) # for reproducibility
    k <- input$clusters
    km <- kmeans(data, centers = k)
    #km = kmeans(df, centers = k, nstart = 25) 
    
    
    plotcluster_plot= function(){
      p <- fviz_cluster(km, data = data, geom = "point", stand = FALSE, 
                        ggtheme = theme_minimal(), main = paste("k-means clustering (k =", k, ")"),
                        label = "none")
      
      # Extract cluster centers
      centers <- as.data.frame(km$centers)
      
      # Add sample names as labels
      p + geom_text(data = centers, aes(x = V1, y = V2, label = rownames(centers)), vjust = -1)
    }
    
    # Plot the clusters
    output$cluster_plot <- renderPlot({
      plotcluster_plot()
    }) 
    
    output$downloadcluster_plot <- downloadHandler(
      filename = "kmean_plot.png",
      content = function(file) {
        png(file)
        plotcluster_plot()
        dev.off()
      }
    )
    
    output$ktable <- renderTable({ 
      # Print cluster assignments
      data_with_clusters <- cbind(gene_names, cluster = km$cluster)
      data_with_clusters
    }) 
    
    
    
  }) 
  
  #pca - cluster
  
  observeEvent(input$pcafile1, {
    req(input$pcafile1)
    
    # Read the uploaded CSV file
    data <- read.csv(input$pcafile1$datapath)
    
    # Update choices for gene_col selectInput
    updateSelectInput(session, "gene_col", choices = colnames(data))
  })
  
  observeEvent(input$pcaanalyze, {
    req(input$pcafile1)
    req(input$gene_col)
    data <- read.csv(input$pcafile1$datapath)
    
    gene_col <- input$gene_col
    # Extract gene names from the selected column
    gene_names <- data[[gene_col]]
    
    # Remove the selected column from the data
    data <- data[, -which(names(data) %in% c(gene_col))]
    
    # Calculate principal components
    results <- prcomp(data, scale = TRUE)
    
    # Reverse the signs
    results$rotation <- -1 * results$rotation 
    rownames(results$rotation) <- colnames(data)
    
    # Display principal components
    results$rotation 
    
    # Reverse the signs of the scores
    results$x <- -1 * results$x
    
    output$biplotpca<- renderPlot({
      biplot(results, scale = 0)
      legend("bottomright", legend = c("Group 1", "Group 2"), col = c("black", "pink"), pch = 1, title = "Groups", cex = 0.8)
    })
    
    output$downloadbiplotpca <- downloadHandler(
      filename = "biplot.png",
      content = function(file) {
        png(file)
        biplot(results, scale = 0)
        legend("bottomright", legend = c("Group 1", "Group 2"), col = c("black", "pink"), pch = 1, title = "Groups", cex = 0.8)
        dev.off()
      }
    ) 
    
    plotscree_plotpca= function(){
      var_explained = results$sdev^2 / sum(results$sdev^2)
      plot(1:length(var_explained), var_explained, type = "b", 
           xlab = "Principal Component", ylab = "Proportion of Variance Explained",
           main = "Scree Plot")
    }
    
    output$scree_plotpca<- renderPlot({
      plotscree_plotpca()
    }) 
    
    output$downloadscree_plotpca <- downloadHandler(
      filename = "screeplot.png",
      content = function(file) {
        png(file)
        plotscree_plotpca()
        dev.off()
      }
    )
    
    plotvariable_importance_plotpca = function(){
      variance_ratio <- results$sdev^2 / sum(results$sdev^2)
      var_importance <- t(results$rotation[,1:2]) * sqrt(variance_ratio[1:2])
      barplot(abs(var_importance), beside = TRUE, col = c("blue", "red"),
              names.arg = colnames(data), 
              main = "Variable Importance Plot",
              xlab = "Principal Component", ylab = "Variable Importance")
      legend("topright", legend = c("PC1", "PC2"), fill = c("blue", "red"), 
             x = "topright", y = 1, bty = "n")
    }
    
    
    output$variable_importance_plotpca<- renderPlot({
      plotvariable_importance_plotpca()
    }) 
    
    output$downloadvariable_importance_plotpca <- downloadHandler(
      filename = "variable.png",
      content = function(file) {
        png(file)
        plotvariable_importance_plotpca()
        dev.off()
      }
    )
    
    plotindividuals_plotpca = function(){
      # Individuals Plots
      plot(results$x[,1], results$x[,2], 
           xlab = "Principal Component 1", ylab = "Principal Component 2",
           main = "Individuals Plot", pch = 19, col = 1:length(gene_names))  
      
      # Add legends for gene names
      legend("topright", legend = gene_names, col = 1:length(gene_names), pch = 19, title = "Gene Names", cex = 0.8)
      
    }
    
    
    output$individuals_plotpca<- renderPlot({
      plotindividuals_plotpca()
    }) 
    
    output$downloadindividuals_plotpca <- downloadHandler(
      filename = "individual.png",
      content = function(file) {
        png(file)
        plotindividuals_plotpca()
        dev.off()
      }
    )
    
  })
  
  # ann
  observeEvent(input$train_model, {
    req(input$annfile)
    req(input$p)
    req(input$hidden_layers)
    
    data <- read.csv(input$annfile$datapath)
    
    last_column_name <- names(data)[ncol(data)]
    
    indexes <- createDataPartition(data[[last_column_name]], p = input$p / 100, list = FALSE)
    train <- data[indexes, ]
    test <- data[-indexes, ]
    
    xtest <- test[, -ncol(test)]
    ytest <- test[[last_column_name]]
    
    hidden_layers <- as.numeric(unlist(strsplit(input$hidden_layers, ",")))
    
    nnet <- neuralnet(formula = as.formula(paste(last_column_name, "~ .")), 
                      data = train, hidden = hidden_layers, linear.output = FALSE)
    
    output$neuralnet_plot <- renderPlot({
      req(nnet)
      plot(nnet, rep = "best")
      #plot(nnet)
    })
    
    output$downloadneuralnet_plot <- downloadHandler(
      filename = "ann_plot.png",
      content = function(file) {
        png(file)
        req(nnet)
        plot(nnet, rep = "best")
        dev.off()
      }
    )
    
    
    ypred <- neuralnet::compute(nnet, xtest)
    
    yhat <- data.frame("yhat" = ifelse(ypred$net.result[, 1] > ypred$net.result[, 2], 
                                       levels(factor(ytest))[1], levels(factor(ytest))[2]))
    
    yhat$yhat <- factor(yhat$yhat, levels = levels(factor(ytest)))
    
    output$confusion_matrix22 <- renderPrint({
      confusionMatrix(factor(ytest), yhat$yhat)
    })
    
    output$downloadconfusion_matrix22 <- downloadHandler(
      filename = "confusion_matrix.txt",
      content = function(file) {
        writeLines(capture.output(confusionMatrix(factor(ytest), yhat$yhat)), file)
      }
    )
    
    
  })
  
  #regression
  observeEvent(input$regButton, {
    data <- reactive({
      req(input$regfile)
      read.csv(input$regfile$datapath)
    })
    
    # Perform LASSO and Ridge regression
    
    plotlasso_plot=function(){
      data_matrix <- as.matrix(data()[, c("logFC", "AveExpr", "t", "P.Value")])
      y <- as.numeric(data()$adj.P.Val)
      lasso_model <- cv.glmnet(data_matrix, y, alpha = 1)
      plot(lasso_model, xvar = "lambda", label = TRUE, col = rainbow(ncol(data_matrix)))
      #legend("topright", legend = colnames(data_matrix), col = rainbow(ncol(data_matrix)), lty = 1)
    }
    
    output$lasso_plot <- renderPlot({
      plotlasso_plot()
    }) 
    
    output$downloadlasso_plot <- downloadHandler(
      filename = "lasso_plot.png",
      content = function(file) {
        png(file)
        plotlasso_plot()
        dev.off()
      }
    )
    
    plotridge_plot = function(){
      data_matrix <- as.matrix(data()[, c("logFC", "AveExpr", "t", "P.Value")])
      y <- as.numeric(data()$adj.P.Val)
      ridge_model <- cv.glmnet(data_matrix, y, alpha = 0)
      plot(ridge_model, xvar = "lambda", label = TRUE, col = rainbow(ncol(data_matrix)))
      #legend("topright", legend = colnames(data_matrix), col = rainbow(ncol(data_matrix)), lty = 1)
    }
    
    output$ridge_plot <- renderPlot({
      plotridge_plot()
    })
    
    output$downloadridge_plot <- downloadHandler(
      filename = "ridge_plot.png",
      content = function(file) {
        png(file)
        plotridge_plot()
        dev.off()
      }
    )
    
  })
  
  #naive
  observeEvent(input$run_naive, {
    req(input$nvfile)
    req(input$pp)
    req(input$features)
    data_nv <- read.csv(input$nvfile$datapath)
    
    # Convert the last column to factor
    data_nv[, ncol(data_nv)] <- as.factor(data_nv[, ncol(data_nv)])
    
    # Split data into train and test sets
    indexes <- createDataPartition(data_nv[[ncol(data_nv)]], p = input$pp / 100, list = FALSE)
    train <- data_nv[indexes, ]
    test <- data_nv[-indexes, ]
    
    # Train Naive Bayes model
    nb_model <- naiveBayes(train[, -ncol(train)], train[, ncol(train)])
    
    # Make predictions on the test set
    nb_pred <- predict(nb_model, test[, -ncol(test)])
    
    # Compute accuracy
    accuracy <- mean(nb_pred == test[, ncol(test)])
    output$accuracyy3 <- renderText(paste("Accuracy:", accuracy)) 
    
    
    plotfeature_plot = function(){
      feature_index <- input$features
      class_colors <- c("red", "blue", "green", "purple", "orange") 
      par(mfrow=c(1,2))
      for (j in 1:length(levels(train[, ncol(train)]))) {
        class_data <- train[train[, ncol(train)] == levels(train[, ncol(train)])[j], feature_index]
        plot(density(class_data), col = class_colors[j], 
             main = paste("Feature", feature_index, "Distribution for", levels(train[, ncol(train)])[j]))
      }
      legend("topright", legend = levels(train[, ncol(train)]), fill = class_colors)
    }
    
    # Display plots for Feature 1
    output$feature_plot <- renderPlot({
      plotfeature_plot()
    })
    
    output$downloadfeature_plot <- downloadHandler(
      filename = "feature_plot.png",
      content = function(file) {
        png(file)
        plotfeature_plot()
        dev.off()
      }
    )
    # Confusion matrix
    conf_matrix <- confusionMatrix(nb_pred, test[, ncol(test)])
    output$confusion_matrixx3 <- renderPrint(conf_matrix)
    
    output$downloadconfusion_matrixx3 <- downloadHandler(
      filename = "confusion_matrix.txt",
      content = function(file) {
        writeLines(capture.output(conf_matrix), file)
      }
    )
    
    # ROC curve
    output$roc_plot <- renderPlot({
      roc_curve <- roc(test[, ncol(test)], as.numeric(nb_pred))
      plot(roc_curve, col = "blue")
    }) 
    
    output$downloadroc_plot <- downloadHandler(
      filename = "roc_plot.png",
      content = function(file) {
        png(file)
        roc_curve <- roc(test[, ncol(test)], as.numeric(nb_pred))
        plot(roc_curve, col = "blue")
        dev.off()
      }
    )
    
    output$ggp_plot <- renderPlot({
      ggpairs(data_nv)
    })
    
    output$downloadggp_plot <- downloadHandler(
      filename = "ggp_plot.png",
      content = function(file) {
        png(file)
        ggpairs(data_nv)
        dev.off()
      }
    )
    
    
    output$miss_plot <- renderPlot({
      missmap(data_nv)
    })
    
    output$downloadmiss_plot <- downloadHandler(
      filename = "miss_plot.png",
      content = function(file) {
        png(file)
        missmap(data_nv)
        dev.off()
      }
    )
    
  })
  
  
}
