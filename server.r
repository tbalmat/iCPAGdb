# Duke University Cross Phenotype Analysis of GWAS Database (iCPAGdb)
# Explore page, November 2020

# Shiny app user interface function

# Information on shiny available at:
# https://shiny.rstudio.com/
# https://github.com/rstudio/shiny
# https://cran.r-project.org/web/packages/shiny/shiny.pdf

options(max.print=1000)      # number of elements, not rows
options(stringsAsFactors=F)
options(scipen=999999)
#options(device="windows")
options(shiny.maxRequestSize=2*10**8) 

library(shiny)
library(DT)
library(RSQLite)
library(ggplot2)
#library(reshape)
library(heatmaply)
library(plotly)
library(future)
library(promises)

# Configure asynchronous execution for future() operations
plan(multisession)

shinyServer(

  function(input, output, session) {

    # Change directories and set some globals
    setwd(c("C:/Projects/Duke/H2P2GenomeWideAssociationStudy/CPAG/iCPAGdb/App/pyCPAG",
            "/srv/shiny-server/CPAG/explore/pyCPAG")[2])
    pyexe <- c("\"C:/Users/Kyung Soon/AppData/Local/Programs/Python/Python37/python.exe\"", "python3")[2]
    threads <- 2
    gwidthMin <- 600
    gheightMin <- 600
    gwidthMax <- 2000
    gheightMax <- 2000
    currentReviewCPAGdata <- data.frame()
    reviewResultsEFOcol <- vector("integer")
    reviewResultsEFOparentCol <- vector("integer")
    reviewResultsEFOparentTerm <- vector("character")
    reviewResultsTableRowFilter <- NULL
    currentExploreCPAGdata <- data.frame()
    exploreResultsEFOcol <- vector("integer")
    exploreResultsEFOparentCol <- vector("integer")
    exploreResultsEFOparentTerm <- vector("character")
    exploreResultsTableRowFilter <- NULL
    largestKnownLogVal <- 347
    inDir <- "input"
    outDir <- "output"
    inDirUserCompute <- "input/userCompute"
    outDirUserCompute <- "output/userCompute"
    exploreTabFirst <- T
    currentUserComputeCPAGdata <- data.frame()
    userComputeResultsEFOcol <- vector("integer")
    userComputeResultsEFOparentCol <- vector("integer")
    userComputeResultsEFOparentTerm <- vector("character")
    userComputeResultsTableRowFilter <- NULL
    notifyDuration <- 5
    nLabelTrunc <- 45
    phenH2P2 <- "limited"
    hideTab("tabsetCPAG", "tabPanelExplore")
    #showTab("tabsetCPAG", "tabPanelExplore")
    

    #######################################################################################################
    # Configure review selection table
    #######################################################################################################

    reviewCPAG <- data.frame("Description"=c("COVID-19",
                                             "Human disease",
                                             "Host-pathogen traits",
                                             "Host-pathogen traits",
                                             "Host-pathogen traits",
                                             "Serum metabolites and xenobiotics",
                                             "Urine metabolites",
                                             "Host-pathogen vs. Human disease",
                                             "Host-pathogen vs. Human disease",
                                             "Host-pathogen vs. Human disease",
                                             "Serum metabolites/xenobiotics vs. Human disease",
                                             "Urine metabolites vs. Human disease",
                                             "All molecular/cellular traits vs. Human disease"),
                             "GWAS_1_source"=c("Ellinghaus",
                                               "NHGRI-EBI GWAS catalog",
                                               "H2P2",
                                               "H2P2",
                                               "H2P2",
                                               "",
                                               "",
                                               "H2P2",
                                               "H2P2",
                                               "H2P2",
                                               "",
                                               "",
                                               ""),
                             "GWAS_2_source"=c("NHGRI, H2P2, Metab, Mol, Clin",
                                               "",
                                               "",
                                               "",
                                               "",
                                               "",
                                               "",
                                               "NHGRI-EBI",
                                               "NHGRI-EBI",
                                               "NHGRI-EBI",
                                               "NHGRI-EBI",
                                               "NHGRI-EBI",
                                               "NHGRI-EBI"),
                             "p_threshold_GWAS_1"=c(HTML("1X10<sup>-5</sup>"),
                                                    HTML("5X10<sup>-8</sup>"),
                                                    HTML("1X10<sup>-5</sup>"),
                                                    HTML("1X10<sup>-5</sup>"),
                                                    HTML("1X10<sup>-5</sup>"),
                                                    HTML("1X10<sup>-5</sup>"),
                                                    HTML("1X10<sup>-5</sup>"),
                                                    HTML("1X10<sup>-5</sup>"),
                                                    HTML("1X10<sup>-5</sup>"),
                                                    HTML("1X10<sup>-5</sup>"),
                                                    HTML("1X10<sup>-5</sup>"),
                                                    HTML("1X10<sup>-5</sup>"),
                                                    HTML("1X10<sup>-5</sup>")),
                             "p_threshold_GWAS_2"=c(HTML("5X10<sup>-8</sup> (NHGRI), 1X10<sup>-5</sup>"),
                                                    "",
                                                    "",
                                                    "",
                                                    "",
                                                    "",
                                                    "",
                                                    HTML("5X10<sup>-8</sup>"),
                                                    HTML("5X10<sup>-8</sup>"),
                                                    HTML("5X10<sup>-8</sup>"),
                                                    HTML("5X10<sup>-8</sup>"),
                                                    HTML("5X10<sup>-8</sup>"),
                                                    HTML("5X10<sup>-8</sup>")),
                             "LD_population"=c("EUR",
                                               "EUR",
                                               "EUR",
                                               "AFR",
                                               "EAS",
                                               "EUR",
                                               "EUR",
                                               "EUR",
                                               "AFR",
                                               "EAS",
                                               "EUR",
                                               "EUR",
                                               "EUR"),
                             "Reference"=c("<a href=https://pubmed.ncbi.nlm.nih.gov/32558485/ target=_blank>Ellinghaus et al., 2020</a>",
                                           "<a href=https://academic.oup.com/nar/article/47/D1/D1005/5184712 target=_blank>Buniello et al., 2019</a>",
                                           "<a href=https://www.cell.com/cell-host-microbe/fulltext/S1931-3128(18)30377-9 target=_blank>Hi-HOST Phenome Project; Wang et al. 2018</a>",
                                           "<a href=https://www.cell.com/cell-host-microbe/fulltext/S1931-3128(18)30377-9 target=_blank>Hi-HOST Phenome Project; Wang et al. 2018</a>",
                                           "<a href=https://www.cell.com/cell-host-microbe/fulltext/S1931-3128(18)30377-9 target=_blank>Hi-HOST Phenome Project; Wang et al. 2018</a>",
                                           "<a href=https://www.nature.com/articles/ng.2982 target=_blank>Shin et al. 2014</a>",
                                           "<a href=https://pubmed.ncbi.nlm.nih.gov/26352407/ target=_blank>Raffler et al. 2015</a>",
                                           "",
                                           "",
                                           "",
                                           "",
                                           "",
                                           ""),
                             "file"=c("top_EllinghausPCs_covid19_pcut1e-5.cpag2_out_20201215_addOntology.csv",
                                      "NHGRI-p5e-08-EUR.csv",
                                      "H2P2-p1e-05-EUR.csv",
                                      "H2P2-p1e-05-AFR.csv",
                                      "H2P2-p1e-05-EAS.csv",
                                      "BloodMetabolites-And-BloodXenobiotic-p1e05-EUR.csv",
                                      "UrineMetabolites-p1e-05-EUR.csv",
                                      "H2P2-p1e-05-NHGRI-p5e-08-EUR.csv",
                                      "H2P2-p1e-05-NHGRI-p5e-08-AFR.csv",
                                      "H2P2-p1e-05-NHGRI-p5e-08-EAS.csv",
                                      "NHGRI-p1e-05-BloodMetabolitesXenobiotic-p1e-05-EUR.csv",
                                      "NHGRI-p1e-05-UrineMetabolites-p1e-05-EUR.csv",
                                      "NHGRI-p1e-05-mol_gwas-p1e-05-EUR.csv")
                            )

    # Render review selection table
    output$reviewSelectionTable <- DT::renderDataTable(
                                     datatable(
                                       reviewCPAG[,c("Description", "GWAS_1_source", "GWAS_2_source", "p_threshold_GWAS_1",
                                                     "p_threshold_GWAS_2", "LD_population", "Reference")],
                                       # Suppress row names
                                       # Do not escape HTML tags
                                       # Automatically hide nav buttons when rows less than pages per row
                                       rownames=F, escape=F, autoHideNavigation=T,
                                       colnames=c("Description",
                                                  HTML("Source (GWAS<sub>1</sub>)"),
                                                  HTML("Source (GWAS<sub>2</sub>)"),
                                                  HTML("P<sub>threshold</sub> (GWAS<sub>1</sub>)"),
                                                  HTML("P<sub>threshold</sub> (GWAS<sub>2</sub>)"),
                                                  "LD Population", "Reference"),
                                       # Eliminate default padding around cells, which causes variable horizontal alignment between rows
                                       # Classes are defined by DT (compact omits padding, cell-border outlines cells, row-border separates
                                       # rows with lines, stripe alternates row background colors, nowrap prohibits wrapping within columns)
                                       # More info at https://datatables.net/manual/styling/classes
                                       class="compact stripe",
                                       # Use cell selection (default is row) to avoid rendering plot when url selected in cell
                                       # Limit selection to single cell (multiple is an option, which causes the _cell_clicked
                                       # list to be appended with another entry each time a cell is clicked)
                                       selection=list(mode="single", target="cell"),
                                       options=list(
                                         # More info at https://datatables.net/reference/option/
                                         # Display table only, no page controls, search fields, etc.
                                         dom="<'row'<'col-sm-12't>>",
                                         # Fixed page length without pagination limits observable rows to the number provided here
                                         pageLength=50,
                                         scrollY="145px",
                                         # Hide column ordering controls
                                         ordering=F,
                                         # Format columns
                                         # Note that col indices are 0-based
                                         # Column widths are proportional to percent of sum of supplied percents provided
                                         # The rendered table fills space allocated in ui.r table def
                                         columnDefs=list(#list(className="dt-right", targets=),
                                                         list(className="dt-center", targets=3:5),
                                                         list(width='15%', targets=0),
                                                         list(width='10%', targets=c(1, 2)),
                                                         list(width='8%', targets=c(3, 4, 5)),
                                                         list(width='20%', targets=6))
                                         # Modify stripes
                                         # Not certain how this functions, but specification overrides or modifies stripe class
                                         #stripeClasses=c("background-color:red", "font-color:red")
                                       )
                                     )
                                   )

    #######################################################################################################
    # Function:  Compose CPAG data table
    # Information on data tables and options available at:
    # https://rstudio.github.io/DT/
    # https://cran.r-project.org/web/packages/DT/DT.pdf
    # https://datatables.net/reference/option/
    #######################################################################################################

    composeCPAGTable <- function(data, truncSNP=T) {

      # Transform values for display
      cdat <- data.frame("Trait1"=data[,"Trait1"],
                         "Trait2"=data[,"Trait2"],
                         #"N_expected"=round(data[,"N_expected"], 4),
                         "Nshare_direct"= data[,"Nshare_direct"],
                         "Nshare_LD"=data[,"Nshare_LDpairComb"],
                         "Nshare_all"=data[,"Nshare_all"],
                         "P_fisher"=-round(log(data[,"P_fisher"])/log(10), 4),
                         "Padj_Bonferroni"=-round(log(data[,"Padj_Bonferroni"])/log(10), 4),
                         "Padj_FDR"=-round(log(data[,"Padj_FDR"])/log(10), 4),
                         "SNPshare_all"=ifelse(rep(truncSNP, nrow(data)),
                                               data[,"SNPshare_all"],
                                               substring(data[,"SNPshare_all"], 1, 50)),
                         "Jaccard"=round(data[,"Jaccard"], 4),
                         "ChaoSorensen"=round(data[,"ChaoSorensen"], 4))
      cdat[which(is.infinite(cdat[,"P_fisher"])),"P_fisher"] <- paste(">", largestKnownLogVal, sep="")
      cdat[which(is.na(cdat[,"P_fisher"])),"P_fisher"] <- "na"
      cdat[which(is.infinite(cdat[,"Padj_Bonferroni"])),"Padj_Bonferroni"] <- paste(">", largestKnownLogVal, sep="")
      cdat[which(is.na(cdat[,"Padj_Bonferroni"])),"Padj_Bonferroni"] <- "na"
      cdat[which(is.na(cdat[,"Padj_Bonferroni"])),"Padj_Bonferroni"] <- "na"
      cdat[which(is.infinite(cdat[,"Padj_FDR"])),"Padj_FDR"] <- paste(">", largestKnownLogVal, sep="")
      cdat[which(is.na(cdat[,"Padj_FDR"])),"Padj_FDR"] <- "na"

      # Configure column headings
      colnames=c(HTML("Trait<sub>1</sub>"), HTML("Trait<sub>2</sub>"), HTML("Nshare<sub>direct</sub>"),
                 HTML("Nshare<sub>LD</sub>"), HTML("Nshare<sub>all</sub>"), HTML("-log<sub>10</sub>(P<sub>Fisher</sub>)"),
                 HTML("-log<sub>10</sub>(P<sub>Bonferroni</sub>)"), HTML("-log<sub>10</sub>(P<sub>FDR</sub>)"),
                 HTML("SNP<sub>shared</sub>"), "Jaccard", "ChaoSorensen")

      # Configure column definitions
      # Note that col indices are 0-based
      columnDefs=list(list(width="20%", targets=c(0, 1)),
                      list(width="5%", targets=c(2, 3, 4, 5, 6, 7, 9, 10)),
                      list(width="5%", targets=8),
                      list(className="dt-right", targets=c(5, 6, 7, 9, 10)))

      # Include hyperlinks for EFO columns, if present
      # Include parent terms
      #efoCol <- which(regexpr("Trait.\\_EFO", colnames(data))>0)
      #efoParentCol <- vector("integer")
      #for(j in efoCol) {
      #  cdat <- cbind(cdat, unlist(lapply(data[,j],
      #                               function(a) {
      #                                 url <- trimws(strsplit(a, ",")[[1]])
      #                                 paste(paste("<a href=", url, " target=_blank>", url, "</a>", sep=""), collapse=" ")
      #                               })))
      #  cdat <- cbind(cdat, data[,j])
      #  colnames <- c(colnames, colnames(data)[j])
      #  # Append corresponding parent column
      #  k <- which(colnames(data)==sub("EFO", "ParentTerm", colnames(data)[j]))
      #  if(length(k)>0) {
      #    efoParentCol <- c(efoParentCol, k[1])
      #    cdat <- cbind(cdat, data[,k[1]])
      #    colnames <- c(colnames, colnames(data)[k[1]])
      #  }
      #}

      # Include EFO parent columns only, if present
      efoParentCol <- which(regexpr("Trait.\\_ParentTerm", colnames(data))>0)
      for(j in efoParentCol) {
        cdat <- cbind(cdat, data[,j])
        colnames <- c(colnames, sub("ParentTerm", "EFO", colnames(data)[j]))
      }

      # Include column width definitions for EFO and parent columns
      # Parent columns only
      #if(length(efoCol)>0)
      #  columnDefs[[length(columnDefs)+1]] <- list(width="5%", targets=10+(1:(length(efoCol)+length(efoParentCol))))
      #if(length(efoParentCol)>0)
      #  columnDefs[[length(columnDefs)+1]] <- list(width="5%", targets=10+(1:(length(efoParentCol))))

      # Include column width definitions for EFO parent columns only
      if(length(efoParentCol)>0)
        columnDefs[[length(columnDefs)+1]] <- list(width="5%", targets=10+(1:(length(efoParentCol))))

      # Compose data table
      datatable(

        cdat,

        # Apply column names
        colnames=colnames,

        # Suppress row names, do not escape HTML tags
        # Automatically hide nav buttons when rows less than pages per row
        rownames=F, escape=F, autoHideNavigation=T,

        # Width does not appear to have any effect (at least when columns have lengthy data)
        # width="1000px",

        # Use cell selection (default is row) to avoid rendering plot when url selected in cell
        # Limit selection to single cell (multiple is an option, which causes the _cell_clicked
        # list to be appended with another entry each time a cell is clicked)
        selection=list(mode="single", target="cell"),

        # Cells can be edited, which modifies the source data frame (what you do with the modified DF is for you to decide)
        #editable="cell",

        # Table appearance can be redefined in CSS options
        class="cell-border stripe", 

        # Include column filters ("bottom" will display them below table)
        filter=c("none", "bottom", "top")[1],

        # Include button extension (for download)
        #extensions="Buttons",

        # Configure other table options
        # Information on data tables options available at https://rstudio.github.io/DT/options.html
        options=list( # Specify what is to appear and in what order
                      # Additional information available at https://datatables.net/reference/option/dom
                      # Info on position objects available at https://stackoverflow.com/questions/49035864/positioning-datatables-elements-with-dom-option
                      # Symbols are: (t)able,
                      #              p(r)ocessing display element (not sure what this is),
                      #              (i)nfo summary (showing rec i of n)
                      #              (p)age control
                      #              (B)uttons (download, print)
                      #              (l)ength of page (number of records) control
                      #              (f)ind
                      # Simple method
                      # dom="tripBlf",
                      # Method to specify row and col layout
                      dom=paste(#"<'row'<'col-sm-3'l><'col-sm-5'f><'col-sm-3'><'col-sm-1'B>>",
                                "<'row'<'col-sm-3'l><'col-sm-8'><'col-sm-1'B>>",
                                "<'row'<'col-sm-12't>>",
                                "<'row'<'col-sm-5'i><'col-sm-7'p>>", sep=""),

                      # Alternative to enabling display of page length (rows per page) and global filter controls
                      #bLengthChange=T, bFilter=T,

                      # Set rows per page
                      pageLength=25, lengthMenu=c(10, 25, 50, 100, 250),

                      # Allow regex style searches
                      #search=list(regex=T, caseInsensitive=T),

                      # Disable auto-width, where col widths a function of max entry
                      autoWidth=F,

                      # Add scroll bar to table, otherwise forced width beyond window causes browser x-scroll-bar to appear
                      scrollX=T,
                      #fixedColumns=T,

                      # Format columns
                      # Note that col indices are 0-based
                      columnDefs=columnDefs

                      # Sort order of columns (ascending for all independent var cols)
                      # Note that col indices are 0-based
                      # Simple method for two fixed columns
                      # order=list(list(0, "asc"), list(1, "asc")),
                      # Method using number of independent vars
                      #order=lapply(1:length(indepVar), function(i) list(i-1, "asc")),

                      # Download button
                      # Although convenient, this download method requires client side table processing (server=F
                      # in renderDataTable()) which significantly increases overhead for large tables (all data must
                      # transferred to the client prior to the first page being viewed)
                      # An alternative to the "buttons" method is to use a download handler (see the
                      #buttons=list(list(extend="collection", buttons=c('csv', 'excel', 'pdf'), text='Download'))
                    )

      )

    }

    #######################################################################################################
    # Function:  Compose observation filter index for review table
    #######################################################################################################

    resultsTableComposeFilterIndex <- function(data, filterTrait, filterSNP, includeCompoundEFO,
                                               efoCol, filterEFOparent, efoParentCol) {

      # Identify observations that satisfy trait filter criteria
      # Omit leading and trailing spaces from filter
      filterTrait <- tolower(trimws(filterTrait))
      if(nchar(filterTrait)>0) {
        k1 <- unique(c(which(regexpr(filterTrait, tolower(data[,"Trait1"]), fixed=T)>0),
                       which(regexpr(filterTrait, tolower(data[,"Trait2"]), fixed=T)>0)))
      } else {
        k1 <- 1:nrow(data)
      }

      # Identify observations that satisfy SNP filter criteria
      filterSNP <- tolower(trimws(filterSNP))
      if(nchar(filterSNP)>0) {
        k2 <- which(regexpr(filterSNP, tolower(data[,"SNPshare_all"]), fixed=T)>0)
      } else {
        k2 <- 1:nrow(data)
      }

      # Display compound EFOs, if requested
      if(includeCompoundEFO & length(efoCol)>0) {
        k3 <- 1:nrow(data)
      } else {
        # Identify oservations with no comma in any phenotype column
        k3 <- setdiff(1:nrow(data),
                      unlist(lapply(efoCol,
                             function(v) which(regexpr(",", data[,v], fixed=T)>0))))
      }

      # Identify observations that satisfy EFO parent filter criteria (using trait 1 or trait 2 EFO parent columns)
      if(length(filterEFOparent)>0  & length(efoParentCol)>0) {
        # Construct regular expression to search multiple terms
        # (?a) executes a case insensitive search
        filterEFO <- paste("(?i)(", paste(filterEFOparent, collapse="|", sep=""), ")", sep="")
        # Index observations with EFO parent column(s) containing terms in filter vector
        k4 <- unique(unlist(lapply(efoParentCol, function(j) which(regexpr(filterEFO, data[,j])>0))))
      } else {
        k4 <- 1:nrow(data)
      }

      # Intersect indices to identify observations satisfying all filters
      k <- intersect(intersect(intersect(k1, k2), k3), k4)
      cat(paste("filter:  ", "n-tr=", length(k1), ", n-SNP=", length(k2), ", n-EFO=", length(k3),
                ", n-EFO-parent=", length(k4), ", n-intersect=", length(k), "\n", sep=""), file=stderr())
      return(k)

    }

    #######################################################################################################
    # Function:  Filter results table, render data table and heatmap
    #######################################################################################################

    renderResults <- function(data, rowFilter, filterTrait, filterSNP, truncSNP, includeCompoundEFO,
                              efoCol, filterEFOparent, efoParentCol, tabset, tableTab, tableName,
                              plotTab, plotName, heatmapMetric, heatmapNphenotype) {

      # Construct filtered observation indices
      resultsTableRowFilter <- resultsTableComposeFilterIndex(data, filterTrait, filterSNP, includeCompoundEFO,
                                                              efoCol, filterEFOparent, efoParentCol)

      # Save filter indices in the environment in which the current function was declared (where variable named in rowFilter
      # was also declared 
      assign(rowFilter, resultsTableRowFilter, envir=environment(renderResults))

      # Display CPAG table
      if(length(resultsTableRowFilter)>0) {
          progress <- shiny::Progress$new()
          progress$set(message="Composing CPAG table", value=0.5)
          output[[tableName]] <- DT::renderDataTable(composeCPAGTable(data[resultsTableRowFilter,], truncSNP),
                                                     # server=F instructs to transfer all rows when rendering, instead of
                                                     # row transfer being synchronized with user page requests
                                                     # This results in considerable additional overhead with large tables,
                                                     # but makes all rows available for download using the data tables
                                                     # download button (otherwise, downloaded rows are limited to those
                                                     # appearing on the page being viewed at time of the download request
                                                     # An alternative to enabling the data tables download button and
                                                     # client-side handling is to implement a 
                                                     server=T)
          progress$close()
          # Compose and render heatmap
          renderHeatmap(data=data[resultsTableRowFilter,],
                        xCol="Trait1", yCol="Trait2", valueCol=heatmapMetric,
                        log10Value=(heatmapMetric %in% c("P_fisher", "Padj_Bonferroni", "Padj_FDR")),
                        nPhenotype=heatmapNphenotype,
                        key.title=ifelse(heatmapMetric %in% c("P_fisher", "Padj_Bonferroni", "Padj_FDR"),
                                         paste("-log10(", heatmapMetric, ")", sep=""),
                                         heatmapMetric),
                        tabset, tableTab, plotTab, plotName)
      } else {
        output[[tableName]] <- DT::renderDataTable(NULL)
        output[[plotName]] <- renderPlotly(NULL)
        #showNotification("No CPAG data exist for specified filter values", type="error", duration=notifyDuration)
        showModal(modalDialog(HTML("No CPAG data exist for specified filter values"),
                              title="iCPAGdb", size="m", easyClose=T, footer=modalButton("OK")))
      }

    }

    ##################################################################################################
    # String trunction function
    # Truncate each element of a to ntrunc characters
    # If multiple distinct values of a have common (ntrunc) leading characters then append -1 to all
    # that contain the first distinct value, -2 to all that contain the second valure, ...
    ##################################################################################################

    truncLabel <- function(a, ntrunc) {
      b <- substring(a, 1, ntrunc)
      for(c in unique(b)) {
        d <- unique(a[which(b==c)])
        if(length(d)>1)
          for(i in 1:length(d)) {
            k <- which(a==d[i])
            b[k] <- paste(b[k], "-", i, sep="")
          }
      }
      b
    }

    #######################################################################################################
    # Alternative cast() function (for when the shape library is not available)
    # Form a matrix with x and y axes, intersecting cells containing z from corresponding x,y rows
    # Note that a single z from those with common x,y is returned
    #######################################################################################################

    xyzMatrix <- function(x, y, z) {
      x0 <- sort(unique(x))
      y0 <- sort(unique(y))
      # Iterate through each unique y
      # Coerce to a matrix in case x or y have a single value
      m <- as.matrix(apply(as.matrix(y0), 1,
                       function(y2) {
                         # Identify rows corresponding to current unique y
                         i <- which(y==y2)
                         # Identify first x in current y rows corresponding to each unique x
                         # Return corresponding z (NA if none)
                         z[i[match(x0, x[i])]]
                       }))
      # Transpose if all of x or y have a single value
      if(any(c(length(x0), length(y0))!=dim(m)))
        m <- t(m)
      rownames(m) <- x0
      colnames(m) <- y0
      m
    }

    #######################################################################################################
    # Function:  Compute the matrix of Euclidian distances between all pairs of rows in a matrix
    # (typically, the Fisher-p matrix in this app)
    # The result is an object of class "dist" similar to that produced by the dist() function
    # This is the lower triangle (less the diagonal) of the matrix of Euclidian row distances, where
    # position i,j corresponds to the difference between rows i and j
    # This function is a substitute for the standard dist() function due to its fixed method of handing NAs
    # NA Fisher-p values indicate that no relationship exists between two phenotypes
    # NA values in the p matrix causes the error:
    # Error in hclustfun(dist) : NA/NaN/Inf in foreign function call (arg 10)
    # during a call to hclust(dist()) by heatmaply()
    # The default dist() method computes Euclidian distances between all rows of a
    # matrix (using the distances between corresponding columnar values)
    # NA elements are omitted, with the result scaled to compensate for omitted values
    # From the dist() doc:  If some columns are excluded in calculating a Euclidean, Manhattan,
    #                       Canberra or Minkowski distance, the sum is scaled up proportionally
    #                       to the number of columns used. If all pairs are excluded when calculating
    #                       a particular distance, the value is NA.
    # It is possible (likely) that some pairs of phenotypes do not have a relationship (in a significant
    # Fisher-p sense), as in the case where AB and AC are significant, but BC is not
    # This results in NA values being returned by dist()
    # To avoid the scaling performed by dist() (if two phenotypes are unrelated, we do not want any
    # adjustment to the distances between any other pairs of phenotypes, while there os no parameter
    # to modify the scaling function of dist()) and to eliminate NAs from distances,
    # a function is defined (naDist) to compute Euclidian distances between rows of the Fisher-p matrix
    # using columns where each are non-NA
    # If all columns contain NA for a particular pair of rows, the the resulting distance is set
    # to the maximum distance observed in all non-NA entries

    naDist <- function(x) {
      if(length(dim(x))==2) {
        # Compute distance between all pairs of rows
        # Return as matrix where row i, col j indicates the difference between rows i and j
        y <- apply(x, 1, function(r1)
                           apply(x, 1, function(r2) {
                                         k <- which(!is.na(r1) & !is.na(r2))
                                         if(length(k)>0) {
                                           sqrt(sum((r1[k]-r2[k])**2))
                                         } else {
                                           NA
                                         }
                                       }))
        y <- y[lower.tri(y, diag=F)]
        # Replace NA values with the maximum observed distance
        k <- which(is.na(y))
        if(length(k)<length(y)) {
          y[k] <- max(y, na.rm=T)
        } else {
          y <- rep(0, length(y))
        }
        attributes(y) <- list("Size"=dim(x)[1], "Diag"=F, "Upper"=F, "method"="euclidian", "Labels"=rownames(x))
        # Establish class after all assignments, in case it was altered
        class(y) <- "dist"
      } else {
        y <- dist(0)
      }
      return(y)
    }

    #######################################################################################################
    # Function:  Compose heatmap
    #######################################################################################################

    composeHeatmap <- function(data, x, y, z, log10z=F, nPhenotype="50", nDir="increasing", key.title="") {

      # Subset to top phenotype pairs, if specified
      if(nPhenotype!="all")
        data <- data[order(data[,z], decreasing=(nDir=="decreasing"))[1:min(as.numeric(nPhenotype), nrow(data))],]

      # heatmaply includes support for a legend and many adjustable visual features
      if(log10z) {
        data[,z] <- -log(data[,z])/log(10)
        # Convert infinite values (attempted to compute log of something near 0) to largest known possible value
        # This avoids numerical problems in the distance function (used for dendrogram associations)
        k <- which(is.infinite(data[,z]))
        data[k,z] <- largestKnownLogVal
      }

      # If all z elements are infinite, convert them to largest known possible value
      #if(all(is.infinite(data[,z])))
      #  data[,z] <- largestKnownLogVal

      # Generate a matrix from specified columns
      # cast() requires the reshape package, which was not available on the server at time of development
      #A <- as.matrix(cast(data, paste(x, "~", y, sep=""), value=z))
      # Truncate x and y labels to nLabelTrunc characters to avoid the observed situation where, with lengthy
      # axis labels, heatmaply does not render all rows and columns (yes, a partial plot can be genereated,
      # with no warning that anything was omitted - ouch!) 
      A <- xyzMatrix(truncLabel(data[,x], nLabelTrunc), truncLabel(data[,y], nLabelTrunc), data[,z])

      # Disable dendrograms when one axis has a single element
      if(all(dim(A)>1)) {
        dendroStyle <- "both"
      } else {
        dendroStyle <- "none"
      }

      # Generate heatmap
      g <- heatmaply(
             A,
             colors=viridis(n=256, alpha=1, begin=0, end=1, option="viridis"),
             # na.value does not appear to have any effect
             # Achieve similar effect with panel border fill in theme
             #na.value="gray50",
             plot_method=c("plotly", "ggplot")[2],
             # Normalize distances used to form dendrograms
             dendrogram=dendroStyle,
             # Parameter to distfun is the input matrix to heatmaply()
             distfun=function(x) {z <- naDist(x); if(max(z)!=0) {z/max(z)} else {z}},
             # Dendrogram line weight
             branches_lwd=0.1,
             # Margins appear to be around the entire plot, including dendrograms
             #margins=c(30, 30, 30, 30),
             # Draw grid to separate cells
             #grid_gap=1,
             # grid_color does not appear to have any effect when plot_method="plotly"
             grid_color="gray95",
             # grid_size does not appear to have any effect
             #grid_size=0.1,
             # Axis label fonts are specified in theme
             #fontsize_row=8,
             #fontsize_col=8,
             key.title=key.title,
             heatmap_layers=list(
                              ggplot2::theme(
                                # Specify a transparent panel background since it is plotted last (yee-gads!)
                                panel.border=element_rect(fill=rgb(0.75, 0.75, 0.75, alpha=0.2), color="black"),
                                # Major grid lines bisect heatmap rectangles, minor grid is not rendered
                                #panel.grid.minor=element_line(color="gray50"),
                                axis.title=element_text(size=10),
                                axis.text.x=element_text(size=8, angle=45, hjust=1, vjust=0.5),
                                axis.text.y=element_text(size=8),
                                legend.text=element_text(size=8),
                                legend.title=element_text(size=10)))
           )

      #print(str(g))

      # Compute plot width and height
      # Trait 1 in rows, trait 2 in columns
      # Allow sufficient pixels for each row and column in the Fisher p-matrix
      # Pad space for dendrogram arcs, labels, and legend
      # Note that x-axis labels are rotated approximately 45 deg
      # Adjust dimensions when an axis has a solitary value, typically from a user upload GWAS
      n1 <- length(unique(data[,"Trait1"]))
      n2 <- length(unique(data[,"Trait2"]))
      if(n2>1) {
        gwidth <- min(20*n2+50+7*max(nchar(data[,"Trait1"]))+80, gwidthMax)
      } else {
        gwidth <- gwidthMin
      }
      if(n1>1) {
        gheight <- min(20*n1+50+5*max(nchar(data[,"Trait2"])), gheightMax)
      } else {
        gheight <- gheightMin
      }
      gwidth <- paste(gwidth, "px", sep="")
      gheight <- paste(gheight, "px", sep="")
      #cat(paste("gwidth: ", gwidth, ", gheight: ", gheight, "\n", sep=""), file=stderr())

      return(list("g"=g, "gwidth"=gwidth, "gheight"=gheight))

    }

    #######################################################################################################
    # Function:  Render heatmap
    #######################################################################################################

    renderHeatmap <- function(data, xCol, yCol, valueCol, log10Value, nPhenotype, key.title,
                              tabset, tableTab, plotTab, plotName) {

      # Plot objects persist when tabs are re-created and a new plot is not rendered (due to error?)
      # Therefore, explicitly nullify the plot
      output[[plotName]] <- renderPlotly(NULL)

      # Compose heatmap from the matrix of values with x, y axes and specified corresponding values
      if(nrow(data)>1) {
        progress <- shiny::Progress$new()
        progress$set(message="Preparing heatmap", value=0.3)
        hmap <- composeHeatmap(data, xCol, yCol, valueCol, log10z=log10Value, nPhenotype=nPhenotype,
                               nDir=ifelse(valueCol %in% c("Jaccard", "ChaoSorensen"), "decreasing", "increasing"),
                               key.title=key.title)

        # Remove and reinsert tab with computed plot region width and height, which are not modifiable
        removeTab(session=session, inputId=tabset, target=plotTab)
        insertTab(session=session, inputId=tabset,
                  tab=tabPanel(title="Heatmap",
                               value=plotTab,
                               fluidRow(width=12,
                                 column(width=10,
                                   HTML("<br><br><br><center>"),
                                   #plotlyOutput(plotName, width=hmap[["gwidth"]], height=hmap[["gheight"]]),
                                   plotlyOutput(plotName, width=hmap[["gwidth"]], hmap[["gheight"]]),
                                   HTML("</center>")
                                 )
                               )
                             ),
                  target=tableTab, position="after", select=F)

        # Display heatmap
        progress$set(message="Rendering heatmap", value=0.6)
        output[[plotName]] <- renderPlotly(hmap[["g"]])
        progress$close()
      } else {
        #showNotification("Insufficient data (cross-GWAS phenotype pairs) to generate heatmap", type="error", duration=notifyDuration)
        showModal(modalDialog(HTML("Insufficient data (cross-GWAS phenotype pairs) to generate heatmap"),
                  "title"="iCPAGdb", size="m", easyClose=T, footer=modalButton("OK")))
      }

    }

    #######################################################################################################
    # Observe event:  Review tab - selection table cell click
    # Note that input$reviewSelectionTable_cells_selected is a variable maintained by the DT package, but
    # is made reactive by Shiny, so that changes to the DT variable cause events here
    #######################################################################################################

    observeEvent(input$reviewSelectionTable_cells_selected, {

        if(length(input$reviewSelectionTable_cell_clicked)>0) {

          # Ignore cell selection when cell contains an HTML anchor tag
          # Note that _cell_clicked returns a list with elements:
          # row = row in data frame corresponding to cell clicked (reordering rows maintains row IDs)
          # col = 0-based col of cell clicked
          # value = data frame cell contents
          # Note that if selection mode of multiple is specified in renderTable(), then multiple cells
          # will be reported by _cell_clicked, all at the top list level
          # Referencing list elements by name returns the first encountered, ignoring multiples
          if(substring(as.character(input$reviewSelectionTable_cell_clicked[["value"]]), 1, 7)!="<a href") {

            # Retrieve file name of specified precomputed results corresponding to selected row number
            fn <- reviewCPAG[input$reviewSelectionTable_cell_clicked[["row"]][1],"file"]
            cat(paste("fileName:  ", fn, "\n", sep=""), file=stderr())

            # Proceed if file exists
            if(length(which(dir(outDir, pattern=paste("^", fn, "$", sep=""))==fn))==1) {

              # Retrieve CPAG results
              # Save in global data frame to be available in other reactive functions
              progress <- shiny::Progress$new()
              progress$set(message="Reading upload data", value=0.5)
              currentReviewCPAGdata <<- read.table(paste(outDir, "/", fn, sep=""), header=T, sep=",", quote="\"")
              progress$close()

              # Index EFO columns
              reviewResultsEFOcol <<- which(regexpr("Trait.\\_EFO", colnames(currentReviewCPAGdata))>0)
              reviewResultsEFOparentCol <<- which(regexpr("Trait.\\_ParentTerm", colnames(currentReviewCPAGdata))>0)

              # Compose list of unique EFO parent terms
              if(length(reviewResultsEFOparentCol)>0) {
                a <- na.omit(unlist(currentReviewCPAGdata[,reviewResultsEFOparentCol]))
                if(length(a)>0) {
                  reviewResultsEFOparentTerm <<- sort(unique(trimws(unlist(strsplit(a, ",")))))
                } else {
                  reviewResultsEFOparentTerm <<- vector("character")
                }
              } else {
                reviewResultsEFOparentTerm <<- vector("character")
              }
              updateSelectInput(session, "reviewSelectionFilterEFOparent", choices=reviewResultsEFOparentTerm)

              # Reset filters
              # Updates are not effective until all output has been pushed to the client
              # This does not guarantee that new values are in effect at the time of the call to renderResults()
              # Therefore, assign values and pass constants as parameter values to renderResults()
              updateTextInput(session, "reviewSelectionFilterTrait", value="")
              updateTextInput(session, "reviewSelectionFilterSNP", value="")
              #updateTextInput(session, "reviewSelectionFilterEFO", value="")
              updateTextInput(session, "reviewSelectionFilterEFOparent", value=vector("character"))
              updateCheckboxInput(session, "reviewSelectionIncludeTableCompoundEFO", value=F)
              renderResults(data=currentReviewCPAGdata,
                            rowFilter="reviewResultsTableRowFilter",
                            filterTrait="",
                            filterSNP="",
                            truncSNP=input$reviewSelectionIncludeTableAllSNPshare,
                            includeCompoundEFO=F,
                            efoCol=reviewResultsEFOcol,
                            filterEFOparent=vector("character"),
                            efoParentCol=reviewResultsEFOparentCol,
                            tabset="tabsetReviewResults",
                            tableTab="tabPanelReviewResultsTable",
                            tableName="reviewResultsTable",
                            plotTab="tabPanelReviewResultsHeatmap",
                            plotName="reviewResultsHeatmap",
                            heatmapMetric=input$reviewSelectionHeatmapMetric,
                            heatmapNphenotype=input$reviewSelectionHeatmapNphenotype)
            } else {
              output$reviewResultsTable <- DT::renderDataTable(NULL)
              output$reviewResultsHeatmap <- renderPlotly(NULL)
              #showNotification("No CPAG data exist for specified GWAS features\n", type="error", duration=notifyDuration)
              showModal(modalDialog(HTML("No CPAG data exist for specified GWAS features"),
                        title="iCPAGdb", size="m", easyClose=T, footer=modalButton("OK")))
            }

          }

        }

      }, ignoreInit=T)

    #######################################################################################################
    # Observe event function:  Review tab - apply filters
    #######################################################################################################

    observeEvent(input$reviewSelectionFilterApply, {

      # Input values are available since they are retrieved as part of the event triggered here 
      renderResults(data=currentReviewCPAGdata,
                    rowFilter="reviewResultsTableRowFilter",
                    filterTrait=input$reviewSelectionFilterTrait,
                    filterSNP=input$reviewSelectionFilterSNP,
                    truncSNP=input$reviewSelectionIncludeTableAllSNPshare,
                    includeCompoundEFO=input$reviewSelectionIncludeTableCompoundEFO,
                    efoCol=reviewResultsEFOcol,
                    filterEFOparent=input$reviewSelectionFilterEFOparent,
                    efoParentCol=reviewResultsEFOparentCol,
                    tabset="tabsetReviewResults",
                    tableTab="tabPanelReviewResultsTable",
                    tableName="reviewResultsTable",
                    plotTab="tabPanelReviewResultsHeatmap",
                    plotName="reviewResultsHeatmap",
                    heatmapMetric=input$reviewSelectionHeatmapMetric,
                    heatmapNphenotype=input$reviewSelectionHeatmapNphenotype)

    }, ignoreInit=T)

    #######################################################################################################
    # Observe event function:  Review tab - clear filters
    #######################################################################################################

    observeEvent(input$reviewSelectionFilterClear, {

      # Updates are not effective until all output has been pushed to the client
      # This does not guarantee that new values are in effect at the time of the call to renderResults()
      # Therefore, assign values and pass constants as parameter values to renderResults()
      #updateTextInput(session, "reviewSelectionFilterTrait", value="")
      #updateTextInput(session, "reviewSelectionFilterSNP", value="")
      #updateTextInput(session, "reviewSelectionFilterEFO", value="")
      #updateTextInput(session, "reviewSelectionFilterEFOparent", value=vector("character"))
      #updateCheckboxInput(session, "reviewSelectionIncludeTableCompoundEFO", value=F)
      renderResults(data=currentReviewCPAGdata,
                    rowFilter="reviewResultsTableRowFilter",
                    filterTrait="",
                    filterSNP="",
                    truncSNP=input$reviewSelectionIncludeTableAllSNPshare,
                    includeCompoundEFO=F,
                    efoCol=reviewResultsEFOcol,
                    filterEFOparent=vector("character"),
                    efoParentCol=reviewResultsEFOparentCol,
                    tabset="tabsetReviewResults",
                    tableTab="tabPanelReviewResultsTable",
                    tableName="reviewResultsTable",
                    plotTab="tabPanelReviewResultsHeatmap",
                    plotName="reviewResultsHeatmap",
                    heatmapMetric=input$reviewSelectionHeatmapMetric,
                    heatmapNphenotype=input$reviewSelectionHeatmapNphenotype)

    }, ignoreInit=T)

    #######################################################################################################
    # Configure a download handler for the results table of the review tab
    # This is an efficient alternative to enabling the download buttons of the review results data table
    # One drawback is that records must be filtered here to agree with those contained in the table
    # (assuming that filters have been applied to the data table)
    #######################################################################################################

    output$reviewResultsDownload <-
      downloadHandler(filename="iCPAGdb-Review-Results.csv",
                      contentType="text/csv",
                      content=function(file) {
                                # Compose filtered observation indices
                                k <- resultsTableComposeFilterIndex(
                                       data=currentReviewCPAGdata,
                                       filterTrait=input$reviewSelectionFilterTrait,
                                       filterSNP=input$reviewSelectionFilterSNP,
                                       includeCompoundEFO=input$reviewSelectionIncludeTableCompoundEFO,
                                       efoCol=reviewResultsEFOcol,
                                       filterEFOparent=input$reviewSelectionFilterEFOparent,
                                       efoParentCol=reviewResultsEFOparentCol)
                                write.table(currentReviewCPAGdata[k,], file, row.names=F, col.names=T, quote=T, sep=",")
                              }
                     )

    #######################################################################################################
    # INACTIVE
    # Observe event function:  Review tab - table filtered
    # Note that input$reviewResultsTable_rows_all is a variable maintained by the DT package, containing
    # the row indices corresponding to the current set of filtered rows
    # It is made reactive by Shiny, so that changes to the DT variable cause events here
    # Note that this event is executed whenever reviewResultsTable is filtered, which may be done
    # repeatedly as the user enters characters into a filter field ("a" followed by "al" followed by
    # "all" ... followed by "allergic"
    # This may cause repeated, and lengthy, rendering of intermediate heatmaps that are not of interest
    # It is also executed when review table is initially rendered
    #######################################################################################################

    # Function is currently inactive
    # Renaming this function disables the event for reviewResultsTable_rows_all
    observeEvent(input$reviewResultsTable_rows_all0, {

      # Update global row index vector
      reviewResultsTableRowFilter <<- input$reviewResultsTable_rows_all
      cat(paste("tfilt, n-row-filter:  ", length(reviewResultsTableRowFilter), "\n", sep=""), file=stderr())

      # Render heatmap based on selected metric
      renderHeatmap(data=currentReviewCPAGdata[reviewResultsTableRowFilter,],
                    xCol="Trait1", yCol="Trait2", valueCol=input$reviewSelectionHeatmapMetric,
                    log10Value=(input$reviewSelectionHeatmapMetric %in% c("P_fisher", "Padj_Bonferroni", "Padj_FDR")),
                    nPhenotype=input$reviewSelectionHeatmapNphenotype,
                    key.title=ifelse(input$reviewSelectionHeatmapMetric %in% c("P_fisher", "Padj_Bonferroni", "Padj_FDR"),
                                     paste("-log10(", input$reviewSelectionHeatmapMetric, ")", sep=""),
                                     input$reviewSelectionHeatmapMetric),
                    tabset="tabsetReviewResults", tableTab="tabPanelReviewResultsTable",
                    plotTab="tabPanelReviewResultsHeatmap", plotName="reviewResultsHeatmap")

    }, ignoreInit=T)

    #######################################################################################################
    # Observe event function:  Review tab - change in include all SNPs in shared SNP option
    # Note that changes are applied to the review rsults table without an reactive function
    # Is this due to input$reviewSelectionIncludeTableAllSNPshare appearing in the observe event for
    # input$reviewSelectionTable_cells_selected, although it is not, itself named in the list of
    # ui variables to be monitored?
    #######################################################################################################

    #observeEvent(input$reviewSelectionIncludeTableAllSNPshare, {}, ignoreInit=T)

    #######################################################################################################
    # Observe event function:  Review tab - change in heatmap metric or top significant phenotype
    # pairs to plot
    #######################################################################################################

    observeEvent(c(input$reviewSelectionHeatmapNphenotype, input$reviewSelectionHeatmapMetric), {

      cat(paste("nph, n-row-filter:  ", length(reviewResultsTableRowFilter), "\n", sep=""), file=stderr())

      # Render heatmap based on selected metric
      renderHeatmap(data=currentReviewCPAGdata[reviewResultsTableRowFilter,],
                    xCol="Trait1", yCol="Trait2", valueCol=input$reviewSelectionHeatmapMetric,
                    log10Value=(input$reviewSelectionHeatmapMetric %in% c("P_fisher", "Padj_Bonferroni", "Padj_FDR")),
                    nPhenotype=input$reviewSelectionHeatmapNphenotype,
                    key.title=ifelse(input$reviewSelectionHeatmapMetric %in% c("P_fisher", "Padj_Bonferroni", "Padj_FDR"),
                                     paste("-log10(", input$reviewSelectionHeatmapMetric, ")", sep=""),
                                     input$reviewSelectionHeatmapMetric),
                    tabset="tabsetReviewResults", tableTab="tabPanelReviewResultsTable",
                    plotTab="tabPanelReviewResultsHeatmap", plotName="reviewResultsHeatmap")

      # Move focus to the plot tab
      updateTabsetPanel(session, "tabsetReviewResults", "tabPanelReviewResultsHeatmap")

    }, ignoreInit=T)

    #######################################################################################################
    # Observe event:  Enab;e/disable "secret" features
    #######################################################################################################

    observeEvent(input$featureEnable, {

      if(tolower(input$featureEnable)=="exp+") {
        showTab("tabsetCPAG", "tabPanelExplore")
      } else if(tolower(input$featureEnable)=="exp-") {
        hideTab("tabsetCPAG", "tabPanelExplore")
      } else if(tolower(input$featureEnable)=="phen+") {
        phenH2P2 <<- "full"
      } else if(tolower(input$featureEnable)=="phen+") {
        phenH2P2 <<- "limited"
      }

    }, ignoreInit=F)

    #######################################################################################################
    # Function:  Explore tab, read CPAG file, render table and heatmap
    #######################################################################################################

    exploreReadCPAGandRender <- function(cpagFile) {

      # Read CPAG results
      currentExploreCPAGdata <<- read.table(cpagFile, header=T, sep=",", quote="\"")

      # Construct observation filter
      exploreResultsTableRowFilter <<- 1:nrow(currentExploreCPAGdata)

      # Index EFO columns
      exploreResultsEFOcol <<- which(regexpr("Trait.\\_EFO", colnames(currentExploreCPAGdata))>0)
      exploreResultsEFOparentCol <<- which(regexpr("Trait.\\_ParentTerm", colnames(currentExploreCPAGdata))>0)

      # Compose list of unique EFO parent terms
      if(length(exploreResultsEFOparentCol)>0) {
        a <- na.omit(unlist(currentExploreCPAGdata[,exploreResultsEFOparentCol]))
        if(length(a)>0) {
          exploreResultsEFOparentTerm <<- sort(unique(trimws(unlist(strsplit(a, ",")))))
        } else {
          exploreResultsEFOparentTerm <<- vector("character")
        }
      } else {
        exploreResultsEFOparentTerm <<- vector("character")
      }
      updateSelectInput(session, "exploreSelectionFilterEFOparent", choices=exploreResultsEFOparentTerm)

      # Render table and heatmap
      renderResults(data=currentExploreCPAGdata,
                    rowFilter="exploreResultsTableRowFilter",
                    filterTrait="",
                    filterSNP="",
                    truncSNP=input$exploreSelectionIncludeTableAllSNPshare,
                    includeCompoundEFO=F,
                    efoCol=exploreResultsEFOcol,
                    filterEFOparent=vector("character"),
                    efoParentCol=exploreResultsEFOparentCol,
                    tabset="tabsetExploreResults",
                    tableTab="tabPanelExploreResultsTable",
                    tableName="exploreResultsTable",
                    plotTab="tabPanelExploreResultsHeatmap",
                    plotName="reviewExploresHeatmap",
                    heatmapMetric=input$exploreSelectionHeatmapMetric,
                    heatmapNphenotype=input$exploreSelectionHeatmapNphenotype)

    }

    #######################################################################################################
    # Observe event function:  Explore tab - compute CPAG button
    #######################################################################################################

    observeEvent(input$exploreCompute, {

      # Compose output file name
      # Format is GWAS1-p1-GWAS2-p2-ldpop.csv
      # Order of appearance (GWAS1, GWAS2) is H2P2, NHGRI, other
      # GWAS2 and p2 are excluded if GWAS1=GWAS2 
      s <- sort(c(input$exploreSource1, input$exploreSource2))
      p <- as.integer(c(input$explorePthresh1, input$explorePthresh2)[order(c(input$exploreSource1, input$exploreSource2))])
      kh <- which(s=="H2P2")
      kn <- which(s=="NHGRI")
      ko <- which(!s %in% c("H2P2", "NHGRI"))
      p[kh] <- max(p[kh], 5, na.rm=T)
      p[kn] <- max(p[kn], 8, na.rm=T)
      p[ko] <- max(p[ko], 3, na.rm=T)
      if(length(kh)>0) {
        cpag <- paste("H2P2-p1e-", sprintf("%02d", p[kh[1]]),
                      ifelse(length(kh)==1,
                             paste("-", s[-kh], ifelse(length(kn)>0, "-p5e-", "-p1e-"), sprintf("%02d", p[-kh]),  sep=""),
                             ""),
                      sep="")
      } else if(length(kn)>0) {
        cpag <- paste("NHGRI-p5e-", sprintf("%02d", p[kn[1]]),
                      ifelse(length(kn)==1, paste("-", s[-kn], "-p5e-", sprintf("%02d", p[-kn]),  sep=""), ""),
                      sep="")
      } else if(s[1]!=s[2]) {
        cpag <- paste(s[1], "-p1e-", sprintf("%02d", p[1]), "-", s[2], "-p1e-", sprintf("%02d", p[2]), sep="")
      } else {
        cpag <- paste(s[1], "-p1e-", sprintf("%02d", p[1]), sep="")
      }
      cpag <- paste(cpag, "-", input$exploreLDpop, ".csv", sep="")
      existFile <- F

      # Test for existense of output file
      if(length(which(dir(outDir, pattern=cpag, ignore.case=T)==cpag))>0) {

        # File exists - read data, render table and heatmap
        exploreReadCPAGandRender(paste(outDir, "/", cpag, sep=""))

      } else {

        # File does not exist
        # Compute results for requested GWAS sets, p-thresholds, and LD population
        # Compose python command to execute CPAG function
        pycmd <- paste(pyexe, " main.py cpagdb --threads ", threads, 
                       " --subtype ", s[1], " --subtype ", s[2],
                       ifelse(s[1]=="H2P2" | s[2]=="H2P2",
                              paste(" --H2P2-Pcut 1e-", p[kh[1]], sep=""), ""),
                       ifelse(s[1]=="NHGRI" | s[2]=="NHGRI",
                              paste(" --NHGRI-Pcut 5e-", p[kn[1]], sep=""), ""),
                       ifelse(!s[1] %in% c("H2P2", "NHGRI") | !s[2] %in% c("H2P2", "NHGRI"),
                              paste(" --Pcut 1e-", p[ko[1]], sep=""), ""),
                       " --lddb-pop ", input$exploreLDpop,
                       " --outfile \"", outDir, "/", cpag, "\"", sep="")
        cat(paste("\npycmd:  ", pycmd, "\n\n", sep=""), file=stderr())
        # Compose commands to include ontology, for NHGRI results only
        if(s[1]=="NHGRI" | s[2]=="NHGRI") {
          pycmd2 <- paste(pyexe, " main.py post_analysis --anno-ontology --anno-cols Trait1",
                          " --infile \"", outDir, "/", cpag, "\" --outfile \"", outDir, "/", cpag, "\"", sep="")
          cat(paste("\npycmd2:  ", pycmd2, "\n\n", sep=""), file=stderr())
        } else {
          pycmd2 <- ""
        }
        if(s[1]=="NHGRI" & s[2]=="NHGRI") {
          pycmd3 <- paste(pyexe, " main.py post_analysis --anno-ontology --anno-cols Trait2",
                          " --infile \"", outDir, "/", cpag, "\" --outfile \"", outDir, "/", cpag, "\"", sep="")
          cat(paste("\npycmd3:  ", pycmd3, "\n\n", sep=""), file=stderr())
        } else {
          pycmd3 <- ""
        }
        # Create progress meter prior to launcing parallel process, since that process is independent of
        # the current Shiny environment
        progress <- shiny::Progress$new()
        progress$set(message=paste("Computing CPAG results. Execution time may be several minutes.",
                                   " Start time is ", format(Sys.time(), "%X %Z"), sep=""), value=0.5)
        # Compute CPAG from within an asynchronous (parallel, independent) process
        # When future() instructions complete, the result (the OS command return value) is piped to a function
        # that is responsible forrendering table and plot results
        # The %...>% is a "promise" pipe that waits on results from the future operation then proceeds
        # The progress object is closed within the promise function so that it is visible while future()
        # instructions are executed
        # Important features:  future() functions asynchronously when plan(multisession) is specified
        # R is, by default, synchronous, executing a single instruction before another
        # With synchronous operation, Shiny waits on completion of any instruction currently being executed
        # before processing any user initiated activity - if user1 initiates a lengthy process in a given Shiny app,
        # user two cannot even load the app until the lengthy operation completes
        # Using asyncgronous operations, lengthy instructions can be executed without delaying other activity
        # The %...>% pipe delays processing of further instructions while the output of lengthy asynchrounous
        # instructions are executed
        future({
          t <- proc.time()
          # Rout stderr to stdout with 2>&1 (but tends to eliminate all output, even when errors occur)
          y <- suppressWarnings(system(pycmd, intern=T))
          if(is.null(attributes(y)) & pycmd2!="")
            y <- suppressWarnings(system(pycmd2, intern=T))
          if(is.null(attributes(y)) & pycmd3!="")
            y <- suppressWarnings(system(pycmd3, intern=T))
          t <- proc.time()-t
          cat(paste(t, "\n", collapse=" "), file=stderr())
          y}) %...>%
            (function(y) {
               progress$close()
               if(is.null(attributes(y))) {
                 exploreReadCPAGandRender(paste(outDir, "/", cpag, sep=""))
               } else {
                 output$exploreResultsTable <- DT::renderDataTable(NULL)
                 output$exploreResultsHeatmap <- renderPlotly(NULL)
                 #showNotification("Error while computing CPAG", type="error", duration=notifyDuration)
                 showModal(modalDialog(HTML("Error while computing CPAG"),
                           title="iCPAGdb", size="m", easyClose=T, footer=modalButton("OK")))

               }
             })

      }

      # Clear filters
      updateTextInput(session, "exploreSelectionFilterTrait", value="")
      updateTextInput(session, "exploreSelectionFilterSNP", value="")
      #updateTextInput(session, "exploreSelectionFilterEFO", value="")
      updateTextInput(session, "exploreSelectionFilterEFOparent", value=vector("character"))
      updateCheckboxInput(session, "exploreSelectionIncludeTableCompoundEFO", value=F)

    }, ignoreInit=T)

    #######################################################################################################
    # Observe event function:  Explore tab - change in heatmap metric or top significant phenotype
    # pairs to plot
    #######################################################################################################

    observeEvent(c(input$exploreSelectionHeatmapNphenotype, input$exploreSelectionHeatmapMetric), {

      cat(paste("nph, n-row-filter:  ", length(exploreResultsTableRowFilter), "\n", sep=""), file=stderr())

      # Render heatmap based on selected metric
      renderHeatmap(data=currentExploreCPAGdata[exploreResultsTableRowFilter,],
                    xCol="Trait1", yCol="Trait2", valueCol=input$exploreSelectionHeatmapMetric,
                    log10Value=(input$exploreSelectionHeatmapMetric %in% c("P_fisher", "Padj_Bonferroni", "Padj_FDR")),
                    nPhenotype=input$exploreSelectionHeatmapNphenotype,
                    key.title=ifelse(input$exploreSelectionHeatmapMetric %in% c("P_fisher", "Padj_Bonferroni", "Padj_FDR"),
                                     paste("-log10(", input$exploreSelectionHeatmapMetric, ")", sep=""),
                                     input$exploreSelectionHeatmapMetric),
                    tabset="tabsetExploreResults", tableTab="tabPanelExploreResultsTable",
                    plotTab="tabPanelExploreResultsHeatmap", plotName="exploreResultsHeatmap")

      # Move focus to the plot tab
      updateTabsetPanel(session, "tabsetExploreResults", "tabPanelExploreResultsHeatmap")

    }, ignoreInit=T)

    #######################################################################################################
    # Observe event function:  Explore tab - apply filters
    #######################################################################################################

    observeEvent(input$exploreSelectionFilterApply, {

      # Input values are available since they are retrieved as part of the event triggered here 
      renderResults(data=currentExploreCPAGdata,
                    rowFilter="exploreResultsTableRowFilter",
                    filterTrait=input$exploreSelectionFilterTrait,
                    filterSNP=input$exploreSelectionFilterSNP,
                    truncSNP=input$exploreSelectionIncludeTableAllSNPshare,
                    includeCompoundEFO=input$exploreSelectionIncludeTableCompoundEFO,
                    efoCol=exploreResultsEFOcol,
                    filterEFOparent=input$exploreSelectionFilterEFOparent,
                    efoParentCol=exploreResultsEFOparentCol,
                    tabset="tabsetExploreResults",
                    tableTab="tabPanelExploreResultsTable",
                    tableName="exploreResultsTable",
                    plotTab="tabPanelExploreResultsHeatmap",
                    plotName="reviewExploresHeatmap",
                    heatmapMetric=input$exploreSelectionHeatmapMetric,
                    heatmapNphenotype=input$exploreSelectionHeatmapNphenotype)

    }, ignoreInit=T)

    #######################################################################################################
    # Observe event function:  Explore tab - clear filters
    #######################################################################################################

    observeEvent(input$exploreSelectionFilterClear, {

      # Updates are not effective until all output has been pushed to the client
      # This does not guarantee that new values are in effect at the time of the call to renderResults()
      # Therefore, assign values and pass constants as parameter values to renderResults()
      #updateTextInput(session, "exploreSelectionFilterTrait", value="")
      #updateTextInput(session, "exploreSelectionFilterSNP", value="")
      #updateTextInput(session, "exploreSelectionFilterEFO", value="")
      #updateTextInput(session, "exploreSelectionFilterEFOparent", value=vector("character"))
      #updateCheckboxInput(session, "exploreSelectionIncludeTableCompoundEFO", value=F)
      renderResults(data=currentExploreCPAGdata,
                    rowFilter="exploreResultsTableRowFilter",
                    filterTrait="",
                    filterSNP="",
                    truncSNP=input$exploreSelectionIncludeTableAllSNPshare,
                    includeCompoundEFO=F,
                    efoCol=exploreResultsEFOcol,
                    filterEFOparent=vector("character"),
                    efoParentCol=exploreResultsEFOparentCol,
                    tabset="tabsetExploreResults",
                    tableTab="tabPanelExploreResultsTable",
                    tableName="exploreResultsTable",
                    plotTab="tabPanelExploreResultsHeatmap",
                    plotName="reviewExploreHeatmap",
                    heatmapMetric=input$exploreSelectionHeatmapMetric,
                    heatmapNphenotype=input$exploreSelectionHeatmapNphenotype)

    }, ignoreInit=T)

    #######################################################################################################
    # Configure a download handler for the results table of the explore tab
    # This is an efficient alternative to enabling the download buttons of the explore results data table
    # One drawback is that records must be filtered here to agree with those contained in the table
    # (assuming that filters have been applied to the data table)
    #######################################################################################################

    output$exploreResultsDownload <-
      downloadHandler(filename="iCPAGdb-Explore-Results.csv",
                      contentType="text/csv",
                      content=function(file) {
                                # Compose filtered observation indices
                                k <- resultsTableComposeFilterIndex(
                                       data=currentExploreCPAGdata,
                                       filterTrait=input$exploreSelectionFilterTrait,
                                       filterSNP=input$exploreSelectionFilterSNP,
                                       includeCompoundEFO=input$exploreSelectionIncludeTableCompoundEFO,
                                       efoCol=exploreResultsEFOcol,
                                       filterEFOparent=input$exploreSelectionFilterEFOparent,
                                       efoParentCol=exploreResultsEFOparentCol)
                                write.table(currentExploreCPAGdata[k,], file, row.names=F, col.names=T, quote=T, sep=",")
                              }
                     )

    #######################################################################################################
    # Observe event function:  User compute tab - browse file selected
    #######################################################################################################

    observeEvent(input$userComputeBrowseFile, {

        # Read first lines of file
        x <- scan(input$userComputeBrowseFile[,"datapath"], "character", n=5, sep="\n", quote="", quiet=T)

        #showNotification(paste("File is ready.  Specify column delimiter and headings then press \"Upload file.\"",
        #                       ifelse(length(grep("\t", x)>0),
        #                              "  NOTE THAT TAB DELIMITERS HAVE BEEN DETECTED IN THE FILE.",
        #                              ""),
        #                       sep=""),
        #                 id="fileInput", duration=NULL, type="message")

        # Display lines and test for tabs
        showModal(modalDialog( 
          #sideBarLayout()
          HTML(paste("UPLOADED FILE STRUCTURE IS:<br><br>",
                     paste(x, collapse="<br>", sep=""),
                     "<br><br>SPECIFY COLUMN DELIMETER AND HEADINGS THEN PROCEED TO &nbsp;&nbsp;\"Compute CPAG\"<br><br>",
                     ifelse(length(grep("\t", x)>0),
                            "<b>NOTE THAT TAB DELIMITERS HAVE BEEN DETECTED IN THE FILE</b>",
                            ""),
                     sep="")),
          title="iCPAGdb", size="l", easyClose=T,
          footer=modalButton("OK")))
 
    }, ignoreInit=T)


    #######################################################################################################
    # Function:  Read uploaded (temp) file and save as a csv
    # Althought the uploaded file name could be passed to the CPAG function, assessing delimiters
    # and saving in a known format takes advantage of file format verification while the current
    # process has control
    # Limited information is available from the CPAG function
    # Further, retaining a copy of user uploaded data permits a review of records if an error is reported 
    #######################################################################################################

    copyValidateUploadFile <- function() {

      uGWASfile <- ""
      uGWAScol <- ""

      # Verify existence of uploaded file
      if(!is.null(input$userComputeBrowseFile[,"name"])) {
        cat(paste("User GWAS file: ", input$userComputeBrowseFile[,"name"], "; ",
                  input$userComputeBrowseFile[,"datapath"], "\n", sep=""), file=stderr())
        progress <- shiny::Progress$new()
        progress$set(message="Reading Data", value=0.33)
        userFile <- read.table(input$userComputeBrowseFile[,"datapath"], header=T,
                    sep=c("Comma"=",", "Tab"="\t")[input$userComputeDelimiter],
                    quote=ifelse(input$userComputeDelimiter=="Comma", "\"", ""))
        progress$close()

        # Locate user specified columns in uploaded file
        # At time of implementation, each upload file is expected to contain a single phenotype
        #k1 <- which(tolower(colnames(userFile))==tolower(input$userComputePhenotypeCol))
        k2 <- which(tolower(colnames(userFile))==tolower(input$userComputeSNPcol))
        k3 <- which(tolower(colnames(userFile))==tolower(input$userComputePcol))
        #if(length(k1)==1) {
          if(length(k2)==1) {
            if(length(k3)==1) {

              # Save colums in global for use in computation function
              # Exclude phenotype in the current version
              #uGWAScol <- colnames(userFile)[c(k1, k2, k3)]
              #names(uGWAScol) <- c("phenotype", "snp", "p")
              uGWAScol <- colnames(userFile)[c(k2, k3)]
              names(uGWAScol) <- c("snp", "p")
              cat(paste("User GWAS columns: ", paste(uGWAScol, collapse=", "), "\n", sep=""), file=stderr())
              cat(paste("User GWAS file rows: ", nrow(userFile), "\n", sep=""), file=stderr())

              # Compose file name to use in making local copy of user's input GWAS and output CPAG results
              # Append random sequence of digits, test for existence, and repeat until
              # unused file name is found
              # Save in global for use in computation function
              uGWASfile <- ""
              while(uGWASfile=="") {
                a <- paste(sample(0:9, 10, replace=T), collapse="", sep="")
                f <- paste(gsub(".txt", "",
                             gsub(".csv", "",
                               gsub(".dat", "", input$userComputeBrowseFile[,"name"], fixed=T), fixed=T), fixed=T),
                           "-", a, sep="")
                if(length(dir(path=inDirUserCompute, pattern=paste(f, ".csv", sep="")))==0 &
                   length(dir(path=outDirUserCompute, pattern=paste(f, ".csv", sep="")))==0)
                  uGWASfile <- f
              }

              # Write local copy of input data and notify user to proceed
              if(uGWASfile!="") {
                progress <- shiny::Progress$new()
                progress$set(message="Reading upload data", value=0.67)
                write.table(userFile[,uGWAScol], paste(inDirUserCompute, "/", uGWASfile, ".csv", sep=""), row.names=F, col.names=T, sep=",", quote=T)
                progress$close()
                result <- ""
                cat(paste("User input GWAS file created: ", uGWASfile, ".csv", "\n", sep=""), file=stderr())
              } else {
                cat(paste("Cannot create user upload file for: ", input$userComputeBrowseFile[,"name"], "\n", sep=""), file=stderr())
                result <- "Cannot process upload file.  Please upload again."
              }

            } else {
              result <- paste("Specified P (significance) column is missing in uploaded GWAS file<br>",
                              "Input columns are: ", paste(colnames(userFile), collapse=", "), sep="")
            }
          } else {
            result <- paste("Specified SNP column is missing in uploaded GWAS file<br>",
                            "Input columns are: ", paste(colnames(userFile), collapse=", "), sep="")
          }
        #} else {
        #  result <- paste("Specified phenotype (trait) column is missing in uploaded GWAS file. ",
        #                  "Input columns are: ", paste(colnames(userFile), collapse=", "), sep="")
        #}
      } else {
        result <- "No file specified.  Please upload a file"
      }

      return(list("result"=result, "uGWASfile"=uGWASfile, "uGWAScol"=uGWAScol))

    }

    #######################################################################################################
    # Observe event function:  User compute tab - compute CPAG button
    #######################################################################################################

    observeEvent(input$userComputeCompute, {

      # Validate and copy input file
      f <- copyValidateUploadFile()
      if(f[["result"]]=="") {

        # Verify existence of copied file
        if(length(which(dir(inDirUserCompute, paste(f[["uGWASfile"]], ".csv", sep=""))==
                        paste(f[["uGWASfile"]], ".csv", sep="")))>0) {

          # Compose output file name
          # Format is UserFile-p1-GWAS2-p2-ldpop.csv
          cpagFile <- paste(f[["uGWASfile"]], "-p1e-", sprintf("%02d", input$userComputePthresh1), "-", input$userComputeSource2,
                            ifelse(input$userComputeSource2=="H2P2",
                                   paste("-p1e-", sprintf("%02d", max(input$userComputePthresh1, 5)), sep=""),
                                   ifelse(input$userComputeSource2=="NHGRI",
                                          paste("-p5e-", sprintf("%02d", max(input$userComputePthresh2, 8)), sep=""),
                                          paste("-p1e-", sprintf("%02d", max(input$userComputePthresh1, 3)), sep=""))),
                            "-", input$userComputeLDpop, ".csv", sep="")

          # Compose python command
          pycmd <- paste(pyexe, " main.py usr-gwas --threads ", threads,
                         " --infile \"", inDirUserCompute, "/", f[["uGWASfile"]], ".csv\"",
                         " --delimitor \",\"",
                         " --usr-pcut 1e-", max(input$userComputePthresh1, 3),
                         #" --usr-pheno-name ", f[["uGWAScol"]]["phenotype"],
                         " --SNPcol \"", f[["uGWAScol"]]["snp"], "\"",
                         " --Pcol \"", f[["uGWAScol"]]["p"], "\"",
                         " --querydb ", input$userComputeSource2,
                         " --cpagdb-pcut",
                         ifelse(input$userComputeSource2=="H2P2",
                                paste(" 1e-", max(input$userComputePthresh2, 5), sep=""),
                                ifelse(input$userComputeSource2=="NHGRI",
                                       paste(" 5e-", max(input$userComputePthresh2, 8), sep=""),
                                       paste(" 1e-", max(input$userComputePthresh1, 3), sep=""))),
                         " --ld-clump 1",
                         " --lddb-pop ", input$userComputeLDpop,
                         " --outfile \"", outDirUserCompute, "/", cpagFile, "\"", sep="")

          cat(paste("\npycmd:  ", pycmd, "\n\n", sep=""), file=stderr())

          # Compose commands to include ontology, for NHGRI results only
          if(input$userComputeSource2=="NHGRI") {
            pycmd2 <- paste(pyexe, " main.py post_analysis --anno-ontology --anno-cols Trait2",
                            " --infile \"", outDirUserCompute, "/", cpagFile, "\" --outfile \"",
                            outDirUserCompute, "/", cpagFile, "\"", sep="")
            cat(paste("\npycmd2:  ", pycmd2, "\n\n", sep=""), file=stderr())
          } else {
            pycmd2 <- ""
          }

          # Create progress meter prior to launcing parallel process, since that process is independent of
          # the current Shiny environment
          progress <- shiny::Progress$new()
          progress$set(message=paste("Computing CPAG results. Estimated execution time is between 30 and 60 seconds.",
                                     " Start time is ", format(Sys.time(), "%X %Z"), sep=""), value=0.5)

          # Compute CPAG from within an asynchronous (parallel, independent) process
          # Wnen future() instructions complete, the result (the OS command return value) is piped to a function
          # that is responsible forrendering table and plot results
          # The %...>% is a "promise" pipe that waits on results from the future operation then proceeds
          # The progress object is closed within the promise function so that it is visible while future()
          # instructions are executed
          # Important features:  future() functions asynchronously when plan(multisession) is specified
          # R is, by default, synchronous, executing a single instruction before another
          # With synchronous operation, Shiny waits on completion of any instruction currently being executed
          # before processing any user initiated activity - if user1 initiates a lengthy process in a given Shiny app,
          # user two cannot even load the app until the lengthy operation completes
          # Using asyncgronous operations, lengthy instructions can be executed without delaying other activity
          # The %...>% pipe delays processing of further instructions while the output of lengthy asynchrounous
          # instructions are executed
          future({
            t <- proc.time()
            # Route stderr to stdout with 2>&1 (but it tends to eliminate all output, even when errors occur)
            y <- suppressWarnings(system(pycmd, intern=T))
            if(is.null(attributes(y)) & pycmd2!="")
              y <- suppressWarnings(system(pycmd2, intern=T))
            t <- proc.time()-t
            cat(paste(t, "\n", collapse=" "), file=stderr())
            y}) %...>%
              (function(y) {
                 progress$close()
                 if(is.null(attributes(y))) {
                   # Read CPAG results
                   # Verify existence of output file
                   if(length(which(dir(outDirUserCompute, pattern=cpagFile)==cpagFile))>0) {
                     currentUserComputeCPAGdata <<- read.table(paste(outDirUserCompute, "/", cpagFile, sep=""),
                                                               header=T, sep=",", quote="\"")
                     # Construct observation filter
                     userComputeResultsTableRowFilter <<- 1:nrow(currentUserComputeCPAGdata)

                     # Index EFO columns
                     userComputeResultsEFOcol <<- which(regexpr("Trait.\\_EFO", colnames(currentUserComputeCPAGdata))>0)
                     userComputeResultsEFOparentCol <<- which(regexpr("Trait.\\_ParentTerm", colnames(currentUserComputeCPAGdata))>0)

                     # Compose list of unique EFO parent terms
                     if(length(userComputeResultsEFOparentCol)>0) {
                       a <- na.omit(unlist(currentUserComputeCPAGdata[,userComputeResultsEFOparentCol]))
                       if(length(a)>0) {
                         userComputeResultsEFOparentTerm <<- sort(unique(trimws(unlist(strsplit(a, ",")))))
                       } else {
                         userComputeResultsEFOparentTerm <<- vector("character")
                       }                     
                     } else {
                       userComputeResultsEFOparentTerm <<- vector("character")
                     }
                     updateSelectInput(session, "userComputeSelectionFilterEFOparent", choices=userComputeResultsEFOparentTerm)

                     # Render table and heatmap
                     renderResults(data=currentUserComputeCPAGdata,
                                   rowFilter="userComputeResultsTableRowFilter",
                                   filterTrait="",
                                   filterSNP="",
                                   truncSNP=input$userComputeSelectionIncludeTableAllSNPshare,
                                   includeCompoundEFO=F,
                                   efoCol=userComputeResultsEFOcol,
                                   filterEFOparent=vector("character"),
                                   efoParentCol=userComputeResultsEFOparentCol,
                                   tabset="tabsetUserComputeResults",
                                   tableTab="tabPanelUserComputeResultsTable",
                                   tableName="userComputeResultsTable",
                                   plotTab="tabPanelUserComputeResultsHeatmap",
                                   plotName="userComputeResultsHeatmap",
                                   heatmapMetric=input$userComputeSelectionHeatmapMetric,
                                   heatmapNphenotype=input$userComputeSelectionHeatmapNphenotype)

                     # Clear filter controls
                     updateTextInput(session, "userComputeSelectionFilterTrait", value="")
                     updateTextInput(session, "userComputeSelectionFilterSNP", value="")
                     #updateTextInput(session, "userComputeSelectionFilterEFO", value="")
                     updateTextInput(session, "userComputeSelectionFilterEFOparent", value=vector("character"))
                     updateCheckboxInput(session, "userComputeSelectionIncludeTableCompoundEFO", value=F)

                   } else {
                     output$userComputeResultsTable <- DT::renderDataTable(NULL)
                     output$userComputeResultsHeatmap <- renderPlotly(NULL)
                     #showNotification(HTML(paste("No results returned for uploaded GWAS requested CPAG parameter values<br>",
                     #                            "Consider relaxing p-thresholds"), sep=""),
                     #                 type="error", duration=2*notifyDuration)
                     showModal(modalDialog(HTML(paste("No results returned for requested CPAG parameter values<br>",
                                                      "Consider relaxing p-thresholds"), sep=""),
                                           title="iCPAGdb", size="m", easyClose=T, footer=modalButton("OK")))
                   }
                 } else {
                   output$userComputeResultsTable <- DT::renderDataTable(NULL)
                   output$userComputeResultsHeatmap <- renderPlotly(NULL)
                   #showNotification(HTML(paste("Error while computing CPAG<br>", paste(y, collapse="<br>", sep=""), sep="")),
                   #                 type="error", duration=3*notifyDuration)
                   showModal(modalDialog(HTML(paste("Error while computing CPAG<br>", paste(y, collapse="<br>", sep=""), sep="")),
                                         title="iCPAGdb", size="m", easyClose=T, footer=modalButton("OK")))
                 }
               })

        } else {
          output$userComputeResultsTable <- DT::renderDataTable(NULL)
          output$userComputeResultsHeatmap <- renderPlotly(NULL)
          #showNotification("Upload file does not exist.  Please upload and try again.", type="error",
          #                 duration=2*notifyDuration)
          showModal(modalDialog(HTML("Upload file does not exist.  Please upload and try again."),
                                title="iCPAGdb", size="m", easyClose=T, footer=modalButton("OK")))
          cat(paste("User compute, upload file does not exist, file:", f, sep=""), file=stderr())
        }

      } else {
        showModal(modalDialog(HTML(f[["result"]]), title="iCPAGdb", size="m", easyClose=T, footer=modalButton("OK")))
      }

    })

    #######################################################################################################
    # Observe event function:  User compute tab - change in heatmap metric or top significant phenotype
    # pairs to plot
    #######################################################################################################

    observeEvent(c(input$userComputeSelectionHeatmapNphenotype, input$userComputeSelectionHeatmapMetric), {

      cat(paste("nph, n-row-filter:  ", length(userComputeResultsTableRowFilter), "\n", sep=""), file=stderr())

      # Render heatmap based on selected metric
      renderHeatmap(data=currentUserComputeCPAGdata[userComputeResultsTableRowFilter,],
                    xCol="Trait1", yCol="Trait2", valueCol=input$userComputeSelectionHeatmapMetric,
                    log10Value=(input$userComputeSelectionHeatmapMetric %in% c("P_fisher", "Padj_Bonferroni", "Padj_FDR")),
                    nPhenotype=input$userComputeSelectionHeatmapNphenotype,
                    key.title=ifelse(input$userComputeSelectionHeatmapMetric %in% c("P_fisher", "Padj_Bonferroni", "Padj_FDR"),
                                     paste("-log10(", input$userComputeSelectionHeatmapMetric, ")", sep=""),
                                     input$userComputeSelectionHeatmapMetric),
                    tabset="tabsetUserComputeResults", tableTab="tabPanelUserComputeResultsTable",
                    plotTab="tabPanelUserComputeResultsHeatmap", plotName="userComputeResultsHeatmap")

      # Move focus to the plot tab
      updateTabsetPanel(session, "tabsetUserComputeResults", "tabPanelUserComputeResultsHeatmap")

    }, ignoreInit=T)

    #######################################################################################################
    # Observe event function:  User Compute tab - apply filters
    #######################################################################################################

    observeEvent(input$userComputeSelectionFilterApply, {

      # Input values are available since they are retrieved as part of the event triggered here 
      renderResults(data=currentUserComputeCPAGdata,
                    rowFilter="userComputeResultsTableRowFilter",
                    filterTrait=input$userComputeSelectionFilterTrait,
                    filterSNP=input$userComputeSelectionFilterSNP,
                    truncSNP=input$userComputeSelectionIncludeTableAllSNPshare,
                    includeCompoundEFO=input$userComputeSelectionIncludeTableCompoundEFO,
                    efoCol=userComputeResultsEFOcol,
                    filterEFOparent=input$userComputeSelectionFilterEFOparent,
                    efoParentCol=userComputeResultsEFOparentCol,
                    tabset="tabsetUserComputeResults",
                    tableTab="tabPanelUserComputeResultsTable",
                    tableName="userComputeResultsTable",
                    plotTab="tabPanelUserComputeResultsHeatmap",
                    plotName="userComputeResultsHeatmap",
                    heatmapMetric=input$userComputeSelectionHeatmapMetric,
                    heatmapNphenotype=input$userComputeSelectionHeatmapNphenotype)

    }, ignoreInit=T)

    #######################################################################################################
    # Observe event function:  User compute tab - clear filters
    #######################################################################################################

    observeEvent(input$userComputeSelectionFilterClear, {

      # Updates are not effective until all output has been pushed to the client
      # This does not guarantee that new values are in effect at the time of the call to renderResults()
      # Therefore, assign values and pass constants as parameter values to renderResults()
      #updateTextInput(session, "userComputeSelectionFilterTrait", value="")
      #updateTextInput(session, "userComputeSelectionFilterSNP", value="")
      #updateTextInput(session, "userComputeSelectionFilterEFO", value="")
      #updateTextInput(session, "userComputeSelectionFilterEFOparent", value=vector("character"))
      #updateCheckboxInput(session, "userComputeSelectionIncludeTableCompoundEFO", value=F)
      renderResults(data=currentUserComputeCPAGdata,
                    rowFilter="userComputeResultsTableRowFilter",
                    filterTrait="",
                    filterSNP="",
                    truncSNP=input$userComputeSelectionIncludeTableAllSNPshare,
                    includeCompoundEFO=F,
                    efoCol=userComputeResultsEFOcol,
                    filterEFOparent=vector("character"),
                    efoParentCol=userComputeResultsEFOparentCol,
                    tabset="tabsetUserComputeResults",
                    tableTab="tabPanelUserComputeResultsTable",
                    tableName="userComputeResultsTable",
                    plotTab="tabPanelUserComputeResultsHeatmap",
                    plotName="userComputeResultsHeatmap",
                    heatmapMetric=input$userComputeSelectionHeatmapMetric,
                    heatmapNphenotype=input$userComputeSelectionHeatmapNphenotype)

    }, ignoreInit=T)

    #######################################################################################################
    # Configure a download handler for the results table of the userCompute tab
    # This is an efficient alternative to enabling the download buttons of the results data table
    # One drawback is that records must be filtered here to agree with those contained in the table
    # (assuming that filters have been applied to the data table)
    #######################################################################################################

    output$userComputeResultsDownload <-
      downloadHandler(filename="iCPAGdb-UserGWAS-Results.csv",
                      contentType="text/csv",
                      content=function(file) {
                                # Compose filtered observation indices
                                k <- resultsTableComposeFilterIndex(
                                       data=currentUserComputeCPAGdata,
                                       filterTrait=input$userComputeSelectionFilterTrait,
                                       filterSNP=input$userComputeSelectionFilterSNP,
                                       includeCompoundEFO=input$userComputeSelectionIncludeTableCompoundEFO,
                                       efoCol=userComputeResultsEFOcol,
                                       filterEFOparent=input$userComputeSelectionFilterEFOparent,
                                       efoParentCol=userComputeResultsEFOparentCol)
                                write.table(currentUserComputeCPAGdata[k,], file, row.names=F, col.names=T, quote=T, sep=",")
                              }
                     )

  }

)
