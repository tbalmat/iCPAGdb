# Duke University Cross Phenotype Analysis of GWAS Database (iCPAGdb)
# Review, explore, and compute page, November 2020

# Shiny app user interface function

# Information on shiny available at:
# https://shiny.rstudio.com/
# https://github.com/rstudio/shiny
# https://cran.r-project.org/web/packages/shiny/shiny.pdf

options(max.print=1000)      # number of elements, not rows
options(stringsAsFactors=F)
options(scipen=999999)
#options(device="windows")
options(shiny.maxRequestSize=10**9) 

library(shiny)
library(shinyjs)
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
    os <- 2
    setwd(c("C:/Projects/Duke/H2P2GenomeWideAssociationStudy/CPAG/iCPAGdb/App-Devel",
            "/srv/shiny-server/CPAG/explore")[os])
    oscmd <- c("cmd /c", "")[os]
    #oscmd <- c("", "")[os]
    oscmdsep <- c("&&", ";")[os]
    pyexe <- c("\"C:/Users/Kyung Soon/AppData/Local/Programs/Python/Python37/python.exe\"", "python3")[os]
    threads <- 2
    gwidth <- 600
    gheightMin <- 600
    gwidthMax <- 2000
    gheightMax <- 2000
    currentReviewCPAGdata <- data.frame()
    reviewResultsTableRowFilter <- NULL
    reviewResultsEFOcol <- vector("integer")
    reviewResultsEFOparentCol <- vector("integer")
    reviewResultsEFOparentTerm <- vector("character")
    currentExploreCPAGdata <- data.frame()
    exploreResultsTableRowFilter <- NULL
    exploreResultsEFOcol <- vector("integer")
    exploreResultsEFOparentCol <- vector("integer")
    exploreResultsEFOparentTerm <- vector("character")
    largestKnownLogVal <- 347
    pyDir <- "pyCPAG"
    appDir <- "app"
    pyInDir <- "input"
    pyOutDir <- "output"
    pyUserComputeInDir <- "userCompute"
    pyUserComputeOutDir <- "userCompute"
    pyInDirUserCompute <- paste(pyInDir, "/", pyUserComputeInDir, sep="")
    pyOutDirUserCompute <- paste(pyOutDir, "/", pyUserComputeOutDir, sep="")
    inDir <- paste(pyDir, "/", pyInDir, sep="")
    outDir <- paste(pyDir, "/", pyOutDir, sep="")
    inDirUserCompute <- paste(inDir, "/", pyUserComputeInDir, sep="")
    outDirUserCompute <- paste(outDir, "/", pyUserComputeOutDir, sep="")
    logDir <- paste(appDir, "/log", sep="")
    logFile0 <- "iCPAGdbLog.dat"
    logFile <- paste(logDir, "/", logFile0, sep="")
    exploreTabFirst <- T
    currentUserComputeCPAGdata <- data.frame()
    userComputeResultsTableRowFilter <- NULL
    userComputeResultsEFOcol <- vector("integer")
    userComputeResultsEFOparentCol <- vector("integer")
    userComputeResultsEFOparentTerm <- vector("character")
    notifyDuration <- 5
    nLabelTrunc <- 45
    phenH2P2 <- "limited"
    pfactorH2P2 <- 1
    pexpH2P2 <- 5
    pfactorNHGRI <- 5
    pexpNHGRI <- 8
    pfactorOther <- 1
    pexpOther <- 3
    pmaxH2P2 <- pfactorH2P2*10**-pexpH2P2
    pmaxNHGRI <- pfactorNHGRI*10**-pexpNHGRI
    pmaxOther <- pfactorOther*10**-pexpOther
    # Hide the Expore tab
    # Note that, although tab contents are hidden in the ui, the tab itself remains in the tabsetPanel
    # Hide the tab here then assign title (using the reactive output variable defined in the ui)
    hideTab("tabsetCPAG", "tabPanelExplore")
    output$exploreTabTitle <- renderText("")
    # Hide the styles tab
    hideTab("tabsetCPAG", "navBarStyles")
    # Clear any previous upload file name, in case app restarted from browser refresh operation
    reset("userComputeBrowseFile")

    #######################################################################################################
    # Function:  Error handler
    #######################################################################################################

    msgWindow <- function(level="ERROR", title="", msg, size="m", buttonText="OK")
      showModal(modalDialog(HTML(paste("<b>", level, "</b><br><br>", paste(msg, collapse="<br><br>", sep=""), sep="")),
                            title=paste("iCPAGdb ", title, sep=""),
                            size=size, easyClose=T, footer=modalButton(buttonText)))

    #######################################################################################################
    # Function:  Write activity entries to log
    #######################################################################################################

    # Create log file, if it does not exist
    # Delimit columns by tabs since some entries have commas and quotes (OS commands, for instance) 
    tryCatch(
      if(length(which(dir(logDir, logFile0)==logFile0))==0)
        write("sessionID\tdateTime\tsection\toperation\tparameter1\tparameter2\tparameter3\tparameter4\tparameter5\tparameter6\tstatus\tnote", logFile),
      warning=function(err) msgWindow("WARNING", "Access or create log", err),
      error=function(err) msgWindow("ERROR", "Access or create log", err)
    )

    logEntryTrim <- function(a)
      gsub("\"", "\\\"", paste(a[which(nchar(gsub(" ", "", a))>0)], collapse="||", sep=""))

    # Create a session ID
    sessionID <- paste(format(Sys.time(), "%Y%m%d%H%M%S"), "-", paste(sample(0:9, 10, replace=T), collapse=""), sep="")

    writeLog <- function(section="", operation="", parameter1="", parameter2="", parameter3="", parameter4="",
                         parameter5="", parameter6="", status="", note="") {

      tryCatch(
        write(paste(sessionID, "\t", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\t",
                    section, "\t", operation, "\t\"",
                    logEntryTrim(parameter1), "\"\t\"", logEntryTrim(parameter2), "\"\t\"",
                    logEntryTrim(parameter3), "\"\t\"", logEntryTrim(parameter4), "\"\t\"",
                    logEntryTrim(parameter5), "\"\t\"", logEntryTrim(parameter6), "\"\t",
                    status, "\t\"", logEntryTrim(note), "\"", sep=""),
              logFile, append=T),
        warning=function(err) msgWindow("WARNING", "Write to log", err),
        error=function(err) msgWindow("ERROR", "Write to log", err)
      )

    }

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

    composeCPAGtable <- function(data, allSNP=F) {

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
                         "SNPshare_all"=ifelse(rep(allSNP, nrow(data)),
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

    ##################################################################################################
    # Function:  Render CPAG table
    ##################################################################################################

    renderCPAGtable <- function(data, allSNP, tableName) {

      progress <- shiny::Progress$new()
      progress$set(message="Composing CPAG table", value=0.5)
      output[[tableName]] <- DT::renderDataTable(composeCPAGtable(data, allSNP),
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

    }

    #######################################################################################################
    # Function:  Compose observation filter index for review table
    #######################################################################################################

    resultsTableComposeFilterIndex <- function(data, filterTrait, filterSNP, includeCompoundEFO,
                                               efoCol, filterEFOparent, efoParentCol) {

      if(nrow(data)>0) {

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

      } else {
        return(vector("integer"))
      }

    }

    #######################################################################################################
    # Function:  Render data table and heatmap
    #######################################################################################################

    renderResults <- function(data, rowFilter, filterTrait, filterSNP, allSNP, includeCompoundEFO,
                              efoCol, filterEFOparent, efoParentCol, tabset, tableTab, tableName,
                              plotTab, plotName, heatmapMetric, heatmapNphenotype) {

      # Construct filtered observation indices
      resultsTableRowFilter <- resultsTableComposeFilterIndex(data, filterTrait, filterSNP, includeCompoundEFO,
                                                              efoCol, filterEFOparent, efoParentCol)

      # Save filter indices in the environment in which the current function was declared (where variable named in rowFilter
      # was also declared 
      assign(rowFilter, resultsTableRowFilter, envir=environment(renderResults))

      if(length(resultsTableRowFilter)>0) {

        # Display CPAG table
        renderCPAGtable(data[resultsTableRowFilter,], allSNP, tableName)

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
        showModal(modalDialog(HTML("No CPAG data exist for specified filter values"),
                              title="iCPAGdb Render Results", size="m", easyClose=T, footer=modalButton("OK")))
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
        showModal(modalDialog(HTML("Insufficient data (cross-GWAS phenotype pairs) to generate heatmap"),
                  "title"="iCPAGdb Render Heatmap", size="m", easyClose=T, footer=modalButton("OK")))
      }

    }

    #######################################################################################################
    # Observe event:  Review tab - selection table cell click
    # Note that input$reviewSelectionTable_cells_selected is a variable maintained by the DT package, but
    # is made reactive by Shiny, so that changes to the DT variable cause events here
    #######################################################################################################

    observeEvent(input$reviewSelectionTable_cells_selected, {

      # This event is executed when a cell in the review table is clicked, but not when a hyperlink
      # within a cell is clicked (hyperlinks are followed by the browser)

      tryCatch(

        if(length(input$reviewSelectionTable_cell_clicked)>0) {

          # Display precomputed results corresponding to selected row of review table
          # Note that _cell_clicked returns a list with elements:
          # row = row in data frame corresponding to cell clicked (reordering rows maintains row IDs)
          # col = 0-based col of cell clicked
          # value = data frame cell contents
          # Note that if selection mode of multiple is specified in renderTable(), then multiple cells
          # will be reported by _cell_clicked, all at the top list level
          # Referencing list elements by name returns the first encountered, ignoring multiples
          if(substring(as.character(input$reviewSelectionTable_cell_clicked[["value"]]), 1, 7)!="<a href") {

            # Get selected row number
            trow <- input$reviewSelectionTable_cell_clicked[["row"]][1]

            # Retrieve file name of specified precomputed results corresponding to selected row number
            fn <- reviewCPAG[trow,"file"]
            cat(paste("fileName:  ", fn, "\n", sep=""), file=stderr())

            # Proceed if file exists
            if(length(which(dir(outDir, pattern=paste("^", fn, "$", sep=""))==fn))==1) {

              # Retrieve CPAG results
              # Save in global data frame to be available in other reactive functions
              progress <- shiny::Progress$new()
              progress$set(message="Reading CPAG data", value=0.5)
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
                            allSNP=input$reviewSelectionIncludeTableAllSNPshare,
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

              # Log GWAS selection
              writeLog(section="Review", operation="selectGWAS",
                       parameter1=reviewCPAG[trow,"Description"],
                       parameter2=reviewCPAG[trow,"GWAS_1_source"],
                       parameter3=reviewCPAG[trow,"GWAS_2_source"],
                       parameter4=reviewCPAG[trow,"p_threshold_GWAS_1"],
                       parameter5=reviewCPAG[trow,"p_threshold_GWAS_2"],
                       parameter6=reviewCPAG[trow,"LD_population"],
                       status="success")

            } else {

              output$reviewResultsTable <- DT::renderDataTable(NULL)
              output$reviewResultsHeatmap <- renderPlotly(NULL)
              showModal(modalDialog(HTML(paste("File <b>", reviewCPAG[trow,"file"], "</b> not found", sep="")),
                        title="iCPAGdb Review", size="m", easyClose=T, footer=modalButton("OK")))
              # Log error
              writeLog(section="Review", operation="selectGWAS", parameter1=reviewCPAG[trow,"file"],
                       status="error", note="File does not exist")

            }
          }
        },

        warning=function(err) msgWindow("WARNING", "Review table row selection", err),
        error=function(err) msgWindow("ERROR", "Review table row selection", err)

      )

    }, ignoreInit=T)

    #######################################################################################################
    # Observe event function:  Review tab - apply filters
    #######################################################################################################

    observeEvent(input$reviewSelectionFilterApply, {

      tryCatch({
        renderResults(data=currentReviewCPAGdata,
                      rowFilter="reviewResultsTableRowFilter",
                      filterTrait=input$reviewSelectionFilterTrait,
                      filterSNP=input$reviewSelectionFilterSNP,
                      allSNP=input$reviewSelectionIncludeTableAllSNPshare,
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
        # Log action
        writeLog(section="Review", operation="filter", status="success")},
        warning=function(err) msgWindow("WARNING", "Review filter", err),
        error=function(err) msgWindow("ERROR", "Review filter", err)
      )

    }, ignoreInit=T)

    #######################################################################################################
    # Observe event function:  Review tab - clear filters
    #######################################################################################################

    observeEvent(input$reviewSelectionFilterClear, {

      # Updates are not effective until all output has been pushed to the client
      # This does not guarantee that new values are in effect at the time of the call to renderResults()
      # Therefore, assign values and pass constants as parameter values to renderResults()
      tryCatch({
        updateTextInput(session, "reviewSelectionFilterTrait", value="")
        updateTextInput(session, "reviewSelectionFilterSNP", value="")
        #updateTextInput(session, "reviewSelectionFilterEFO", value="")
        updateTextInput(session, "reviewSelectionFilterEFOparent", value=vector("character"))
        updateCheckboxInput(session, "reviewSelectionIncludeTableCompoundEFO", value=F)
        renderResults(data=currentReviewCPAGdata,
                      rowFilter="reviewResultsTableRowFilter",
                      filterTrait="",
                      filterSNP="",
                      allSNP=input$reviewSelectionIncludeTableAllSNPshare,
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
                      heatmapNphenotype=input$reviewSelectionHeatmapNphenotype)},
        warning=function(err) msgWindow("WARNING", "Review filter clear", err),
        error=function(err) msgWindow("ERROR", "Review filter clear", err)
      )

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
                                tryCatch({
                                  # Compose filtered observation indices
                                  k <- resultsTableComposeFilterIndex(
                                         data=currentReviewCPAGdata,
                                         filterTrait=input$reviewSelectionFilterTrait,
                                         filterSNP=input$reviewSelectionFilterSNP,
                                         includeCompoundEFO=input$reviewSelectionIncludeTableCompoundEFO,
                                         efoCol=reviewResultsEFOcol,
                                         filterEFOparent=input$reviewSelectionFilterEFOparent,
                                         efoParentCol=reviewResultsEFOparentCol)
                                  # Log action
                                  writeLog(section="Review", operation="download", parameter1="results", status="success")
                                  write.table(currentReviewCPAGdata[k,], file, row.names=F, col.names=T, quote=T, sep=",")},
                                  warning=function(err) msgWindow("WARNING", "Review download results", err),
                                  error=function(err) msgWindow("ERROR", "Review download results", err)
                                )
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
    #######################################################################################################

    observeEvent(input$reviewSelectionIncludeTableAllSNPshare, {

      tryCatch({
        renderCPAGtable(currentReviewCPAGdata[reviewResultsTableRowFilter,],
                        input$reviewSelectionIncludeTableAllSNPshare,
                        "reviewResultsTable")
        # Log action
        writeLog(section="Review", operation="allSNP", status="")},
        warning=function(err) msgWindow("WARNING", "Review all SNPs", err),
        error=function(err) msgWindow("ERROR", "Review all SNPs", err)
      )

    }, ignoreInit=T)

    #######################################################################################################
    # Observe event function:  Review tab - change in heatmap metric or top significant phenotype
    # pairs to plot
    #######################################################################################################

    observeEvent(c(input$reviewSelectionHeatmapNphenotype, input$reviewSelectionHeatmapMetric), {

      tryCatch({
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
        # Log action
        writeLog(section="Review", operation="heatmapMetricTopAdjust", status="success")},
        warning=function(err) msgWindow("WARNING", "Review heatmap n, metric adjust", err),
        error=function(err) msgWindow("ERROR", "Review heatmap n, metric adjust", err)
      )

    }, ignoreInit=T)

    #######################################################################################################
    # Observe event:  Enable/disable "secret" features
    #######################################################################################################

    observeEvent(input$featureEnable, {

      tryCatch({
        if(tolower(input$featureEnable)=="exp+") {
          showTab("tabsetCPAG", "tabPanelExplore")
          # Assign visible title
          output$exploreTabTitle <- renderText("Explore iCPAGdb associations")
        } else if(tolower(input$featureEnable)=="exp-") {
          hideTab("tabsetCPAG", "tabPanelExplore")
        } else if(tolower(input$featureEnable)=="phen+") {
          phenH2P2 <<- "full"
        } else if(tolower(input$featureEnable)=="phen-") {
          phenH2P2 <<- "limited"
        }
        # Log action
        if(input$featureEnable!="")
          writeLog(section="", operation="featureEnable", parameter1=input$featureEnable, status="success")},
        warning=function(err) msgWindow("WARNING", "Feature enable", err),
        error=function(err) msgWindow("ERROR", "Feature enable", err)
      )

    # Do not ignore, so that initial environment consistent with any values in secret field
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
                    allSNP=input$exploreSelectionIncludeTableAllSNPshare,
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
    # Function:  Verify parameter values on Explore tab
    # Present a dialog window when:
    # Sources are identical but p-thresholds are different
    # Any p-threshold is above the max allowed for the corresponding source
    # Prompt for auto-correction or cancel (for manual correction)
    # Auto-correction update then re-attempts computation (follow reactive events)
    # Use parameter variables for evaluation, instead of direct reference to reactive values, since
    # some are updated from within the script and the new values are not reflected in the
    # eactive variables until data are sent, or flushed, to the browser, which may be after
    # further script instructions reference the variables
    #######################################################################################################

    exploreVerifyParameterValues <- function(source1, source2, pfactor1, pexp1, pfactor2, pexp2) {

      msg <- ""

      # Compare p-thresholds with max allowed values
      # Evaluate equivalence of sources and p-threshold values
      p1 <- as.numeric(pfactor1)*10**-pexp1
      p2 <- as.numeric(pfactor2)*10**-pexp2
      if(source1=="H2P2" & p1>pmaxH2P2 |
         source1=="NHGRI" & p1>pmaxNHGRI |
         source1!="H2P2" & source1!="NHGRI" & p1>pmaxOther) {
        msg <- "p-threshold<sub>1</sub> above maximum allowed value<br><br>"
        buttonID <- "exploreComputeSetPMax1"
        buttonText <- "Use maximum allowed value"
      } else if(source2=="H2P2" & p2>pmaxH2P2 |
                source2=="NHGRI" & p2>pmaxNHGRI |
                source2!="H2P2" & source2!="NHGRI" & p2>pmaxOther) {
        msg <- "p-threshold<sub>2</sub> above maximum allowed value<br><br>"
        buttonID <- "exploreComputeSetPMax2"
        buttonText <- "Use maximum allowed value"
      } else if(source1==source2 & p1!=p2) {
        msg <- "GWAS sources are identical, but p-thresholds are different<br><br>"
        buttonID <- "exploreComputeSetPsmall"
        buttonText <- "Use lesser p-threshold for both sets"
      }

      if(msg=="") {
        return(list("success"=T))
      } else {
        return(list("success"=F, "msg"=msg, "buttonID"=buttonID, "buttonText"=buttonText))
      }

      return(success)
    }

    #######################################################################################################
    # Function:  Explore tab compute CPAG
    # Use parameter variables for evaluation, instead of direct reference to reactive values, since
    # some are updated from within the script and the new values are not reflected in the
    # reactive variables until data are sent, or flushed, to the browser, which may be after
    # further script instructions reference the variables
    #######################################################################################################

    exploreComputeCPAG <- function(source1, source2, pfactor1, pexp1, pfactor2, pexp2, ldpop) {

      q <- tryCatch({

        pycmd <- ""
        pycmd2 <- ""
        pycmd3 <- ""

        # Verify CPAG computation parameter values
        parV <- exploreVerifyParameterValues(source1, source2, pfactor1, pexp1, pfactor2, pexp2)
        if(parV[["success"]]) {

          # Identify which GWAS sources are in request 
          s <- c(source1, source2)
          kh <- which(s=="H2P2")
          kn <- which(s=="NHGRI")
          ko <- which(!s %in% c("H2P2", "NHGRI"))
          # Use lesser of p values when sources are identical or when neither set is H2P2 nor NHGRI
          # Note that the interface allows specification of different p-thresholds with two non(H2P2 or NHRI)
          # sets, but the CPAG function has a single parameter for non(H2P2 or NHRI) data sets
          pfac <- c(pfactor1, pfactor2)
          px <- c(pexp1, pexp2)
          if(s[1]==s[2] | !s[1] %in% c("H2P2", "NHGRI") & !s[2] %in% c("H2P2", "NHGRI"))
            if(as.integer(pfac[2])*10**-pexp2<as.integer(pfac[1])*10**-pexp1) {
              pfac[1] <- pfac[2]
              px <- c(pexp2, pexp2)
            } else {
              pfac[2] <- pfac[1]
              px <- c(pexp1, pexp1)
            }

          # Compose output file name
          # Format is GWAS1-p1-GWAS2-p2-ldpop.csv
          # Order of appearance (GWAS1, GWAS2) is H2P2, NHGRI, other
          # GWAS2 and p2 are excluded if GWAS1=GWAS2
          if(length(kh)>0) {
            cpag <- paste("H2P2-p", pfac[kh[1]], "e-", sprintf("%02d", px[kh[1]]),
                          # Include second GWAS only if different from H2P2
                          ifelse(length(kh)==1,
                                 paste("-", s[-kh], "-p", pfac[-kh], "e-", sprintf("%02d", px[-kh]), sep=""),
                                 ""),
                          sep="")
          } else if(length(kn)>0) {
            cpag <- paste("NHGRI-p", pfac[kn[1]], "e-", sprintf("%02d", px[kn[1]]),
                          # Include second GWAS only if different from NHGRI
                          ifelse(length(kn)==1,
                                 paste("-", s[-kn], "-p", pfac[-kn], "e-", sprintf("%02d", px[-kn]), sep=""),
                                 ""),
                          sep="")
          } else if(s[1]!=s[2]) {
            k <- order(s)
            cpag <- paste(s[k[1]], "-p", pfac[k[1]], "e-", sprintf("%02d", px[k[1]]), "-",
                          s[k[2]], "-p", pfac[k[2]], "e-", sprintf("%02d", px[k[2]]), sep="")
          } else {
            cpag <- paste(s[1], "-p", pfac[1], "e-", sprintf("%02d", px[1]), sep="")
          }
          cpag <- paste(cpag, "-", ldpop, ".csv", sep="")
          existFile <- F

          # Test for existence of output file
          if(length(which(dir(outDir, pattern=cpag, ignore.case=T)==cpag))>0) {

            # File exists - read data, render table and heatmap
            exploreReadCPAGandRender(paste(outDir, "/", cpag, sep=""))
            # Log action
            writeLog(section="Explore", operation="computeCPAG",
            parameter1=cpag, note="CPAG results exist - no computation required", status="success")

          } else {

            # File does not exist
            # Compute results for requested GWAS sets, p-thresholds, and LD population
            # Compose python command to execute CPAG function
            pycmd <- paste(oscmd, " cd ", pyDir, " ", oscmdsep, " ",
                           pyexe, " main.py cpagdb --threads ", threads, 
                           paste(paste(" --subtype ", unique(sort(c(s[kh], s[kn], s[ko]))), sep=""), collapse=""),
                           ifelse(s[1]=="H2P2" | s[2]=="H2P2",
                                  paste(" --H2P2-Pcut ", pfac[kh[1]], "e-", px[kh[1]], sep=""), ""),
                           ifelse(s[1]=="NHGRI" | s[2]=="NHGRI",
                                  paste(" --NHGRI-Pcut ", pfac[kn[1]], "e-", px[kn[1]], sep=""), ""),
                           # Note that the lesser of p-thresholds is used when both sets non-H2P2-NHGRI
                           ifelse(!s[1] %in% c("H2P2", "NHGRI") | !s[2] %in% c("H2P2", "NHGRI"),
                                  paste(" --Pcut ", pfac[ko[1]], "e-", px[ko[1]], sep=""), ""),
                           " --lddb-pop ", ldpop,
                           " --outfile \"", pyOutDir, "/", cpag, "\"", sep="")

            # Compose commands to include ontology, for NHGRI results only
            if(length(kh)==0 & length(kn)>0) {
              pycmd2 <- paste(oscmd, " cd ", pyDir, " ", oscmdsep, " ",
                              pyexe, " main.py post_analysis --anno-ontology --anno-cols Trait1", " --infile \"",
                              pyOutDir, "/", cpag, "\" --outfile \"", pyOutDir, "/", cpag, "\"", sep="")
            } else {
              pycmd2 <- ""
            }
            if(length(kh)==1 & length(kn)==1 | length(kn)==2) {
              pycmd3 <- paste(oscmd, " cd ", pyDir, " ", oscmdsep, " ",
                              pyexe, " main.py post_analysis --anno-ontology --anno-cols Trait2", " --infile \"",
                              pyOutDir, "/", cpag, "\" --outfile \"", pyOutDir, "/", cpag, "\"", sep="")
            } else {
              pycmd3 <- ""
            }

            cat(paste("\npycmd:  ", pycmd, "\n\n", sep=""), file=stderr())
            cat(paste("\npycmd2:  ", pycmd2, "\n\n", sep=""), file=stderr())
            cat(paste("\npycmd3:  ", pycmd3, "\n\n", sep=""), file=stderr())

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
            # Using asynchronous operations, lengthy instructions can be executed without delaying other activity
            # The %...>% pipe delays processing of further instructions while the output of lengthy asynchrounous
            # instructions are executed
            # Note that errors during the system() call raise errors here, with the system message(s) as text
            # If output reaches the pipe, noe system() error occurred
            future({
              t0 <- proc.time()
              cpagPhase <- 1
              q0 <- suppressWarnings(system(pycmd, intern=T))
              if(is.null(attr(q0, "status")))
                attr(q0, "status") <- 0
              if(attr(q0, "status")==0  & pycmd2!="") {
                cpagPhase <- 2
                q0 <- suppressWarnings(system(pycmd2, intern=T))
                if(is.null(attr(q0, "status")))
                  attr(q0, "status") <- 0
              }
              if(attr(q0, "status")==0  & pycmd3!="") {
                cpagPhase <- 3
                q0 <- suppressWarnings(system(pycmd3, intern=T))
                if(is.null(attr(q0, "status")))
                  attr(q0, "status") <- 0
              }
              t0 <- (proc.time()-t0)["elapsed"]
              list("cpagPhase"=cpagPhase, "q0"=q0, "t0"=t0)
            }) %...>%
            (function(q) {
               progress$close()
               cat(paste("Execution time: ", q[["t0"]], "\n", collapse=" "), file=stderr())
               if(attr(q[["q0"]], "status")==0) {

                 # Verify existence of output file
                 if(length(which(dir(outDir, pattern=cpag)==cpag))>0) {
                   exploreReadCPAGandRender(paste(outDir, "/", cpag, sep=""))

                   # Clear filters
                   updateTextInput(session, "exploreSelectionFilterTrait", value="")
                   updateTextInput(session, "exploreSelectionFilterSNP", value="")
                   #updateTextInput(session, "exploreSelectionFilterEFO", value="")
                   updateTextInput(session, "exploreSelectionFilterEFOparent", value=vector("character"))
                   updateCheckboxInput(session, "exploreSelectionIncludeTableCompoundEFO", value=F)

                   # Log success
                   writeLog(section="Explore", operation="computeCPAG",
                            parameter1=pycmd, parameter2=pycmd2, parameter3=pycmd3, status="success")

                 } else {
                   # Treat error here since the outer function (where future environment spawned)
                   # continues processing while the future/promise sequence executed independently
                   writeLog(section="Explore", operation="computeCPAG", parameter1=q[["q0"]],
                            parameter2=ifelse(q[["cpagPhase"]]==1, pycmd, ifelse(q[["cpagPhase"]]==2, pycmd2, pycmd3)),
                            parameter3=q[["q0"]], parameter4=q[["t0"]], note="CPAG output file not found", status="error")
                   msgWindow("ERROR", "Explore compute CPAG",
                             c(paste("Computed CPAG file ", outDirUserCompute, "/", cpagFile, " does not exist", sep=""),
                               "Occurred during execution of:",
                               parameter2=ifelse(q[["cpagPhase"]]==1, pycmd, ifelse(q[["cpagPhase"]]==2, pycmd2, pycmd3)),
                               "<b>CPAG output follows</b>", q[["q0"]]))
                 }

               } else {

                 # Treat error here since the outer function continues processing while future/promise executes independently
                 writeLog(section="Explore", operation="computeCPAG", parameter1=q[["q0"]],
                          parameter2=ifelse(q[["cpagPhase"]]==1, pycmd, ifelse(q[["cpagPhase"]]==2, pycmd2, pycmd3)),
                          parameter3=q[["t0"]], note="Error during CPAG computation", status="error")
                 msgWindow("ERROR", "Explore compute CPAG", c(q[["q0"]], "Occurred during execution of:",
                           parameter2=ifelse(q[["cpagPhase"]]==1, pycmd, ifelse(q[["cpagPhase"]]==2, pycmd2, pycmd3))))
               }

            })

          }

        } else {

          # Open a dialog, window if any rules violated
          # Present options for correction or canceling
          # Follow the reactive events triggered from within the dialog for path to compute function
          showModal(modalDialog(HTML(parV[["msg"]]),
                                HTML("Choose from the following:"),
                                div(
                                  actionButton(parV[["buttonID"]], parV[["buttonText"]]),
                                  style="display:inline-block; vertical-align:top; margin-top:-5px; margin-left:20px"
                                ),
                                div(
                                  modalButton("Cancel"),
                                  style="display:inline-block; vertical-align:top; margin-top:-5px; margin-left:10px"
                                ),
                                title="iCPAGdb Explore", size="m", easyClose=T,
                                footer=NULL))

            # Return a non-error class object
            ""

        }},
        warning=function(err) err,
        error=function(err) err
      )

      # Return message and CPAG commands, if any errors
      if(!any(class(q) %in% c("warning", "error"))) {
        "success"
      } else {
        errorCondition(c(q[["message"]], pycmd, pycmd2, pycmd3))
      }

    }

    #######################################################################################################
    # Observe event function:  Explore tab - set both p-thresholds to the lesser of the two
    # This event is the result of pressing a "modify values" button in a modal dialog window
    # Reactive values (screen controls) are programmatically updated from within this function, but
    # the corresponding variables (input$exploreX) do not reflect the new values until a flush to
    # the browser is accomplished, which is prompted by user action (press a button, change an on-screen
    # value)
    # Therefore pass updated values to compute and data verification in function parameters
    #######################################################################################################

    observeEvent(input$exploreComputeSetPsmall, {

      p1 <- as.numeric(input$explorePfactor1)*10**-input$explorePexp1
      p2 <- as.numeric(input$explorePfactor2)*10**-input$explorePexp2
      if(p1<p2) {
        pfactor1 <- input$explorePfactor1
        pexp1 <- input$explorePexp1
        pfactor2 <- pfactor1
        pexp2 <- pexp1
        updateRadioButtons(session, "explorePfactor2", selected=pfactor2)
        updateSliderInput(session, "explorePexp2", value=pexp2)
      } else {
        pfactor2 <- input$explorePfactor2
        pexp2 <- input$explorePexp2
        pfactor1 <- pfactor2
        pexp1 <- pexp2
        updateRadioButtons(session, "explorePfactor1", selected=pfactor1)
        updateSliderInput(session, "explorePexp1", value=pexp1)
      }
      removeModal()
      q <- tryCatch(
        exploreComputeCPAG(input$exploreSource1, input$exploreSource2, pfactor1, pexp1,
                           pfactor2, pexp2, input$exploreLDpop),
        warning=function(err) err,
        error=function(err) err
      )
      if(any(class(q) %in% c("warning", "error"))) {
        output$exploreResultsTable <- DT::renderDataTable(NULL)
        output$exploreResultsHeatmap <- renderPlotly(NULL)
        writeLog(section="Explore", operation="computeCPAG", note=q[["message"]], status="error")
        msgWindow("ERROR", "Explore compute CPAG setPsmall", q[["message"]])
      }

    }, ignoreInit=T)

    #######################################################################################################
    # Observe event function:  Explore tab - set p-threshold 1 to max allowable value
    # Use local variables and pass values as parameters due to reactive variables not reflecting
    # updated values until flush
    #######################################################################################################

    observeEvent(input$exploreComputeSetPMax1, {

      if(input$exploreSource1=="H2P2") {
        pfactor1 <- pfactorH2P2
        pexp1 <- pexpH2P2
      } else if(input$exploreSource1=="NHGRI") {
        pfactor1 <- pfactorNHGRI
        pexp1 <- pexpNHGRI
      } else {
        pfactor1 <- pfactorOther
        pexp1 <- pexpOther
      }
      updateRadioButtons(session, "explorePfactor1", selected=pfactor1)
      updateSliderInput(session, "explorePexp1", value=pexp1)
      removeModal()
      q <- tryCatch(
        exploreComputeCPAG(input$exploreSource1, input$exploreSource2, pfactor1=pfactor1, pexp1=pexp1,
                           input$explorePfactor2, input$explorePexp2, input$exploreLDpop),
        warning=function(err) err,
        error=function(err) err
      )
      if(any(class(q) %in% c("warning", "error"))) {
        output$exploreResultsTable <- DT::renderDataTable(NULL)
        output$exploreResultsHeatmap <- renderPlotly(NULL)
        writeLog(section="Explore", operation="computeCPAG", note=q[["message"]], status="error")
        msgWindow("ERROR", "Explore compute CPAG setMaxP1", q[["message"]])
      }

    }, ignoreInit=T)

    #######################################################################################################
    # Observe event function:  Explore tab - set p-threshold 2 to max allowable value
    # Use local variables and pass values as parameters due to reactive variables not reflecting
    # updated values until flush
    #######################################################################################################

    observeEvent(input$exploreComputeSetPMax2, {

      if(input$exploreSource2=="H2P2") {
        pfactor2 <- pfactorH2P2
        pexp2 <- pexpH2P2
      } else if(input$exploreSource2=="NHGRI") {
        pfactor2 <- pfactorNHGRI
        pexp2 <- pexpNHGRI
      } else {
        pfactor2 <- pfactorOther
        pexp2 <- pexpOther
      }
      updateRadioButtons(session, "explorePfactor2", selected=pfactor2)
      updateSliderInput(session, "explorePexp2", value=pexp2)
      removeModal()
      q <- tryCatch(
        exploreComputeCPAG(input$exploreSource1, input$exploreSource2, input$explorePfactor1, input$explorePexp1,
                           pfactor2, pexp2, input$exploreLDpop),
        warning=function(err) err,
        error=function(err) err
      )
      if(any(class(q) %in% c("warning", "error"))) {
        output$exploreResultsTable <- DT::renderDataTable(NULL)
        output$exploreResultsHeatmap <- renderPlotly(NULL)
        writeLog(section="Explore", operation="computeCPAG", note=q[["message"]], status="error")
        msgWindow("ERROR", "Explore compute CPAG setMaxP2", q[["message"]])
      }

    }, ignoreInit=T)

    #######################################################################################################
    # Observe event function:  Explore tab - compute CPAG button
    #######################################################################################################

    observeEvent(input$exploreCompute, {

      # Log request
      writeLog(section="Explore", operation="computeReq",
               parameter1=input$exploreSource1, parameter2=input$exploreSource2,
               parameter3=paste(input$explorePfactor1, "e-", sprintf("%02d", input$explorePexp1), sep=""),
               parameter4=paste(input$explorePfactor2, "e-", sprintf("%02d", input$explorePexp2), sep=""),
               parameter5=input$exploreLDpop, status="")

      q <- tryCatch(
        exploreComputeCPAG(input$exploreSource1, input$exploreSource2, input$explorePfactor1,
                           input$explorePexp1, input$explorePfactor2, input$explorePexp2,
                           input$exploreLDpop),
        warning=function(err) err,
        error=function(err) err
      )
      if(any(class(q) %in% c("warning", "error"))) {
        output$exploreResultsTable <- DT::renderDataTable(NULL)
        output$exploreResultsHeatmap <- renderPlotly(NULL)
        writeLog(section="Explore", operation="computeCPAG", note=q[["message"]], status="error")
        msgWindow("ERROR", "Explore compute CPAG", q[["message"]])
      }

    }, ignoreInit=T)

    #######################################################################################################
    # Observe event function:  Explore tab - change in include all SNPs in shared SNP option
    #######################################################################################################

    observeEvent(input$exploreSelectionIncludeTableAllSNPshare, {

      tryCatch({
        renderCPAGtable(currentExploreCPAGdata[exploreResultsTableRowFilter,],
                        input$exploreSelectionIncludeTableAllSNPshare,
                        "exploreResultsTable")
        # Log action
        writeLog(section="Explore", operation="allSNP", status="")},
        warning=function(err) msgWindow("WARNING", "Explore all SNPs", err),
        error=function(err) msgWindow("ERROR", "Explore all SNPs", err)
      )

    }, ignoreInit=T)

    #######################################################################################################
    # Observe event function:  Explore tab - change in heatmap metric or top significant phenotype
    # pairs to plot
    #######################################################################################################

    observeEvent(c(input$exploreSelectionHeatmapNphenotype, input$exploreSelectionHeatmapMetric), {

      tryCatch({
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
        # Log action
        writeLog(section="Explore", operation="heatmapMetricTopAdjust", status="success")},
        warning=function(err) msgWindow("WARNING", "Explore heatmap n, metric adjust", err),
        error=function(err) msgWindow("ERROR", "Explore heatmap n, metric adjust", err)
      )

    }, ignoreInit=T)

    #######################################################################################################
    # Observe event function:  Explore tab - apply filters
    #######################################################################################################

    observeEvent(input$exploreSelectionFilterApply, {

      tryCatch({
        renderResults(data=currentExploreCPAGdata,
                      rowFilter="exploreResultsTableRowFilter",
                      filterTrait=input$exploreSelectionFilterTrait,
                      filterSNP=input$exploreSelectionFilterSNP,
                      allSNP=input$exploreSelectionIncludeTableAllSNPshare,
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
        # Log action
        writeLog(section="Explore", operation="filter", status="success")},
        warning=function(err) msgWindow("WARNING", "Explore filter", err),
        error=function(err) msgWindow("ERROR", "Explore filter", err)
      )
 
    }, ignoreInit=T)

    #######################################################################################################
    # Observe event function:  Explore tab - clear filters
    #######################################################################################################

    observeEvent(input$exploreSelectionFilterClear, {

      # Updates are not effective until all output has been pushed to the client
      # This does not guarantee that new values are in effect at the time of the call to renderResults()
      # Therefore, assign values and pass constants as parameter values to renderResults()
      tryCatch({
        updateTextInput(session, "exploreSelectionFilterTrait", value="")
        updateTextInput(session, "exploreSelectionFilterSNP", value="")
        #updateTextInput(session, "exploreSelectionFilterEFO", value="")
        updateTextInput(session, "exploreSelectionFilterEFOparent", value=vector("character"))
        updateCheckboxInput(session, "exploreSelectionIncludeTableCompoundEFO", value=F)
        renderResults(data=currentExploreCPAGdata,
                      rowFilter="exploreResultsTableRowFilter",
                      filterTrait="",
                      filterSNP="",
                      allSNP=input$exploreSelectionIncludeTableAllSNPshare,
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
                      heatmapNphenotype=input$exploreSelectionHeatmapNphenotype)},
        warning=function(err) msgWindow("WARNING", "Explore filter clear", err),
        error=function(err) msgWindow("ERROR", "Explore filter clear", err)
      )

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
                                tryCatch({
                                  # Compose filtered observation indices
                                  k <- resultsTableComposeFilterIndex(
                                         data=currentExploreCPAGdata,
                                         filterTrait=input$exploreSelectionFilterTrait,
                                         filterSNP=input$exploreSelectionFilterSNP,
                                         includeCompoundEFO=input$exploreSelectionIncludeTableCompoundEFO,
                                         efoCol=exploreResultsEFOcol,
                                         filterEFOparent=input$exploreSelectionFilterEFOparent,
                                         efoParentCol=exploreResultsEFOparentCol)
                                  # Log action
                                  writeLog(section="Explore", operation="download", parameter1="results", status="success")
                                  write.table(currentExploreCPAGdata[k,], file, row.names=F, col.names=T, quote=T, sep=",")},
                                  warning=function(err) msgWindow("WARNING", "Explore download results", err),
                                  error=function(err) msgWindow("ERROR", "Explore download results", err)
                                )
                              }
                     )

    #######################################################################################################
    # Observe event function:  User Compute tab - browse file selected
    #######################################################################################################

    observeEvent(input$userComputeBrowseFile, {

      tryCatch({
        # Read first lines of file
        x <- scan(input$userComputeBrowseFile[,"datapath"], "character", n=5, sep="\n", quote="", quiet=T)
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
          title="iCPAGdb Compute", size="l", easyClose=T,
          footer=modalButton("OK")))
        # Log action
        writeLog(section="Compute", operation="uploadFile",
                 parameter1=input$userComputeBrowseFile[,"name"],
                 parameter2=input$userComputeBrowseFile[,"size"],
                 parameter3=input$userComputeBrowseFile[,"datapath"],
                 status="success")},
        warning=function(err) msgWindow("WARNING", "User Compute browse files", err),
        error=function(err) msgWindow("ERROR", "User Compute browse files", err)
      )

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
                # Log action
                writeLog(section="Compute", operation="uploadCopy",
                         parameter1=input$userComputeBrowseFile[,"name"],
                         parameter2=paste(uGWASfile, ".csv", sep=""),
                         parameter3=nrow(userFile),
                         parameter4=input$userComputeSNPcol,
                         parameter5=input$userComputePcol,
                         parameter6=input$userComputeDelimiter,
                         status="success")
              } else {
                cat(paste("Cannot create user upload file for: ", input$userComputeBrowseFile[,"name"], "\n", sep=""), file=stderr())
                result <- "Cannot process upload file.  Please upload again."
                writeLog(section="Compute", operation="uploadCopy",
                         parameter1=input$userComputeBrowseFile[,"name"],
                         parameter2="Cannot copy upload file",
                         status="error")
              }

            } else {
              result <- paste("Specified P (significance) column is missing in uploaded GWAS file<br>",
                              "Input columns are: ", paste(colnames(userFile), collapse=", "),
                              "<br>Be sure to verify delimiter style", sep="")
              writeLog(section="Compute", operation="uploadCopy",
                       parameter1=input$userComputeBrowseFile[,"name"],
                       parameter2=input$userComputePcol,
                       parameter3=input$userComputeDelimiter,
                       note="Specified P column missing",
                       status="error")
            }
          } else {
            result <- paste("Specified SNP column is missing in uploaded GWAS file<br>",
                            "Input columns are: ", paste(colnames(userFile), collapse=", "),
                            "<br>Be sure to verify delimiter style", sep="")
            writeLog(section="Compute", operation="uploadCopy",
                     parameter1=input$userComputeBrowseFile[,"name"],
                     parameter2=input$userComputeSNPcol,
                     parameter3=input$userComputeDelimiter,
                     note="Specified SNP column missing",
                     status="error")
          }
        #} else {
        #  result <- paste("Specified phenotype (trait) column is missing in uploaded GWAS file<br> ",
        #                  "Input columns are: ", paste(colnames(userFile), collapse=", "),
        #                  "<br>Be sure to verify delimiter style", sep="")
        #  writeLog(section="Compute", operation="uploadCopy",
        #           parameter1=input$userComputeBrowseFile[,"name"],
        #           parameter2=input$userComputePhenotypeCol,
        #           parameter=input$userComputeDelimiter,
        #           note="Specified phenotype column missing",
        #           status="error")
        #}
      } else {
        result <- "No file specified.  Please upload a file."
        writeLog(section="Compute", operation="uploadCopy",
                 parameter1=input$userComputeBrowseFile[,"name"],
                 note="File name missing in upload specification",
                 status="error")
      }

      return(list("result"=result, "uGWASfile"=uGWASfile, "uGWAScol"=uGWAScol))

    }

    #######################################################################################################
    # Function:  Verify parameter values on User Compute tab
    # Present a dialog window when p-threshold is above the max allowed for the corresponding source
    # Prompt for auto-correction or cancel (for manual correction)
    # Auto-correction update then re-attempts computation (follow reactive events)
    # Use parameter variables for evaluation, instead of direct reference to reactive values, since
    # some are updated from within the script and the new values are not reflected in the
    # eactive variables until data are sent, or flushed, to the browser, which may be after
    # further script instructions reference the variables
    #######################################################################################################

    userComputeVerifyParameterValues <- function(source2, pfactor2, pexp2) {

      msg <- ""

      # Compare p-thresholds with max allowed values
      p2 <- as.numeric(pfactor2)*10**-pexp2
      if(source2=="H2P2" & p2>pmaxH2P2 |
         source2=="NHGRI" & p2>pmaxNHGRI |
         source2!="H2P2" & source2!="NHGRI" & p2>pmaxOther) {
        msg <- "p-threshold<sub>2</sub> above maximum allowed value<br><br>"
        buttonID <- "userComputeComputeSetPMax2"
        buttonText <- "Use maximum allowed value"
      }

      if(msg=="") {
        return(list("success"=T))
      } else {
        return(list("success"=F, "msg"=msg, "buttonID"=buttonID, "buttonText"=buttonText))
      }

    }

    #######################################################################################################
    # Function:  User Compute tab - compute CPAG
    # Use parameter variables for evaluation, instead of direct reference to reactive values, since
    # some are updated from within the script and the new values are not reflected in the
    # reactive variables until data are sent, or flushed, to the browser, which may be after
    # further script instructions reference the variables
    #######################################################################################################

    userComputeComputeCPAG <- function(source2, pfactor1, pexp1, pfactor2, pexp2, ldpop) {

      q <- tryCatch({

        pycmd <- ""
        pycmd2 <- ""
        pycmd3 <- ""

        # Validate and copy input file
        f <- copyValidateUploadFile()
        if(f[["result"]]=="") {

          # Verify existence of copied file
          if(length(which(dir(inDirUserCompute, paste(f[["uGWASfile"]], ".csv", sep=""))==
                          paste(f[["uGWASfile"]], ".csv", sep="")))>0) {

            # Verify CPAG computation parameter values
            parV <- userComputeVerifyParameterValues(source2, pfactor2, pexp2)
            if(parV[["success"]]) {

              # Compose output file name
              # Format is UserFile-p1-GWAS2-p2-ldpop.csv
              cpagFile <- paste(f[["uGWASfile"]], "-p", pfactor1, "e-", sprintf("%02d", pexp1), "-",
                                source2, "-p", pfactor2, "e-", sprintf("%02d", pexp2), "-",
                                ldpop, ".csv", sep="")

              # Compose python command
              pycmd <- paste(oscmd, " cd ", pyDir, " ", oscmdsep, " ",
                             pyexe, " main.py usr-gwas --threads ", threads,
                             " --infile \"", pyInDirUserCompute, "/", f[["uGWASfile"]], ".csv\"",
                             " --delimitor \",\"",
                             " --usr-pcut ", pfactor1, "e-", pexp1,
                             #" --usr-pheno-name ", f[["uGWAScol"]]["phenotype"],
                             " --SNPcol \"", f[["uGWAScol"]]["snp"], "\"",
                             " --Pcol \"", f[["uGWAScol"]]["p"], "\"",
                             " --querydb ", source2,
                             " --cpagdb-pcut ", pfactor2, "e-", pexp2,
                             " --ld-clump 1",
                             " --lddb-pop ", ldpop,
                             " --outfile \"", pyOutDirUserCompute, "/", cpagFile, "\" 2>&1", sep="")

              # Compose command to include ontology, for NHGRI results only
              if(source2=="NHGRI") {
                pycmd2 <- paste(oscmd, " cd ", pyDir, " ", oscmdsep, " ",
                                pyexe, " main.py post_analysis --anno-ontology --anno-cols Trait2", " --infile \"",
                                pyOutDirUserCompute, "/", cpagFile, "\" --outfile \"", pyOutDirUserCompute,
                                "/", cpagFile, "\" 2>&1", sep="")
              } else {
                pycmd2 <- ""
              }

              cat(paste("\npycmd:  ", pycmd, "\n\n", sep=""), file=stderr())
              cat(paste("\npycmd2:  ", pycmd2, "\n\n", sep=""), file=stderr())

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
              # Note that errors during the system() call raise errors here, with the system message(s) as text
              # If output reaches the pipe, noe system() error occurred
              future({
                t0 <- proc.time()
                cpagPhase <- 1
                q0 <- suppressWarnings(system(pycmd, intern=T))
                if(is.null(attr(q0, "status")))
                  attr(q0, "status") <- 0
                if(attr(q0, "status")==0  & pycmd2!="") {
                  cpagPhase <- 2
                  q0 <- suppressWarnings(system(pycmd2, intern=T))
                  if(is.null(attr(q0, "status")))
                    attr(q0, "status") <- 0
                }
                t0 <- (proc.time()-t0)["elapsed"]
                list("cpagPhase"=cpagPhase, "q0"=q0, "t0"=t0)
              }) %...>%
              (function(q) {
                 progress$close()

                 cat(paste("Execution time: ", q[["t0"]], "\n", collapse=" "), file=stderr())
                 if(attr(q[["q0"]], "status")==0) {

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
                                   allSNP=input$userComputeSelectionIncludeTableAllSNPshare,
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

                     # Log success
                     writeLog(section="Compute", operation="computeCPAG", parameter1=pycmd,
                              parameter2=pycmd2, parameter3=q[["t0"]], status="success")

                   } else {
                     # Treat error here since the outer function (where future environment spawned)
                     # continues processing while the future/promise sequence executed independently
                     # Return contents of CPAG function call, since it may contain useful information on why
                     # a file was not generated (a common cause is the simple absence of an intersection
                     # of SNPs and phenotypes between the two GWAS sets at the specified p-threshold levels)
                     writeLog(section="Compute", operation="computeCPAG", parameter1=q[["q0"]],
                              parameter2=ifelse(q[["cpagPhase"]]==1, pycmd, pycmd2),
                              parameter3=q[["q0"]],  parameter4=q[["t0"]],
                              note="CPAG output file not found", status="error")
                     msgWindow("ERROR", "User Compute compute CPAG",
                               c(paste("Computed CPAG file ", outDirUserCompute, "/", cpagFile, " does not exist", sep=""),
                                 "Occurred during execution of:", ifelse(q[["cpagPhase"]]==1, pycmd, pycmd2),
                                 "<b>CPAG output follows</b>", q[["q0"]]))
                   }

                 } else {
                   # Treat error here since the outer function continues processing while future/promise executes independently
                   writeLog(section="Compute", operation="computeCPAG", parameter1=q[["q0"]],
                            parameter2=ifelse(q[["cpagPhase"]]==1, pycmd, pycmd2), parameter3=q[["t0"]],
                            note="Error during CPAG computation", status="error")
                   msgWindow("ERROR", "User Compute compute CPAG", c(q[["q0"]], "Occurred during execution of:",
                             ifelse(q[["cpagPhase"]]==1, pycmd, pycmd2)))
                 }

              })

            } else {

              # Open a dialog, window if any rules violated
              # Present options for correction or canceling
              # Follow the reactive events triggered from within the dialog for path back to compute function
              showModal(modalDialog(HTML(parV[["msg"]]),
                                    HTML("Choose from the following:"),
                                    div(
                                      actionButton(parV[["buttonID"]], parV[["buttonText"]]),
                                      style="display:inline-block; vertical-align:top; margin-top:-5px; margin-left:20px"
                                    ),
                                    div(
                                      modalButton("Cancel"),
                                      style="display:inline-block; vertical-align:top; margin-top:-5px; margin-left:10px"
                                    ),
                                    title="iCPAGdb Compute", size="m", easyClose=T,
                                    footer=NULL))

              # Return a non-error class object
              ""

            }
          } else {
            errorCondition(paste("Upload file ", f[["uGWASfile"]], " does not exist.  Please upload and try again.", sep=""))
          }
        } else {
          errorCondition(paste(f[["result"]], sep=""))
        }},
        warning=function(err) err,
        error=function(err) err
      )

      # Return message and CPAG commands, if any errors
      if(!any(class(q) %in% c("warning", "error"))) {
        ""
      } else {
        errorCondition(c(q[["message"]], pycmd, pycmd2, pycmd3))
      }

    }

    #######################################################################################################
    # Observe event function:  User Compute tab - set p-threshold 2 to max allowable value
    # Use local variables and pass values as parameters due to reactive variables not reflecting
    # updated values until flush
    #######################################################################################################

    observeEvent(input$userComputeComputeSetPMax2, {

      if(input$userComputeSource2=="H2P2") {
        pfactor2 <- pfactorH2P2
        pexp2 <- pexpH2P2
      } else if(input$userComputeSource2=="NHGRI") {
        pfactor2 <- pfactorNHGRI
        pexp2 <- pexpNHGRI
      } else {
        pfactor2 <- pfactorOther
        pexp2 <- pexpOther
      }
      updateRadioButtons(session, "userComputePfactor2", selected=pfactor2)
      updateSliderInput(session, "userComputePexp2", value=pexp2)
      removeModal()
      tryCatch(
        userComputeComputeCPAG(input$userComputeSource2, input$userComputePfactor1,
                               input$userComputePexp1, pfactor2, pexp2, input$userComputeLDpop),
        warning=function(err) err,
        error=function(err) err
      )
      if(any(class(q) %in% c("warning", "error"))) {
        output$userComputeResultsTable <- DT::renderDataTable(NULL)
        output$userComputeResultsHeatmap <- renderPlotly(NULL)
        writeLog(section="Compute", operation="computeCPAG", note=q[["message"]], status="error")
        msgWindow("ERROR", "User Compute compute CPAG setPsmall", q[["message"]])
      }

    }, ignoreInit=T)

    #######################################################################################################
    # Observe event function:  User Compute tab - compute CPAG button
    #######################################################################################################

    observeEvent(input$userComputeCompute, {

      # Log request
      writeLog(section="Compute", operation="computeReq",
               parameter1=input$userComputeSource1, parameter2=input$userComputeSource2,
               parameter3=paste(input$userComputePfactor1, "e-", sprintf("%02d", input$userComputePexp1), sep=""),
               parameter4=paste(input$userComputePfactor2, "e-", sprintf("%02d", input$userComputePexp2), sep=""),
               parameter5=input$userComputeLDpop, status="")

      q <- tryCatch(
        userComputeComputeCPAG(input$userComputeSource2, input$userComputePfactor1, input$userComputePexp1,
                               input$userComputePfactor2, input$userComputePexp2, input$userComputeLDpop),
        warning=function(err) err,
        error=function(err) err
      )
      if(any(class(q) %in% c("warning", "error"))) {
        output$userComputeResultsTable <- DT::renderDataTable(NULL)
        output$userComputeResultsHeatmap <- renderPlotly(NULL)
        writeLog(section="Compute", operation="computeCPAG", note=q[["message"]], status="error")
        msgWindow("ERROR", "User Compute compute CPAG", q[["message"]])
      }

    }, ignoreInit=T)

    #######################################################################################################
    # Observe event function:  User Compute tab - change in include all SNPs in shared SNP option
    #######################################################################################################

    observeEvent(input$userComputeSelectionIncludeTableAllSNPshare, {

      tryCatch({
        renderCPAGtable(currentUserComputeCPAGdata[userComputeResultsTableRowFilter,],
                        input$userComputeSelectionIncludeTableAllSNPshare,
                        "userComputeResultsTable")
        # Log action
        writeLog(section="Compute", operation="allSNP", status="")},
        warning=function(err) msgWindow("WARNING", "Usere Compute all SNPs", err),
        error=function(err) msgWindow("ERROR", "UserCompute all SNPs", err)
      )

    }, ignoreInit=T)

    #######################################################################################################
    # Observe event function:  User Compute tab - change in heatmap metric or top significant phenotype
    # pairs to plot
    #######################################################################################################

    observeEvent(c(input$userComputeSelectionHeatmapNphenotype, input$userComputeSelectionHeatmapMetric), {

      tryCatch({
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
        # Log action
        writeLog(section="Compute", operation="heatmapMetricTopAdjust", status="success")},
        warning=function(err) msgWindow("WARNING", "User Compute heatmap n, metric adjust", err),
        error=function(err) msgWindow("ERROR", "User Compute heatmap n, metric adjust", err)
      )

    }, ignoreInit=T)

    #######################################################################################################
    # Observe event function:  User Compute tab - apply filters
    #######################################################################################################

    observeEvent(input$userComputeSelectionFilterApply, {

      tryCatch({
        renderResults(data=currentUserComputeCPAGdata,
                      rowFilter="userComputeResultsTableRowFilter",
                      filterTrait=input$userComputeSelectionFilterTrait,
                      filterSNP=input$userComputeSelectionFilterSNP,
                      allSNP=input$userComputeSelectionIncludeTableAllSNPshare,
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
        # Log action
        writeLog(section="Compute", operation="filter", status="success")},
        warning=function(err) msgWindow("WARNING", "UserCompute filter", err),
        error=function(err) msgWindow("ERROR", "User Compute filter", err)
      )

    }, ignoreInit=T)

    #######################################################################################################
    # Observe event function:  User Compute tab - clear filters
    #######################################################################################################

    observeEvent(input$userComputeSelectionFilterClear, {

      # Updates are not effective until all output has been pushed to the client
      # This does not guarantee that new values are in effect at the time of the call to renderResults()
      # Therefore, assign values and pass constants as parameter values to renderResults()
      tryCatch({
        updateTextInput(session, "userComputeSelectionFilterTrait", value="")
        updateTextInput(session, "userComputeSelectionFilterSNP", value="")
        #updateTextInput(session, "userComputeSelectionFilterEFO", value="")
        updateTextInput(session, "userComputeSelectionFilterEFOparent", value=vector("character"))
        updateCheckboxInput(session, "userComputeSelectionIncludeTableCompoundEFO", value=F)
        renderResults(data=currentUserComputeCPAGdata,
                      rowFilter="userComputeResultsTableRowFilter",
                      filterTrait="",
                      filterSNP="",
                      allSNP=input$userComputeSelectionIncludeTableAllSNPshare,
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
                      heatmapNphenotype=input$userComputeSelectionHeatmapNphenotype)},
        warning=function(err) msgWindow("WARNING", "User Compute filter clear", err),
        error=function(err) msgWindow("ERROR", "User Compute filter clear", err)
      )

    }, ignoreInit=T)

    #######################################################################################################
    # Configure a download handler for the sample GWAS file on the userCompute tab
    #######################################################################################################

    output$userComputeSampleGWASdownload <-
      downloadHandler(filename="iCPAGdb-Sample-GWAS-top_EllinghausPCs_covid19.csv",
                      contentType="txt/csv",
                      content=function(file) {
                                tryCatch({
                                  x <- read.table(paste(inDir, "/top_EllinghausPCs_covid19-Sample.csv", sep=""), header=T, sep=",")
                                  # Log action
                                  writeLog(section="Compute", operation="download", parameter1="sampleInput", status="success")
                                  write.table(x, file, row.names=F, col.names=T, quote=F, sep=",")},
                                  warning=function(err) msgWindow("WARNING", "User Compute download sample GWAS", err),
                                  error=function(err) msgWindow("ERROR", "User Compute download sample GWAS", err)
                                )
                              }
                     )

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
                                tryCatch({
                                  # Compose filtered observation indices
                                  k <- resultsTableComposeFilterIndex(
                                         data=currentUserComputeCPAGdata,
                                         filterTrait=input$userComputeSelectionFilterTrait,
                                         filterSNP=input$userComputeSelectionFilterSNP,
                                         includeCompoundEFO=input$userComputeSelectionIncludeTableCompoundEFO,
                                         efoCol=userComputeResultsEFOcol,
                                         filterEFOparent=input$userComputeSelectionFilterEFOparent,
                                         efoParentCol=userComputeResultsEFOparentCol)
                                  # Log action
                                  writeLog(section="Compute", operation="download", parameter1="results", status="success")
                                  write.table(currentUserComputeCPAGdata[k,], file, row.names=F, col.names=T, quote=T, sep=",")},
                                  warning=function(err) msgWindow("WARNING", "Usr Compute download results", err),
                                  error=function(err) msgWindow("ERROR", "User Compute download results", err)
                                )
                              }
                     )

  }

)
