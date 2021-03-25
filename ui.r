# Duke University Cross Phenotype Analysis of GWAS Database (iCPAGdb)
# November 2020

# Shiny app user interface function

# Information on shiny available at:
# https://shiny.rstudio.com/
# https://github.com/rstudio/shiny
# https://cran.r-project.org/web/packages/shiny/shiny.pdf

library(shiny)
library(shinyjs)
library(DT)
library(plotly)
library(RSQLite)

# Set current working directory to project python and data resources
os <- 2
setwd(c("C:/Projects/Duke/H2P2GenomeWideAssociationStudy/CPAG/iCPAGdb/App-Devel", "/srv/shiny-server/CPAG/explore")[os])

# Specify database location
dbloc <- "pyCPAG/db/cpag_gwasumstat_v1.1.db"

# Connect to GWAS database
db <- dbConnect(SQLite(), dbloc)

# Compose GWAS sources for input selection
# Assign names for selectInput use
# Force order of select studies
gwasSource <- c(unique(c("NHGRI", "H2P2", sort(dbGetQuery(db, "select distinct source from GWAStable order by source")[,1]))),
                "mol_gwas", "clin_gwas")
names(gwasSource) <- c(gwasSource[1:(length(gwasSource)-2)], "Molecular", "Clinical")

dbDisconnect(db)

shinyUI(

  fluidPage(

    useShinyjs(),
    #includeCSS("../app/www/iCPAGdb-spacelab.css"),
    #shinythemes::themeSelector(),

    title="iCPAGdb",

    # Configure some styles
    # It is better to apply these to individual objects when possible, to guarantee application
    # For instance, the navbarPage, below, specifies a theme, which causes styles established here to be ignored
    # Defining styles as the final step of the navbarPage definition overrides styles defined in the theme
    # Use the "inspection" feature of your browser (right-button over object of interest) to reveal classes
    tags$head(
      # Reposition and alter appearance of notification window
      tags$style(HTML(".shiny-notification {font-size:14px; color:black; font-weight:bold; width:50%;
                      height=200px; position:fixed; top:calc(50%); left:calc(25%)}")),
      # Make text and background of fileInput progress bar transparent
      # Otherwise, "File uploaded" message appears once file upoaded, but prior to a read.table operation accomplished
      # Note that background:transparent causes the progress bar to disappear, which is useful for fileInput(), but also
      # causes the bar to disappear when using shiny::Progress$new(), the standard progress object
      # Progress$new() has a style parameter, but accepts only two values: notification or old
      # Notification uses the .progress-bar style below, old seems to behave  similarly
      # Unfortunately the unmodifiable text of the fileInput() progress bar and the fixed style of Progress bars
      # does not give much flexibility for tailoring progress bar appearance
      #tags$style(HTML(".progress-bar {color: transparent!important; background:transparent!important}")),
      # Change color of text that appears inside of the bar
      # This affects the bar for file input objects, but not for progress#new() objects, since their text appears
      # below the bar
      #tags$style(HTML(".progress-bar {color: transparent!important}")),
      # Change the color of the progress bar (file input and progress meter)
      #tags$style(HTML(".progress-bar {background-color: gray;}"))
      #HTML(".shiny-file-input-progress {color: transparent!important}")
      # Hide file input progress bar
      #tag$style(HTML(".shiny-file-input-progress {display: none}")),
      # Customize the modal window
      #tags$style(".modal-body {padding: 10px}
      #            .modal-content  {-webkit-border-radius: 6px !important;-moz-border-radius: 6px !important;border-radius: 6px !important;}
      #            .modal-dialog { width: 240px; display: inline-block; text-align: left; vertical-align: top;}
      #            .modal-header {background-color: #339FFF; border-top-left-radius: 6px; border-top-right-radius: 6px}
      #            .modal { text-align: right; padding-right:10px; padding-top: 24px;}
      #            .close { font-size: 16px}"))
      #tags$style(HTML(".modal-lg {position: relative; display: flex; flex-direction: column; margin-top: 50%}"))
      #tags$style(HTML(".modal-lg {width: 50%; margin-top: 10%}")),
      #tags$style(HTML(".modal {margin-top: 10%}")),
      #tags$style(HTML(".modal-header {color: #ffffff; background-color: #0066cc; border-top-left-radius: 6px; border-top-right-radius: 6px}"))
      # Adjust style of action buttons
      # Button appearance instructions are in iCPAGdb-spacelab.css (.btn-default class)
      # An attempt was made to include corresponding css tags here, but they are over-ridden when a css file is used
      # Note that the active Shiny process may have to be reloaded in order for changes to the css file to take effect
      #tags$style(HTML(".btn {color:white; background:linear-gradient(#54b4eb, #2fa4e7 60%, #0088dd)}")),
      # Adjust style of file input button
      #tags$style(HTML(".btn-file {color:white; background:linear-gradient(#54b4eb, #2fa4e7 60%, #0088dd)}")),
    ),

    div(

      div(
        HTML("<H3>iCPAGdb - A hypothesis engine for cross-phenotype genetic associations connecting molecular, cellular, and human disease phenotypes</H3><br>"),
        style="display:inline-block; vertical-align:top; margin-top:0px; width:97%"
      ),
      # Feature enable element
      div(
        tags$style(HTML("#featureEnable {font-size: 10px; border-style: none; box-shadow: none}")),
        textInput("featureEnable", ""),
        style="display:inline-block; vertical-align:top; margin-top:0px; width:50px"
      ),

      navbarPage(
        id="tabsetCPAG", title="", inverse=T,
        # Themes have fixed search locations:
        # theme=shinytheme("xyz") searches for xyz.min.css in the shinythemes/css directory
        # of the shinythemes package directory in the R installation directory
        # theme="xyz.css" searches for file xyz.css in the www directory within the directory that the
        # shiny app was executed from within, not the current working directory, as set by setwd() 
        # The following method does not report an error, but the theme file is not loaded
        # theme="app/www/iCPAGdb-spacelab.css",
        # The following method does not report an error, but the theme file is not loaded
        # theme="C:/Users/Kyung Soon/Documents/R/win-library/4.0/shinythemes/shinythemes/css/spacelab.min.css",
        # The following methods load specified themes
        # Use the shinytheme() function (that returns a character vector)
        # This will locate the css file in the shinythemes package directory
        # Note that font urls within the css reference relative directory locations, which are
        # interpreted with respect to the css file location
        theme=shinythemes::shinytheme("spacelab"),
        # Use a character vector indicating the desired theme (file name)
        # This method is identical to theme=shinythemes::shinytheme("spacelab")
        # theme="shinythemes/css/spacelab.min.css",
        # Grab custom theme css from www directory
        # Note that font file references must be altered and/or fonts must be copied to
        # appropriate directories as referenced in css @font-face {} instructions
        # This method does not function (for me, anyway) in a Shiny Server environment
        # theme="iCPAGdb-spacelab.css",
        # theme="spacelab.min.css",

        # Review iCPAGdb
        tabPanel(title="Review iCPAGdb", value="tabPanelReview",
          div(
            # Prompts
            div(
              sidebarPanel(width=12,
                # Select a data set
                div(
                  div(
                    HTML("<b>1. Select a data set to review</b>"),
                    style="vertical-align:top; margin-top:0px;"
                  ),
                  div(
                    DT::dataTableOutput(outputId="reviewSelectionTable", width="98%"),
                    style="width:100%; vertical-align:top; margin-top:15px"
                  ),
                  style="margin-top:0px"
                ),
                # Filter
                div(
                  HTML("<br><b>2. Filter</b><br>"),
                  div(
                    checkboxInput("reviewSelectionIncludeTableAllSNPshare", "Include all SNPs in table"),
                    style="display:inline-block; vertical-align:top; margin-top:20px; width:12%"
                  ),
                  div(
                    textInput("reviewSelectionFilterTrait", "Trait filter", width="85%"),
                    style="display:inline-block; vertical-align:top; margin-top:0px; width:12%"
                  ),
                  div(
                    textInput("reviewSelectionFilterSNP", "SNP filter", width="85%"),
                    style="display:inline-block; vertical-align:top; margin-top:0px; width:12%"
                  ),
                  #div(
                  #  textInput("reviewSelectionFilterEFO", HTML("EFO filter <i>comma separated</i>"), width="85%"),
                  #  style="display:inline-block; vertical-align:top; margin-top:0px; width:12%"
                  #),
                  div(
                    selectInput("reviewSelectionFilterEFOparent", HTML("EFO filter <i>select multiple</i>"),
                                choices=vector("character"), multiple=T, width="85%"),
                    style="width:250px; display:inline-block; vertical-align:top; margin-top:0px; width:12%"
                  ),
                  div(
                    checkboxInput("reviewSelectionIncludeTableCompoundEFO", "Include compound EFOs"),
                    style="display:inline-block; vertical-align:top; margin-top:20px; width:12%"
                  ),
                  #div(
                  #  HTML("<a href=\"https://www.ebi.ac.uk/efo/\" target=\"_blank\">Experimental Factor Ontology (EFO)</a>"),
                  #  style="display:inline-block; vertical-align:top; margin-top:50px; width:12%"
                  #),
                  div(
                    actionButton("reviewSelectionFilterApply", "Apply filters", width="100px"),
                    style="display:inline-block; vertical-align:top; margin-top:25px; margin-left:0px"
                  ),
                  div(
                    actionButton("reviewSelectionFilterClear", "Clear filters", width="100px"),
                    style="display:inline-block; vertical-align:top; margin-top:25px; margin-left:10px"
                  ),
                  div(
                    downloadButton("reviewResultsDownload", "Download filtered records"),
                    style="display:inline-block; vertical-align:top; margin-top:25px; margin-left:10px"
                  ),
                  div(
                    div(
                      radioButtons("reviewSelectionHeatmapMetric", "Heatmap metric",
                                   choices=c("Fisher"="P_fisher", "Bonferroni"="Padj_Bonferroni", "FDR"="Padj_FDR",
                                             "Jaccard"="Jaccard", "Chao-Sorensen"="ChaoSorensen"),
                                   inline=T, selected="ChaoSorensen"),
                      style="display:inline-block; vertical-align:top; margin-top:30px; margin-left:0px"
                    ),
                    div(
                      radioButtons("reviewSelectionHeatmapNphenotype", "Display top significant phenotype pairs in heatmap",
                                   choices=c("10"="10", "25"="25", "50"="50", "100"="100", "250"="250", "500"="500", "1,000"="1000", "all"="all"),
                                   inline=T, selected="25"),
                      style="display:inline-block; vertical-align:top; margin-top:30px; margin-left:75px"
                    ),
                    style="margin-top:-20px"
                  ),
                  style="margin-top:0px; margin-bottom:-10px"
                )
              ),
              style="margin-left:-15px; margin-top:0px; margin-right:-30px"
            ),
            # Review CPAG results
            div(
              tabsetPanel(id="tabsetReviewResults", type="pills",
                # Table
                tabPanel(title="Table", value="tabPanelReviewResultsTable",
                  HTML("<br>"),
                  DT::dataTableOutput(outputId="reviewResultsTable", width="98%")
                ),
                # Heatmap
                tabPanel(title="Heatmap", value="tabPanelReviewResultsHeatmap",
                  HTML("<br><center>"),
                  plotlyOutput("reviewResultsHeatmap"),
                  HTML("</center>")
                )
              ),
              style="width=100%; margin-top:20px;"
            ),
            style="margin-top:-40px; margin-left:-15px"
          ),
        ),

        # Explore CPAG associations
        # Note that "hidden" hides tab contents, but not the tab itself
        # Begin with a blank title and assign it after tab is hidden (top of server script)
        # Prblem with above is that shinyjs::show() does not present contents of tab
        # Title is "Explore iCPAGdb associations"
        # uiOutput() creates a reactive output variable that appears as HTML on the page
        # Since title is assigned to this variable, the contents will appear in the tab title
        # The initial value is empty and is modified by the server script
        tabPanel(title=uiOutput("exploreTabTitle"), value="tabPanelExplore",
          conditionalPanel("tabPanelExploreVisible",
            div(
              # Prompts
              div(
                sidebarPanel(width=12,
                  # Compute
                  HTML("<b>1. Compute</b>"),
                  div(
                    div(
                      selectInput(inputId="exploreSource1", label="GWAS source one", choices=gwasSource, width="85%"),
                      style="display:inline-block; vertical-align:top; margin-top:0px; width:12%"
                    ),
                    div(
                      selectInput(inputId="exploreSource2", label="GWAS source two", choices=gwasSource, width="85%"),
                      style="display:inline-block; vertical-align:top; margin-top:0px; width:12%"
                    ),

                    # P threshold 1
                    div(
                      div(
                        div(
                          HTML("<b>p-threshold<sub>1</sub> (factor<sub>1</sub> X 10<sup>-x<sub>1</sub></sup>)</b>"),
                          style="margin-top:-30px; width:100%"
                        ),
                        div(
                          radioButtons(inputId="explorePfactor1", HTML("factor<sub>1</sub>"),
                                       choices=c(1, 5), inline=T, selected=5),
                          style="display:inline-block; vertical-align:top; margin-top:20px; width:30%"
                        ),
                        div(
                          sliderInput(inputId="explorePexp1", HTML("x<sub>1</sub>"),
                                      min=3, max=20, step=1, value=8, width="85%"),
                          style="display:inline-block; vertical-align:top; margin-top:5px; width:65%"
                        ),
                        style="display:inline-block; vertical-align:top; margin-top:0px; width:48%"
                      ),
                      # P threshold 2
                      div(
                        div(
                          HTML("<b>p-threshold<sub>2</sub> (factor<sub>2</sub> X 10<sup>-x<sub>2</sub></sup>)</b>"),
                          style="margin-top:-30px; width:100%"
                        ),
                        div(
                          radioButtons(inputId="explorePfactor2", HTML("factor<sub>2</sub>"),
                                       choices=c(1, 5), inline=T, selected=5),
                          style="display:inline-block; vertical-align:top; margin-top:20px; width:30%"
                        ),
                        div(
                          sliderInput(inputId="explorePexp2", HTML("x<sub>2</sub>"),
                                      min=3, max=20, step=1, value=8, width="85%"),
                          style="display:inline-block; vertical-align:top; margin-top:5px; width:65%"
                        ),
                        style="display:inline-block; vertical-align:top; margin-top:0px; width:48%"
                      ),
                      div(
                        HTML("Note: p-threshold maximums are H2P2 = 1X10<sup>-5</sup>, NHGRI = 5X10<sup>-8</sup>, all others = 1X10<sup>-3</sup>"),
                        style="margin-top:-5px; width:90%"
                      ),
                      style="display:inline-block; vertical-align:top; margin-top:-5px; width:36%"
                    ),
                    div(
                      radioButtons(inputId="exploreLDpop", "LD 1000 Genomes population",
                                   choices=c("European"="EUR", "African"="AFR", "Asian"="EAS"), inline=T),
                      style="display:inline-block; vertical-align:top; margin-top:5px; width:15%"
                    ),
                    # Exec button
                    div(
                      actionButton(inputId="exploreCompute", "Compute CPAG"),
                      style="display:inline-block; vertical-align:top; margin-top:25px; width:12%"
                    ),
                    style="margin-top:20px"
                  ),
                  HTML("<hr style='height:1px;color:black;background-color:black'>"),
                  # Filter
                  HTML("<b>2. Filter</b>"),
                  div(
                    div(
                      checkboxInput("exploreSelectionIncludeTableAllSNPshare", "Include all SNPs in table"),
                      style="display:inline-block; vertical-align:top; margin-top:20px; width:12%"
                    ),
                    div(
                      textInput("exploreSelectionFilterTrait", "Trait filter", width="85%"),
                      style="display:inline-block; vertical-align:top; margin-top:0px; width:12%"
                    ),
                    div(
                      textInput("exploreSelectionFilterSNP", "SNP filter", width="85%"),
                      style="display:inline-block; vertical-align:top; margin-top:0px; width:12%"
                    ),
                    #div(
                    #  textInput("exploreSelectionFilterEFO", HTML("EFO filter <i>comma separated</i>")),
                    #  style="display:inline-block; vertical-align:top; margin-top:0px; width:250px"
                    #),
                    div(
                      selectInput("exploreSelectionFilterEFOparent", HTML("EFO filter <i>select multiple</i>"),
                                  choices=vector("character"), multiple=T, width="85%"),
                      style="display:inline-block; vertical-align:top; margin-top:0px; width:12%"
                    ),
                    div(
                      checkboxInput("exploreSelectionIncludeTableCompoundEFO", "Include compound EFOs"),
                      style="display:inline-block; vertical-align:top; margin-top:20px; width:12%"
                    ),
                    #div(
                    #  HTML("<a href=\"https://www.ebi.ac.uk/efo/\" target=\"_blank\">Experimental Factor Ontology (EFO)</a>"),
                    #  style="display:inline-block; vertical-align:top; margin-top:0px; width:12%"
                    #),
                    div(
                      actionButton("exploreSelectionFilterApply", "Apply filters", width="100px"),
                      style="display:inline-block; vertical-align:top; margin-top:25px"
                    ),
                    div(
                      actionButton("exploreSelectionFilterClear", "Clear filters", width="100px"),
                      style="display:inline-block; vertical-align:top; margin-top:25px; margin-left:10px"
                    ),
                    div(
                      downloadButton("exploreResultsDownload", "Download filtered records"),
                      style="display:inline-block; vertical-align:top; margin-top:25px; margin-left:10px"
                    ),
                    style="margin-top:-5px"
                  ),
                  div(
                    div(
                      radioButtons("exploreSelectionHeatmapMetric", "Heatmap metric",
                                   choices=c("Fisher"="P_fisher", "Bonferroni"="Padj_Bonferroni", "FDR"="Padj_FDR",
                                             "Jaccard"="Jaccard", "Chao-Sorensen"="ChaoSorensen"),
                                   inline=T, selected="ChaoSorensen"),
                      style="display:inline-block; vertical-align:top; margin-top:20px"
                    ),
                    div(
                      radioButtons("exploreSelectionHeatmapNphenotype", "Display top significant phenotype pairs in heatmap",
                                   choices=c("10"="10", "25"="25", "50"="50", "100"="100", "250"="250", "500"="500", "1,000"="1000", "all"="all"),
                                   inline=T, selected="25"),
                      style="display:inline-block; vertical-align:top; margin-top:20px; margin-left:75px"
                    ),
                    style="margin-top:-5px"
                  )
                ),
                style="margin-left:-15px; margin-top:0px; margin-right:-30px"
              ),
              # Explore CPAG results
              div(
                tabsetPanel(id="tabsetExploreResults", type="pills",
                  # Table
                  tabPanel(title="Table", value="tabPanelExploreResultsTable",
                    HTML("<br>"),
                    DT::dataTableOutput(outputId="exploreResultsTable", width="98%")
                  ),
                  # Heatmap
                  tabPanel(title="Heatmap", value="tabPanelExploreResultsHeatmap",
                    HTML("<br><center>"),
                    plotlyOutput("exploreResultsHeatmap"),
                    HTML("</center>")
                  )
                ),
                style="width=100%; margin-top:20px"
              ),
              style="margin-top:-40px; margin-left:-15px;"
            )
          )
        ),

        # Upload and compute CPAG
        tabPanel(title="Upload GWAS and compute CPAG", value="userComputeGWAS",
          div(
            # Prompts
            div(
              sidebarPanel(width=12,
                # Upload
                div(
                  HTML("<b>1. Upload a GWAS file</b>"),
                  div(
                    HTML("Note:  Maximum file size is 1GB.  Expected upload time is appproxiamtaely 30 seconds per 100MB.  For faster upload, reduce your input file to two columns (SNP and P-value) and/or pre-clump for only lead SNPs at your desired threshold. Upload progress is indicated in the bar below the \"Browse\" button.  To download a sample GWAS file for review, click here: &nbsp;&nbsp;"),
                    downloadLink(outputId="userComputeSampleGWASdownload", label="sample GWAS file (severe COVID-19, Ellinghaus et al. 2020)"),
                    style="margin-top:10px; margin-bottom:10px; margin-right:25%"
                  ),
                  # Controls
                  div(
                    div(
                      fileInput("userComputeBrowseFile", "Choose file", multiple=FALSE,
                                accept=c("text/csv", "text/comma-separated-values,text/plain", ".csv"), width="85%"),
                      style="display:inline-block; vertical-align:top; margin-top:0px; width:24%"
                    ),
                    div(
                      radioButtons("userComputeDelimiter", "Delimiter", choices=c("Comma", "Tab"), inline=T, selected="Comma"),
                      style="display:inline-block; vertical-align:top; margin-top:0px; margin-left:5px; width:12%"
                    ),
                    div(
                      disabled(textInput(inputId="userComputePhenotypeCol", label="Trait (phenotype) column",
                               value="One phenotype per file", width="85%")),
                      style="display:inline-block; vertical-align:top; margin-top:0px; width:12%"
                    ),
                    div(
                      textInput(inputId="userComputeSNPcol", label="SNP column", width="85%"),
                      style="display:inline-block; vertical-align:top; margin-top:0px; width:12%"
                    ),
                    div(
                      textInput(inputId="userComputePcol", label="P (significance) column", width="68%"),
                      style="display:inline-block; vertical-align:top; margin-top:0px; width:15%"
                    ),
                    # Upload button
                    #div(
                    #  actionButton(inputId="userComputeUploadFile", "Upload file", width="120px"),
                    #  style="display:inline-block; vertical-align:top; margin-top:40px; width:12%"
                    #),
                    style="margin-top:20px; margin-bottom:0px"
                  ),
                  style="margin-top:0px; margin-bottom:-30px"
                ),
                HTML("<hr style='height:1px;color:black;background-color:black'>"),
                # Compute
                div(
                  HTML("<b>2. Compute</b>"),
                  div(
                    div(
                      selectInput(inputId="userComputeSource1", label="GWAS source one", choices="User Supplied GWAS", width="85%"),
                      style="display:inline-block; vertical-align:top; margin-top:0px; width:12%"
                    ),
                    div(
                      selectInput(inputId="userComputeSource2", label="GWAS source two", choices=gwasSource, width="85%"),
                      style="display:inline-block; vertical-align:top; margin-top:0px; width:12%"
                    ),
                    # P threshold 1
                    div(
                      div(
                        div(
                          HTML("<b>p-threshold<sub>1</sub> (factor<sub>1</sub> X 10<sup>-x<sub>1</sub></sup>)</b>"),
                          style="margin-top:-30px; width:100%"
                        ),
                        div(
                          radioButtons(inputId="userComputePfactor1", HTML("factor<sub>1</sub>"),
                                       choices=c(1, 5), inline=T, selected=5),
                          style="display:inline-block; vertical-align:top; margin-top:20px; width:30%"
                        ),
                        div(
                          sliderInput(inputId="userComputePexp1", HTML("x<sub>1</sub>"),
                                      min=3, max=20, step=1, value=8, width="85%"),
                          style="display:inline-block; vertical-align:top; margin-top:5px; width:65%"
                        ),
                        style="display:inline-block; vertical-align:top; margin-top:0px; width:48%"
                      ),
                      # P threshold 2
                      div(
                        div(
                          HTML("<b>p-threshold<sub>2</sub> (factor<sub>2</sub> X 10<sup>-x<sub>2</sub></sup>)</b>"),
                          style="margin-top:-30px; width:100%"
                        ),
                        div(
                          radioButtons(inputId="userComputePfactor2", HTML("factor<sub>2</sub>"),
                                       choices=c(1, 5), inline=T, selected=5),
                          style="display:inline-block; vertical-align:top; margin-top:20px; width:30%"
                        ),
                        div(
                          sliderInput(inputId="userComputePexp2", HTML("x<sub>2</sub>"),
                                      min=3, max=20, step=1, value=8, width="85%"),
                          style="display:inline-block; vertical-align:top; margin-top:5px; width:65%"
                        ),
                        style="display:inline-block; vertical-align:top; margin-top:0px; width:48%"
                      ),
                      div(
                        HTML("Note: p-threshold maximums are H2P2 = 1X10<sup>-5</sup>, NHGRI = 5X10<sup>-8</sup>, all others = 1X10<sup>-3</sup>"),
                        style="margin-top:-5px; width:90%"
                      ),
                      style="display:inline-block; vertical-align:top; margin-top:-5px; width:36%"
                    ),
                    div(
                      radioButtons(inputId="userComputeLDpop", "LD 1000 Genomes population",
                                   choices=c("European"="EUR", "African"="AFR", "Asian"="EAS"), inline=T),
                      style="display:inline-block; vertical-align:top; margin-top:0px; margin-left:5px; width:15%"
                    ),
                    # Exec button
                    div(
                      actionButton(inputId="userComputeCompute", "Compute CPAG", width="130px"),
                      style="display:inline-block; vertical-align:top; margin-top:20px; width:12%"
                    ),
                    style="margin-top:20px; margin-bottom:0px"
                  ),
                  style="margin-top:0px; margin-bottom:0px"
                ),
                HTML("<hr style='height:1px;color:black;background-color:black'>"),
                # Filter
                div(
                  HTML("<b>3. Filter</b>"),
                  div(
                    div(
                      checkboxInput("userComputeSelectionIncludeTableAllSNPshare", "Include all SNPs in table"),
                      style="display:inline-block; vertical-align:top; margin-top:20px; width:12%"
                    ),
                    div(
                      textInput("userComputeSelectionFilterTrait", "Trait filter", width="85%"),
                      style="display:inline-block; vertical-align:top; margin-top:0px; width:12%"
                    ),
                    div(
                      textInput("userComputeSelectionFilterSNP", "SNP filter", width="85%"),
                      style="display:inline-block; vertical-align:top; margin-top:0px; width:12%"
                    ),
                    #div(
                    #  textInput("userComputeSelectionFilterEFO", HTML("EFO filter <i>comma separated</i>")),
                    #  style="display:inline-block; vertical-align:top; margin-top:0px; width:250px"
                    #),
                    div(
                      selectInput("userComputeSelectionFilterEFOparent", HTML("EFO filter <i>select multiple</i>"),
                                  choices=vector("character"), multiple=T, width="85%"),
                      style="display:inline-block; vertical-align:top; margin-top:0px; width:12%"
                    ),
                    div(
                      checkboxInput("userComputeSelectionIncludeTableCompoundEFO", "Include compound EFOs"),
                      style="display:inline-block; vertical-align:top; margin-top:20px; width:12%"
                    ),
                    #div(
                    #  HTML("<a href=\"https://www.ebi.ac.uk/efo/\" target=\"_blank\">Experimental Factor Ontology (EFO)</a>"),
                    #  style="display:inline-block; vertical-align:top; margin-top:0px; width:12%"
                    #),
                    div(
                      actionButton("userComputeSelectionFilterApply", "Apply filters", width="100px"),
                      style="display:inline-block; vertical-align:top; margin-top:25px"
                    ),
                    div(
                      actionButton("userComputeSelectionFilterClear", "Clear filters", width="100px"),
                      style="display:inline-block; vertical-align:top; margin-top:25px; margin-left:10px"
                    ),
                    div(
                      downloadButton("userComputeResultsDownload", "Download filtered records"),
                      style="display:inline-block; vertical-align:top; margin-top:25px; margin-left:10px"
                    ),
                    style="margin-top:0px; margin-bottom:0px"
                  ),
                  div(
                    div(
                      radioButtons("userComputeSelectionHeatmapMetric", "Heatmap metric",
                                   choices=c("Fisher"="P_fisher", "Bonferroni"="Padj_Bonferroni", "FDR"="Padj_FDR",
                                             "Jaccard"="Jaccard", "Chao-Sorensen"="ChaoSorensen"),
                                   inline=T, selected="ChaoSorensen"),
                      style="display:inline-block; vertical-align:top; margin-top:20px"
                    ),
                    div(
                      radioButtons("userComputeSelectionHeatmapNphenotype", "Display top significant phenotype pairs in heatmap",
                                   choices=c("10"="10", "25"="25", "50"="50", "100"="100", "250"="250", "500"="500", "1,000"="1000", "all"="all"),
                                   inline=T, selected="25"),
                      style="display:inline-block; vertical-align:top; margin-top:20px; margin-left:75px"
                    ),
                    style="margin-top:0px; margin-bottom:0px"
                  ),
                  style="margin-top:0px; margin-bottom:0px"
                )
              ),
              style="margin-top:0px"
            ),
            style="margin-top:-40px; margin-left:-15px; margin-right:-30px"
          ),
          # User compute tab CPAG results
          div(
            tabsetPanel(id="tabsetUserComputeResults", type="pills",
              # Table
              tabPanel(title="Table", value="tabPanelUserComputeResultsTable",
                HTML("<br>"),
                DT::dataTableOutput(outputId="userComputeResultsTable", width="98%")
              ),
              # Heatmap
              tabPanel(title="Heatmap", value="tabPanelUserComputeResultsHeatmap",
                HTML("<br><center>"),
                plotlyOutput("userComputeResultsHeatmap"),
                HTML("</center>")
              )
            ),
            style="width=100%; margin-top:20px"
          ),
          style="margin-left:-15px; margin-top:-25px"
        ),

        # Quick-start
        tabPanel(title="Quick-start", value="tabPanelQuickStart",
          div(includeHTML("app/QuickStart.html"), style="margin-top:15px; margin-left:20%; margin-right:20%")
        ),

        # About
        tabPanel(title="About iCPAGdb", value="tabPanelAbout",
          div(includeHTML("app/About.html"), style="margin-top:15px; margin-left:20%; margin-right:20%")
        ),

        # Bibliography
        tabPanel(title="Bibliography", value="tabPanelBibliography",
          div(includeHTML("app/Bibliography.html"), style="margin-top:15px; margin-left:20%; margin-right:20%")
        ),

        # Define styles within (as final step of) navbarPage to overide specific theme styles
        # Hide panel, otherwise clickable, yet invisible elements appear in the bar
        # with contents equal to the contents of the style tag, below
        shinyjs::hidden(
          tabPanel(title="", value="navBarStyles",
            # Use the "inspection" feature of your browser (right-button over object of interest) to reveal classes
            # Grab colors with the browser color sampler (dropper)
            # Buttons, default and hover
            tags$style(HTML(".btn-default{color:#ffffff;background:linear-gradient(#6c93be, #446e9b 60%, #3e648d);border-color:#446e9b}")),
            tags$style(HTML(".btn-default:hover{color:#ffffff;background:linear-gradient(#6e96c2, #45709e 80%, #416994); border-color:#447e9b}")),
            # Hide title, round corners of bar
            tags$style(HTML(".navbar-brand {display: none;}")),
            tags$style(HTML(".navbar {border-top-left-radius:7px; border-top-right-radius:7px}")),
            # Modal box appearance
            tags$style(HTML(".modal {margin-top:10%}")),
            tags$style(HTML(".modal-header {background:linear-gradient(#6c93be, #446e9b 60%, #3e648d); border-radius:5px;}")),
            # Opacity of main window
            tags$style(HTML(".modal-backdrop.in{opacity:0.15;filter:alpha(opacity=50)}")),
            # Title
            tags$style(HTML(".modal-title{color:#ffffff; margin:0;line-height:1.42857143}")),
            # tabsetPanel buttons
            # In the rendered HTML, the buttons appear as anchor tags within in unordered lists
            # View (toggle) pseudo class of object using the browser inspector to review pseudo classes (hover,
            # active, focus, focus-within) and related css properties
            # Pill button height (all pseudo classes)
            tags$style(HTML(".nav-pills{line-height:0.75}")),
            #tags$style(HTML(".nav-pills>li>a:hover{color:green;}")),
            #tags$style(HTML(".nav-pills>li>a:active{color:red;}")),
            # Create a hidden variable to instruct on tabPanelExplore visibility
            checkboxInput("tabPanelExploreVisible", "", F)
          )
        )

      ),

      # Feature enable element
      #div(
      #  textInput("featureEnable", ""),
      #  style="width:50px"
      #),

      style="margin-left: 10px"
    )
  )
)
