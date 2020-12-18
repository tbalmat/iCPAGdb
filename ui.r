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
setwd(c("C:/Projects/Duke/H2P2GenomeWideAssociationStudy/CPAG/iCPAGdb/App/pyCPAG",
        "/srv/shiny-server/CPAG/explore/pyCPAG")[1])

# Specify database location
dbloc <- "db/cpag_gwasumstat_v1.1.db"

# Connect to GWAS database
db <- dbConnect(RSQLite::SQLite(), dbloc)

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
    #includeCSS("App/V2/style.css"),

    title="iCPAGdb",

    tags$head(
      # Reposition and alter appearance of notification window
      tags$style(
        HTML(".shiny-notification {font-size:14px; color:black; font-weight:bold; width:50%; height=200px; position:fixed;
              top:calc(50%); left:calc(25%)}")
      ),
      # Make text and background of fileInput progress bar transparent
      # Otherwise, "File uploaded" message appears once file upoaded, but prior to a read.table operation accomplished
      # Note that background:transparent causes the progress bar to disappear, which is useful for fileInput(), but also
      # causes the bar to disappear when using shiny::Progress$new(), the standard progress object
      # Progress$new() has a style parameter, but accepts on;y two values: notification or old
      # Notification uses the .progress-bar style below, old seems to behave  similarly
      # Unfortunately the unmodifiable text of the fileInput() progress bar and the fixed style of Progress bars
      # does not give much flexibility for tailoring progress bar appearance
      tags$style(
        #HTML(".progress-bar {color: transparent!important; background:transparent!important}")
        #HTML(".progress-bar {color: transparent!important}")
        # Here's the real answer!
        HTML(".shiny-file-input-progress {display: none}")
      )
    ),

    div(

      div(
        HTML("<H3>iCPAGdb - A hypothesis engine for cross-phenotype genetic associations connecting molecular, cellular, and human disease phenotypes</H3><br>"),
        style="display:inline-block; vertical-align:top; margin-top:0px; width:97%"
      ),
      # Feature enable element
      div(
        textInput("featureEnable", ""),
        style="display:inline-block; vertical-align:top; margin-top:0px; width:50px"
      ),

      tabsetPanel(id="tabsetCPAG",

        # Review iCPAGdb
        tabPanel(title="Review iCPAGdb", value="tabPanelReview",
          div(
            # Prompts
            HTML("<br>"),
            sidebarPanel(width=12,
              div(
                HTML("<b>1. Select a data set to review</b>"),
                style="vertical-align:top; margin-top:0px;"
              ),
              div(
                DT::dataTableOutput(outputId="reviewSelectionTable", width="98%"),
                style="width:100%; vertical-align:top; margin-top:20px"
              ),
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
                actionButton("reviewSelectionFilterApply", "Apply filters", width="100px",
                             style="color:white; background:linear-gradient(#54b4eb, #2fa4e7 60%, #0088dd)"),
                style="display:inline-block; vertical-align:top; margin-top:25px; margin-left:0px"
              ),
              div(
                actionButton("reviewSelectionFilterClear", "Clear filters", width="100px",
                             style="color:white; background:linear-gradient(#54b4eb, #2fa4e7 60%, #0088dd)"),
                style="display:inline-block; vertical-align:top; margin-top:25px; margin-left:10px"
              ),
              div(
                downloadButton("reviewResultsDownload", "Download filtered records",
                                style="color:white; background:linear-gradient(#54b4eb, #2fa4e7 60%, #0088dd)"),
                style="display:inline-block; vertical-align:top; margin-top:25px; margin-left:10px"
              ),
              HTML("<br>"),
              div(
                radioButtons("reviewSelectionHeatmapMetric", "Heatmap metric",
                             choices=c("Fisher"="P_fisher", "Bonferroni"="Padj_Bonferroni", "FDR"="Padj_FDR",
                                       "Jaccard"="Jaccard", "Chao-Sorensen"="ChaoSorensen"),
                             inline=T, selected="P_fisher"),
                style="display:inline-block; vertical-align:top; margin-top:30px; margin-left:0px"
              ),
              div(
                radioButtons("reviewSelectionHeatmapNphenotype", "Display top significant phenotype pairs in heatmap",
                             choices=c("10"="10", "25"="25", "50"="50", "100"="100", "250"="250", "500"="500", "1,000"="1000", "all"="all"),
                             inline=T, selected="25"),
                style="display:inline-block; vertical-align:top; margin-top:30px; margin-left:75px"
              )
            ),
            style="height:90px; margin-left:-15px"
          ),
          # Review CPAG results
          div(
            tabsetPanel(id="tabsetReviewResults",
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
            style="width=100%; margin-top:20px"
          )
        ),

        # Explore CPAG associations
        tabPanel(title="Explore iCPAGdb associations", value="tabPanelExplore",
          div(
            # Prompts
            HTML("<br>"),
            sidebarPanel(width=12,
              HTML("<b>1. Compute</b><br><br>"),
              div(
                selectInput(inputId="exploreSource1", label="GWAS source one", choices=gwasSource, width="85%"),
                style="display:inline-block; vertical-align:top; margin-top:10px; width:12%"
              ),
              div(
                selectInput(inputId="exploreSource2", label="GWAS source two", choices=gwasSource, width="85%"),
                style="display:inline-block; vertical-align:top; margin-top:10px; width:12%"
              ),
              div(
                sliderInput(inputId="explorePthresh1", HTML("p-threshold<sub>1</sub> (1X10<sup>-<it>x</it></sup>)"),
                            min=3, max=20, step=1, value=7, width="85%"),
                style="display:inline-block; vertical-align:top; margin-top:0px; width:12%"
              ),
              div(
                sliderInput(inputId="explorePthresh2", HTML("p-threshold<sub>2</sub> (1X10<sup>-<it>x</it></sup>)"),
                            min=3, max=20, step=1, value=7, width="85%"),
                style="display:inline-block; vertical-align:top; margin-top:0px; width:12%"
              ),
              div(
                HTML("<b>p-threshold adjustments</b><br>H2P2 min = 1X10<sup>-5</sup><br>NHGRI min = 5X10<sup>-8</sup><br>All others, min = 1X10<sup>-3</sup>"),
                style="display:inline-block; vertical-align:top; margin-top:-5px; width:12%"
              ),
              div(
                radioButtons(inputId="exploreLDpop", "LD 1000 Genomes population",
                             choices=c("European"="EUR", "African"="AFR", "Asian"="EAS"), inline=T),
                style="display:inline-block; vertical-align:top; margin-top:20px; width:15%"
              ),
              # Exec button
              div(
                actionButton(inputId="exploreCompute", "Compute CPAG",
                             style="color:white; background:linear-gradient(#54b4eb, #2fa4e7 60%, #0088dd)"),
                style="display:inline-block; vertical-align:top; margin-top:40px; width:12%"
              ),
              HTML("<hr style='height:1px;color:black;background-color:black'>"),
              HTML("<b>2. Filter</b><br>"),
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
                actionButton("exploreSelectionFilterApply", "Apply filters", width="100px",
                             style="color:white; background:linear-gradient(#54b4eb, #2fa4e7 60%, #0088dd)"),
                style="display:inline-block; vertical-align:top; margin-top:25px"
              ),
              div(
                actionButton("exploreSelectionFilterClear", "Clear filters", width="100px",
                             style="color:white; background:linear-gradient(#54b4eb, #2fa4e7 60%, #0088dd)"),
                style="display:inline-block; vertical-align:top; margin-top:25px; margin-left:10px"
              ),
              div(
                downloadButton("exploreResultsDownload", "Download filtered records",
                                style="color:white; background:linear-gradient(#54b4eb, #2fa4e7 60%, #0088dd)"),
                style="display:inline-block; vertical-align:top; margin-top:25px; margin-left:10px"
              ),
              HTML("<br>"),
              div(
                radioButtons("exploreSelectionHeatmapMetric", "Heatmap metric",
                             choices=c("Fisher"="P_fisher", "Bonferroni"="Padj_Bonferroni", "FDR"="Padj_FDR",
                                       "Jaccard"="Jaccard", "Chao-Sorensen"="ChaoSorensen"),
                             inline=T, selected="P_fisher"),
                style="display:inline-block; vertical-align:top; margin-top:20px; width:24%"
              ),
              div(
                radioButtons("exploreSelectionHeatmapNphenotype", "Display top significant phenotype pairs in heatmap",
                             choices=c("10"="10", "25"="25", "50"="50", "100"="100", "250"="250", "500"="500", "1,000"="1000", "all"="all"),
                             inline=T, selected="25"),
                style="display:inline-block; vertical-align:top; margin-top:20px; width:24%"
              )
            ),
            style="margin-left:-15px"
          ),
          # Explore CPAG results
          div(
            tabsetPanel(id="tabsetExploreResults",
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
          )
        ),

        # Upload and compute CPAG
        tabPanel(title="Upload GWAS and compute CPAG", value="userComputeGWAS",
          div(
            # Prompts
            HTML("<br>"),
            sidebarPanel(width=12,
              HTML("<b>1. Upload a GWAS file</b><br><br>"),
              div(
                fileInput("userComputeBrowseFile", "Choose File", multiple=FALSE,
                          accept=c("text/csv", "text/comma-separated-values,text/plain", ".csv"), width="85%"),
                style="display:inline-block; vertical-align:top; margin-top:10px; width:24%"
              ),
              div(
                radioButtons("userComputeDelimiter", "Delimiter", choices=c("Comma", "Tab"), inline=T, selected="Comma"),
                style="display:inline-block; vertical-align:top; margin-top:10px; margin-left:5px; width:12%"
              ),
              div(
                disabled(textInput(inputId="userComputePhenotypeCol", label="Trait (phenotype) column",
                         value="One phenotype per file", width="85%")),
                style="display:inline-block; vertical-align:top; margin-top:10px; width:12%"
              ),
              div(
                textInput(inputId="userComputeSNPcol", label="SNP column", width="85%"),
                style="display:inline-block; vertical-align:top; margin-top:10px; width:12%"
              ),
              div(
                textInput(inputId="userComputePcol", label="P (significance) column", width="68%"),
                style="display:inline-block; vertical-align:top; margin-top:10px; width:15%"
              ),
              # Upload button
              div(
                actionButton(inputId="userComputeUploadFile", "Upload file", width="120px",
                             style="color:white; background:linear-gradient(#54b4eb, #2fa4e7 60%, #0088dd)"),
                style="display:inline-block; vertical-align:top; margin-top:40px; width:12%"
              ),
              HTML("<hr style='height:1px;color:black;background-color:black'>"),
              HTML("<b>2. Compute</b><br><br>"),
              div(
                selectInput(inputId="userComputeSource1", label="GWAS source one", choices="User Supplied GWAS", width="85%"),
                style="display:inline-block; vertical-align:top; margin-top:10px; width:12%"
              ),
              div(
                selectInput(inputId="userComputeSource2", label="GWAS source two", choices=gwasSource, width="85%"),
                style="display:inline-block; vertical-align:top; margin-top:10px; width:12%"
              ),
              div(
                sliderInput(inputId="userComputePthresh1", HTML("p-threshold<sub>1</sub> (1X10<sup>-<it>x</it></sup>)"),
                            min=3, max=20, step=1, value=7, width="85%"),
                style="display:inline-block; vertical-align:top; margin-top:0px; width:12%"
              ),
              div(
                sliderInput(inputId="userComputePthresh2", HTML("p-threshold<sub>2</sub> (1X10<sup>-<it>x</it></sup>)"),
                            min=3, max=20, step=1, value=7, width="85%"),
                style="display:inline-block; vertical-align:top; margin-top:0px; width:12%"
              ),
              div(
                HTML("<b>p-threshold adjustments</b><br>H2P2 min = 1X10<sup>-5</sup><br>NHGRI min = 5X10<sup>-8</sup><br>All others, min = 1X10<sup>-3</sup>"),
                style="display:inline-block; vertical-align:top; margin-top:-5px; width:12%"
              ),
              div(
                radioButtons(inputId="userComputeLDpop", "LD 1000 Genomes population",
                             choices=c("European"="EUR", "African"="AFR", "Asian"="EAS"), inline=T),
                style="display:inline-block; vertical-align:top; margin-top:20px; width:15%"
              ),
              # Exec button
              div(
                actionButton(inputId="userComputeCompute", "Compute CPAG", width="120px",
                             style="color:white; background:linear-gradient(#54b4eb, #2fa4e7 60%, #0088dd)"),
                style="display:inline-block; vertical-align:top; margin-top:40px; width:12%"
              ),
              HTML("<hr style='height:1px;color:black;background-color:black'>"),
              HTML("<b>3. Filter</b><br>"),
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
                actionButton("userComputeSelectionFilterApply", "Apply filters", width="100px",
                             style="color:white; background:linear-gradient(#54b4eb, #2fa4e7 60%, #0088dd)"),
                style="display:inline-block; vertical-align:top; margin-top:25px"
              ),
              div(
                actionButton("userComputeSelectionFilterClear", "Clear filters", width="100px",
                             style="color:white; background:linear-gradient(#54b4eb, #2fa4e7 60%, #0088dd)"),
                style="display:inline-block; vertical-align:top; margin-top:25px; margin-left:10px"
              ),
              div(
                downloadButton("userComputeResultsDownload", "Download filtered records",
                                style="color:white; background:linear-gradient(#54b4eb, #2fa4e7 60%, #0088dd)"),
                style="display:inline-block; vertical-align:top; margin-top:25px; margin-left:10px"
              ),
              HTML("<br>"),
              div(
                radioButtons("userComputeSelectionHeatmapMetric", "Heatmap metric",
                             choices=c("Fisher"="P_fisher", "Bonferroni"="Padj_Bonferroni", "FDR"="Padj_FDR",
                                       "Jaccard"="Jaccard", "Chao-Sorensen"="ChaoSorensen"),
                             inline=T, selected="P_fisher"),
                style="display:inline-block; vertical-align:top; margin-top:20px; width:24%"
              ),
              div(
                radioButtons("userComputeSelectionHeatmapNphenotype", "Display top significant phenotype pairs in heatmap",
                             choices=c("10"="10", "25"="25", "50"="50", "100"="100", "250"="250", "500"="500", "1,000"="1000", "all"="all"),
                             inline=T, selected="25"),
                style="display:inline-block; vertical-align:top; margin-top:20px; width:24%"
              )
            ),
            style="margin-left:-15px"
          ),
          # User compute tab CPAG results
          div(
            tabsetPanel(id="tabsetUserComputeResults",
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
          )
        ),

        # Bibliography
        tabPanel(title="Bibliography", value="tabPanelBibliography",
          div(
            HTML("<ul>
                    <br>
                    <li>
                      Buniello, MacArthur et al, 2019 <a href=https://academic.oup.com/nar/article/47/D1/D1005/5184712 target=_blank>The NHGRI-EBI GWAS Catalog of published genome-wide association studies, targeted arrays and summary statistics 2019</a>. Nucleic Acids Research, 2019, Vol. 47 (Database issue): D1005-D1012.
                    </li>
                    <br>
                    <li>
                      Severe Covid-19 GWAS Group, Ellinghaus D, Degenhardt F, et al. <a href=https://pubmed.ncbi.nlm.nih.gov/32558485/ target=_blank>Genomewide Association Study of Severe Covid-19 with Respiratory Failure</a>. N Engl J Med. 2020;383(16):1522-1534. doi:10.1056/NEJMoa2020283
                    </li>
                    <br>
                    <li>
                      Raffler, J., Friedrich, N., Arnold, M., Kacprowski, T., Rueedi, R., Altmaier, E., Bergmann, S., Budde, K., Gieger, C., Homuth, G., et al. (2015). <a href=https://pubmed.ncbi.nlm.nih.gov/26352407/ target=_blank>Genome-Wide Association Study with Targeted and Non-targeted NMR Metabolomics Identifies 15 Novel Loci of Urinary Human Metabolic Individuality</a>. PLoS Genet 11, e1005487
                    </li>
                    <br>
                    <li>
                      Shin, S.Y., Fauman, E.B., Petersen, A.K., Krumsiek, J., Santos, R., Huang, J., Arnold, M., Erte, I., Forgetta, V., Yang, T.P., et al. (2014). <a href=https://www.nature.com/articles/ng.2982 target=_blank>An atlas of genetic influences on human blood metabolites</a>. Nat Genet 46, 543-550
                    </li>
                    <br>
                    <li>
                      Wang, L., Oehlers, S.H., Espenschied, S.T., Rawls, J.F., Tobin, D.M., and Ko, D.C. (2015). <a href=https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0722-1 target=_blank>CPAG: software for leveraging pleiotropy in GWAS to reveal similarity between human traits links plasma fatty acids and intestinal inflammation</a>. Genome Biol 16, 190
                    </li>
                    <br>
                    <li>
                      Wang, L., Pittman, K.J., Barker, J.R., Salinas, R.E., Stanaway, I.B., Williams, G.D., Carroll, R.J., Balmat, T., Ingham, A., Gopalakrishnan, A.M., et al. (2018). <a href=https://www.cell.com/cell-host-microbe/fulltext/S1931-3128(18)30377-9 target=_blank>An Atlas of Genetic Variation Linking Pathogen-Induced Cellular Traits to Human Disease</a>. Cell Host Microbe 24, 308-323 e306
                    </li>
                    <br>
                    <li>
                      Welter, D., MacArthur, J., Morales, J., Burdett, T., Hall, P., Junkins, H., Klemm, A., Flicek, P., Manolio, T., Hindorff, L., et al. (2014). <a href=https://pubmed.ncbi.nlm.nih.gov/24316577/ target=_blank>The NHGRI GWAS Catalog, a curated resource of SNP-trait associations</a>. Nucleic Acids Res 42, D1001-1006
                    </li>
                  </ul>"
                ),
            style="margin-top:15px; margin-left:20px"
          )
        )

      ),

      # Feature enable element
      #div(
      #  textInput("featureEnable", ""),
      #  style="width:50px"
      #),

      style="margin-left: 20px"
    )
  )
)
