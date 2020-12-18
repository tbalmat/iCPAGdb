# Duke University Cross Phenotype Analysis of GWAS
# Explore page, Novenber 2020

# Launch app ui.r and server.r from specified appDir
# Note the specification of a tcp port that the process will listen on for http requests

library("shiny")

# Specify directory containing ui.r and server.r

# Local dir
ad <- "C:/Projects/Duke/H2P2GenomeWideAssociationStudy/CPAG/iCPAGdb/App/App-Devel"

# Execute 
runApp(appDir=ad,
       launch.browser=T,
       host = getOption("shiny.host", "127.0.0.1"),
       port=4295,
       display.mode = c("auto", "normal", "showcase"))