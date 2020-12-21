options(max.print=1000)      # number of elements, not rows
options(stringsAsFactors=F)
options(scipen=999999)
options(device="windows")

##################################################################################################
# Test python execution
##################################################################################################

setwd("C:\\Projects\\Duke\\H2P2GenomeWideAssociationStudy\\CPAG\\iCPAGdb\\App-Devel\\R-Devel")
dir()

#pyexe <- "\"C:\\Users\\Kyung Soon\\AppData\\Local\\Microsoft\\WindowsApps\\python.exe\""
pyexe <- "\"C:\\Users\\Kyung Soon\\AppData\\Local\\Programs\\Python\\Python37\\python.exe\""

script <- "hello.py"
system(paste(pyexe, " ", script, sep=""), intern=T)

##################################################################################################
# Install python packages
##################################################################################################

# Upgrade pip
system(paste(pyexe, " -m pip install --upgrade pip", sep=""), intern=T)

# joblib
system(paste(pyexe, " -m pip install joblib", sep=""), intern=T)

# tqdm
system(paste(pyexe, " -m pip install tqdm", sep=""), intern=T)

# Install past version of scipy that avoids WinError 126 ("specified module cannot be found" in _dlopen()) 
#system(paste(pyexe, " -m pip install scipy==1.4.1", sep=""), intern=T)
system(paste(pyexe, " -m pip install scipy", sep=""), intern=T)

# Downgrade numpy to a version known to function with python 3.7 on Windows (also necessary to avoid WinError 126)
system(paste(pyexe, " -m pip install numpy==1.17", sep=""), intern=T)

#system(paste(pyexe, " -m pip install .qhull", sep=""), intern=T)
#system(paste(pyexe, " -m pip install _distributor_init", sep=""), intern=T)

# pandas
system(paste(pyexe, " -m pip install pandas", sep=""), intern=T)

# statsmodels
system(paste(pyexe, " -m pip install statsmodels", sep=""), intern=T)

# List installed package versions
system(paste(pyexe, " -m pip list", sep=""), intern=T)

##################################################################################################
# Test python imports
##################################################################################################

write("import numpy", "import.py")
write("import scipy", "import.py")
write("import scipy.stats", "import.py")
write("from scipy.spatial.distance import cdist", "R\\import.py")
system(paste(pyexe, " import.py", sep=""), intern=T)

##################################################################################################
# SQLite queries
##################################################################################################

#install.packages("RSQLite")
library(RSQLite)

db <- dbConnect(RSQLite::SQLite(), "../pyCPAG/db/cpag_gwasumstat_v1.1.db")
dbGetQuery(db, "select * from sqlite_master")
dbGetQuery(db, "select source, count(1) as n from GWAStable group by source")
dbGetQuery(db, "select * from GWAStable where SNP in(select distinct SNP from GWASTable limit 20)")
dbGetQuery(db, "select * from GWAStable where source='NHGRI' limit 20")
hist(log(dbGetQuery(db, "select * from GWAStable where source='NHGRI' and pval<1e-100")[,"pval"])/log(10))
hist(log(dbGetQuery(db, "select * from GWAStable where source='BloodMetabolites' and pval<1e-5")[,"pval"])/log(10))
hist(log(dbGetQuery(db, "select * from GWAStable where source='BloodXenobiotic' and pval<1e-5")[,"pval"])/log(10))
dbDisconnect(db)

db <- dbConnect(RSQLite::SQLite(), "../../pyCPAG/db/cpag_gwasumstat_v1.1_EUR.ld0.4.db")
dbGetQuery(db, "select * from sqlite_master")
dbGetQuery(db, "select * from LDSNP limit 100")
dbDisconnect(db)

##################################################################################################
# CPAG
##################################################################################################

setwd("C:\\Projects\\Duke\\H2P2GenomeWideAssociationStudy\\CPAG\\iCPAGdb\\App-Devel\\pyCPAG")

fout <- c("output/H2P2-p1e-05-EUR.csv", "output/H2P2-p1e-05-AFR.csv", "output/H2P2-p1e-05-EAS.csv",
          "output/NHGRI-p5e-08-EUR.csv", "output/H2P2-p1e-05-NHGRI-p5e-08-EUR.csv", "output/cpag.csv")[4]

# Within a single subtype
#pycmd <- paste(pyexe, " main.py cpagdb --threads 2 --subtype H2P2 --H2P2-Pcut 1e-5 --lddb-pop EUR --outfile ", fout, sep="")
subtype <- c("BloodMetabolites", "BloodXenobiotic", "H2P2", "NHGRI", "UrineMetabolites", "mol_gwas", "clin_gwas", "x")[7]
pycmd <- paste(pyexe, " main.py cpagdb --threads 2 --subtype ", subtype, " --Pcut 1e-100 --outfile ", fout, sep="")

# Between two subtypes
pycmd <- paste(pyexe, " main.py cpagdb --threads 2 --subtype H2P2 --subtype NHGRI --H2P2-Pcut 1e-7 --Pcut 1e-100 --outfile ", fout, sep="")
pycmd <- paste(pyexe, " main.py cpagdb --threads 2 --subtype H2P2 --subtype BloodCytokine --H2P2-Pcut 1e-7 --Pcut 1e-10 --outfile ", fout, sep="")

# User supplied GWAS and CPAG DB traits
fin <- c("Okbay-SWB.csv", "Okbay-Neuroticism.csv", "GPC-2-Extraversion.csv", "LowAtrialFibrillation.csv",
         "VanderKalkBirthLength.csv", "top_EllinghausPCs_covid19.csv")[6]
x <- read.table(paste("C:\\Projects\\Duke\\H2P2GenomeWideAssociationStudy\\CPAG\\ExtGWAS\\GRASP\\", fin, sep=""), header=T, sep=",")
x[,"p"] <- 10**(-x[,"log_10_p"])
hist(x[,"log_10_p"])
write.table(x, paste("input\\", fin, sep=""), row.names=F, col.names=T, sep=",", quote=F)
snpcol <- c("rsID", "avsnp150")[2]
pcol <- c("p", "p_value")[2]

pycmd <- paste(pyexe,
               " main.py usr-gwas --threads 2",
               " --infile \"input\\", fin, "\" --delimitor \",\" --usr-pcut 1e-5 --SNPcol \"", snpcol, "\" --Pcol \"", pcol, "\"",
               #" --querydb mol_gwas",
               " --cpagdb-pcut 5e-8 --H2P2-Pcut 1e-5 --ld-clump 1",
               " --outfile ", fout, sep="")

#dir("output")
system(paste("rm ", fout, sep=""), intern=T)
t <- proc.time()
y <- suppressWarnings(system(pycmd, intern=T))
proc.time()-t
attributes(y)

k <- length(y)
y[k]
k <- k-1

z <- read.table(fout, header=T, sep=",", quote="\"")

##################################################################################################
# CPAG help
##################################################################################################

pycmd <- paste(pyexe, " main.py -h", sep="")
system(pycmd, intern=T)

##################################################################################################
# Examine uniqueness of shared SNPs
##################################################################################################

setwd("C:\\Projects\\Duke\\H2P2GenomeWideAssociationStudy\\CPAG\\iCPAGdb\\App-Devel\\pyCPAG\\output")
dir()
fout <- c("H2P2-p1e-05-EUR.csv", "H2P2-p1e-05-AFR.csv", "H2P2-p1e-05-EAS.csv",
          "NHGRI-p5e-08-EUR.csv", "H2P2-p1e-05-NHGRI-p5e-08-EUR.csv")[4]
x <- read.table(fout, header=T, sep=",", quote="\"")
nrow(x)

length(unique(x[,"SNPshare_all"]))
y <- strsplit(x[,"SNPshare_all"], "|", fixed=T)

# Identify phenotypes appearing in both traits 1 and 2
z <- unlist(lapply(unique(x[,"Trait1"]), function(t) length(which(x[,"Trait2"]==t))))
table(z)

k <- which(x[,"Trait1"]=="body height")
k <- k[order(x[k,"Jaccard"])]
k <- order(x[,"P_fisher"])[1:50]
z <- data.frame(x[k,c("Trait1", "Trait2")], -log(x[k,"P_fisher"])/log(10), x[k,c("Jaccard", "ChaoSorensen")])

##################################################################################################
# Examine relationships between Fisher-P, Jaccard, and Chao-Sorensen statistics
##################################################################################################

setwd("C:\\Projects\\Duke\\H2P2GenomeWideAssociationStudy\\CPAG\\iCPAGdb\\App-Devel\\pyCPAG\\output")
dir()
fout <- c("H2P2-p1e-05-EUR.csv", "H2P2-p1e-05-AFR.csv", "H2P2-p1e-05-EAS.csv",
          "NHGRI-p5e-08-EUR.csv", "H2P2-p1e-05-NHGRI-p5e-08-EUR.csv")[4]
x <- read.table(fout, header=T, sep=",", quote="\"")
nrow(x)

k <- sample(which(x[,"P_fisher"]<1e-50), 1000, replace=F)
#dev.new()
plot(x=x[k,"P_fisher"], y=x[k,"Jaccard"])

##################################################################################################
# Chao-Sorensen, Jaccard similarity indices (have read references relating these to Kruskal-Wallis tests)
##################################################################################################

install.packages("CommEcol")
library(CommEcol)

# From docs:
aa <- c(1, 2, 4, 5)
bb <- c(1, 2, 0, 5)
cc <- c(0, 2, 3, 3)
dat3 <- rbind(aa, bb, cc) 
colnames(dat3) <- c("sp1", "sp2", "sp3", "sp4")
dat3

vegdist(dat3, method='jaccard', binary=FALSE) 
# Notice dissimilarity between the pair aa-bb is the same of the pair aa-cc. 
# In fact, bb and cc differ from aa in the same way (one species, and 
# 4 exclusive individuals).

dis.chao(dat3, index=c("jaccard", "sorensen")[1], version="prob")
# The probability version of the Chao index, however, produce different 
# dissimilarities for the pairs aa-bb and aa-cc 
# (aa-cc is less dissimilar than aa-bb).

# Reduced values indicate increased similarity (for both Jaccard and Sorensen methods
aa <- c(1, 12, 3, 4)
bb <- c(1, 0, 3, 5)
cc <- c(0, 2, 3, 7)
dat3 <- rbind(aa, bb, cc) 
colnames(dat3) <- c("sp1", "sp2", "sp3", "sp4")
dat3

remove.packages("CommEcol")

##################################################################################################
# Heatmap
##################################################################################################

install.packages("heatmaply")
install.packages("reshape")
library(heatmaply)
library(reshape)

setwd("C:/Projects/Duke/H2P2GenomeWideAssociationStudy/CPAG/iCPAGdb/App-Devel/R-Devel")
dir("../pyCPAG/output")
z <- read.table("../pyCPAG/output/H2P2-p1e-05-AFR.csv", header=T, sep=",", quote="\"")
z <- z[which(regexpr("rs105", z[,"SNPshare_all"])>0),]
z <- read.table("../pyCPAG/output/top_EllinghausPCs_covid19-5286078011-p1e-07-NHGRI-p5e-08-EUR.csv", header=T, sep=",", quote="\"")

min(z[,"P_fisher"])
length(which(z[,"P_fisher"]<1e-100))

z2 <- as.matrix(cast(z, Trait1~Trait2, value="P_fisher"))
class(z2)
dim(z2)
rownames(z2)
colnames(z2)
z2[1,1]
z2 <- -log(z2)/log(10)
k <- which(is.infinite(z2))
z2[k] <- max(z2[-k], 100, na.rm=T)

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

#z2[which(is.na(z2))] <- 0
z3 <- as.matrix(dist(z2))
z3[which(is.na(z3))] <- 0
hclust(z3)
x <- dist(z2[8:10, 8:10])
attributes(x)
as.matrix(y)

x <- matrix(c(1, 2, 3, 1, 2, 3, 1, NA, 3), nrow=3)
y <- dist(x)
attributes(y)
str(y)

x <- z2

# Alternative to dist() that omits NA values
naDist <- function(x) {
  if(dim(x)[1]>0) {
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
    # rownames(y) <- rownames(x)
    # y <- as.dist(y)
    y <- y[lower.tri(y, diag=F)]
    attributes(y) <- list("Size"=dim(x)[1], "Diag"=F, "Upper"=F, "method"="euclidian", "Labels"=rownames(x))
  } else {
    y <- dist(0)
  }
  # Replace NA values with the maximum observed distance
  y[which(is.na(y))] <- max(y, na.rm=T)
  # Establish class after all assignments, since it may have been altered
  class(y) <- "dist"
  return(y)
}

y <- naDist(z2)
class(y) 
str(y)
length(y)
hclust(y)

rownames(y) <- rownames(z2)
class(as.dist(y))
str(as.dist(y))

#install.packages("usedist")
#library(usedist)
#edit(dist_make)

heatmaply(z2, file="heatmaply_plot.html", plot_method="plotly", distfun=naDist, na.value="gray50")
browseURL("heatmaply_plot.html")

#library(ggplot2)
#ggplot() + geom_line(aes(x=1:10, y=10:10))

x <- matrix(1:100, nrow=10)
rownames(x) <- 1:10
colnames(x) <- 1:10
print(heatmaply(x))

g3 <- heatmaply(z2,
          file="heatmaply_plot.html",
          colors = viridis(n=256, alpha=1, begin=1, end=0, option="viridis"),
          # na.value does not appear to have any effect
          na.value="gray50",
          plot_method=c("plotly", "ggplot")[2],
          dendrogram=c("both", "none")[2],
          #distfun=naDist,
          branches_lwd=0.2,
          #grid_gap=1,
          # grid_color does not appear to have any effect when plot_method="plotly"
          grid_color="gray95",
          # Axis label fonts configured in theme
          #fontsize_row=8,
          #fontsize_col=8,
          key.title="-log10(Fisher-p)",
          heatmap_layers=list(
                           ggplot2::theme(
                             #plot.title=element_text(size=12, hjust=0.5),
                             #plot.caption=element_text(size=12, hjust=0.5),
                             #panel.background=element_rect(fill="gray75"),
                             #panel.grid.major.x=element_blank(),
                             #panel.grid.major.y=element_blank(),
                             #panel.grid.minor=element_blank(),
                             panel.border=element_rect(fill=rgb(0.8, 0.8, 0.8, alpha=0.2), color="black"),
                             axis.title.x=element_text(size=10),
                             axis.title.y=element_text(size=10),
                             axis.text.x=element_text(size=8, angle=45, hjust=1, vjust=0.5),
                             axis.text.y=element_text(size=8),
                             #legend.position="bottom",
                             #legend.background=element_rect(color="gray"),
                             #legend.key=element_rect(fill="white"),
                             #legend.box="horizontal",
                             legend.text=element_text(size=8),
                             legend.title=element_text(size=10))),
         )

browseURL("heatmaply_plot.html")


##################################################################################################
# heatmaply examples from docs
##################################################################################################

# different colors
g <- heatmaply(mtcars, colors = heat.colors(200), width="800px")
# using special scale_fill_gradient_fun colors
heatmaply(mtcars, scale_fill_gradient_fun = scale_color_gradient())
dev.off()
heatmaply(mtcars, file="heatmaply_plot.html")
browseURL("heatmaply_plot.html")

##################################################################################################
# plot_ly()
##################################################################################################

g2 <- plot_ly(z = ~volcano, type = "surface", width=800, height=400)
print(g2)

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

##################################################################################################
# Heatmap that will not render
# Truncarion of long axis labels appears to correct the problem
# Also, on some platforms, heatmaps with lengthy labels will render, but with some x or y
# coordinates missing (yikes!)
##################################################################################################

library(heatmaply)
library(reshape)

setwd("C:/Projects/Duke/H2P2GenomeWideAssociationStudy/CPAG/iCPAGdb/App-Devel/R-Devel")
dir("../pyCPAG/output")
z <- read.table("../pyCPAG/output/top_EllinghausPCs_covid19-5286078011-p1e-07-NHGRI-p5e-08-EUR.csv", header=T, sep=",", quote="\"")

#q <- data.frame(z[,"Trait2"], "w"=0)

nLabTrunc <- 20
#z2 <- as.matrix(cast(data.frame("Trait1"=truncLabel(z[,"Trait1"], nLabTrunc),
#                                "Trait2"=truncLabel(z[,"Trait2"], nLabTrunc),
#                                "z"=z[,c("P_fisher", "Padj_FDR", "Padj_Bonferroni", "Jaccard", "ChaoSorensen")[2]]),
#                     Trait1~Trait2, value="z"))

z2 <- xyzMatrix(truncLabel(z[,"Trait1"], nLabTrunc),
                truncLabel(z[,"Trait2"], nLabTrunc),
                z[,c("P_fisher", "Padj_FDR", "Padj_Bonferroni", "Jaccard", "ChaoSorensen")[2]])

z2 <- -log(z2)/log(10)
k <- which(is.infinite(z2))
z2[k] <- 350
dim(z2)
heatmaply(z2,
          file="heatmaply_plot.html",
          plot_method=c("plotly", "ggplot")[1],
          dendrogram=c("both", "none")[2])
browseURL("heatmaply_plot.html")

##################################################################################################
# Misc
##################################################################################################

dir(pattern="^stats.py$")
dir("output", pattern="H2P2-H2P2-p1e-05-p1e-05-EUR.csv")

##################################################################################################
# Review post-analysis ontology
# Empirical use of the post_analysis option of the CPAG function indicates that experimental
# factor ontology is available for NHGRI phenotypes only
##################################################################################################

setwd("C:\\Projects\\Duke\\H2P2GenomeWideAssociationStudy\\CPAG\\iCPAGdb\\App-Devel\\pyCPAG")

fin <- c("output/NHGRI-p5e-08-EUR.csv", "output/H2P2-p1e-05-EUR.csv", "output/cpag.csv")[1]
fout <- c(fin, sub(".csv", "-Ontology.csv", fin, fixed=T))[1]

pycmd <- paste(pyexe, " main.py post_analysis --anno-ontology --anno-cols Trait2 --infile ", fin,
               " --outfile ", fout, sep="")

dir("output")
system(paste("rm ", fout, sep=""), intern=T)
t <- proc.time()
y <- suppressWarnings(system(pycmd, intern=T))
proc.time()-t
attributes(y)

k <- length(y)
y[k]
k <- k-1

z1 <- read.table(fout, header=T, sep=",", quote="\"")
z2 <- read.table(paste("db/gwas-efo-trait-mappings_download_", as.character(Sys.Date()), ".txt", sep=""), header=T, sep="\t", quote="\"")

z3 <- unique(z2[,"Parent.term"])
z4 <- sort(unique(z2[,"EFO.term"]))

z1[,"Trait1"] <- tolower(z1[,"Trait1"])
z1[,"Trait2"] <- tolower(z1[,"Trait2"])
z2[,"Disease.trait"] <- tolower(z2[,"Disease.trait"])
z2[,"EFO.term"] <- tolower(z2[,"EFO.term"])

z2[which(z2[,"Disease.trait"]==z1[2,"Trait1"]),]

# Join composed results (with ontology) to downloaded ontology dictionary
# Traits exist in both EFO.term and Disease.trait columns of the dictionary
# However, more matches are made (unique traits) with EFO.term than with Disease.trait
# This indicates that the post_analysis algorithm matches the specified data column to EFO.term
z5 <- merge(data.frame("i"=1:nrow(z1), z1[,c("Trait1", "Trait2", "Trait2_EFO")]),
            data.frame(z2[,c("Disease.trait", "EFO.term", "EFO.URI", "Parent.term")]),
            by.x="Trait2", by.y=c("EFO.term", "Disease.trait")[1])
length(unique(z5[,"Trait2"]))
length(unique(z1[,"Trait2"]))

# Do all traits have a corresponding EFO?
table(nchar(z1[,"Trait2_EFO"]), useNA="always")
table(nchar(z1[,"Trait1_EFO"]), useNA="always")

# Inspect phenotypes with misssing EFO
k <- which(is.na(z1[,"Trait2_EFO"]))
x1 <- z1[k,]
sort(unique(z1[k,"Trait2"]))
unique(z1[which(z1[,"Trait2"] %in% c("coronary artery disease", "kidney stone")),"Trait2_EFO"])

# Identify observations with compound EFOs
k <- which(regexpr(",", z1[,"Trait2_EFO"], fixed=T)>0)
length(k)
head(z1[k,"Trait2_EFO"])

# Evaluate EFO sources (leading characters of URLs)
# Parse compund EFOs
x1 <- trimws(unlist(strsplit(z1[,"Trait2_EFO"], ",")))
head(x1)
# Omit text following final slash in each URL
i <- gregexpr("/", x1, fixed=T)
unique(unlist(lapply(1:length(x1), function(j) substring(x1[j], 1, i[[j]][length(i[[j]])]))))
unique(x1[which(substring(x1, 1, 4)!="http")])

# Inspect NA EFOs
x1[9988:9999]
z1[which(regexpr("http://www.ebi.ac.uk/efo/EFO_0006925", z1[,"Trait2_EFO"], fixed=T)>1),"Trait2_EFO"]

# Attempt to identify source of empty URLs
which(substring(gsub(" ", "", z1[,"Trait2_EFO"]), 1, 1)==",")
unique(substring(z1[,"Trait2_EFO"], nchar(z1[,"Trait2_EFO"]), nchar(z1[,"Trait2_EFO"])))
# NA is the culprit
z1[which(is.na(substring(z1[,"Trait2_EFO"], nchar(z1[,"Trait2_EFO"]), nchar(z1[,"Trait2_EFO"])))), "Trait2_EFO"]
which(regexpr(",,", gsub(" ", "", z1[,"Trait2_EFO"]), fixed=T)>0)

# Inspect EFOs from the three non-empty URL sources
for(url in c("http://www.ebi.ac.uk/efo/", "http://purl.obolibrary.org/obo/", "http://www.orpha.net/ORDO/"))
  print(z1[which(substring(z1[,"Trait2_EFO"], 1, nchar(url))==url)[1:10],"Trait2_EFO"])

# Retrieve EFOs from GWAS table
# At time of development, experimental factor ontology data are available for NHGRI phenotypes only
library(RSQLite)
db <- dbConnect(RSQLite::SQLite(), "../pyCPAG/db/cpag_gwasumstat_v1.1.db")
z7 <- dbGetQuery(db, "select distinct trait, efo from GWAStable")
dbGetQuery(db, "select distinct source from GWAStable where efo is not null")
dbDisconnect(db)

dbGetQuery(db, "select * from GWAStable limit 1")

##################################################################################################
# Compose CPAG results with EFO columns for NHGRI phenotypes
# Specify Trait1 and Trait2 for within NHGRI phenotypes
# Specify trait orresponding to NHGRI phenotypes for other combinations
# Original CPAG files are modified and replaced, with EFO columns appended
##################################################################################################

setwd("C:\\Projects\\Duke\\H2P2GenomeWideAssociationStudy\\CPAG\\iCPAGdb\\App-Devel\\pyCPAG")
dir("output")

f <- c("output/NHGRI-p5e-08-EUR.csv", "output/H2P2-p1e-05-NHGRI-p5e-08-EUR.csv",
       "output/H2P2-p1e-05-NHGRI-p5e-08-AFR.csv", "output/H2P2-p1e-05-NHGRI-p5e-08-EAS.csv")[4]

pycmd <- paste(pyexe, " main.py post_analysis --anno-ontology --anno-cols Trait2 --infile ", f, " --outfile ", f, sep="")
y <- suppressWarnings(system(pycmd, intern=T))
attributes(y)

x <- read.table(f, header=T, sep=",", quote="\"")

##################################################################################################
# Search EFOs
##################################################################################################

setwd("C:\\Projects\\Duke\\H2P2GenomeWideAssociationStudy\\CPAG\\iCPAGdb\\App-Devel\\pyCPAG")
f <- c("output/NHGRI-p5e-08-EUR.csv", "output/H2P2-p1e-05-EUR.csv", "output/cpag.csv")[1]
x <- read.table(f, header=T, sep=",", quote="\"")
grep("(EFO_0004842)", x[,"Trait1_EFO"])

##################################################################################################
# Enumerate phenotypes and EFOs when compounded
##################################################################################################

setwd("C:\\Projects\\Duke\\H2P2GenomeWideAssociationStudy\\CPAG\\iCPAGdb\\App-Devel\\pyCPAG")
f <- c("output/NHGRI-p5e-08-EUR.csv", "output/H2P2-p1e-05-EUR.csv", "output/cpag.csv")[1]
x <- read.table(f, header=T, sep=",", quote="\"")

# Compare the frequency of phenotypes to their corresponding EFO
y <- rbind(setNames(x[,c("Trait1", "Trait1_ParentTerm")], c("ph", "efo")),
           setNames(x[,c("Trait2", "Trait2_ParentTerm")], c("ph", "efo")))
nrow(y)
y <- y[!duplicated(y),]
nrow(y)

nph <- unlist(lapply(strsplit(y[,"ph"], ","), length))
nefo <- unlist(lapply(strsplit(y[,"efo"], ","), length))
length(which(nph==1 & nefo==1))
k <- which(nph>1 & nefo==1)
length(k)
y[k,]
length(which(nph==1 & nefo>1))

k <- which(y[,"ph"] %in% c("coronary artery disease", "kidney stone"))
length(k)

k <- which(nph>1 & nefo==1)
z <- data.frame(nph[k], nefo[k], y[k,c("ph","efo")], 0)

write.table(z, "CompoundPhenotypes-SingleEFO.csv", row.names=F, col.names=T, sep=",", quote=T)

y[which(regexpr("(?i)(hepatitis c)", y[,"ph"])>0),]
y[which(regexpr("(?i)(hepatocellular carcinoma)", y[,"ph"])>0),]

z[which(regexpr("(?i)(stroke, coronary heart disease)", z[,"ph"])>0),]

z[which(regexpr("(?i)(asthma)", z[,"ph"])>0),]

##################################################################################################
# Evaluate compound EFOs
##################################################################################################

setwd("C:\\Projects\\Duke\\H2P2GenomeWideAssociationStudy\\CPAG\\iCPAGdb\\App-Devel\\pyCPAG")
f <- c("output/NHGRI-p5e-08-EUR.csv", "output/H2P2-p1e-05-EUR.csv", "output/cpag.csv")[1]
x <- read.table(f, header=T, sep=",", quote="\"")

k <- unique(unlist(lapply(c("Trait1_EFO", "Trait2_EFO"), function(v) which(regexpr(",", x[,v], fixed=T)<0))))
head(x[k,c("Trait1_EFO", "Trait2_EFO")])
nrow(x)-length(k)

length(which(regexpr(",", x[,"Trait1_EFO"], fixed=T)<0 & regexpr(",", x[,"Trait2_EFO"], fixed=T)<0))

length(do.call(intersect,
               lapply(c("Trait1_EFO", "Trait2_EFO"), function(v) which(regexpr(",", x[,v], fixed=T)<0))))

##################################################################################################
# Review COVID results
##################################################################################################

setwd("C:\\Projects\\Duke\\H2P2GenomeWideAssociationStudy\\CPAG\\iCPAGdb\\App-Devel\\pyCPAG")
x <- read.table("input/top_EllinghausPCs_covid19_pcut1e-5.cpag2_out_20201022_addOntology.csv", header=T, sep=",", quote="\"")

# Duke GWAS phenotypes in COVID CPAG are from multiple studies (H2P2, NHGRI, etc.)
library(RSQLite)
db <- dbConnect(RSQLite::SQLite(), "../pyCPAG/db/cpag_gwasumstat_v1.1.db")
y <- dbGetQuery(db, "select distinct trait from GWAStable where source='NHGRI'")

unique(x[,"Trait2"])
unique(y)
x[which(!unique(x[,"Trait2"]) %in% y),"Trait2"]

##################################################################################################
# Construct sample H2P2 GWAS file
##################################################################################################

library(RSQLite)
setwd("C:\\Projects\\Duke\\H2P2GenomeWideAssociationStudy\\CPAG\\iCPAGdb\\App-Devel\\pyCPAG")
db <- dbConnect(RSQLite::SQLite(), "../pyCPAG/db/cpag_gwasumstat_v1.1.db")
x <- dbGetQuery(db, "select trait as Phenotype, SNP, pval as p from GWAStable where source='H2P2' and pval<1e-07")
write.table(x, "input/UserGWAS.csv", row.names=F, col.names=T, quote=T, sep=",")

##################################################################################################
# Evaluate stdout of python script
##################################################################################################

setwd("C:\\Projects\\Duke\\H2P2GenomeWideAssociationStudy\\CPAG\\iCPAGdb\\App-Devevl\\pyCPAG")

pycmd <- paste(pyexe, " main.py cpagdb --threads 2 --subtype H2P2 --H2P2-Pcut 1e-7 --lddb-pop EUR --outfile output/cpag.csv", sep="")

# intern=T instructs to capture output (stdout, stderr)
y <- system(pycmd, intern=T)
attributes(y)
z <- data.frame("out"=y)

# wait=F instructs to execute system call asynchronously
# Must be used with intern=F
y <- system(pycmd, intern=F, wait=F)
attributes(y)
z <- data.frame("out"=y)

##################################################################################################
# Plot execution time vs. phenotype pairs
##################################################################################################

# Pairs
x <- c(17214, 117855, 1299078, 147015, 1841697)

# Seconds to complete
y <- c(60, 404, 3600+15*60+49, 457, 3600+28*60+41)

dev.new()
plot(x[order(x)], y[order(y)], type="l")

m <- lm(y~x)
abline(m[["coefficients"]][1], m[["coefficients"]][2], lty="dotted")

##################################################################################################
# insert Padj_FDR column into Ellinghaus results
##################################################################################################

setwd("C:/Projects/Duke/H2P2GenomeWideAssociationStudy/CPAG/iCPAGdb/App-Devel/pyCPAG/output/")
x <- read.table("top_EllinghausPCs_covid19_pcut1e-5.cpag2_out_20201022_addOntology.csv", header=T, sep=",", quote="\"")
x[,"Padj_FDR"] <- 1
write.table(x, "top_EllinghausPCs_covid19_pcut1e-5.cpag2_out_20201022_addOntology.csv", row.names=F, col.names=T, sep=",", quote=T)

