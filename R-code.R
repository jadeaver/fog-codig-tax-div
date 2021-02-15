#R Code used for data analyses reported in Deaver et al., 2021

library(phyloseq)
library(vegan)
library(indicspecies)

library(ape)
library(readxl)

######################################################################################
#Creating a phyloseq object

tree.unrooted <- read_tree("FastTree_otu.nhx")

#Rooting tree
#Function from: https://github.com/joey711/phyloseq/issues/597
pick_new_outgroup <- function(tree.unrooted){
  require("magrittr")
  require("data.table")
  require("ape") # ape::Ntip
  # tablify parts of tree that we need.
  treeDT <- 
    cbind(
      data.table(tree.unrooted$edge),
      data.table(length = tree.unrooted$edge.length)
    )[1:Ntip(tree.unrooted)] %>% 
    cbind(data.table(id = tree.unrooted$tip.label))
  # Take the longest terminal branch as outgroup
  new.outgroup <- treeDT[which.max(length)]$id
  return(new.outgroup)
}

new.outgroup <- pick_new_outgroup(tree.unrooted)
rootedTree <- ape::root(tree.unrooted, outgroup=new.outgroup, resolve.root=TRUE)

otudf <- read_excel("otu_T.xlsx")
taxdf <- read_excel("taxa_T.xlsx")
meta <- read_excel("metadata.xlsx")

otumat <- data.matrix(otudf)
taxmat <- as.matrix.data.frame(taxdf)

#Add otu names (taxa and abundance table in order otu1 to otu1417)
rownames(otumat) <- paste0("otu", 1:nrow(otumat))
rownames(taxmat) <- paste0("otu", 1:nrow(otumat))

OTU <- otu_table(otumat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)

physeq <- phyloseq(OTU, TAX, rootedTree)
physeq

sampledata <- sample_data(meta)
row.names(sampledata)<-sample_names(physeq)

physeq <- merge_phyloseq(physeq, sampledata)
physeq

uni <- UniFrac(physeq, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)

######################################################################################
#Creating NMDS plot

wUF.dist <- UniFrac(physeq, weighted=TRUE, normalized=TRUE)
wUF.nmds <- metaMDS(wUF.dist,method="NMDS",k=2)
wUF.nmds

wUF.nmds$points
wUF.nmds$data

plot(wUF.nmds, display = "sites", type = "n")

group = c(rep("Group1", 8), rep("Group2", 8))
colors = c(rep("blue", 8), rep("orange", 8))

for(i in unique(group)) {
  ordihull(wUF.nmds$point[grep(i, group),], draw="polygon",
           groups = group[group == i],col = colors[grep(i,group)],label=F)
}

orditorp(wUF.nmds, display = "sites", col = c(rep("black",8),                          
                                                 rep("black", 8)), air = 1, cex = 0.9) 

######################################################################################
#ANOSIM and Mantel tests

anosim(wUF.dist, meta$FOG, permutations = 500000)
anosim(wUF.dist, meta$Overload, permutations = 500000)

pH <- meta$pH
tCOD <- meta$tCOD
NH3 <- meta$NH3
Phosp <- meta$Phosp

dist.pH <- dist(pH, method = "euclidean")
dist.tCOD <- dist(tCOD, method = "euclidean")
dist.NH3 <- dist(NH3, method = "euclidean")
dist.Phosp <- dist(Phosp, method = "euclidean")

dist.abund <- uni

mantel(dist.abund, dist.pH, method = "spearman", permutations = 500000, na.rm = TRUE)
mantel(dist.abund, dist.tCOD, method = "spearman", permutations = 500000, na.rm = TRUE)
mantel(dist.abund, dist.NH3, method = "spearman", permutations = 500000, na.rm = TRUE)
mantel(dist.abund, dist.Phosp, method = "spearman", permutations = 500000, na.rm = TRUE)

######################################################################################
#Inidcator Species analysis

genus.group <- tax_glom(physeq, "genus" )

#Converting OTU table in phyloseq object to data object compatible with vegan/indicspecies
#Function from: https://jacobrprice.github.io/2017/08/26/phyloseq-to-vegan-and-back.html
psotu2veg <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

otu.table<- psotu2veg(genus.group)

status <- meta$Overload

inv <- multipatt(otu.table, status, func = "r.g", control = how(nperm=500000))
summary(inv)

######################################################################################