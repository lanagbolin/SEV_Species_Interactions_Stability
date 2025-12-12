# Extracting plant phylogenies for Catherine and Steve to use in pollinator effectivness meta-analyses 3/2019
# this code not quite working perfectly.  TPL finds synonyms, but those have to be replaced in the original species list (.csv file)
# by hand, by consulting the "unmatchedTPL" file

# Start with Tim phylogeny
# this contains some Sev species added to the Qian & Ji 2016 phylogeny of 31,383 species
# which in turn is an updated and fixed version of the Zanne et al. 2014 phylogeny of 32,223 land plant species, branch lengths in MY

library(Taxonstand)
library(pez)
library(vegan)
library(picante)
library(lattice)
library(Hmisc)
library(readxl)
library(phytools)
library(visreg)
library(tidyverse)
library(ggtree)
library(phylotools)
setwd("C:/Users/Ken/Documents/Advisees & Collaborators/Catherine Cumberland/R analyses")

#input Catherine species lists
data<-read.csv("Catherine species list all v2.csv") #this is the full list; phylogeny is built for this, then pruned based on the below
SVGplantlist<-as.vector(read.csv("SVGplantlist.csv"))
Apisplantlist<-as.vector(read.csv("Apisplantlist.csv"))

# input megaphylogeny 
# For TPL, accepted spacer is " " (blank space): Genus species 
# For most R phylo tools, accepted spacer is underscore: Genus_species
tim_tree<-read.tree("PhytoPhylo_manadd5.tre")

# make list of genera in tim_tree
tim_list<-tim_tree$tip.label #extract a species list from phylogeny
#write.csv(tim_list, file = "tim species list.csv") # make a species list in excel so that a file of genera can be constructed by hand in excel
tim_genera<-read.csv("tim_genera.csv") # import file of genera

# Make lists of my taxa that do and don't match tim_tree
# Note: setdiff returns items in the first set, but not in the second set.  Thus set order matters.
matched<-intersect(data$Gen_sp,tim_tree$tip.label) 
unmatched<-setdiff(data$Gen_sp,tim_tree$tip.label) 

#write.csv(matched, file = "matched.csv")
#write.csv(unmatched, file = "unmatched.csv")


# Now use TPL (function that accesses The Plant List) to find synonyms for unmatched spp
# Note: TPL Will give error messages if there are multiple synonyms for a name, 
#     but congeneric.merge (below) will still add one of the synonyms to the tree 
unmatched2<-gsub('_', ' ', unmatched) #replace Genus_species with Genus species for purposes of accessing TPL
unmatchedTPL<-TPL(splist=unmatched2, infra = TRUE,  
                  corr = TRUE, diffchar = 2, max.distance = 1, version = "1.1")
unmatchedTPL$taxon<-paste(unmatchedTPL$New.Genus,unmatchedTPL$New.Species, sep='_') # Make a new column with the official TPL taxon name, spacer = underscore
unmatchedTPL$taxonold<-paste(unmatchedTPL$Genus,unmatchedTPL$Species, sep='_') # Make a new column with the input (old) taxon name, spacer = underscore
write.csv(unmatchedTPL, file ="unmatchedTPL.csv")

# Update the list of matches
matched2<-union(matched,intersect(unmatchedTPL$taxon,tim_tree$tip.label))

# Update list of TPL names of spp STILL not matched
unmatched3<-setdiff(unmatchedTPL$taxon,tim_tree$tip.label) #species still not matched: TPL names

# Some spp in unmatched3 are in TPL genera not present in Zanne tree, so they will not be added in congeneric.merge (below)
# So, recover the original (non-TPL) generic names of these species to try to merge them as well
temp1<-unmatchedTPL[unmatchedTPL$taxon %in% unmatched3,]  #use vector of unmatched species names to recover the dataframe 
temp1genera<-temp1$New.Genus                              #extract TPL genus names for unmatched3 vector
unmatched3gen<-setdiff(temp1genera,tim_genera)          #vector of TPL genera in unmatched3 that are NOT present in Zanne tree
temp2<-unmatchedTPL[unmatchedTPL$New.Genus %in% unmatched3gen,] #use above vector to recover the dataframe
unmatched4<-temp2$taxonold                                # list of original (non-TPL) names of species in TPL genera NOT present in Zanne tree
unmatched5<-union(unmatched3, unmatched4)                 # combined list of TPL and non-TPL names

# Try to merge each of these spp into tim_tree by "replacing all members of the clade it belongs to with a polytomy" 
# Behavior of congeneric.merge: it places the polytomy halfway between tips and the common ancestor of the genus and its sister group
# Technically, it drops every member of the genus except one, places a node halfway along that remaining branch, then adds back in every 
# member of the genus (plus the new one(s) you are merging) to that node
tim_modified<-congeneric.merge(as.character(unmatched5), tree = tim_tree,split = "_")

# Prune the tree: drop the species that are in tim_tree, but not in my set, where my set includes the exact matches plus the ones stuck in as polytomies
# Note that the pruned tree will be a combo of my names and TPL names
pruned.tree<-drop.tip(tim_modified, setdiff(tim_modified$tip.label, union(matched2,unmatched5)))


# SUMMARY OF TAXON MATCHING
cat(length(unmatched)+length(matched), "taxa in my list") 
cat(length(matched)," of my original taxa matched to tim_tree")
cat(length(intersect(unmatchedTPL$taxon,tim_tree$tip.label))," additional taxa matched to tim_tree after standardizing to TPL name")
cat(length(unmatched3)," additional taxa ATTEMPTED to add via polytomies")
cat(length(tim_modified$tip.label)-length(tim_tree$tip.label), " taxa actually added via polytomies")
cat(length(pruned.tree$tip.label), " total taxa in pruned tree")
xx<-length(unmatched3)-(length(tim_modified$tip.label)-length(tim_tree$tip.label))
cat("(",xx, " taxa missing b/c genus not present in tim_tree)")
cat("(",length(unmatched)+length(matched)-xx-length(pruned.tree$tip.label), " additional taxa missing -- doubles removed b/c synonymous??)")

plot(pruned.tree)

#replace the working tip names (Gen_sp) with the original tip names (Plant_sp.)
sublist <- data.frame(data$Gen_sp,data$Plant_sp.) #2 columns: old names, new names
pruned.tree2 <- sub.taxa.label(pruned.tree, sublist) 
plot(pruned.tree2)

#subset the phylogenies for specific analyses
plantsSVG<-keep.tip(pruned.tree2,SVGplantlist$Plant_sp.)
plot(plantsSVG)
plantsApis<-keep.tip(pruned.tree2,Apisplantlist$Plant_sp.)
plot(plantsApis)

#write tree files
write.tree(pruned.tree2, file="C:/Users/Ken/Documents/Advisees & Collaborators/Catherine Cumberland/R analyses/tree_all.tre") #Newick format.  Default is to put it in My Documents
write.nexus(pruned.tree2, file="tree_all.nex") #Nexus format
write.tree(plantsSVG, file="C:/Users/Ken/Documents/Advisees & Collaborators/Catherine Cumberland/R analyses/tree_SVG.tre") #Newick format.  Default is to put it in My Documents
write.nexus(plantsSVG, file="tree_SVG.nex") #Nexus format
write.tree(plantsApis, file="C:/Users/Ken/Documents/Advisees & Collaborators/Catherine Cumberland/R analyses/tree_Apis.tre") #Newick format.  Default is to put it in My Documents
write.nexus(plantsApis, file="tree_Apis.nex") #Nexus format

