# Lake Mendota Winter Tag Time Series Analysis for Kali Denis BIO 152 Project
# Authors: Elizabeth McDaniel and Kali Denis
# emails: emcdaniel@wisc.edu and kdenis@wisc.edu 

# Packages 
library(ape)
library(dplyr)
library(ggplot2)
library(gplots)
library(lme4)
library(phangorn)
library(plotly)
library(tidyr)
library(vegan)
library(VennDiagram)
library(phyloseq)
library(stringr)
library(OTUtable)

# Import data
# OTU table
OTU = read.table("Data/WinterTags.final.shared", header=TRUE, sep="\t")
# Taxonomy of each OTU
tax = read.table("Data/WinterTags.OTU.taxonomy", header=TRUE, sep="\t")

# Row names match
row.names(OTU) = OTU$Group
# Remove label, numOTUs, Group columns, not OTU counts
OTU.clean = OTU[,-which(names(OTU) %in% c("label", "numOtus", "Group", "X"))]
# Taxonomy table
row.names(tax) =tax$OTU
tax.clean = tax[row.names(tax) %in% colnames(OTU.clean),]
tax.clean = separate(tax.clean, Taxonomy, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"), sep=";")
# Get rid of underscores
tax.clean = as.data.frame(sapply(tax.clean, gsub, pattern="k__", replacement=""))
tax.clean = as.data.frame(sapply(tax.clean, gsub, pattern="p__", replacement=""))
tax.clean = as.data.frame(sapply(tax.clean, gsub, pattern="c__", replacement=""))
tax.clean = as.data.frame(sapply(tax.clean,gsub, pattern="o__", replacement=""))
tax.clean = as.data.frame(sapply(tax.clean, gsub, pattern = "f__", replacement=""))
tax.clean = as.data.frame(sapply(tax.clean, gsub, pattern = "g__", replacement=""))
tax.clean = as.data.frame(sapply(tax.clean, gsub, pattern = "s__", replacement=""))
row.names(tax.clean) = tax.clean$OTU
# make phyloseq object to combine taxonomy and OTU tables
OTU.UF = otu_table(as.matrix(OTU.clean), taxa_are_rows=FALSE)
tax.UF = tax_table(as.matrix(tax.clean))
# Calculate a rooted phylogenetic tree
seqs <- read.dna("Data/WinterTag-repseqs.fasta", format="fasta")
d <- dist.dna(seqs, model="raw")
tree <- nj(d)
physeq = phyloseq(OTU.UF, tax.UF, tree)

# Visualizations
plot_bar(physeq, fill="Phylum") + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")
