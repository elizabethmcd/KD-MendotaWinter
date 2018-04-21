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

# Import data
# OTU table
OTU = read.table("Data/WinterTags.final.an.unique_list.shared", header=TRUE, sep="\t")
# Taxonomy of each OTU
tax = read.table("Data/WinterTags.taxonomy.reformatted", header=FALSE, sep=";")

# Row names match
row.names(OTU) = OTU$Group
# Remove label, numOTUs, Group columns, not OTU counts
OTU.clean = OTU[,-which(names(OTU) %in% c("lablel", "numOtus", "Group"))]
# Taxonomy table
row.names(tax) =tax$OTU
