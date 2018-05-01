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
library(tidyr)
library(magrittr)

# Import data
# OTU table
OTU = read.table("Data/WinterTags.final.shared", header=TRUE, sep="\t")
# Taxonomy of each OTU
tax = read.table("Data/WinterTags.OTU.taxonomy", header=TRUE, sep="\t")
# Sample data
meta = read.table("Data/sample-metadata.txt", header=TRUE, row.names=1, sep="\t")
meta.UF = sample_data(meta)


# Row names match
row.names(OTU) = OTU$Group
sample_id <- OTU[,2]
# Remove label, numOTUs, Group columns, not OTU counts
OTU.clean = OTU[,-which(names(OTU) %in% c("label", "Group", "numOtus", "X"))]
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

# Phyloseq object
physeq = phyloseq(OTU.UF, tax.UF, meta.UF, tree)

# Visualizations
# All Phyla
p1 <- plot_bar(physeq, fill="Phylum", x="Date") + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")
# Subset by Actinobacteria
subset <- subset_taxa(physeq, Phylum=="Actinobacteria")
p2 <- plot_bar(subset, fill="Family", x="Date")

# Plot the 10 most abundant taxa in heatmap by Family
fam <- prune_taxa(names(sort(taxa_sums(physeq), TRUE)[1:10]), physeq)
p3 <- plot_heatmap(fam, sample.label="Date", taxa.label="Family", sample.order="Date")

# Ordinations with abundant phyla
phylum.sum = tapply(taxa_sums(physeq), tax_table(physeq)[, "Phylum"], sum, na.rm=TRUE)
top5phyla = names(sort(phylum.sum, TRUE))[1:5]
phys = prune_taxa((tax_table(physeq)[, "Phylum"] %in% top5phyla), physeq)
phy.ord <- ordinate(phys, method="NMDS", distance="bray")
p4 = plot_ordination(phys, phy.ord, type="taxa", color="Phylum", title="taxa")
p4
p4 <- p4 + facet_wrap(~Phylum, 3)

# Plot by sample of abundant phyla
# Mapping colors onto variables doesn't work with only one variable, add a dummy column to the sample_data
# NMDS for correspondance/variance between objects in an ordination, represents the pairwise dissimilarity between objects in low-dimensional space, where points represent objects/samples, more similar are closer together, axes are arbitrary as also is the ordination
no = plot_ordination(phys, phy.ord, type="samples", color="Date")
no
# split graphic where samples/OTUs separated
# Plot of ordination of most abundant phyla and give a visualization of beta-diversity at the phylum level
p5 = plot_ordination(phys, phy.ord, type="split", color="Phylum", label="Date")
p5

# Statistically test beta=diversity
BC.dist = vegdist(OTU.clean, distance="bray")
anosim(BC.dist, meta$Date, permutations=6000)
# no significant difference between OTUs of samples

# Save plots
# All phyla
ggsave("WinterTags-all-phyla.png", plot=p1)
ggsave("WinterTags-Actino.png", plot=p2)
ggsave("WinterTags-Actino-abundance.png", plot=p3)
ggsave("WinterTags-ord-phyla.png", plot=p4)
ggsave("WinterTags-sample-ord-phyla.png", plot=p5)
