# Authors: Alexandra Linz from Microbial Observatory Scripts, repurposed by Elizabeth McDaniel 

# Goal: make an OTU fasta file

# I need: a key that says what seqIDs went into which OTU
# something that says what the representative sequence for each OTU is
# the seqID fasta file
# The taxonomy file reports a consensus of all seqs in an OTU. No rep seqs were chosen
# The list file gives what sequences ended up in which OTU, separated by certain cutoffs. Make sure to know what cutoff/label was used, for this amount of sequences will match the OTU and taxonomy files. The list file doesn't give names of the OTUs, which can be changed when select for the correct label 

library(dplyr)
library(stringr)

# Don't read in header for list, will include the first sequence for each OTU
list <- read.table("Data/WinterTags.final.list", header = F, row.names = 1, fill = T, colClasses = c("character"))
# Only 3rd row with 0.03 label to correspond to 3325 OTUs like tax and shared files
list.lab <- slice(list, 3)
# Rid of sequence count
OTU.list <- list.lab[,-c(1)] 
# Get rid of empty columns
new.list <- OTU.list[!sapply(OTU.list, function (x) all(is.na(x) | x ==""))]
# Count OTU + # and rename columns
count <- seq(0001, 3325, by=1)
zcount <- str_pad(count, 4, pad="0")
OTU <- "Otu"
OTUnames <- paste(OTU, zcount, sep="")
names(new.list) <- OTUnames

# read in raw fasta file with dashes, because needs to be the same length for distance calculations
fasta <- read.table("Data/wintertags.good.good.unique.good.filter.unique.precluster.pick.fasta", header = F, colClasses = c("character"))


# Approach:
# - loop through 0.02 row
# - split up mulitple entries and choose one randomly if necessary
# - once seq id is chosen, add > and match to entry in fasta file
# - get sequence and replace seqID with otu number
# - create and export a fasta file with otus and rep seqs

# Remove OTUs that didn't make it through subsampling - about 100
# Load taxonomy file to see which OTUs to keep

taxonomy <- read.table("Data/WinterTags.OTU.taxonomy", header = T, sep="\t")
otus <- as.character(taxonomy[, 1])
keep <- match(otus, colnames(new.list))
cand <- new.list[1, keep]

new_fasta <-c()
for(i in 1:length(cand)){
  entry <- cand[i]
  if(nchar(entry) > 6){
    all <- strsplit(as.character(entry), split = ",")
    entry <- sample(all[[1]], 1)
  }
  search <- paste(">", entry, sep = "")
  seqid <- match(search, fasta$V1)
  
  new_fasta <- append(new_fasta, paste(">", colnames(new.list)[keep[i]], sep = ""), length(new_fasta))
  new_fasta <- append(new_fasta, fasta$V1[seqid+1], length(new_fasta))
}

write.table(new_fasta, file="WinterTag-repseqs.fasta", quote = F, row.names = F, col.names = F)

# This gives a FASTA file with gaps for calculating distances matrix with OTUs are the representative sequences, and get the tree file for phyloseq


