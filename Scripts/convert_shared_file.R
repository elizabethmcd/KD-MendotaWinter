# Convert shared file for 5M project

# File requirements: shared file must be under dataEdited/ folder.

# This script will read in a .shared file and convert it to a relatived
# abundance table. This particular script has been edited to work on
# CHTC.
# Command line script for SDUE project:
# R CMD BATCH '--args shared.file=*.shared' convert_shared_file.R


# Set up user input and read in file ####
########################
userprefs <- commandArgs(trailingOnly = TRUE)
print(userprefs)
if(length(userprefs) != 1){
  stop("Wrong number of arguments, dummy.")
}

OTU.table <- read.table(userprefs[1],
                        header = TRUE,
                        sep = "\t",
                        stringsAsFactors = FALSE)
#####



# Edit OTU table before transformation ####
########################
# Save sample names in separate vector
sampleID <- OTU.table[, 2]

# Remove the first three columns
# First column is a label of the clustering cutoff
# Second column is sample name.
# Never fear, these will be added in as row names!
# Third column is the total number of OTUs
OTU.table <- OTU.table[, -c(1:3)]

# Add sample names as row names
row.names(OTU.table) <- sampleID
#####


# Calculate relative abundances for each OTU. ####
########################
# Calculate the row sums for each
sums.of.rows <- rowSums(OTU.table)

# Apply function to get the relative abundance table.
REL.OTU.table <- apply(OTU.table, # We're applying this function the OTU.table variable
                         2, # function is applied over columns (cols = 2)
                         function(x) {
                           (x/sums.of.rows) # Want to divide the number of hits of a particular OTU in a sample by the total number of hits in that sample.
                           })
t.REL.OTU.table <- t(REL.OTU.table)

# Make the sampleIDs a part of the tab-delimited file
# Want to make sure there isn't a gap in top-left corner.
REL.OTU.table.OTUs <- cbind(colnames(REL.OTU.table), t.REL.OTU.table)
colnames(REL.OTU.table.OTUs) <- c("OtuID", sampleID)
#####


# Write-out OTU table
########################
# Write out the OTU table to a tab-delimited file
write.table(REL.OTU.table.OTUs,
            file = paste(strsplit(userprefs[1],
                                  ".shared")[[1]][1],
                         ".abund",
                         sep = ""),
            # Don't use quotes in the output file.
            quote = FALSE,
            # The OTU names are stored as the column names.
            # The sampleID column also has the sampleID header.
            col.names = TRUE,
            row.names = FALSE,
            # Separate with a tab.
            sep = "\t")
#####
