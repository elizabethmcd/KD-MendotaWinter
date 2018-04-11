#!/bash/sh

# Directory: /Volumes/mcmahonlab/home/bpeterson26/extractionTest/dataRaw/samplesToCompare/

# Remove old mothur log files before starting
rm mothur.*.logfile

# Blast OTUs against FW blast database
ncbi-blast-2.6.0+/bin/blastn -query 5M.fasta -task megablast -db FWonly_11Feb2016_1452_ready.fasta.db -out OTU.custom.blast -outfmt 11 -max_target_seqs 5

# Reformat the blast results
ncbi-blast-2.6.0+/bin/blast_formatter -archive OTU.custom.blast -outfmt "6 qseqid pident length qlen qstart qend" -out OTU.custom.blast.table

# Use Robin's script to calculate the full length pident of each blast result.
Rscript calc_full_length_pident.R OTU.custom.blast.table OTU.custom.blast.table.modified

# Use Robin's script to pull out sequence IDs that have greater than 98%
# identity to FW database entry
Rscript filter_seqIDs_by_pident.R OTU.custom.blast.table.modified ids.above.98 98 TRUE

# Use Robin's script to pull out sequence IDs that have less than 98%
# identity to FW database entry
Rscript filter_seqIDs_by_pident.R OTU.custom.blast.table.modified ids.below.98 98 FALSE

# Generate a sanity check plot that can be seen in R
# Run Robin's R script to generate the plots.
#Rscript plot_blast_hit_stats.R OTU.custom.blast.table.modified 98 plots

# Python code to pull out OTUs that hit nothing on the blast
python find_seqIDs_blast_removed.py 5M.fasta OTU.custom.blast.table.modified ids.missing

# Combine the list of sequence IDs that didn't hit anything on the blast with
# the sequence IDs that had a pident below 98.
cat ids.below.98 ids.missing > ids.below.98.all

# Create a fasta file with the sequences that correspond to each group (above
# or below 98%). This also uses a python code written by Robin.
python create_fastas_given_seqIDs.py ids.above.98 5M.fasta otus.above.98.fasta
python create_fastas_given_seqIDs.py ids.below.98.all 5M.fasta otus.below.98.fasta

# Classify each group of OTUs using mothur with the two different databases.
# Sequences with pident above 98 are classified using the FW database.
# Sequences with pident below 98 are classified using the Gg database.
mothur/mothur "#classify.seqs(fasta=otus.above.98.fasta, template=FWonly_11Feb2016_1452_ready.fasta, taxonomy=FWonly_11Feb2016_1452.taxonomy, method=wang, probs=T, processors=2)"
mothur/mothur "#classify.seqs(fasta=otus.below.98.fasta, template=Gg.fasta, taxonomy=Gg.taxonomy, method=wang, probs=T, processors=2)"

# Concatonate the two separate classifications together.
cat otus.above.98.FWonly_11Feb2016_1452.wang.taxonomy otus.below.98.Gg.wang.taxonomy > 5M.taxonomy

# Reformat the OTU classification files so that it is delimited by semicolons,
# rather than spaces or tabs.
sed 's/[[:blank:]]/\;/' <5M.taxonomy > 5M.taxonomy.reformatted
mv 5M.taxonomy.reformatted 5M.taxonomy
