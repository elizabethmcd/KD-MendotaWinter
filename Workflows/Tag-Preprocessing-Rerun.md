# Preprocessing Winter Tag Samples with Mothur Re-run

This workflow references several sources: 

- [Mothur MiSeq SOP](https://www.mothur.org/wiki/MiSeq_SOP) from the Schloss Lab
- [Microbiota Analysis with Mothur Materials](http://rpubs.com/dillmcfarlan/mothurSOP) from Dr. Kim Dill-McFarland's Biotech Center Workshop
- [Microbiota Analysis in R Materials](http://rpubs.com/dillmcfarlan/R_microbiotaSOP) from Dr. Kim Dill-McFarland's Biotech Center Workshop
- [Breakaway Package](https://adw96.github.io/breakaway/) for species richness
- [Denef Lab Verucco Repo](https://github.com/DenefLab/Verruco) including Edna Chiang's scripts for her recent Verucco paper

As of Saturday 2018-04-21, there is something wrong with our initially generated count table, so we are rerunning the preprocessing steps to troubleshoot getting representative OTUs for making phylogenetic trees and objects for all downstream usage in R. 

First run `fastq.info` on each fastq file to generate pairs of fasta and quality files. Concatenate all generated fasta files to WinterTags.fasta Then correlate the file with the sample: 

```
make.group(fasta=ERR1547084_1.fasta-ERR1547102_1.fasta-ERR1547116_1.fasta-ERR1547121_1.fasta-ERR1547062_1.fasta-ERR1547087_1.fasta-ERR1547108_1.fasta-ERR1547117_1.fasta-ERR1547126_1.fasta, groups=A-B-C-D-E-F-G-H-J)
```

Screen the sequences: 
```
screen.seqs(fasta=wintertags.fasta, group=current, maxambig=0, maxlength=300, maxhomop=8)
screen.seqs(fasta=current, group=current, minlength=100, maxlength=100)
```

Resulting summary: 

```
		Start	End	NBases	Ambigs	Polymer	NumSeqs
Minimum:	1	100	100	0	2	1
2.5%-tile:	1	100	100	0	3	107091
25%-tile:	1	100	100	0	3	1070902
Median: 	1	100	100	0	4	2141804
75%-tile:	1	100	100	0	4	3212706
97.5%-tile:	1	100	100	0	5	4176517
Maximum:	1	100	100	0	8	4283607
Mean:	1	100	100	0	3.76369
# of Seqs:	4283607
```

Get unique sequences with `unique.seqs()`. Use the flag `fasta=current` for example since mothur keeps track of the most recent fasta, count, name, list files etc. so we don't have to manually type the filenames, which leaves room for a lot of human error at multiple steps. Then `count.seqs` for simplify the names and groups. 

Resulting summary: 

```
		Start	End	NBases	Ambigs	Polymer	NumSeqs
Minimum:	1	76	76	0	2	1
2.5%-tile:	1	77	77	0	3	10671
25%-tile:	1	92	92	0	3	106709
Median: 	1	100	100	0	4	213418
75%-tile:	1	100	100	0	4	320127
97.5%-tile:	1	100	100	0	6	416165
Maximum:	1	100	100	0	8	426835
Mean:	1	95.79	95.79	0	3.9049
# of Seqs:	426835
```

Then align to the V4 region of the Silva database: 
```
align.seqs(fasta=wintertags.good.good.unique.fasta, reference=silva.v4.fasta, flip=T)

# Summary: 

		Start	End	NBases	Ambigs	Polymer	NumSeqs
Minimum:	0	0	0	0	1	1
2.5%-tile:	8	2445	99	0	3	6235
25%-tile:	8	2445	99	0	3	62350
Median: 	8	2445	99	0	4	124700
75%-tile:	8	2445	99	0	4	187049
97.5%-tile:	8	2448	99	0	5	243164
Maximum:	9582	9582	100	0	8	249398
Mean:	99.1701	2508.88	97.6488	0	3.85655
# of Seqs:	249398
```

Then trim the alignment for those that start at 8 and end at 2447. Filter the sequences and get the unique ones to cut down on redundancy again: 
```
screen.seqs(fasta=wintertags.good.good.unique.align, count=wintertags.good.good.count_table, summary=current, start=8, end=2447)
filter.seqs(fasta=current, vertical=T, trump=.)
unique.seqs(fasta=current, count=current)
pre.cluster(fasta=current, count=current, diffs=2)

# summary

		Start	End	NBases	Ambigs	Polymer	NumSeqs
Minimum:	1	137	89	0	3	1
2.5%-tile:	1	138	97	0	3	127
25%-tile:	1	138	98	0	4	1270
Median: 	1	138	98	0	4	2539
75%-tile:	1	140	99	0	5	3808
97.5%-tile:	1	140	99	0	6	4950
Maximum:	1	140	99	0	8	5076
Mean:	1	138.936	98.407	0	4.17671
# of Seqs:	5076

```
Then remove chimeras. This might be where some difficulty came before, so make sure to just use the current flags since mothur keeps track of what is the most current fasta/count files and there is no user error at this point, especially which count table to use. 

```
chimera.uchime(fasta=current, count=current, dereplicate=t)
remove.seqs(fasta=current, accnos=current)
```

This removed 255 sequences from the fasta file. Then cluster the OTUs by performing distance calculations. Use a cutoff of 0.20. The final sequence number is 4821 after removing chimeras. Then cluster them based on those calculated distances. 

```
dist.seqs(fasta=current, cutoff=0.20)
cluster(column=current, count=current)
```

This then outputs a .list file. Then get representative OTUs for classification and making phylogenetic trees. This is the step where previous errors have been encountered about missing sequences in the count table. However, Alex has a script that gets representative OTUs from the list and sequence files and then correlates it to the OTU. Then make the shared table to show how many OTUs for each sample.  

```
make.shared(list=current, count=current, label=0.03)
```

For downstream analyses and classification, the fasta, shared, and list files are needed, which have been moved to a `final` directory. 

```
accnos=wintertags.good.good.unique.good.filter.unique.precluster.uchime.accnos
column=wintertags.good.good.unique.good.filter.unique.precluster.pick.dist
fasta=wintertags.good.good.unique.good.filter.unique.precluster.pick.fasta
list=wintertags.good.good.unique.good.filter.unique.precluster.pick.an.unique_list.list
rabund=wintertags.good.good.unique.good.filter.unique.precluster.pick.an.unique_list.A.rabund
shared=wintertags.good.good.unique.good.filter.unique.precluster.pick.an.unique_list.shared
count=wintertags.good.good.unique.good.filter.unique.precluster.uchime.pick.count_table
summary=wintertags.good.good.unique.good.filter.unique.precluster.summary
```

Now classifying the OTUs using a combination of greengenes and blast to identify with the Freshwater Training Set using TaxASS and Robin's scripts. The first thing is to make a BLAST database of the Freshwater Training set. I pulled the fasta file off of Ben's database folder of the fileshare, and I will have a copy myself. The command: 

```
makeblastdb -dbtype nucl -in FWonly_11Feb2016_1452_ready.fasta -input_type fasta -parse_seqids -out FWonly_11Feb2016_1452_ready.fasta.db
```

Also reformat the fasta file to get rid of hyphens: `sed -e 's/-//g' WinterTags.final.fasta > WinterTags.final.fasta.reformatted`. 

Then query the FW database and reformat the blast results: 

```
blastn -query WinterTags.final.fasta.reformatted -task megablast -db FWonly_11Feb2016_1452_ready.fasta.db -out WinterOTU.custom.blast -outfmt 11 -max_target_seqs 5
# reformat
blast_formatter -archive WinterOTU.custom.blast -outfmt "6 qseqid pident length qlen qstart qend" -out WinterOTU.custom.blast.table
```

Then using Robin's set of scripts to calculate full length PIDent, pull out sequence IDs greater than and less than 98%, some sanity plots in R. An important note is that the step for getting representative OTUs creates really weird header lines, so keep the fasta file without representative sequences for classification purposes, and use the representative OTU later for creating trees, or creating a tree from the raw fasta and calculating a tree in R. Steps for using Robin's Rscripts: 

```
# calculate full length 
Rscript ../Scripts/calc_full_length_pident.R WinterOTU.custom.blast.table WinterOTU.custom.blast.table.modified

# above 98 match
Rscript ../Scripts/filter_seqIDs_by_pident.R WinterOTU.custom.blast.table.modified ids.above.98 98 TRUE

# below 98
Rscript ../Scripts/filter_seqIDs_by_pident.R WinterOTU.custom.blast.table.modified ids.below.98 98 FALSE

# pull out OTUs that hit nothing on blast
python ../Scripts/find_seqIDs_blast_removed.py WinterTags.final.fasta.reformatted WinterOTU.custom.blast.table.modified ids.missing

# combine list of seequence IDs that didn't hit anythign on blast with those that had pident below 98
cat ids.below.98 ids.missing > ids.below.98.all

# create fasta file with sequences that correspond to each group
python ../Scripts/create_fastas_given_seqIDs.py ids.above.98 WinterTags.final.fasta.reformatted above.98.fasta

python ../Scripts/create_fastas_given_seqIDs.py ids.below.98.all WinterTags.final.fasta.reformatted below.98.fasta
```

Now the different fasta files are setup to classify the groups of OTUs using mothur with the two different databases. Sequences with a pident above 98 are classified using the FW database. And sequences below 98 are classified using the Greengenes database within mothur. The fasta and taxonomy files for both the freshwater database and Greengenes. 

Within mothur for both classifications: 

```
classify.seqs(fasta=above.98.fasta, template=FWonly_11Feb2016_1452_ready.fasta, taxonomy=FWonly_11Feb2016_1452.taxonomy, method=wang, probs=T, processors=14)

classify.seqs(fasta=below.98.fasta, template=Gg.fasta, taxonomy=Gg.taxonomy, method=wang, probs=T, processors=14)

# concatenate 
cat below.98.Gg.wang.taxonomy above.98.FWonly_11Feb2016_1452.wang.taxonomy > WinterTags.final.cat.taxonomy
```

Then get OTUs for the classifications 
```
classify.otu(list=wintertags.good.good.unique.good.filter.unique.precluster.pick.an.unique_list.list, taxonomy=WinterTags.final.cat.taxonomy, count=wintertags.good.good.unique.good.filter.unique.precluster.uchime.pick.count_table, label=0.03, cutoff=80, basis=otu, probs=F)

# rename the final taxonomy file

cp wintertags.good.good.unique.good.filter.unique.precluster.pick.an.unique_list.0.03.cons.taxonomy WinterTags.OTU.taxonomy
```

The important files are now named:
```
fasta: WinterTags.final.fasta.reformatted
list: WinterTags.final.list
taxonomy: WinterTags.OTU.taxonomy
shared: WinterTags.final.shared
```
