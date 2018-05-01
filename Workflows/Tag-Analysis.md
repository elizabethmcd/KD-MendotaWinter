# Preprocessing Winter Tag Samples with Mothur

This workflow references several sources: 

- [Mothur MiSeq SOP](https://www.mothur.org/wiki/MiSeq_SOP) from the Schloss Lab
- [Microbiota Analysis with Mothur Materials](http://rpubs.com/dillmcfarlan/mothurSOP) from Dr. Kim Dill-McFarland's Biotech Center Workshop
- [Microbiota Analysis in R Materials](http://rpubs.com/dillmcfarlan/R_microbiotaSOP) from Dr. Kim Dill-McFarland's Biotech Center Workshop
- [Breakaway Package](https://adw96.github.io/breakaway/) for species richness
- [Denef Lab Verucco Repo](https://github.com/DenefLab/Verruco) including Edna Chiang's scripts for her recent Verucco paper

## Tuesday 2018-04-10

We are analyzing 10 16S samples from the Lake Mendota Tag Timeseries dataset with a focus on winter samples (Dec, Jan, Feb) of various years. These were sequenced by the Earth Microbiome Project (EMP) at UCSF, and are not interleaved paired reads. This deviates from the Mothur SOP for the first couple of steps, mostly with the `make.conitgs` step. 

### Making FASTQ files

We ran the command `fastq.info` on each sample to create pairs of .fasta and .qual files that mothur can read. Now we are combining all 10 samples into one big group for further anlaysis with the following command: 

```
make.group(fasta=ERR1547040_1.fasta-ERR1547084_1.fasta-ERR1547102_1.fasta-ERR1547116_1.fasta-ERR1547121_1.fasta-ERR1547062_1.fasta-ERR1547087_1.fasta-ERR1547108_1.fasta-ERR1547117_1.fasta-ERR1547126_1.fasta, groups=A-B-C-D-E-F-G-H-J-K)

```

We have intentionally skipped naming a group the letter I because it looks like the number 1. 

In order to get one giant fasta file from all 10 samples, we will merge them together to get summary statistics and run filtering on everything downstream. 

```
merge.files(input=ERR1547040_1.fasta-ERR1547084_1.fasta-ERR1547102_1.fasta-ERR1547116_1.fasta-ERR1547121_1.fasta-ERR1547062_1.fasta-ERR1547087_1.fasta-ERR1547108_1.fasta-ERR1547117_1.fasta-ERR1547126_1.fasta, output=WinterTags.fasta)
```

Note, we are using the most current version of mothur, and the `merge.files` command doesn't seem to be working for whatever reason, so I concatenated the files and it seems to be the same thing. Whereas a different version of mothur I have installed performs `merge.files` just fine. Moving on with the summary statistics.
 
 Preliminary summary statistics: 

 ```

		Start	End	NBases	Ambigs	Polymer	NumSeqs
Minimum:	1	76	76	0	2	1
2.5%-tile:	1	88	88	0	3	126288
25%-tile:	1	100	100	0	3	1262878
Median: 	1	100	100	0	4	2525756
75%-tile:	1	100	100	0	4	3788634
97.5%-tile:	1	100	100	0	5	4925224
Maximum:	1	100	100	0	25	5051511
Mean:	1	98	98	0	3
# of Seqs:	5051511

It took 27 secs to summarize 5051511 sequences.
```

Now to filter the sequences by length and ambiguous bases: 

```
screen.seqs(fasta=wintertags.fasta, group=merge.groups, summary=wintertags.summary, maxambig=0, minlength=100, maxlength=100, maxhomop=8)
```

Summary of filtered sequences: 

```

		Start	End	NBases	Ambigs	Polymer	NumSeqs
Minimum:	1	100	100	0	2	1
2.5%-tile:	1	100	100	0	3	107091
25%-tile:	1	100	100	0	3	1070903
Median: 	1	100	100	0	4	2141806
75%-tile:	1	100	100	0	4	3212709
97.5%-tile:	1	100	100	0	5	4176521
Maximum:	1	100	100	0	8	4283611
Mean:	1	100	100	0	3
# of Seqs:	4283611

It took 23 secs to summarize 4283611 sequences.
```

Get unique sequences: 
```
unique.seqs(fasta=wintertags.good.good.fasta)
```

After getting unique sequences:
```
		Start	End	NBases	Ambigs	Polymer	NumSeqs
Minimum:	1	100	100	0	2	1
2.5%-tile:	1	100	100	0	3	6236
25%-tile:	1	100	100	0	3	62351
Median: 	1	100	100	0	4	124701
75%-tile:	1	100	100	0	4	187051
97.5%-tile:	1	100	100	0	5	243166
Maximum:	1	100	100	0	8	249401
Mean:	1	100	100	0	3
# of Seqs:	249401

It took 2 secs to summarize 249401 sequences.
```

Simplify names and groups with `counts`: 
```
count.seqs(name=wintertags.good.good.names, group=merge.good.groups)
```

End of the day on 2018-04-10, we have formatted the fastq files to individual fasta and quality files, manually concatenated those because `merge.files` isn't working, then filtered the sequences by length, got rid of duplicated sequences and performed counts. Next will be aligning the specific 16S to the Silva database. 

Running a summary on the count table will also show the unique sequences and the total number of sequences for simplification purposes: 

```
		Start	End	NBases	Ambigs	Polymer	NumSeqs
Minimum:	1	100	100	0	2	1
2.5%-tile:	1	100	100	0	3	107091
25%-tile:	1	100	100	0	3	1070903
Median: 	1	100	100	0	4	2141806
75%-tile:	1	100	100	0	4	3212709
97.5%-tile:	1	100	100	0	5	4176521
Maximum:	1	100	100	0	8	4283611
Mean:	1	100	100	0	3
# of unique seqs:	249401
total # of seqs:	4283611

It took 4 secs to summarize 4283611 sequences.
```

Now aligning sequences to a reference alignment with `align.seqs`. We will be using Silva wiht the full-length sequences. I will align to the full database just because I'm not too familiar with this dataset, and this might be a different version of Silva. Will probably end up with the normal V4 region anyways, but again because of version issues not 100% sure. The full Silva NR database is about 10GB and that has been moved to `/home/emcdaniel/Databases/Silva.nr_v132` on the fileshare for easy access in the future. This will be moved to the current working directory where we are running mothur, and then use the reduced database with our region of interest for future runs. The Silva NR database used for this anlayses is the 132 release. 

```
align.seqs(fasta=wintertags.good.good.unique.fasta, reference=silva.nr_v132.align, flip=T)
```

The `flip=T` flag improves the alignment so both the forward and reverse complement of the sequence is checked. After aligning to the full database: 

```
It took 1499 secs to align 249401 sequences.

[WARNING]: 3541 of your sequences generated alignments that eliminated too many bases, a list is provided in wintertags.good.good.unique.flip.accnos.
[NOTE]: 1638 of your sequences were reversed to produce a better alignment.
```

And the summary report: 
```
		Start	End	NBases	Ambigs	Polymer	NumSeqs
Minimum:	0	0	0	0	1	1
2.5%-tile:	13862	16307	100	0	3	107091
25%-tile:	13862	16307	100	0	3	1070903
Median: 	13862	16307	100	0	4	2141806
75%-tile:	13862	16307	100	0	4	3212709
97.5%-tile:	13862	16309	100	0	5	4176521
Maximum:	43116	43116	100	0	8	4283611
Mean:	13883	16325	99	0	3
# of unique seqs:	249401
total # of seqs:	4283611

It took 26 secs to summarize 4283611 sequences.
```

Most of the sequences are aligning at 13862 and 16307/9, not sure what is up with the 7/9 base pairing. There is definitely one that is start and end at the same site, 43116. Now using the fasta aligned file `wintertags.good.good.unique.align` and will further screen the alignments before moving on to clustering and classifications. 

To fit Robin's scripts we might have to reformat the fasta file names so the sequence name is >OtuXXXX. She might have a script that does that already. After clustering and getting representative OTUs and a shared table, we will deviate from the mothur SOP to classifying with BLAST/TaxASS with the Freshwater Training Set. The pipeline and scripts will use the 98% cutoff to use TaxASS for classification, and the rest with Mothur using Greengenes or something comparable. Then we will have classified OTUs and be able to do some analyses on our winter samples. 

## Wednesday 2018-04-11

Rescreening the aligned sequences to the full Silva NR database. Most of the sequences are within 13862 start site and 16309 end site. The mean start/stop site gives a different range, but I think we can go with this for now which is where most of them are ending. 

```
screen.seqs(fasta=wintertags.good.good.unique.align, count=wintertags.good.good.count_table, summary=wintertags.good.good.unique.summary, start=13862, end=16309)

## screen again to still get rid of bad sequences

screen.seqs(fasta=wintertags.good.good.unique.good.align, count=wintertags.good.good.good.count_table, summary=wintertags.good.good.unique.good.summary, start=13862, end=16309)
```

Weird sequence that is much longer than the overlapping region, hopefully can get rid of it with the `filter` function to get rid of overhangs. 

```
filter.seqs(fasta=wintertags.good.good.unique.good.good.align, vertical=T, trump=.)
# unique sequences
unique.seqs(fasta=wintertags.good.good.unique.good.good.filter.fasta, count=wintertags.good.good.good.good.count_table)
```

Now to pre-cluster, current files are: 
```
wintertags.good.good.unique.good.good.filter.count_table
wintertags.good.good.unique.good.good.filter.unique.fasta
pre.cluster(fasta=wintertags.good.good.unique.good.good.filter.unique.fasta, count=wintertags.good.good.unique.good.good.filter.count_table, diffs=2)
# summary
summary.seqs(fasta=wintertags.good.good.unique.good.good.filter.unique.precluster.fasta, count=wintertags.good.good.unique.good.good.filter.unique.precluster.count_table)

		Start	End	NBases	Ambigs	Polymer	NumSeqs
Minimum:	1	141	60	0	3	1
2.5%-tile:	1	161	99	0	3	2654
25%-tile:	1	163	100	0	4	26537
Median: 	1	163	100	0	4	53074
75%-tile:	1	163	100	0	4	79611
97.5%-tile:	1	163	100	0	5	103494
Maximum:	3	163	100	0	8	106147
Mean:	1	162	99	0	4
# of unique seqs:	4274
total # of seqs:	106147
```

Now we are down to around ~4000 unique sequences to work with - I still haven't figured out what the start and end sites really mean, but most of them go from 1 to 163. 

Getting rid of chimeras: 
```
# Identifying chimeras
chimera.uchime(fasta=wintertags.good.good.unique.good.good.filter.unique.precluster.fasta, count=wintertags.good.good.unique.good.good.filter.unique.precluster.count_table, dereplicate=t)
# Getting rid of chimeras
remove.seqs(fasta=wintertags.good.good.unique.good.good.filter.unique.precluster.fasta, count=wintertags.good.good.unique.good.good.filter.unique.precluster.count_table, accnos=wintertags.good.good.unique.good.good.filter.unique.precluster.denovo.uchime.accnos)
# Summary
summary.seqs(fasta=wintertags.good.good.unique.good.good.filter.unique.precluster.pick.fasta, count=wintertags.good.good.unique.good.good.filter.unique.precluster.pick.count_table)

		Start	End	NBases	Ambigs	Polymer	NumSeqs
Minimum:	1	141	60	0	3	1
2.5%-tile:	1	161	99	0	3	2644
25%-tile:	1	163	100	0	4	26436
Median: 	1	163	100	0	4	52871
75%-tile:	1	163	100	0	4	79306
97.5%-tile:	1	163	100	0	5	103098
Maximum:	3	163	100	0	8	105741
Mean:	1	162	99	0	4
# of unique seqs:	4031
total # of seqs:	105741
```

Now to classify sequences preliminarily to get rid of undesirable sequences such as archaea, chloroplasts, and mitochondria. I never made the smaller database so we will just use the full Silva NR database and the NR taxonomy files. That should still work it will just take longer. 
```
classify.seqs(fasta=wintertags.good.good.unique.good.good.filter.unique.precluster.pick.fasta, count=wintertags.good.good.unique.good.good.filter.unique.precluster.pick.count_table, reference=silva.nr_v132.align, taxonomy=silva.nr_v132.tax, cutoff=80)

# Remove the lineages of archaea, eukarya, but keep unclassified for now
remove.lineage(fasta=wintertags.good.good.unique.good.good.filter.unique.precluster.pick.fasta, count=wintertags.good.good.unique.good.good.filter.unique.precluster.pick.count_table, taxonomy=wintertags.good.good.unique.good.good.filter.unique.precluster.pick.nr_v132.wang.taxonomy, taxon=unknown;-Archaea;-Eukaryota;)

# Summary
summary.seqs(fasta=wintertags.good.good.unique.good.good.filter.unique.precluster.pick.pick.fasta, count=wintertags.good.good.unique.good.good.filter.unique.precluster.pick.pick.count_table)

		Start	End	NBases	Ambigs	Polymer	NumSeqs
Minimum:	1	141	60	0	3	1
2.5%-tile:	1	161	99	0	3	2644
25%-tile:	1	163	100	0	4	26436
Median: 	1	163	100	0	4	52871
75%-tile:	1	163	100	0	4	79306
97.5%-tile:	1	163	100	0	5	103098
Maximum:	3	163	100	0	8	105741
Mean:	1	162	99	0	4
# of unique seqs:	4031
total # of seqs:	105741
```

This didn't seem to remove any sequences, especially the unknown ones. Ben doesn't remove sequences in using TaxASS, so I think things will get removed in the end anyways. 

Rename the filenames for ease, which is done inside of mothur, filenames are now: 
```
fasta: WinterTags.final.fasta
count table: WinterTags.final.count_table
```

After this we will cluster OTUs inside of mothur and get representative OTU sequences for final classification with TaxASS first and then GreenGenes. This will probably take a lot of time for the clustering and distance calculations.

### Defining OTUs and performing distance calculations 

Defining Operational Taxonomic Units (OTUs) = sequence proxy for a microbial species by calculating distance between sequences/how much they diverge, then cluster the distance calculations with a certain cutoff point. We will be using the following cutoffs as proxies for each delimiter: 

- Species: 97% similar
- Genus: 95% similar
- Family: 90% similar

```
dist.seqs(fasta=WinterTags.final.fasta)
```

In the future when thinking about larger datasets, this function can have cutoffs set to decrease the distance matrix. 

We will use the average neighbor clustering method. There are other ways to cluster OTUs, which Ananke does it temporally for time series analyses. We only have a couple of samples for winter months that aren't really "time-series," but this is something we can think about in the future. 

```
cluster.split(column=WinterTags.final.dist, count=WinterTags.final.count_table, method=average)
``` 

Important note - there is a bug in the mothur version on the VM, which is the newest version, so I did the clustering on my personal computer and that seemed to work fine. Surprisingly this didn't take forever, probably because we only have 9 samples and a really reduced amount of sequences. We can now continue with classification using TaxASS and then GreenGenes for whatever TaxASS doesn't hit. 

The current files are now: 

- `WinterTags.final.an.unique_list.list`
- `WinterTags.final.count_table`
- `WinterTags.final.dist` 
- `WinterTags.final.fasta`

Which are now on my personal computer and need to be moved back to the virtual machine for further analyses. 

### Friday 2018-04-20

Now classifying the OTUs using a combination of greengenes and blast to identify with the Freshwater Training Set using TaxASS and Robin's scripts. The first thing is to make a BLAST database of the Freshwater Training set. I pulled the fasta file off of Ben's database folder of the fileshare, and I will have a copy myself. The command: 

```
makeblastdb -dbtype nucl -in FWonly_11Feb2016_1452_ready.fasta -input_type fasta -parse_seqids -out FWonly_11Feb2016_1452_ready.fasta.db
```

This creates a lot of database files for BLAST to use when searching our sets against the freshwater training set, and then parsing the percent identities to see if it's a good hit with the freshwater training set, or it should be classified with another reference set. Then blast against the freshwater database. First hyphens should be completely removed  and smooshed together with `sed -e 's/-//g' WinterTags.final.fasta > WinterTags.nohyph.fasta`

Then query the FW database and reformat the blast results: 

```
blastn -query WinterTags.final.fasta -task megablast -db FW_only11Feb2016_1452_ready.fasta.db -out WinterOTU.custom.blast -outfmt 11 -max_target_seqs 5
# reformat
blast_formatter -archive WinterOTU.custom.blast -outfmt "6 qseqid pident length qlen qstart qend" -out WinterOTU.custom.blast.table
```

Then using Robin's set of scripts to calculate full length PIDent, pull out sequence IDs greater than and less than 98%, some sanity plots in R. Steps for using Robin's Rscripts: 

```
# calculate full length 
Rscript calc_full_length_pident.R OTU.custom.blast.table OTU.custom.blast.table.modified
# above 98 match
Rscript filter_seqIDs_by_pident.R OTU.custom.blast.table.modified ids.above.98 98 TRUE
# below 98
Rscript filter_seqIDs_by_pident.R OTU.custom.blast.table.modified ids.below.98 98 FALSE
# pull out OTUs that hit nothing on blast
python find_seqIDs_blast_removed.py 5M.fasta OTU.custom.blast.table.modified ids.missing
# combine list of seequence IDs that didn't hit anythign on blast with those that had pident below 98
cat ids.below.98 ids.missing > ids.below.98.all
# create fasta file with sequences that correspond to each group
python create_fastas_given_seqIDs.py ids.above.98 fasta above.fasta
python create_fastas_given_seqIDs.py ids.below.98.all fasta below.fasta
```

Now the different fasta files are setup to classify the groups of OTUs using mothur with the two different databases. Sequences with a pident above 98 are classified using the FW database. And sequences below 98 are classified using the Greengenes database within mothur. The fasta and taxonomy files for both the freshwater database and Greengenes. 

Within mothur for both classifications: 

```
classify.seqs(fasta=WinterTags.above.98.fasta, template=FWonly_11Feb2016_1452_ready.fasta, taxonomy=FWonly_11Feb2016_1452.taxonomy, method=wang, probs=T, processors=10)
classify.seqs(fasta=WinterTags.below.98.fasta, template=../Gg.fasta, taxonomy=Gg.taxonomy, method=wang, probs=T. processors=10)
```

From a first glance, it doesn't take too terribly long to generate/search the FW database, but it takes a while to go through the GreenGenes database and search it, because it's a larger database and probably takes a bit more time to sift through. GreenGenes contains incomplete 16S sequences and therefore isn't good for alignment purposes, but in terms of classification purposes is more promiscuous for identifying uncultured taxa. 

Thinking about the input files for downstream R analyses and data exploration, you give it a unique list or shared file and the taxonomy file. So once greengenes is done running, concatenate the Greengenes and FW classifications, and give it the unique list table of OTUs. 

Also reformat the taxonomy file so delimited by semicolons: `sed 's/[[:blank:]]/\;/' <5M.taxonomy > 5M.taxonomy.reformatted` or `sed 's/[[:blank:]]/\;/' <WinterTags.taxonomy> > WinterTags.taxonomy.reformatted`. 

### Saturday 2018-04-21

Problem I have now run into that I didn't realize is that I should've normalized the list table and that is what is used to put into R. I either didn't save the shared table or didn't make it because it was near the normal classification steps. But I still have the final count table and the final list file, so that should be easy to make the shared table, which is what they use for rarefaction and normalization steps. To make the shared file within mothur: 

```
mothur > make.shared(list=WinterTags.final.an.unique_list.list, count=WinterTags.final.count_table, label=0.03)
```

This will give us the species level OTUs, although can you really identify species level at OTU I think not. However something is wrong with the current list file, which is generated from the `cluster.split` command. Can go back and regenerate that with: 

```
mothur > cluster.split(column=WinterTags.final.dist, count=WinterTags.final.count_table, method=average)
```

This is where I had issues with the newest version of mothur and had to do things on my local computer and that version of mothur. So try regenerating that file. 

Just realized that we forgot to run the classify.otu script that will turn the taxonomy file into a list of OTUs and not the sequence headers, and that's what is similar between the shared file and the taxonomy file. That command is: 

```
mothur > classify.otu(list=WinterTags.final.an.unique_list.list, taxonomy=WinterTags.taxonomy, count=WinterTags.final.count_table, label=0.03, cutoff=80, basis=otu, probs=F)
# reformat
sed 's/[[:blank:]]/\;/' WinterTags.final.an.unique_list.0.03.cons.taxonomy > WinterTags.taxonomyOTUs.reformatted
```

To make a phylogenetic tree of OTUs, need to pull out representative sequences for each OTU. This tree is also needed for further statistics and plotting in R. For getting representative sequences, you can either calculate the best possible representative sequence fo reach OTU with the smallest maximum distance to the other sequences in an OTU group, or choose the most abundant unique sequence within each OTU group as the representative - this runs much faster and with less RAM than the other method, not the best for all uniques within an OTU. So I will run the most abundant one for now: 

```
mothur > get.oturep(list=WinterTags.final.an.unique_list.list, fasta=WinterTags.final.fasta, count=WinterTags.final.count_table, label=0.03, method=abundance)
```

Now there is the problem of sequences not being in the list file and the count file, so regenerate the list file by redoing the distance matrix and whatnot. Which I thought I already did this with the distance file, but it might be because I did some of this with different versions of mothur. 

Ok so it looks like I have to re-run everything because something is very wrong with the generated count table and I don't know which step that happened at. I am starting from the very beginning of the preprocessing steps and starting with the raw fastq files. So before the classification steps, make sure the representative OTUs can be pulled out correctly before going through all that mess. 


