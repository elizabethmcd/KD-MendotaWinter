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

Now aligning sequences to a reference alignment with `align.seqs`. We will be using Silva wiht the full-length sequences. I will align to the full database just because I'm not too familiar with this dataset, and this might be a different version of Silva. Will probably end up with the normal V4 region anyways, but again because of version issues not 100% sure. The full Silva NR database is about 10GB and that has been moved to `/home/emcdaniel/Databases/Silva.nr_v132` on the fileshare for easy access in the future. This will be moved to the current working directory where we are running mothur, and then use the reduced database with our region of interest for future runs.  

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