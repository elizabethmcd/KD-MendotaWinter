#! /bin/bash

# RRR 7-30-16
# download fastq from multiple SRR accession numbers at once
# file "SraAccList.txt" is downloaded from NCBI Sequence Read Archive, top right, doesn't look like button but it is.

# syntax in terminal:
# ./download_SRRs.sh SraAccList.txt


AccFile="$1"

echo "$AccFile"

while read "SRR"; do
   echo "$SRR"
   fastq-dump --split-files --skip-technical --outdir ./ "$SRR"
done < "$AccFile"

echo "ta da"

exit

# How the script works:
# $1 is a built in bash variable that is the first argument following the script name when sourced.
# SRR is a variable that I defined. I named it SRR b/c each one is a SRA Run Accession Number. ex: SRR1531662
# the file is not called until the end of the while loop.  apparently this "slurps" it in.
# cannot have $ when variables are defines, only use $ when they're called.
# but why can use quotes when define SRR but not when define AccFile ?
# read is a bash function.  It reads one line of a file at a time.