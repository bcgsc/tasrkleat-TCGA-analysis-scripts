#!/bin/bash

wget --no-clobber http://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa.gz
wget --no-clobber http://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz

gunzip -v Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa.gz
gunzip -v Homo_sapiens.GRCh37.75.gtf.gz

if which samtools > /dev/null; then
    echo "indexing fasta file..."
    samtools faidx Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa

    echo "indexing gtf file..."
    (\grep ^"#" Homo_sapiens.GRCh37.75.gtf; \grep -v ^"#" Homo_sapiens.GRCh37.75.gtf | sort -k1,1 -k4,4n) | bgzip > Homo_sapiens.GRCh37.75.sorted.gtf.gz;
    tabix Homo_sapiens.GRCh37.75.sorted.gtf.gz 
else
    echo 'please install samtools and run `samtools faidx Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa`'
fi
