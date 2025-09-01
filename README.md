## Centrostephanus rodgersii Genome Assembly Workflow


#### HiFi-ASM with RAFT
```bash
#First run of hifiasm with 4 threads to obtain error corrected reads and coverage estimate
hifiasm -o errorcorrect -t 32 --write-ec HQ_assembly_input.fastq.gz 2> errorcorrect.log
COVERAGE=$(grep "homozygous" errorcorrect.log | tail -1 | awk '{print $6}')

# Second run of hifiasm to obtain all-vs-all read overlaps as a paf file
hifiasm -o getOverlaps -t 32 --dbg-ovec errorcorrect.ec.fa 2> getOverlaps.log

# Merge cis and trans overlaps
cat getOverlaps.0.ovlp.paf getOverlaps.1.ovlp.paf > overlaps.paf
COVERAGE=40

# RAFT fragments the error corrected reads
raft -e ${COVERAGE} -o fragmented errorcorrect.ec.fa overlaps.paf

# Final hifiasm run to obtain assembly of fragmented reads, single round of errore correction (-r1)
hifiasm -o manual --dual-scaf --ul fastq.gz  --h1 R1.fq.gz --h2 R2.fq.gz   -t 32 -r1 fragmented.reads.fasta 2> finalasm.log
```

#### Mitochondrial genome
```bash
# make blast database
makeblastdb -in genome.fa -dbtype nucl -parse_seqids
# use genbank Strongylocentrotus purpuratus mtgenome as query  
blastn -query sc_mito.fa -db genome.fa -out mito_search -outfmt 6

# use samtools to see single hit contig properties
samtools faidx genome.fa
grep "contig_name" genome.fa.fai
 
# pull this contig out of the genome fasta (file.list contains the tg identifier) 
nano file.list
seqtk subseq genome.fa file.list > mt_contig.fa

# can now use to include one single version in haplotype genome assemblies
# start by combining haplotypes
cat hap1.fa hap2.fa > haps.fa

# make new list of contigs removing contig 6 and 4 that are the mito genome
seqtk subseq haps.fa haps_wo_mito.txt > haps_nomg.fa

# then cat no mg with the mito genome nomg/wo_mg = no mito genome, wmg = with mito genome
cat haps_nomg.fa mito_gnome.fa > haps_wmg.fa
```

#### Haplotypes purge_dups
```bash
# run first two lines of purge_dups then pull apart and run each hap independently
bwa index haps.fa
bwa mem -t 32 haps.fa AG0956_R1_trimmed.fastq.gz AG0956_R2_trimmed.fastq.gz > haps_align_Ill.sam

# index haplotypes, create bed to split files to be used for each purge dups on each hap
awk '{OFS="\t"; print $1, "0", $2}' hap1.fa.fai > hap1.bed
awk '{OFS="\t"; print $1, "0", $2}' hap2.fa.fai > hap2.bed

# first output from purgedups step: align indexed haps to a sam file 
samtools view -h -o hap1.sam -L hap1.bed haps_align_Ill.sam 
samtools view -h -o hap2.sam -L hap2.bed haps_align_Ill.sam 

# to check if this step was done correctly
samtools flagstat hap1.sam 
samtools flagstat hap2.sam

# filter before putting into rest of purge_dups pipeline
samtools view -h -F 256  hap1.sam > hap1_filt.sam
samtools view -h -F 256  hap2.sam > hap2_filt.sam

## check again
samtools flagstat hap1_filt.sam
samtools flagstat hap2_filt.sam

### hap1 purge_dups
# need to load BWA minimap2 SAMtools purge_dups 

# move to next step after aligning, filtering
samtools view -bh hap1_filt.sam > hap1_filt.bam
ngscstat hap1_filt.bam
calcuts TX.stat > cutoffs 2>calcults.log
split_fa hap1.fa > hap1.split
minimap2 -xasm5 -DP hap1.split hap1.split | gzip -c - > hap1.split.self.paf.gz
purge_dups -2 -T cutoffs -c TX.base.cov hap1.split.self.paf.gz > dups.bed 2> purge_dups.log
get_seqs dups.bed -e hap1.fa

### repeat for hap2
samtools view -bh hap2_filt.sam > hap2_filt.bam
ngscstat hap2_filt.bam
calcuts TX.stat > cutoffs 2>calcults.log
split_fa hap2.fa > hap2.split
minimap2 -xasm5 -DP hap2.split hap2.split | gzip -c - > hap2.split.self.paf.gz
purge_dups -2 -T cutoffs -c TX.base.cov hap2.split.self.paf.gz > dups.bed 2> purge_dups.log
get_seqs dups.bed -e hap2.fa
```
#### NextPolish2 
```bash
# combine haplotypes
cat hap1_wmg hap1_wmg > haps_purged_wmg.fa

# need to load Merqury Winnowmap SAMtools
# NextPolish2 step 1 to create input files

# create k15 first
meryl count k=15 output merylDB haps_purged_wmg.fa 
meryl print greater-than distinct=0.9998 merylDB > repetitive_haps_purg_wmg_k15.txt
winnowmap -k 15 -W repetitive_haps_purg_wmg_k15.txt -ax map-ont haps_purged_wmg.fa HQ_assembly_input.fastq.gz UL_20kb_Q10.fastq.gz | samtools view -hb | samtools sort > haps_purg_wmg_k15.bam

#yak for k21
./yak/yak count -k 21 -b 37 -t 16 -o k21.yak <(zcat AG0956_R1_trimmed.fastq.gz) <(zcat AG0956_R2_trimmed.fastq.gz)
#yak for k31
./yak/yak count -k 31 -b 37 -t 16 -o k31.yak <(zcat AG0956_R1_trimmed.fastq.gz) <(zcat AG0956_R2_trimmed.fastq.gz)

# NextPolish2 step 2
# need to load Perl Rust yak OpenSSL CMake
nextPolish2 -t 8 haps_purg_wmg_k15.bam haps_purged_wmg.fa k21.yak k31.yak > haps_purged_wmg.np2.fa

## QC step
# check using merqury qv, need to split by haplotypes
seqtk subseq haps_purged_wmg.np2.fa hap1_wmg_splitnp2.txt > hap1_purged_wmg.np2.fa
seqtk subseq haps_purged_wmg.np2.fa hap2_wmg_splitnp2.txt > hap2_purged_wmg.np2.fa

# need meryl database first
meryl k=19 count Combined_AG0956_trimmed.fastq.gz output Combined_AG0956_trimmed.meryl
merqury.sh Combined_AG0956_trimmed.meryl  hap1_purged_wmg.np2.fa hap2_purged_wmg.np2.fa haps_19mer

# use this file to see sub 50 qv contigs 
hap1_purged_wmg.np2.fa.qv hap2_purged_wmg.np2.fa.qv
# created txt file with resulting contigs using above file - check sub 50 QV small contigs 
# make new list of contigs removing all small contigs and put into hap1_purged_wmg_np2.txt file 
seqtk subseq hap1_purged_wmg.np2.fa haps_purged_wmg_np2.txt > hap1_purged_nomg.np2.fa
# repeat for hap 2
### merqury qv looks good, continue on using these to assemblies for HiC
haps_purged_wmg.np2.fa, hap1_purged_wmg.np2.fa, hap2_purged_wmg.np2.fa

# run compleasm on both haps
compleasm.py run -a hap1_purged_wmg.np2.fa -o out_directory -l metazoa -t 8
# looks good, continue to HiC
```
#### Hi-C using YAHS pipeline

``` bash
# running first steps for YAHS
# cat haplotypes together haps_purged_wmg.np2.fa

# load SAMtools BWA samblaster
# align reads with BWA using -5sp for Hi-C reads
samtools faidx haps_purged_wmg.np2.fa
bwa mem -5SP -t12 haps_purged_wmg.np2.fa Urchin_HiC_combined_filtered_R1.fq.gz Urchin_HiC_combined_filtered_R2.fq.gz -o haps_hic.sam


