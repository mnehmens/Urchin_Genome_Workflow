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

### Flag PCR Duplicates with SAMBLASTER 
samblaster -i haps_hic.sam -o haps_hic_marked_byread.sam

### Remove unmmaped and non-primary aligned reads, sort and index bam files
samtools view -S -b -h -@ 12 -F 2316 haps_hic_marked_byread.sam > haps_hic_presort_marked.bam
samtools sort -@ 12 haps_hic_presort_marked.bam -o haps_hic.bam

# Create .bed file containing each contig and start and end position - use .fai file to do this.
awk '{print $1, 1, $2}' haps_purged_wmg.np2.fa.fai > hap1.bed
# duplicated bed file, removed so appropriate haplotype in each, so have hap1.bed hap2.bed

# Create bam with hap1 and hap 2 bed files
# load SAMtools
samtools view -L hap1.bed -bh haps_hic.bam > hap1_sorted.bam
samtools view -L hap2.bed -bh haps_hic.bam > hap2_sorted.bam

# in order to write the fastq files, you need to sort by read name 
samtools sort -@ 16  -n -o hap1_sorted.name.bam hap1_sorted.bam
samtools sort  -@ 16  -n -o hap2_sorted.name.bam hap2_sorted.bam

# Having sorted, we then get samtools fastq to write out the reads
samtools fastq  -@16 -s hap1_singletons.fastq -1 hap1_paired_1.fastq -2 hap1_paired_2.fastq  hap1_sorted.name.bam
samtools fastq  -@16 -s hap2_singletons.fastq -1 hap2_paired_1.fastq -2 hap2_paired_2.fastq  hap2_sorted.name.bam

# Repeat steps, but one hap at a time
# laod SAMtools BWA
bwa mem -5SP -t 16 hap1_purged_wmg.np2.fa hap1_paired_1.fastq hap1_paired_2.fastq > hap1_remapped.sam
samtools view -S -b -h -@ 16 -F 2316 hap1_remapped.sam > hap1_presort.bam
samtools sort -@ 16 hap1_presort.bam -o hap1_sorted.bam

bwa mem -5SP -t 16 hap2_purged_wmg.np2.fa hap2_paired_1.fastq hap2_paired_2.fastq > hap2_remapped.sam
samtools view -S -b -h -@ 16 -F 2316 hap2_remapped.sam > hap2_presort.bam
samtools sort -@ 16 hap2_presort.bam -o hap2_sorted.bam

# Now can go into normal YAHS pipeline
yahs --no-mem-check hap1_purged_wmg.np2.fa hap1_sorted.bam -o hap1_yahs4
yahs --no-mem-check hap2_purged_wmg.np2.fa hap2_sorted.bam -o hap2_yahs4

# Haps YAHS juicer for JBAT
yahs/juicer pre -a -o hap1_JBAT hap1_yahs4.bin hap1_yahs4_scaffolds_final.agp hap1_purged_wmg.np2.fa.fai >hap1_JBAT.log 2>&1
cat hap1_JBAT.log  | grep PRE_C_SIZE | awk '{print $2" "$3}' >hap1_JBAT.log.chrom.size
yahs/juicer pre -a -o hap2_JBAT hap2_yahs4.bin hap2_yahs4_scaffolds_final.agp hap2_purged_wmg.np2.fa.fai >hap2_JBAT.log 2>&1
cat hap2_JBAT.log  | grep PRE_C_SIZE | awk '{print $2" "$3}' >hap2_JBAT.log.chrom.size

### Juicer tools, for HiC file ####
java -Xmx36G -jar juicer/CPU/common/juicer_tools.jar pre hap1_JBAT.txt hap1_JBAT.hic hap1_JBAT.log.chrom.size
java -Xmx36G -jar juicer/CPU/common/juicer_tools.jar pre hap2_JBAT.txt hap2_JBAT.hic hap2_JBAT.log.chrom.size

# produce the final step of juicer using the juicebox output 
yahs/juicer post -o test_hap2_JBAT TEST_hap2_JBAT.review.assembly hap2_JBAT.liftover.agp hap2_purged_wmg.np2.fa
yahs/juicer post -o hap1_JBAT hap1_JBAT.review.assembly hap1_JBAT.liftover.agp hap1_purged_wmg.np2.fa

# Now in scaffolds by haplotype, change name to easier identify
grep ">" hap1_JBAT.FINAL.fa > original.names.hap1_JABT.FINAL.txt
grep ">" hap2_JBAT.FINAL.fa > original.names.hap2_JABT.FINAL.txt

# Resave as new.names.hap*_JBAT.FINAL.txt and open in editor to put new names removing the ">" before the new name BE CAREFUL TO KEEP ORDER SAME
seqkit fx2tab hap1_JBAT.FINAL.fa | cut -f 2 | paste new.names.hap1_JBAT.FINAL.txt - | seqkit tab2fx > hap1_JBAT.renamed.FINAL.fa
seqkit fx2tab hap2_JBAT.FINAL.fa | cut -f 2 | paste new.names.hap2_JBAT.FINAL.txt - | seqkit tab2fx > hap2_JBAT.renamed.FINAL.fa

# Now
hap1_JBAT.FINAL.agp test_hap2_JBAT.FINAL.agp
agptools rename rename_agp_hap1.txt hap1_JBAT.FINAL.agp > hap1_JBAT.renamed.FINAL.agp
agptools rename rename_agp_hap2.txt hap2_JBAT.FINAL.agp > hap2_JBAT.renamed.FINAL.agp

# QC check
# merqury and compleasm - looks good
# qualimap, map ONT first on concatenated haps
minimap2 -t 24 -ax lr:hqae haps_JBAT.renamed.FINAL.fa assembly_input.fastq.gz assembly_input2.fastq.gz > haps_JBAT_ONTmap.sam
samtools sort -@ 32 -T ali.tmp haps_JBAT_ONTmap.sam > haps_JBAT_ONTmap.bam

#in qualimap script needed to change MaxPermSize=1024m to MaxMetaspaceSize re:suggestion online to get newer java to run
./qualimap bamqc -bam haps_JBAT_ONTmap.bam -outdir results --java-mem-size=32G
```
#### quarTeT

```bash
# Use quarTeT to find telomeres and centromeres
# Telominer and Centrominer (repeat for hap2)
# load Python minimap2 MUMmer trf BLAST gnuplot R
python3 ./quartet.py TeloExplorer -i hap1_JBAT.renamed.FINAL.fa -c animal -p hap1_renamed
python3 ./quartet.py CentroMiner -i hap1_JBAT.renamed.FINAL.fa -p hap1_renamed
```
#### Manual Curation

```bash
# Trim excessive N's due to introduction during dual scaffold option assembly
# Use combined haplotypes
# Remove the second mt genome in the combined haps fasta
seqtk subseq haps_JBAT.renamed.FINAL.fa haps_rm_2ndmtg.txt > haps_1mtg_FINAL.fa
# Trim N's using custom python script, replacing with 10,000 bp
python trim_ns.py haps_1mtg_FINAL.fa haps_1mtg_trimN.fa 10000

# Rename scaffolds into pesudo-chromosomes
# use cut to isolate the columns you want, generated an alternative list of scaffold IDs to replace the old ones, using the same order that they are in the fasta file and then seqkit tab2fx to bring it back together
seqkit fx2tab hap1_1mtg_FINAL_trimN.fa | cut -f 2 | paste hap1_newChromName.txt - | seqkit tab2fx > Crod1.0_chrom_hap1.fa
seqkit fx2tab hap2_1mtg_FINAL_trimN.fa | cut -f 2 | paste hap2_newChromName.txt - | seqkit tab2fx > Crod1.0_chrom_hap2.fa


# Checking it worked correctly 
grep ">" Crod1.0_chrom_hap1.fa > hap1_NEWTEST.txt
grep ">" Crod1.0_chrom_hap2.fa > hap2_NEWTEST.txt
# check in DGENIES -- all looks good, use to figure out orientation, and reverse compliment any not in same direction

# Need to reverse compliment some scaffolds to match both haps. Generate a single fasta per contig using seqkit split -i
seqkit split -i Crod1.0_chrom_hap2.fa
seqkit split -i Crod1.0_chrom_hap1.fa

# Doing reverse compliment, then will combine 
seqtk seq -r $i > ${i%.fa}.rc.fa

# Running the revesre compliment for the correct files e.g. for hap2 chromosome 3 to reverse compliment (rc)
seqtk seq -r Crod1.0_chrom_hap2.part_Crod_chr03_h2.fa > Crod1.0_chrom_hap2.part_Crod_chr03_h2.rc.fa

# Cat back together all scaffolds in order e.g. hap1_crhom1, chrom2, chrom3, etc...
# Check in DGENIES summary out file, run seqkit stats between:
# hap1_1mtg_FINAL_trimN.fa and Crod1.0_chrom_hap1.fa
# hap2_1mtg_FINAL_trimN.fa, Crod1.0_chrom_hap2.fa and with reverse compliment Crod1.0_chrom_hap2_wRC.fa
# Check order to make sure everything was put in correct order
grep ">" Crod1.0_chr_mt_hap1.fa > testorderhap1.txt
grep ">" Crod1.0_chr_mt_wrc_hap2.fa > testorderhap2.txt
# all looks good, proceed

# Add back in mt genome per haplotype
e.g. cat Crod1.0_chrom_hap1.fa mito_gnome.fa > Crod1.0_chr_mt_hap1.fa

# Final quality check at this step before annotation 
quast Crod1.0_chr_mt_hap1.fa -o quasthap1
quast Crod1.0_chr_mt_wrc_hap2.fa -o quasthap2
# looks good, go to EarlGrey
```

#### EarlGrey
```bash
# Run EarlGrey, complete for both haplotypes
#define the path to container image and export it to a variable
export CONTAINER_IMG="PATH/TO/earlg-working/earlgrey-4.2.4.aimg"

#simplify the apptainer exec command by exporting it to a variable. 
#This will help with reduce the length of earlGrey functional command
export CMD="apptainer exec ${CONTAINER_IMG}"
 ${CMD} earlGrey -g Crod1.0_chr_mt_hap1.fa -s centrostephanusRodgersii -o ./EarlyGrey_Hap1 -t 16

# Using the output *-families.fa.strained and combing by haplotypes
cat hap1_-families.fa.strained hap2_-families.fa.strained > haps_-families.fa.strained

# Use CD-HIT est - tested different parameters, this was best
cd-hit-est -i haps_families.fa.strained -o est_haps.fa -aS 0.8 -c 0.95 -G 0 -n 10 -M 24000 -T 8








