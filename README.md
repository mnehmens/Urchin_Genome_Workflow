## Urchin_Genome_Workflow

### Steps for genome assembly of C. rodgersii
```bash
# Raw ONT data trimmed adapters, ready for use from AW: folder Recalled_114_Apr2023
# Used NanoPlot to visualize data: 
NanoPlot --fastq Combined_pass.pc.fastq.gz   --huge  --plots kde --color lightgray --outdir /nesi/nobackup/ga03714/Melissa/Fastq/Recalled_114_Apr2023/read_qc$ -t 2
NanoPlot --fastq Combined_fail.pc.fastq.gz   --huge  --plots kde --color lightgray --prefix failed_reads_test --outdir /nesi/nobackup/ga03714/Melissa/Fastq/Recalled_114_Apr2023/read_qc$ -t 2
```
Evaluated desired cutoffs for quality and size, decided to combine pass/fail files 

```bash
# Combined_pass.pc.fastq.gz, Combined_fail.pc.fastq.gz into one file Combined_pass_and_fail.pc.fastq.gz
cat Combined_pass.pc.fastq.gz Combined_fail.pc.fastq.gz > Combined_pass_and_fail.pc.fastq.gz
```

Used combined file to filter into binned sizes first 
```bash
gunzip -c Combined_pass_and_fail.pc.fastq.gz | chopper --minlength 50000 | gzip > All_50kbplus.fastq.gz
gunzip -c Combined_pass_and_fail.pc.fastq.gz | chopper --minlength 20000 --maxlength 49999 | gzip > All_20kb_to_50kb.fastq.gz
gunzip -c Combined_pass_and_fail.pc.fastq.gz | chopper --minlength 10000 --maxlength 19999 | gzip > All_10kb_to_20kb.fastq.gz
gunzip -c Combined_pass_and_fail.pc.fastq.gz | chopper --minlength 5000 --maxlength 9999 | gzip > All_5kb_to_10kb.fastq.gz
gunzip -c Combined_pass_and_fail.pc.fastq.gz | chopper --minlength 1000 --maxlength 4999 | gzip > All_1kb_to_5kb.fastq.gz
gunzip -c Combined_pass_and_fail.pc.fastq.gz | chopper  --maxlength 999 | gzip > All_sub1kb.fastq.gz
```

Next used size files to filter by quality 
```bash
gunzip -c All_50kbplus.fastq.gz | chopper --quality 15 | gzip > Q15_50kbplus.fastq.gz
gunzip -c All_20kb_to_50kb.fastq.gz | chopper --quality 15 | gzip > Q15_20kb_to_50kb.fastq.gz
gunzip -c All_10kb_to_20kb.fastq.gz | chopper --quality 15 | gzip > Q15_10kb_to_20kb.fastq.gz
gunzip -c All_5kb_to_10kb.fastq.gz | chopper --quality 15 | gzip > Q15_5kb_to_10kb.fastq.gz
gunzip -c All_50kbplus.fastq.gz | chopper --quality 7 | gzip > Q7_50kbplus.fastq.gz
gunzip -c All_20kb_to_50kb.fastq.gz | chopper --quality 7 | gzip > Q7_20kb_to_50kb.fastq.gz
```

### Summary stats for each file:
```bash
All_50kbplus.fastq.gz
All_20kb_to_50kb.fastq.gz
All_10kb_to_20kb.fastq.gz
All_5kb_to_10kb.fastq.gz
All_1kb_to_5kb.fastq.gz
All_sub1kb.fastq.gz
Q15_50kbplus.fastq.gz
Q15_20kb_to_50kb.fastq.gz
Q15_10kb_to_20kb.fastq.gz
Q15_5kb_to_10kb.fastq.gz
Q7_50kbplus.fastq.gz
Q7_20kb_to_50kb.fastq.gz
```

We use seqkit
```bash
ml SeqKit 
for i in All*.fastq.gz
do
seqkit stats $i >> allStats.txt
done

for i in Q*.fastq.gz
do 
seqkit stats $i
done
```

Used quality and size filtered files to feed into assmeblers Shasta, FLYE and HiFiASM, example of each:

### Shasta

```bash
#!/bin/bash -e
#SBATCH --job-name=shasta_Q15_5Kbgreater
#SBATCH --mail-user=m.nehmens@massey.ac.nz
#SBATCH --mail-type=END
#SBATCH --time=10:00:00
#SBATCH --mem=960G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --profile=task
#SBATCH --account=ga03714

module purge
cd /nesi/nobackup/ga03714/Melissa/Fastq/Recalled_114_Apr2023/

/nesi/nobackup/ga03714/Melissa/Software/Shasta/shasta-Linux-0.11.1 --input Q15_50kbplus.fastq Q15_20kb_to_50kb.fastq Q15_10kb_to_20kb.fastq Q15_5kb_to_10kb.fastq --config Nanopore-R10-Fast-Nov2022 --assemblyDirectory /nesi/nobackup/ga03714/Melissa/Assemblies/Shasta/Shasta_Q15_5kbplus_120423
```

### FLYE
```bash
#!/bin/bash -e
#SBATCH --job-name=FLYEq15
#SBATCH --mail-user=m.nehmens@massey.ac.nz
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#SBATCH --mail-type=END
#SBATCH --time=48:00:00
#SBATCH --mem=400G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --profile=task
#SBATCH --account=ga03714

module purge
module load Flye/2.9-gimkl-2020a-Python-3.8.2
cd /nesi/nobackup/ga03714/Melissa/Fastq/Recalled_114_Apr2023/

flye --nano-hq Q15_50kbplus.fastq.gz Q15_20kb_to_50kb.fastq.gz Q15_10kb_to_20kb.fastq.gz Q15_5kb_to_10kb.fastq.gz --out-dir /nesi/nobackup/ga03714/Melissa/Assemblies/FLYE/flye_Q15_5Kb_11April  -t 32
```

###HiFiASM
```bash
#!/bin/bash -e
#SBATCH --job-name=hifi_Q15_l0
#SBATCH --output=MCN_%j.out
#SBATCH --error=MCN_%j.err
#SBATCH --mail-user=m.nehmens@massey.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --time=4:00:00
#SBATCH --mem=75G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --profile=task
#SBATCH --account=ga03714

module purge
module load hifiasm/0.15.5-GCC-9.2.0
cd /nesi/nobackup/ga03714/Melissa/Assemblies/Hifiasm

hifiasm -o Q15_5kbPlus_l0_01May2023 -t 32 Q15_10kb_to_20kb.fastq.gz  Q15_20kb_to_50kb.fastq.gz  Q15_50kbplus.fastq.gz  Q15_5kb_to_10kb.fastq.gz  
```

###BUSCO used to get metrics on assemblies, example:
```bash
#!/bin/bash -e
#SBATCH --job-name=busco_FLYEQ155Kb
#SBATCH --output=MCN_%j.out
#SBATCH --error=MCN_%j.err
#SBATCH --mail-user=m.nehmens@massey.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --time=24:00:00
#SBATCH --mem=18G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --account=ga03714

module purge
module load BUSCO/5.3.2-gimkl-2020a
cd /nesi/nobackup/ga03714/Melissa/Assembly_QC

busco -i /nesi/nobackup/ga03714/Melissa/Assemblies/FLYE/flye_Q15_5Kb_11April/assembly.fasta -c 8 -o flye_Q15_5Kb_11April_busco5.3.2 -m genome -l metazoa
```

###Illumina data
Arrived, need to trim adaptors using fastp, needed both full pathway and new output file directory:
```bash
#!/bin/bash -e
#SBATCH --job-name=fastp_Illumina
#SBATCH --output=MCN_%j.out
#SBATCH --error=MCN_%j.err
#SBATCH --mail-user=m.nehmens@massey.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --time=3:00:00
#SBATCH --mem=12G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --account=ga03714

module purge

module load fastp/0.23.2-GCC-11.3.0

cd /nesi/nobackup/ga03714/Melissa/Illumina
fastp -w 4 -i /nesi/nobackup/ga03714/Melissa/Illumina/AG0956_001_S469_R1_001.fastq.gz -I /nesi/nobackup/ga03714/Melissa/Illumina/AG0956_001_S469_R2_001.fastq.gz -o /nesi/nobackup/ga03714/Melissa/Illumina/trimmed/AG0956_R1_trimmed.fastq.gz -O /nesi/nobackup/ga03714/Melissa/Illumina/trimmed/AG0956_R2_trimmed.fastq.gz
```

###Filtered trimmed Illumina

jellyfish:
```bash
#!/bin/bash -e
#SBATCH --job-name=jellyfish
#SBATCH --output=MCN_%j.out
#SBATCH --error=MCN_%j.err
#SBATCH --mail-user=m.nehmens@massey.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --time=8:00:00
#SBATCH --mem=50G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --profile=task
#SBATCH --account=ga03714

module purge
module load Jellyfish/2.3.0-gimkl-2020a Trimmomatic/0.39-Java-1.8.0_144
cd /nesi/nobackup/ga03714/Melissa/Illumina/trimmed

# need to concatenate read files (this is post fastp trimming of adapters)
cat AG0956_R1_trimmed.fastq.gz  AG0956_R2_trimmed.fastq.gz > All_kina_trimmed_reads.fq.gz
# use trimmomatic to obtain only 140bp reads, and write uncompressed output
trimmomatic SE -threads 10 All_kina_trimmed_reads.fq.gz All_kina_140bptrimmed_reads.fq CROP:140 MINLEN:140
# Count kmers
jellyfish count -C -m 21 -s 1000000000 -t 10  -o All_kina_140bptrimmed_reads.jf All_kina_140bptrimmed_reads.fq
# Create histogram
jellyfish histo -t 10 All_kina_140bptrimmed_reads.jf > All_kina_140bptrimmed_reads.histo

# Seqkit on new trimmed Illumina file 
seqkit stats All_kina_140bptrimmed_reads.fq >> IlluminaStats.txt
```

###[Genomescope]( http://genomescope.org/analysis.php?code=KmHYaCKuuK8wLs21dAt3 "website for Illumina data run with genomescope")
Change read length to 140, load All_kina_140bptrimmed_reads.histo Illumina file and run.<br>
Started mapping with purge_haplotigs (code below), but realized need to use purge_dups
```bash
#!/bin/bash -e
#SBATCH --job-name=mapping_BWA_samtools
#SBATCH --output=MCN_%j.out
#SBATCH --error=MCN_%j.err
#SBATCH --mail-user=m.nehmens@massey.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --time=8:00:00
#SBATCH --mem=20G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --profile=task
#SBATCH --account=ga03714

module purge
module load BWA/0.7.17-gimkl-2017a purge_haplotigs/1.1.2-gimkl-2022a-Perl-5.34.1
cd /nesi/nobackup/ga03714/Melissa/Mapping/

bwa index /nesi/nobackup/ga03714/Melissa/Mapping/Q15_5kb_130423.bp.p_ctg.fa

bwa mem -t 12 /nesi/nobackup/ga03714/Melissa/Mapping/Q15_5kb_130423.bp.p_ctg.fa AG0956_R1_trimmed.fastq.gz AG0956_R2_trimmed.fastq.gz > illumina_mapped_to_Q15_5kb_130423.sam

samtools sort -@ 12 -T ali.tmp illumina_mapped_to_Q15_5kb_130423.sam > purgehap_27April.bam

samtools flagstat purgehap_27April.bam > purgehap_27April.flagstat.txt

purge_haplotigs readhist -b purgehap_27April.bam -g /nesi/nobackup/ga03714/Melissa/Mapping/Q15_5kb_130423.bp.p_ctg.fa -t 24
```

### KAT
compare to see what the assembly looks like 
```bash
#!/bin/bash -e
#SBATCH --job-name=KAT_test_asmb10
#SBATCH --output=MCN_%j.out
#SBATCH --error=MCN_%j.err
#SBATCH --mail-user=m.nehmens@massey.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --time=10:00:00
#SBATCH --mem=200G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --profile=task
#SBATCH --account=ga03714

module purge
module load KAT/2.4.2-gimkl-2018b-Python-3.7.3
cd /nesi/nobackup/ga03714/Melissa/Compare

kat comp -t 24 -o KAT_asmb10_04May 'Illumina_R1_trimmed.fastq Illumina_R2_trimmed.fastq' assembly10.fasta

# Created symbolic links to make life easier when running Merqury 
ln -s /nesi/nobackup/ga03714/Melissa/Assemblies/FLYE/assembly3_18April/assembly.fasta assembly3.fasta
```

### Merqury and meryl
#!/bin/bash -e
#SBATCH --job-name=merqury_asmb3
#SBATCH --mail-user=m.nehmens@massey.ac.nz
#SBATCH --output=MCN_%j.out
#SBATCH --error=MCN_%j.err
#SBATCH --mail-type=ALL
#SBATCH --time=5:00:00
#SBATCH --mem=20G
#SBATCH --ntasks=1
#SBATCH --profile=task
#SBATCH --account=ga03714


module purge
export PATH=/nesi/nobackup/ga03714/Melissa/Software/Merqury/meryl-1.4/bin:/nesi/nobackup/ga03714/Melissa/Software/Merqury/merqury:$PATH
export MERQURY=/nesi/nobackup/ga03714/Melissa/Software/Merqury/merqury
ml SAMtools BEDTools R

cd /nesi/nobackup/ga03714/Melissa/Software/Merqury/merqury

merqury.sh AG0956_trimmed_test1.4.meryl assembly3.fasta Flye_asmb3_test

# Purge_dups workflow, sorted out, just waiting to get to assembly we want 

#!/bin/bash -e
#SBATCH --job-name=purgeDups_getseq
#SBATCH --output=MCN_%j.out
#SBATCH --error=MCN_%j.err
#SBATCH --mail-user=m.nehmens@massey.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --time=8:00:00
#SBATCH --mem=24G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --profile=task
#SBATCH --account=ga03714

module purge
module load purge_dups/1.2.6-gimkl-2022a-Python-3.10.5 
cd /nesi/nobackup/ga03714/Melissa/Mapping

# Step 1 aligns primary assembly and ONT data
minimap2 -x map-ont Q15_5kb_130423.bp.p_ctg.fa All_Q15_fastq.gz > HiFi130423_Q15.paf.gz
# creates .base.cov and .stat files
pbcstat HiFi130423_Q15.paf.gz
# calculates read-depth and cutoffs
calcuts PB.stat > cutoffs 2> calcults.log
# split the primary assembly 
split_fa Q15_5kb_130423.bp.p_ctg.fa > HiFi130423.split
# Uses the split primary assembly and does a self assignment
minimap2 -xasm5 -DP HiFi130423.split HiFi130423.split > HiFi130423.split.self.paf.gz
# purges the haplotigs, gives a duplication bed file and log of purged haplotigs - do we need a -e flag here?
purge_dups -2 -T cutoffs -c PB.base.cov HiFi130423.split.self.paf.gz > dups.bed 2> purge_dups.log
# taking the duplication bed file and the original assembly and getting the purged (primary and haplotig) sequences - results in purged.fa (use this) and hap.fa (alternatives, might be useful)
get_seqs dups.bed Q15_5kb_130423.bp.p_ctg.fa

 # Below gives same results for the output purged.fa, using the -c flag, just wanted to see if a difference
#!/bin/bash -e
#SBATCH --job-name=purgeDups_align-c
#SBATCH --output=MCN_%j.out
#SBATCH --error=MCN_%j.err
#SBATCH --mail-user=m.nehmens@massey.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --time=1:00:00
#SBATCH --mem=12G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --profile=task
#SBATCH --account=ga03714

module purge
module load purge_dups/1.2.6-gimkl-2022a-Python-3.10.5 
cd /nesi/nobackup/ga03714/Melissa/Mapping

# Step 1 aligns primary assembly and ONT data
#minimap2 -x map-ont Q15_5kb_130423.bp.p_ctg.fa All_Q15_fastq.gz | gzip -c -> HiFi130423_Q15-c.paf.gz
# creates .base.cov and .stat files
#pbcstat HiFi130423_Q15-c.paf.gz
# calculates read-depth and cutoffs? (or whatever .stat output is)
calcuts PB.stat > cutoffs-c 2> calcults-c.log
# split the primary assembly 
split_fa Q15_5kb_130423.bp.p_ctg.fa > HiFi130423.split
# Uses the split primary assembly and does a self assignment
minimap2 -xasm5 -DP HiFi130423.split HiFi130423.split | gzip -c -> HiFi130423-c.split.self.paf.gz
# purges the haplotigs, gives a duplication bed file and log of purged haplotigs - do we need a -e flag here?
purge_dups -2 -T cutoffs-c -c PB.base.cov HiFi130423-c.split.self.paf.gz > dups.bed 2> purge_dups.log
# taking the duplication bed file and the original assembly and getting the purged (primary and haplotig) sequences - results in purged.fa (use this) and hap.fa (alternatives, might be useful)
get_seqs dups.bed Q15_5kb_130423.bp.p_ctg.fa


# Trying ModEst for genome size - see how it stacks up to other tools for same purpose 
# Annabel figured out best to just use aligned and sorted .bam files to skip a few steps in backmap.pl
# used bwa for Illumina, and minimap2 for ONT for .sam files, then created .bam using samtools sort
# Will run this code when files are done (23 May)

#!/bin/bash -e
#SBATCH --job-name=modest_test
#SBATCH --output=MCN_%j.out
#SBATCH --error=MCN_%j.err
#SBATCH --mail-user=m.nehmens@massey.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --time=18:00:00
#SBATCH --mem=20G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --profile=task 
#SBATCH --account=brins03581

module purge
cd  /nesi/nobackup/ga03714/Melissa/Software/backmap/backmap
module load SAMtools BWA minimap2 BEDTools MultiQC R Perl

perl backmap.pl -b Purged_AllQ15_ONT.bam  -b Purged_Illumina.bam -o ModEst_test -t 16 -nq
