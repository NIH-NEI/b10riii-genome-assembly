# b10riii-genome-assembly

## 1. Short reads - kmer estimation
`sinteractive --cpus-per-task=60 --mem=1507g --gres=lscratch:300`\
`module load kmergenie`\
`cd /b10riii/RawData/illumina/Flowcell_HWWNGDRXY/illumina_reads`\
`ls *.fastq.gz > inputs`\
`kmergenie inputs -t60`

## 2. Short reads - Genome size estimation
`cd /b10riii/RawData/illumina/Flowcell_HWWNGDRXY/illumina_reads`\
`jellyfish count -t 60 -C -m 117 -s 30G -o jelly.jf -F 2 <(zcat F1_S1_R1_001.fastq.gz) <(zcat F1_S1_R2_001.fastq.gz)`\
`jellyfish stats -o jellystats.txt -v jelly.jf`\
`jellyfish histo -o jellyhisto.txt jelly.jf`

### 2.1. Generate plots and analysis in R using jellyhisto.txt
`sum(as.numeric(jel[5:9989,1]*jel[5:9989,2]))`\
#Peak at 18\
`sum(as.numeric(jel[5:9989,1]*jel[5:9989,2]))/18`\
#R commands in jellyanalysis.R in Tools folder

### 2.2. GenomeScope - analyze global genome properties
`Rscript genomescope.R jellyhisto.txt 21 150 /data/CaspiWGSData/b10riii/Results/genomescope`

## 3. Short reads - QA and QC 
`module load fastqc`\
`module load multiqc`\
`cd /b10riii/RawData/illumina/Flowcell_HWWNGDRXY/illumina_reads/`\
`fastqc F1_S1_R1_001.fastq.gz -t 20`\
`fastqc F1_S1_R2_001.fastq.gz -t 10`\
`multiqc *fastqc.zip`

### 3.1. Adapter identification, removal and trimming
`module load bbtools`\
`bbtools bbduk in1=F1_S1_R1_001.fastq.gz in2=F1_S1_R2_001.fastq.gz k=23 ref=adapters.fasta stats=bbdukstats.txt` \
`bbtools bbduk in1=F1_S1_R1_001.fastq.gz in2=F1_S1_R2_001.fastq.gz out1=F1_S1_R1_001.trimmed.fastq.gz out2=F1_S1_R2_001.trimmed.fastq.gz ktrim=r k=21 mink=11 hdist=1 ref=adapters.fasta stats=bbduk-trim-stats.txt 2> log.txt` \
#Check after trimming \
`bbtools bbduk in1=F1_S1_R1_001.trimmed.fastq.gz in2=F1_S1_R2_001.trimmed.fastq.gz k=23 ref=adapters.fasta stats=bbdukstatsafter.txt`

## 4. Long reads - Adapter identification and removal
`module load bamtools`\
`module load blast`\
`export PATH=$PATH:/data/CaspiWGSData/b10riii/Tools/HiFiAdapterFilt`\
`export PATH=$PATH:/data/CaspiWGSData/b10riii/Tools/HiFiAdapterFilt/DB`\
`bash pbadapterfilt.sh -t 30`

## 5. Hicanu assembly
`module load canu`\
`canu -p hicanudefault2 -d hicanudefault2 genomeSize=2.7g -pacbio-hifi *.fastq.gz`\
`canu -p hicanudefault3 -d hicanudefault3 genomeSize=2.7g -pacbio-hifi *.fastq.gz batThreads=24 batMemory=247g`

### 5.1. hicanu assembly quality
`module load quast`\
#Loading quast, version 5.0.2\
`quast.py -o quast/hicanu pacbio/hicanu/hicanudefault3/hicanudefault3.contigs.fasta pacbio/hicanu/hicanudefault4/hicanudefault4.contigs.fasta ../RawData/GCA_000001635.9_GRCm39_genomic.fna.gz`\
`quast.py -o quast/hicanu/vsmm39 -r ../RawData/GCA_000001635.9_GRCm39_genomic.fna.gz -g ../RawData/GCA_000001635.9_GRCm39_genomic.gff.gz pacbio/hicanu/hicanudefault3/hicanudefault3.contigs.fasta pacbio/hicanu/hicanudefault4/hicanudefault4.contigs.fasta`

## 6. pbipa assembly
`ipa local --nthreads 20 --njobs 1 --run-dir /data/CaspiWGSData/b10riii/Results/pacbio/pbipa/default -i /b10riii/RawData/pacbio/pacbio_ccs/F1_1.ccs.fastq -i /b10riii/RawData/pacbio/pacbio_ccs/F1_2.ccs.fastq -i /b10riii/RawData/pacbio/pacbio_ccs/F1_3.ccs.fastq -i /b10riii/RawData/pacbio/pacbio_ccs/F1_4.ccs.fastq`

### 6.1. pbipa assembly filtered
`ipa local --nthreads 20 --njobs 1 --run-dir /b10riii/Results/pacbio/pbipa/defaultfiltered -i /b10riii/Results/pacbio/hifiadapterfilt/F1_1.ccs.filt.fastq.gz -i /b10riii/Results/pacbio/hifiadapterfilt/F1_2.ccs.filt.fastq.gz -i /b10riii/Results/pacbio/hifiadapterfilt/F1_3.ccs.filt.fastq.gz -i /b10riii/Results/pacbio/hifiadapterfilt/F1_4.ccs.filt.fastq.gz`\
#to check snakemake command\
`ipa local -i ../Results/pacbio/pbipa/default/input.fofn --only-print --nthreads 20 --njobs 1`

### 6.2. pbipa assembly quality
`quast.py -o ../Results/quast/pbipa ../Results/pacbio/pbipa/default/assembly-results/final.p_ctg.fasta ../Results/pacbio/pbipa/default/assembly-results/final.a_ctg.fasta ../Results/pacbio/pbipa/default2/assembly-results/final.p_ctg.fasta ../Results/pacbio/pbipa/default2/assembly-results/final.a_ctg.fasta ../RawData/GCA_000001635.9_GRCm39_genomic.fna.gz --threads 60`

## 7. hifiasm assembly
`hifiasm -o ../Results/hifiasm/default.asm -t 600 ../RawData/pacbio/pacbio_ccs/F1_1.ccs.fastq ../RawData/pacbio/pacbio_ccs/F1_2.ccs.fastq ../RawData/pacbio/pacbio_ccs/F1_3.ccs.fastq ../RawData/pacbio/pacbio_ccs/F1_4.ccs.fastq`\
`sbatch --cpus-per-task=60 --mem=1507g --time=10-00:00:00 --partition=largemem /b10riii/Tools/hifiasm60filtered.sbatch`
#convert gfa to fa\
`awk '/^S/{print ">"$2;print $3}' default-60-filtered.asm.bp.p_ctg.gfa > default-60-filtered.asm.bp.p_ctg.fa`\
`awk '/^S/{print ">"$2;print $3}' default.asm.bp.p_ctg.gfa > default.asm.bp.p_ctg.fa`

### 7.1 quast with hicanu, pbipa, hifiasm
`quast.py -o /b10riii/Results/quast/hicanu-pbipa-hifiasm /b10riii/Results/pacbio/hicanu/hicanudefault3/hicanudefault3.contigs.fasta /b10riii/Results/pacbio/hicanu/hicanudefault4/hicanudefault4.contigs.fasta /b10riii/Results/pacbio/pbipa/default/assembly-results/final.p_ctg.fasta /b10riii/Results/pacbio/pbipa/default2/assembly-results/final.p_ctg.fasta /b10riii/Results/pacbio/hifiasm/default-60/default.asm.bp.p_ctg.fa /b10riii/Results/pacbio/hifiasm/default-60-filtered/default-60-filtered.asm.bp.p_ctg.fa /b10riii/Results/chromosome/test/GCF_000001635.27_GRCm39_genomic.fna.asmdefaultfiltered.bp.p_ctg.fa.split.reconciled.fa /b10riii/RawData/ncbi_dataset/data/GCF_000001635.27/GCF_000001635.27_GRCm39_genomic.fna --threads 120`

## 8. Polishing using long reads
`source /b10riii/Tools/conda/etc/profile.d/conda.sh`\
`conda activate pbmm2`\
`cd /b10riii/Results/polished/pacbiopolished`\
`echo "/b10riii/RawData/pacbio/pacbio_subreads/F1_1.subreads.bam" > subreadbams.fofn`\
`echo "/b10riii/RawData/pacbio/pacbio_subreads/F1_2.subreads.bam" >> subreadbams.fofn`\
`echo "/data/CaspiWGSData/b10riii/RawData/pacbio/pacbio_subreads/F1_3.subreads.bam" >> subreadbams.fofn`\
`echo "/data/CaspiWGSData/b10riii/RawData/pacbio/pacbio_subreads/F1_4.subreads.bam" >> subreadbams.fofn`\
`pbmm2 align /b10riii/Results/pacbio/assemblies/asmdefault.bp.p_ctg.fa subreadbams.fofn asmaligned.bam --sort -j 80 -J 60 -m 32G --preset SUBREAD`\
#-j, -J number of threads for alignment and sorting\

#parallelize alignment step\
#below generates the swarm file\
`for i in `cat sequence_heads.txt`; do echo "source /b10riii/Tools/conda/etc/profile.d/conda.sh ; conda activate pbmm2 ; export TMPDIR=/lscratch/\$SLURM_JOB_ID ; cd /b10riii/Results/polished/pacbiopolished/gcpp-parallel/alignments ; pbmm2 align /b10riii/Results/polished/pacbiopolished/gcpp-parallel/"$i".fa /b10riii/Results/polished/pacbiopolished/subreadbams.fofn asmdefaultfiltered-norm-aligned-"$i".bam --sort -j 56 -J 30 -m 8G --preset SUBREAD --log-level DEBUG --log-file asmdefaultfiltered-norm-aligned-"$i".log"; done`\
`swarm -f /b10riii/Tools/pbmm2-alignment-norm-parallel.swarm -g 247 -t 56 --gres=lscratch:800`\
`ls -lh | grep "samtools" | cut -d'.' -f4 | sort | uniq -c`

## 9. Hybrid assembly
`source /b10riii/Tools/conda/etc/profile.d/conda.sh`
`conda create -n python37 python=3.7.7 lxml=4.5.0 pandas=1.0.3 natsort=7.0.1 drmaa=0.7.9 numpy=1.18.4 numba=0.42.0 pandasql=0.7.3 scikit-learn=0.22.1 pyyaml=5.3.1 intervaltree=3.0.2 scipy=1.4.1 imbalanced-learn=0.6.2 matplotlib=3.1.3 lightgbm=2.3.0 xlrd=1.2.0 pytest=5.4.2 pytest-forked pytest-xdist pytest-cov coverage=5.4 hmmlearn xgboost=0.90 joblib=0.13.2 xlrd`\
`conda activate python37`\
`cd perl-bionano/`\
`pip install tools/pipeline/1.0/bionano_packages/pyBionano`\
`pip install tools/pipeline/1.0/bionano_packages/SVConfModels`\
`pip install tools/pipeline/1.0/bionano_packages/lohdetection`\
`cd /b10riii/Tools/`\
#edited align_molecules.pl to fix python2.7 to python \
#run bionano-solve-hybrid-batch.sh\
#/b10riii/Tools/perl-bionano/tools/pipeline/Solve3.7_03302022_283/HybridScaffold/03302022/scripts/align_molecules.pl\
#bionano uses 2.7 in code, but install requires 3.7.7\
`ln -s /b10riii/Tools/conda/envs/python37/bin/python /b10riii/Tools/conda/envs/python37/bin/python2.7`

`perl /b10riii/Tools/perl-bionano/tools/pipeline/1.0/HybridScaffold/1.0/hybridScaffold.pl -n /b10riii/Results/pacbio/assemblies/defaultfiltered.asm.bp.p_ctg.fa -b /b10riii/RawData/bionano/Assembly_data_delivery/output/contigs/exp_refineFinal1/EXP_REFINEFINAL1.cmap -c /b10riii/Tools/perl-bionano/tools/pipeline/1.0/HybridScaffold/1.0/hybridScaffold_DLE1_config.xml -r /b10riii/Tools/perl-bionano/tools/pipeline/1.0/RefAligner/1.0/RefAligner -o /b10riii/Results/hybridscaffold/hifiasm-filtered -f -g -B 2 -N 2`

## 10. Polishing using short reads
`cd /b10riii/Results/polished/pacbiopolished/gcpp-parallel/polished_seqs`\
`cut -f1 /b10riii/Results/pacbio/assemblies/asmdefaultfiltered.bp.p_ctg.fa.fai > sequence_heads.txt`\
`cd /b10riii/Results/polished/pacbiopolished/gcpp-parallel`\

#below generates the swarm file\
`for i in `cat sequence_heads.txt`; do echo "source /b10riii/Tools/conda/etc/profile.d/conda.sh ; module load bamtools ; TMPDIR=/b10riii/ ; cd /b10riii/Results/polished/pacbiopolished/gcpp-parallel/polished_seqs ; bamtools merge -in asmfilteredaligned-norm.bam -out "$i".bam -region "$i" ; conda activate pbindex ; pbindex "$i".bam"; done > bam-split-index.swarm`\

`swarm -f /b10riii/Tools/bam-split-index.swarm`\ 

`conda activate gcpp`\
`gcpp -r "$i".fasta -o polished_seqs/"$i".polished.fastq -o polished_seqs/"$i".polished.vcf -o polished_seqs/"$i".polished.gff -n 5 --annotateGFF --reportEffectiveCoverage "$i".bam"'
gcpp -r /b10riii/Results/polished/pacbiopolished/gcpp-parallel/ptg000111l.fa -o /b10riii/Results/polished/pacbiopolished/gcpp-parallel/polished_seqs/ptg000111l.polished.fastq --reportEffectiveCoverage /b10riii/Results/polished/pacbiopolished/gcpp-parallel/polished_seqs/ptg000111l.bam

Create index
cd /data/CaspiWGSData/b10riii/Results/polished/pacbiopolished/gcpp-parallel/polished_seqs
samtools faidx /data/CaspiWGSData/b10riii/Results/pacbio/assemblies/asmdefaultfiltered.bp.p_ctg.fa
cut -f1 /data/CaspiWGSData/b10riii/Results/pacbio/assemblies/asmdefaultfiltered.bp.p_ctg.fa.fai > sequence_heads.txt
# also faidx index split fasta files - not used
cd /data/CaspiWGSData/b10riii/Results/polished/pacbiopolished/gcpp-parallel
for i in `cat sequence_heads.txt`; do samtools faidx "$i".fa ; done

Test slice
bamtools merge -in asmfilteredaligned-norm.bam -out ptg000002l.bam -region ptg000002l

#below generates the swarm file
for i in `cat sequence_heads.txt`; do echo "source /data/CaspiWGSData/b10riii/Tools/conda/etc/profile.d/conda.sh ; module load bamtools ; TMPDIR=/data/CaspiWGSData/b10riii/ ; cd /data/CaspiWGSData/b10riii/Results/polished/pacbiopolished/gcpp-parallel/polished_seqs ; bamtools merge -in asmfilteredaligned-norm.bam -out "$i".bam -region "$i" ; conda activate pbindex ; pbindex "$i".bam"; done > bam-split-index.swarm

swarm -f /data/CaspiWGSData/b10riii/Tools/bam-split-index.swarm 

Test polishing (didnt work)
source /data/CaspiWGSData/b10riii/Tools/conda/etc/profile.d/conda.sh
conda activate gcpp
gcpp -r "$i".fasta -o polished_seqs/"$i".polished.fastq -o polished_seqs/"$i".polished.vcf -o polished_seqs/"$i".polished.gff -n 5 --annotateGFF --reportEffectiveCoverage "$i".bam"'
gcpp -r /data/CaspiWGSData/b10riii/Results/polished/pacbiopolished/gcpp-parallel/ptg000111l.fa -o /data/CaspiWGSData/b10riii/Results/polished/pacbiopolished/gcpp-parallel/polished_seqs/ptg000111l.polished.fastq --reportEffectiveCoverage /data/CaspiWGSData/b10riii/Results/polished/pacbiopolished/gcpp-parallel/polished_seqs/ptg000111l.bam

Test polishing using gcpp -w windows option
https://github.com/PacificBiosciences/pbbioconda/issues/285
cut -f1-2 /data/CaspiWGSData/b10riii/Results/pacbio/assemblies/asmdefaultfiltered.bp.p_ctg.fa.fai | awk ' { print $1 ":0-" $2 } ' > contigs.txt

#Test command

gcpp -w ptg000111l:0-16631 -r /data/CaspiWGSData/b10riii/Results/pacbio/assemblies/asmdefaultfiltered.bp.p_ctg.fa -o /data/CaspiWGSData/b10riii/Results/polished/pacbiopolished/gcpp-parallel/polished_seqs/ptg000111l.polished.fasta /data/CaspiWGSData/b10riii/Results/polished/pacbiopolished/gcpp-parallel/polished_seqs/ptg000111l.bam --log-level TRACE --log-file ptg000111l.log

#generate swarm file
bash /data/CaspiWGSData/b10riii/Tools/swarm-gcpp-windows-generator.sh > gcpp-swarm-with-windows.swarm

swarm -f /data/CaspiWGSData/b10riii/Tools/gcpp-swarm-with-windows.swarm -g 121 -t 20 --gres=lscratch:300

# merge polished fasta files
/data/CaspiWGSData/b10riii/Results/polished/pacbiopolished/gcpp-parallel/polished_seqs
ls *polished.fasta | wc -l
wc -l sequence_heads.txt 
cat *polished.fasta > arrow-polished-1.fasta
grep ">" arrow-polished-1.fasta
grep ">" arrow-polished-1.fasta | wc -l


# quast after and before polishing
module load quast
quast.py -o /data/CaspiWGSData/b10riii/Results/quast/polished-hifiasm-filtered /data/CaspiWGSData/b10riii/Results/pacbio/hifiasm/default-60-filtered/default-60-filtered.asm.bp.p_ctg.fa /data/CaspiWGSData/b10riii/Results/polished/pacbiopolished/gcpp-parallel/polished_seqs/arrow-polished-1/arrow-polished-1.fasta /data/CaspiWGSData/b10riii/RawData/ncbi_dataset/data/GCF_000001635.27/GCF_000001635.27_GRCm39_genomic.fna --threads 120

=========================
Arrow - polish 2

#pbmm2 align
source /data/CaspiWGSData/b10riii/Tools/conda/etc/profile.d/conda.sh ; TMPDIR="/data/CaspiWGSData/b10riii/" ; conda activate pbmm2 ; cd /data/CaspiWGSData/b10riii/Results/polished/pacbiopolished ; pbmm2 align /data/CaspiWGSData/b10riii/Results/polished/pacbiopolished/gcpp-parallel/polished_seqs/arrow-polished-1/arrow-polished-1.fasta /data/CaspiWGSData/b10riii/Results/polished/pacbiopolished/subreadbams.fofn asmfilteredaligned-norm-polished-1.bam --sort -j 160 -J 160 -m 8G --preset SUBREAD --log-level INFO --log-file asmdefaultfiltered-norm-polished-1.log

swarm -f /data/CaspiWGSData/b10riii/Tools/pbmm2-alignment-norm-parallel-polished-1.swarm -g 1507 -t 60 --gres=lscratch:800 --partition largemem



=======================
# Polishing with Pilon
# https://timkahlke.github.io/LongRead_tutorials/ECR_P.html
# trimmed fastq illumina location
/data/CaspiWGSData/b10riii/Results/illumina/bbduk

### parallelize pilon polishing - 1 day for 1 contigs
#https://raw.githubusercontent.com/harish0201/General_Scripts/master/Fasta_splitter.py
cd /data/CaspiWGSData/b10riii/Results/polished/illuminapolished/pilon-parallel/
python /data/CaspiWGSData/b10riii/Tools/Fasta_splitter.py arrow-polished-1.fasta > sequence_heads.txt

module load pilon
#pilon  1.24

# index genome
module load bwa
bwa 0.7.17
cd /data/CaspiWGSData/b10riii/Results/polished/illuminapolished
ln -s /data/CaspiWGSData/b10riii/Results/polished/pacbiopolished/gcpp-parallel/polished_seqs/arrow-polished-1/arrow-polished-1.fasta .
bwa index arrow-polished-1.fasta
#About 1 hr

# map reads
bwa mem -t 260 arrow-polished-1.fasta F1_S1_R1_001.trimmed.fastq.gz F1_S1_R2_001.trimmed.fastq.gz -o bwa_mapping_illumina_for_pilon_on_arrow_polished-1.sam
swarm -f /data/CaspiWGSData/b10riii/Tools/pilon-bwa.swarm -g 1507 -t 60 --gres=lscratch:600 --partition largemem
pilon alignment - 60 core machine - over night. 747 gb

# samtools view/sort/index
module load samtools ; cd /data/CaspiWGSData/b10riii/Results/polished/illuminapolished ; samtools view -@ 120 -Sb bwa_mapping_illumina_for_pilon_on_arrow_polished-1.sam > bwa_mapping_illumina_for_pilon_on_arrow_polished-1.bam ; samtools sort -@ 120 -o bwa_mapping_illumina_for_pilon_on_arrow_polished-1_sorted.bam bwa_mapping_illumina_for_pilon_on_arrow_polished-1.bam ; samtools index -@ 120 bwa_mapping_illumina_for_pilon_on_arrow_polished-1_sorted.bam
swarm -f /data/CaspiWGSData/b10riii/Tools/pilon-samtools.swarm -g 1507 -t 60 --gres=lscratch:600 --partition largemem
# about 2 hours

# pilon polish
module load pilon ; cd /data/CaspiWGSData/b10riii/Results/polished/illuminapolished ; java -Xmx${SLURM_MEM_PER_NODE}m -jar $PILON_JAR --genome arrow-polished-1.fasta --bam bwa_mapping_illumina_for_pilon_on_arrow_polished-1_sorted.bam --outdir testpilon
swarm -f /data/CaspiWGSData/b10riii/Tools/pilon.swarm -g 1507 -t 60 --gres=lscratch:600 --partition largemem


-------
Pilon altogether didnt work - 1 contig more than 1 day.
Parallel using individual contigs.

#below generates the swarm file to split pilon bam, sort and index per contig and run pilon
for i in `cat sequence_heads.txt`; do echo "module load bamtools ; module load samtools ; module load pilon ; TMPDIR=/data/CaspiWGSData/b10riii/ ; cd /data/CaspiWGSData/b10riii/Results/polished/illuminapolished/pilon-parallel ; bamtools merge -in bwa_mapping_illumina_for_pilon_on_arrow_polished-1_sorted.bam -out "$i".bam -region "$i" ; samtools sort -o "$i"_sorted.bam "$i".bam ; samtools index "$i"_sorted.bam ; java -Xmx${SLURM_MEM_PER_NODE}m -jar $PILON_JAR --genome "$i".fa --bam "$i"_sorted.bam --output "$i" --outdir testpilon" ; done > pilon_parallel.swarm

# Find/replace ptg000001l|arrow to ptg000001l\|arrow in the swarm file

# example swarm command
module load bamtools ; module load samtools ; module load pilon ; TMPDIR=/data/CaspiWGSData/b10riii/ ; cd /data/CaspiWGSData/b10riii/Results/polished/illuminapolished/pilon-parallel ; bamtools merge -in bwa_mapping_illumina_for_pilon_on_arrow_polished-1_sorted.bam -out ptg000001l\|arrow.bam -region ptg000001l\|arrow ; samtools sort -o ptg000001l\|arrow_sorted.bam ptg000001l\|arrow.bam ; samtools index ptg000001l\|arrow_sorted.bam ; java -Xmx1543168m -jar /usr/local/apps/pilon/1.24/pilon-1.24.jar --genome ptg000001l\|arrow.fa --bam ptg000001l\|arrow_sorted.bam --output ptg000001l\|arrow --outdir testpilon

cd /data/CaspiWGSData/b10riii/Results/polished/illuminapolished/pilon-parallel

# run normal swarm
swarm -f /data/CaspiWGSData/b10riii/Tools/pilon_parallel.swarm -g 247 -t 28 --gres=lscratch:400
#58l takes 3 days
# run 58l alone in a different swarm job with largemem
swarm -f /data/CaspiWGSData/b10riii/Tools/pilon_parallel-58.swarm -g 3000 -t 60 --gres=lscratch:700 --partition largemem



# combine fastas
/data/CaspiWGSData/b10riii/Results/polished/illuminapolished/pilon-parallel/testpilon

cat *arrow.fasta > ../../pilon-parallel/arrow-1-pilon-1-polished-1.fasta
grep ">" arrow-1-pilon-1-polished-1.fasta
grep ">" arrow-1-pilon-1-polished-1.fasta | wc -l

/data/CaspiWGSData/b10riii/Results/polished/illuminapolished/pilon-parallel/arrow-1-pilon-1-polished-1.fasta

# quast after and before pilon polishing
module load quast
quast.py -o /data/CaspiWGSData/b10riii/Results/quast/pilon-polished-hifiasm-filtered /data/CaspiWGSData/b10riii/Results/pacbio/hifiasm/default-60-filtered/default-60-filtered.asm.bp.p_ctg.fa /data/CaspiWGSData/b10riii/Results/polished/pacbiopolished/gcpp-parallel/polished_seqs/arrow-polished-1/arrow-polished-1.fasta /data/CaspiWGSData/b10riii/Results/polished/illuminapolished/pilon-parallel/arrow-1-pilon-1-polished-1.fasta /data/CaspiWGSData/b10riii/RawData/ncbi_dataset/data/GCF_000001635.27/GCF_000001635.27_GRCm39_genomic.fna --threads 120
=======================

Bionano



#bionano solve install local
#tried conda, docker - errors emailed helix
# installed singularity in aws ubuntu linux, created perl container, installed conda and python, used container in biowulf


singularity build --remote --sandbox singularity-bionano.img docker://ubuntu:20.04
singularity shell --writable singularity-bionano.img/

module load singularity
singularity shell -B /usr/lib/locale/:/usr/lib/locale/ --bind /data/CaspiWGSData/ singularity-bionano-perl-python.sif 
source /data/CaspiWGSData/b10riii/Tools/conda/etc/profile.d/conda.sh
conda create -n python37 python=3.7.7 lxml=4.5.0 pandas=1.0.3 natsort=7.0.1 drmaa=0.7.9 numpy=1.18.4 numba=0.42.0 pandasql=0.7.3 scikit-learn=0.22.1 pyyaml=5.3.1 intervaltree=3.0.2 scipy=1.4.1 imbalanced-learn=0.6.2 matplotlib=3.1.3 lightgbm=2.3.0 xlrd=1.2.0 pytest=5.4.2 pytest-forked pytest-xdist pytest-cov coverage=5.4 hmmlearn xgboost=0.90 joblib=0.13.2 xlrd
conda activate python37
cd perl-bionano/
pip install tools/pipeline/1.0/bionano_packages/pyBionano
pip install tools/pipeline/1.0/bionano_packages/SVConfModels
pip install tools/pipeline/1.0/bionano_packages/lohdetection
cd /data/CaspiWGSData/b10riii/Tools/

#edited align_molecules.pl to fix python2.7 to python 
#run bionano-solve-hybrid-batch.sh
#/home/nagarajanv/Desktop/CaspiWGSData/b10riii/Tools/perl-bionano/tools/pipeline/Solve3.7_03302022_283/HybridScaffold/03302022/scripts/align_molecules.pl

singularity shell -B /usr/lib/locale/:/usr/lib/locale/ --bind /data/CaspiWGSData/ singularity-bionano-perl-python-R.sif
source /data/CaspiWGSData/b10riii/Tools/conda/etc/profile.d/conda.sh
conda activate python37
cd perl-bionano
#bionano uses 2.7 in code, but install requires 3.7.7
ln -s /data/CaspiWGSData/b10riii/Tools/conda/envs/python37/bin/python /data/CaspiWGSData/b10riii/Tools/conda/envs/python37/bin/python2.7


/data/CaspiWGSData/b10riii/Tools/perl-bionano/tools/pipeline/1.0/HybridScaffold/1.0/hybridScaffold.pl
-n /data/CaspiWGSData/b10riii/Results/pacbio/assemblies/defaultfiltered.asm.bp.p_ctg.fa
/data/CaspiWGSData/b10riii/RawData/bionano/Assembly_data_delivery/output/contigs/exp_refineFinal1/EXP_REFINEFINAL1.cmap
/data/CaspiWGSData/b10riii/Tools/perl-bionano/tools/pipeline/1.0/HybridScaffold/1.0/hybridScaffold_DLE1_config.xml
/data/CaspiWGSData/b10riii/Tools/perl-bionano/tools/pipeline/1.0/RefAligner/1.0/RefAligner
/data/CaspiWGSData/b10riii/Results/hybridscaffold/hifiasm-filtered

-x
-y
-m /data/CaspiWGSData/b10riii/RawData/bionano/Molecule_files/Merged_2Runs_splenocytes_-_Molecule_Merge_RawMolecules.bnx
-p /data/CaspiWGSData/b10riii/Tools/perl-bionano/tools/pipeline/1.0/Pipeline/1.0
-q /data/CaspiWGSData/b10riii/Tools/perl-bionano/tools/pipeline/1.0/Pipeline/1.0/optArguments.xml
-e /data/CaspiWGSData/b10riii/RawData/bionano/Assembly_data_delivery/output/contigs/auto_noise/autoNoise1.errbin





perl /data/CaspiWGSData/b10riii/Tools/perl-bionano/tools/pipeline/1.0/HybridScaffold/1.0/hybridScaffold.pl -n /data/CaspiWGSData/b10riii/Results/pacbio/assemblies/defaultfiltered.asm.bp.p_ctg.fa -b /data/CaspiWGSData/b10riii/RawData/bionano/Assembly_data_delivery/output/contigs/exp_refineFinal1/EXP_REFINEFINAL1.cmap -c /data/CaspiWGSData/b10riii/Tools/perl-bionano/tools/pipeline/1.0/HybridScaffold/1.0/hybridScaffold_DLE1_config.xml -r /data/CaspiWGSData/b10riii/Tools/perl-bionano/tools/pipeline/1.0/RefAligner/1.0/RefAligner -o /data/CaspiWGSData/b10riii/Results/hybridscaffold/hifiasm-filtered -f -g -B 2 -N 2

Mon Jun 13 22:48:14 EDT 2022, Mon Jun 13 23:28:40 EDT 2022
# memory usage 1092, cores 64 100% in merging step
# about 1 hr

# load singularity

cd /data/CaspiWGSData/b10riii/Tools/
module load singularity
singularity shell -B /usr/lib/locale/:/usr/lib/locale/ --bind /data/CaspiWGSData/ singularity-bionano-perl-python-R-time-zip.sif
source /data/CaspiWGSData/b10riii/Tools/conda/etc/profile.d/conda.sh
conda activate python37
cd perl-bionano


/data/CaspiWGSData/b10riii/Tools/perl-bionano/tools/pipeline/1.0/HybridScaffold/1.0/hybridScaffold.pl
-n /data/CaspiWGSData/b10riii/Results/polished/illuminapolished/pilon-parallel/arrow-1-pilon-1-polished-1.fasta
/data/CaspiWGSData/b10riii/RawData/bionano/Assembly_data_delivery/output/contigs/exp_refineFinal1/EXP_REFINEFINAL1.cmap
/data/CaspiWGSData/b10riii/Tools/perl-bionano/tools/pipeline/1.0/HybridScaffold/1.0/hybridScaffold_DLE1_config.xml
/data/CaspiWGSData/b10riii/Tools/perl-bionano/tools/pipeline/1.0/RefAligner/1.0/RefAligner
/data/CaspiWGSData/b10riii/Results/hybridscaffold/arrow-pilon-polished

-x
-y
-m /data/CaspiWGSData/b10riii/RawData/bionano/Molecule_files/Merged_2Runs_splenocytes_-_Molecule_Merge_RawMolecules.bnx
-p /data/CaspiWGSData/b10riii/Tools/perl-bionano/tools/pipeline/1.0/Pipeline/1.0
-q /data/CaspiWGSData/b10riii/Tools/perl-bionano/tools/pipeline/1.0/Pipeline/1.0/optArguments.xml
-e /data/CaspiWGSData/b10riii/RawData/bionano/Assembly_data_delivery/output/contigs/auto_noise/autoNoise1.errbin

arrow and pilon polished scaffold
perl /data/CaspiWGSData/b10riii/Tools/perl-bionano/tools/pipeline/1.0/HybridScaffold/1.0/hybridScaffold.pl -n /data/CaspiWGSData/b10riii/Results/polished/illuminapolished/pilon-parallel/arrow-1-pilon-1-polished-1.fasta -b /data/CaspiWGSData/b10riii/RawData/bionano/Assembly_data_delivery/output/contigs/exp_refineFinal1/EXP_REFINEFINAL1.cmap -c /data/CaspiWGSData/b10riii/Tools/perl-bionano/tools/pipeline/1.0/HybridScaffold/1.0/hybridScaffold_DLE1_config.xml -r /data/CaspiWGSData/b10riii/Tools/perl-bionano/tools/pipeline/1.0/RefAligner/1.0/RefAligner -o /data/CaspiWGSData/b10riii/Results/hybridscaffold/arrow-pilon-polished -f -g -B 2 -N 2

arrow-ony polished scaffold
perl /data/CaspiWGSData/b10riii/Tools/perl-bionano/tools/pipeline/1.0/HybridScaffold/1.0/hybridScaffold.pl -n /data/CaspiWGSData/b10riii/Results/polished/pacbiopolished/gcpp-parallel/polished_seqs/arrow-polished-1/arrow-polished-1.fasta -b /data/CaspiWGSData/b10riii/RawData/bionano/Assembly_data_delivery/output/contigs/exp_refineFinal1/EXP_REFINEFINAL1.cmap -c /data/CaspiWGSData/b10riii/Tools/perl-bionano/tools/pipeline/1.0/HybridScaffold/1.0/hybridScaffold_DLE1_config.xml -r /data/CaspiWGSData/b10riii/Tools/perl-bionano/tools/pipeline/1.0/RefAligner/1.0/RefAligner -o /data/CaspiWGSData/b10riii/Results/hybridscaffold/arrow-only-polished -f -g -B 2 -N 2


=========================

QUAST with unpolished and polished scaffolds
unpolished scaffold
/data/CaspiWGSData/b10riii/Results/hybridscaffold/hifiasm-60-filtered/hybrid_scaffolds/EXP_REFINEFINAL1_bppAdjust_cmap_default-60-filtered_asm_bp_p_ctg_fa_NGScontigs_HYBRID_SCAFFOLD.fasta
arrow only polished scaffold
/data/CaspiWGSData/b10riii/Results/hybridscaffold/arrow-only-polished/hybrid_scaffolds/EXP_REFINEFINAL1_bppAdjust_cmap_arrow-polished-1_fasta_NGScontigs_HYBRID_SCAFFOLD.fasta
polished scaffold
/data/CaspiWGSData/b10riii/Results/hybridscaffold/arrow-pilon-polished/hybrid_scaffolds/EXP_REFINEFINAL1_bppAdjust_cmap_arrow-1-pilon-1-polished-1_fasta_NGScontigs_HYBRID_SCAFFOLD.fasta


module load quast
quast.py -o /data/CaspiWGSData/b10riii/Results/quast/polished-scaffolded /data/CaspiWGSData/b10riii/Results/hybridscaffold/hifiasm-60-filtered/hybrid_scaffolds/EXP_REFINEFINAL1_bppAdjust_cmap_default-60-filtered_asm_bp_p_ctg_fa_NGScontigs_HYBRID_SCAFFOLD.fasta /data/CaspiWGSData/b10riii/Results/hybridscaffold/arrow-only-polished/hybrid_scaffolds/EXP_REFINEFINAL1_bppAdjust_cmap_arrow-polished-1_fasta_NGScontigs_HYBRID_SCAFFOLD.fasta /data/CaspiWGSData/b10riii/Results/hybridscaffold/arrow-pilon-polished/hybrid_scaffolds/EXP_REFINEFINAL1_bppAdjust_cmap_arrow-1-pilon-1-polished-1_fasta_NGScontigs_HYBRID_SCAFFOLD.fasta /data/CaspiWGSData/b10riii/RawData/ncbi_dataset/data/GCF_000001635.27/GCF_000001635.27_GRCm39_genomic.fna --threads 120

=================================

Reference based chromosome scaffolder
/data/CaspiWGSData/b10riii/Tools/MaSuRCA-4.1.0/bin

./../../../Tools/MaSuRCA-4.1.0/bin/chromosome_scaffolder.sh -r /data/CaspiWGSData/b10riii/RawData/ncbi_dataset/data/GCF_000001635.27/GCF_000001635.27_GRCm39_genomic.fna -q /data/CaspiWGSData/b10riii/Results/pacbio/assemblies/asmdefaultfiltered.bp.p_ctg.fa -t 120 -nb

# Ragtag based chromosome scaffolder
source /data/CaspiWGSData/b10riii/Tools/conda/etc/profile.d/conda.sh ; TMPDIR="/data/CaspiWGSData/b10riii/" ; conda activate ragtag ; /data/CaspiWGSData/b10riii/Results/chromosome/testragtag/ ; ragtag.py scaffold /data/CaspiWGSData/b10riii/RawData/ncbi_dataset/data/GCF_000001635.27/GCF_000001635.27_GRCm39_genomic.fna /data/CaspiWGSData/b10riii/Results/pacbio/assemblies/asmdefaultfiltered.bp.p_ctg.fa -t 120 -o chromosomeragtagtest

swarm -f /data/CaspiWGSData/b10riii/Tools/chromosome-ragtag.swarm -g 247 -t 28 --gres=lscratch:300
# an hour

# compare masurca and ragtag quast
module load quast
quast.py -o /data/CaspiWGSData/b10riii/Results/quast/chromosome/test-masurca-ragtag /data/CaspiWGSData/b10riii/Results/chromosome/test/GCF_000001635.27_GRCm39_genomic.fna.asmdefaultfiltered.bp.p_ctg.fa.split.reconciled.fa /data/CaspiWGSData/b10riii/Results/chromosome/testragtag/chromosomeragtagtest/ragtag.scaffold.fasta /data/CaspiWGSData/b10riii/RawData/ncbi_dataset/data/GCF_000001635.27/GCF_000001635.27_GRCm39_genomic.fna --threads 120

masurca is the best based on n50.

/data/CaspiWGSData/b10riii/Results/hybridscaffold/hifiasm-60-filtered/hybrid_scaffolds/EXP_REFINEFINAL1_bppAdjust_cmap_default-60-filtered_asm_bp_p_ctg_fa_NGScontigs_HYBRID_SCAFFOLD.fasta
/data/CaspiWGSData/b10riii/Results/hybridscaffold/arrow-only-polished/hybrid_scaffolds/EXP_REFINEFINAL1_bppAdjust_cmap_arrow-polished-1_fasta_NGScontigs_HYBRID_SCAFFOLD.fasta 
/data/CaspiWGSData/b10riii/Results/hybridscaffold/arrow-pilon-polished/hybrid_scaffolds/EXP_REFINEFINAL1_bppAdjust_cmap_arrow-1-pilon-1-polished-1_fasta_NGScontigs_HYBRID_SCAFFOLD.fasta

/data/CaspiWGSData/b10riii/Results/chromosome/final

ln -s /data/CaspiWGSData/b10riii/Results/hybridscaffold/hifiasm-60-filtered/hybrid_scaffolds/EXP_REFINEFINAL1_bppAdjust_cmap_default-60-filtered_asm_bp_p_ctg_fa_NGScontigs_HYBRID_SCAFFOLD.fasta unpolished-scaffold.fasta
ln -s /data/CaspiWGSData/b10riii/Results/hybridscaffold/arrow-only-polished/hybrid_scaffolds/EXP_REFINEFINAL1_bppAdjust_cmap_arrow-polished-1_fasta_NGScontigs_HYBRID_SCAFFOLD.fasta arrow-polished-scaffold.fasta
ln -s /data/CaspiWGSData/b10riii/Results/hybridscaffold/arrow-pilon-polished/hybrid_scaffolds/EXP_REFINEFINAL1_bppAdjust_cmap_arrow-1-pilon-1-polished-1_fasta_NGScontigs_HYBRID_SCAFFOLD.fasta arrow-pilon-polished-scaffold.fasta

unpolished-scaffold.fasta
arrow-polished-scaffold.fasta
arrow-pilon-polished-scaffold.fasta

cd /data/CaspiWGSData/b10riii/Results/chromosome/final/unpolished/masurca/ ; bash /data/CaspiWGSData/b10riii/Tools/MaSuRCA-4.1.0/bin/chromosome_scaffolder.sh -r /data/CaspiWGSData/b10riii/RawData/ncbi_dataset/data/GCF_000001635.27/GCF_000001635.27_GRCm39_genomic.fna -q /data/CaspiWGSData/b10riii/Results/chromosome/final/unpolished-scaffold.fasta -t 120 -nb

cd /data/CaspiWGSData/b10riii/Results/chromosome/final/arrow/masurca/ ; bash /data/CaspiWGSData/b10riii/Tools/MaSuRCA-4.1.0/bin/chromosome_scaffolder.sh -r /data/CaspiWGSData/b10riii/RawData/ncbi_dataset/data/GCF_000001635.27/GCF_000001635.27_GRCm39_genomic.fna -q /data/CaspiWGSData/b10riii/Results/chromosome/final/arrow-polished-scaffold.fasta -t 120 -nb

cd /data/CaspiWGSData/b10riii/Results/chromosome/final/arrow-pilon/masurca/ ; bash /data/CaspiWGSData/b10riii/Tools/MaSuRCA-4.1.0/bin/chromosome_scaffolder.sh -r /data/CaspiWGSData/b10riii/RawData/ncbi_dataset/data/GCF_000001635.27/GCF_000001635.27_GRCm39_genomic.fna -q /data/CaspiWGSData/b10riii/Results/chromosome/final/arrow-pilon-polished-scaffold.fasta -t 120 -nb

source /data/CaspiWGSData/b10riii/Tools/conda/etc/profile.d/conda.sh ; TMPDIR="/data/CaspiWGSData/b10riii/" ; conda activate ragtag ; cd /data/CaspiWGSData/b10riii/Results/chromosome/final/unpolished/ragtag/ ; ragtag.py scaffold /data/CaspiWGSData/b10riii/RawData/ncbi_dataset/data/GCF_000001635.27/GCF_000001635.27_GRCm39_genomic.fna /data/CaspiWGSData/b10riii/Results/chromosome/final/unpolished-scaffold.fasta -t 120 -o unpolished

source /data/CaspiWGSData/b10riii/Tools/conda/etc/profile.d/conda.sh ; TMPDIR="/data/CaspiWGSData/b10riii/" ; conda activate ragtag ; cd /data/CaspiWGSData/b10riii/Results/chromosome/final/arrow/ragtag/ ; ragtag.py scaffold /data/CaspiWGSData/b10riii/RawData/ncbi_dataset/data/GCF_000001635.27/GCF_000001635.27_GRCm39_genomic.fna /data/CaspiWGSData/b10riii/Results/chromosome/final/arrow-polished-scaffold.fasta -t 120 -o arrow-polished

source /data/CaspiWGSData/b10riii/Tools/conda/etc/profile.d/conda.sh ; TMPDIR="/data/CaspiWGSData/b10riii/" ; conda activate ragtag ; cd /data/CaspiWGSData/b10riii/Results/chromosome/final/arrow-pilon/ragtag/ ; ragtag.py scaffold /data/CaspiWGSData/b10riii/RawData/ncbi_dataset/data/GCF_000001635.27/GCF_000001635.27_GRCm39_genomic.fna /data/CaspiWGSData/b10riii/Results/chromosome/final/arrow-pilon-polished-scaffold.fasta -t 120 -o arrow-pilon-polished

swarm -f /data/CaspiWGSData/b10riii/Tools/chromosome-ragtag.swarm -g 247 -t 28 --gres=lscratch:300

/data/CaspiWGSData/b10riii/Results/chromosome/final/unpolished/masurca/GCF_000001635.27_GRCm39_genomic.fna.unpolished-scaffold.fasta.split.reconciled.fa
/data/CaspiWGSData/b10riii/Results/chromosome/final/unpolished/ragtag/ragtag.scaffold.fasta
/data/CaspiWGSData/b10riii/Results/chromosome/final/arrow/masurca/GCF_000001635.27_GRCm39_genomic.fna.arrow-polished-scaffold.fasta.split.reconciled.fa
/data/CaspiWGSData/b10riii/Results/chromosome/final/arrow/ragtag/ragtag.scaffold.fasta
/data/CaspiWGSData/b10riii/Results/chromosome/final/arrow-pilon/masurca/GCF_000001635.27_GRCm39_genomic.fna.arrow-pilon-polished-scaffold.fasta.split.reconciled.fa
/data/CaspiWGSData/b10riii/Results/chromosome/final/arrow-pilon/ragtag/ragtag.scaffold.fasta

chromosome-masurca-ragtag.swarm

#quast for chromosomes unpolished, arrow polished, arrow-pilon polished - both masurca and ragtag
quast.py -o /data/CaspiWGSData/b10riii/Results/quast/chromosome/final /data/CaspiWGSData/b10riii/Results/chromosome/final/unpolished/masurca/GCF_000001635.27_GRCm39_genomic.fna.unpolished-scaffold.fasta.split.reconciled.fa /data/CaspiWGSData/b10riii/Results/chromosome/final/unpolished/ragtag/ragtag.scaffold.fasta /data/CaspiWGSData/b10riii/Results/chromosome/final/arrow/masurca/GCF_000001635.27_GRCm39_genomic.fna.arrow-polished-scaffold.fasta.split.reconciled.fa /data/CaspiWGSData/b10riii/Results/chromosome/final/arrow/ragtag/ragtag.scaffold.fasta /data/CaspiWGSData/b10riii/Results/chromosome/final/arrow-pilon/masurca/GCF_000001635.27_GRCm39_genomic.fna.arrow-pilon-polished-scaffold.fasta.split.reconciled.fa /data/CaspiWGSData/b10riii/Results/chromosome/final/arrow-pilon/ragtag/ragtag.scaffold.fasta /data/CaspiWGSData/b10riii/RawData/ncbi_dataset/data/GCF_000001635.27/GCF_000001635.27_GRCm39_genomic.fna --threads 200

======================
#Busco
#Use custom conda install busco5, not the biowulf

busco -c 120 -m genome -l glires_odb10 -i /data/CaspiWGSData/b10riii/Results/chromosome/test/GCF_000001635.27_GRCm39_genomic.fna.asmdefaultfiltered.bp.p_ctg.fa.split.reconciled.fa -o test

source /data/CaspiWGSData/b10riii/Tools/conda/etc/profile.d/conda.sh ; TMPDIR="/data/CaspiWGSData/b10riii/" ; conda activate busco5 ; cd /data/CaspiWGSData/b10riii/Results/busco ; busco -c 120 -m genome -l glires_odb10 -i /data/CaspiWGSData/b10riii/Results/chromosome/test/GCF_000001635.27_GRCm39_genomic.fna.asmdefaultfiltered.bp.p_ctg.fa.split.reconciled.fa -o test

source /data/CaspiWGSData/b10riii/Tools/conda/etc/profile.d/conda.sh ; TMPDIR="/data/CaspiWGSData/b10riii/" ; conda activate busco5 ; cd /data/CaspiWGSData/b10riii/Results/busco ; busco -c 120 -m genome --auto-lineage-euk -i /data/CaspiWGSData/b10riii/Results/chromosome/test/GCF_000001635.27_GRCm39_genomic.fna.asmdefaultfiltered.bp.p_ctg.fa.split.reconciled.fa -o testautoeuk

swarm -f /data/CaspiWGSData/b10riii/Tools/busco-hifiasm-60-filtered.sbatch -g 247 -t 28 --gres=lscratch:300
# less than 12 hours 96% complete test chromosome

source /data/CaspiWGSData/b10riii/Tools/conda/etc/profile.d/conda.sh ; TMPDIR="/data/CaspiWGSData/b10riii/" ; conda activate busco5 ; cd /data/CaspiWGSData/b10riii/Results/busco ; busco -c 120 -m genome --auto-lineage-euk -i /data/CaspiWGSData/b10riii/RawData/ncbi_dataset/data/GCF_000001635.27/GCF_000001635.27_GRCm39_genomic.fna -o testautoeuk-mm39

# RUN EUK busco and gliers both

busco-final.swarm

ln -s /data/CaspiWGSData/b10riii/Results/busco/arrow-masurca-autoeuk-mm39/auto_lineage/run_eukaryota_odb10/short_summary.txt short_summary_arrow-ma-auto.txt
ln -s /data/CaspiWGSData/b10riii/Results/busco/arrow-masurca-gliers-mm39/run_glires_odb10/short_summary.txt short_summary_arrow-ma-gliers.txt
ln -s /data/CaspiWGSData/b10riii/Results/busco/arrow-pilon-masurca-autoeuk-mm39/auto_lineage/run_eukaryota_odb10/short_summary.txt short_summary_arrow-pilon-ma-auto.txt
ln -s /data/CaspiWGSData/b10riii/Results/busco/arrow-pilon-masurca-gliers-mm39/run_glires_odb10/short_summary.txt short_summary_arrow-pilon-ma-gliers.txt
ln -s /data/CaspiWGSData/b10riii/Results/busco/arrow-pilon-ragtag-autoeuk-mm39/auto_lineage/run_eukaryota_odb10/short_summary.txt short_summary_arrow-pilon-rag-auto.txt
ln -s /data/CaspiWGSData/b10riii/Results/busco/arrow-pilon-ragtag-gliers-mm39/run_glires_odb10/short_summary.txt short_summary_arrow-pilon-rag-gliers.txt
ln -s /data/CaspiWGSData/b10riii/Results/busco/arrow-ragtag-autoeuk-mm39/auto_lineage/run_eukaryota_odb10/short_summary.txt short_summary_arrow-rag-auto.txt
ln -s /data/CaspiWGSData/b10riii/Results/busco/arrow-ragtag-gliers-mm39/run_glires_odb10/short_summary.txt short_summary_arrow-rag-gliers.txt
ln -s /data/CaspiWGSData/b10riii/Results/busco/ref-autoeuk-mm39/auto_lineage/run_eukaryota_odb10/short_summary.txt short_summary_ref-auto.txt
ln -s /data/CaspiWGSData/b10riii/Results/busco/ref-gliers-mm39/run_glires_odb10/short_summary.txt short_summary_ref-gliers.txt
ln -s /data/CaspiWGSData/b10riii/Results/busco/unpolished-masurca-autoeuk-mm39/auto_lineage/run_eukaryota_odb10/short_summary.txt short_summary_unpolished-ma-auto.txt
ln -s /data/CaspiWGSData/b10riii/Results/busco/unpolished-masurca-gliers-mm39/run_glires_odb10/short_summary.txt short_summary_unpolished-ma-gliers.txt
ln -s /data/CaspiWGSData/b10riii/Results/busco/unpolished-ragtag-autoeuk-mm39/auto_lineage/run_eukaryota_odb10/short_summary.txt short_summary_unpolished-rag-auto.txt
ln -s /data/CaspiWGSData/b10riii/Results/busco/unpolished-ragtag-gliers-mm39/run_glires_odb10/short_summary.txt short_summary_unpolished-rag-gliers.txt

mutliqc results for above;
/data/CaspiWGSData/b10riii/Results/busco/allsummaries/multiqc_report.html

Final, best N50, best busco - gliers - arrow only
/data/CaspiWGSData/b10riii/Results/chromosome/final/arrow/ragtag/ragtag.scaffold.fasta

-------------
quast on arrow-ragtag and reference
quast.py -o /data/CaspiWGSData/b10riii/Results/quast/chromosome/final-one /data/CaspiWGSData/b10riii/Results/chromosome/final/arrow/ragtag/ragtag.scaffold.fasta /data/CaspiWGSData/b10riii/RawData/ncbi_dataset/data/GCF_000001635.27/GCF_000001635.27_GRCm39_genomic.fna --threads 200

busco on arrow-ragtag only
source /data/CaspiWGSData/b10riii/Tools/conda/etc/profile.d/conda.sh ; TMPDIR="/data/CaspiWGSData/b10riii/" ; conda activate busco5 ; cd /data/CaspiWGSData/b10riii/Results/busco ; busco -c 20 -m genome -l glires_odb10 -i /data/CaspiWGSData/b10riii/Results/chromosome/final/arrow/ragtag/ragtag.scaffold.fasta -o arrow-ragtag-gliers-only-mm39

--------------
Rename chromosome numbers

/data/CaspiWGSData/b10riii/Results/b10riii.fa

sed -i 's/NW_023337853.1_RagTag/bNW_023337853.1 Mus musculus strain B10RIII unplaced unlocalized genomic scaffold/g' /data/CaspiWGSData/b10riii/Results/b10riii.fa

sed -i 's/Super-Scaffold_100310/Super-Scaffold_100310 Mus musculus strain B10RIII unplaced unlocalized genomic scaffold/g' /data/CaspiWGSData/b10riii/Results/b10riii.fa
sed -i 's/Super-Scaffold_100970/Super-Scaffold_100970 Mus musculus strain B10RIII unplaced unlocalized genomic scaffold/g' /data/CaspiWGSData/b10riii/Results/b10riii.fa
sed -i 's/Super-Scaffold_17/Super-Scaffold_17 Mus musculus strain B10RIII unplaced unlocalized genomic scaffold/g' /data/CaspiWGSData/b10riii/Results/b10riii.fa
sed -i 's/Super-Scaffold_187/Super-Scaffold_187 Mus musculus strain B10RIII unplaced unlocalized genomic scaffold/g' /data/CaspiWGSData/b10riii/Results/b10riii.fa
sed -i 's/Super-Scaffold_191/Super-Scaffold_191 Mus musculus strain B10RIII unplaced unlocalized genomic scaffold/g' /data/CaspiWGSData/b10riii/Results/b10riii.fa
sed -i 's/Super-Scaffold_202/Super-Scaffold_202 Mus musculus strain B10RIII unplaced unlocalized genomic scaffold/g' /data/CaspiWGSData/b10riii/Results/b10riii.fa
sed -i 's/Super-Scaffold_295/Super-Scaffold_295 Mus musculus strain B10RIII unplaced unlocalized genomic scaffold/g' /data/CaspiWGSData/b10riii/Results/b10riii.fa
sed -i 's/Super-Scaffold_329/Super-Scaffold_329 Mus musculus strain B10RIII unplaced unlocalized genomic scaffold/g' /data/CaspiWGSData/b10riii/Results/b10riii.fa
sed -i 's/Super-Scaffold_40/Super-Scaffold_40 Mus musculus strain B10RIII unplaced unlocalized genomic scaffold/g' /data/CaspiWGSData/b10riii/Results/b10riii.fa
sed -i 's/Super-Scaffold_45/Super-Scaffold_45 Mus musculus strain B10RIII unplaced unlocalized genomic scaffold/g' /data/CaspiWGSData/b10riii/Results/b10riii.fa
sed -i 's/Super-Scaffold_587/Super-Scaffold_587 Mus musculus strain B10RIII unplaced unlocalized genomic scaffold/g' /data/CaspiWGSData/b10riii/Results/b10riii.fa
sed -i 's/Super-Scaffold_61/Super-Scaffold_61 Mus musculus strain B10RIII unplaced unlocalized genomic scaffold/g' /data/CaspiWGSData/b10riii/Results/b10riii.fa
sed -i 's/Super-Scaffold_71/Super-Scaffold_71 Mus musculus strain B10RIII unplaced unlocalized genomic scaffold/g' /data/CaspiWGSData/b10riii/Results/b10riii.fa
sed -i 's/Super-Scaffold_74/Super-Scaffold_74 Mus musculus strain B10RIII unplaced unlocalized genomic scaffold/g' /data/CaspiWGSData/b10riii/Results/b10riii.fa
sed -i 's/Super-Scaffold_77/Super-Scaffold_77 Mus musculus strain B10RIII unplaced unlocalized genomic scaffold/g' /data/CaspiWGSData/b10riii/Results/b10riii.fa
sed -i 's/Super-Scaffold_78/Super-Scaffold_78 Mus musculus strain B10RIII unplaced unlocalized genomic scaffold/g' /data/CaspiWGSData/b10riii/Results/b10riii.fa
sed -i 's/Super-Scaffold_79/Super-Scaffold_79 Mus musculus strain B10RIII unplaced unlocalized genomic scaffold/g' /data/CaspiWGSData/b10riii/Results/b10riii.fa


>bNC_000067.7 Mus musculus strain B10.RIII chromosome 1
>bNC_000068.8 Mus musculus strain B10.RIII chromosome 2
>bNC_000069.7 Mus musculus strain B10.RIII chromosome 3
>bNC_000070.7 Mus musculus strain B10.RIII chromosome 4
>bNC_000071.7 Mus musculus strain B10.RIII chromosome 5
>bNC_000072.7 Mus musculus strain B10.RIII chromosome 6
>bNC_000073.7 Mus musculus strain B10.RIII chromosome 7
>bNC_000074.7 Mus musculus strain B10.RIII chromosome 8
>bNC_000075.7 Mus musculus strain B10.RIII chromosome 9
>bNC_000076.7 Mus musculus strain B10.RIII chromosome 10
>bNC_000077.7 Mus musculus strain B10.RIII chromosome 11
>bNC_000078.7 Mus musculus strain B10.RIII chromosome 12
>bNC_000079.7 Mus musculus strain B10.RIII chromosome 13
>bNC_000080.7 Mus musculus strain B10.RIII chromosome 14
>bNC_000081.7 Mus musculus strain B10.RIII chromosome 15
>bNC_000082.7 Mus musculus strain B10.RIII chromosome 16
>bNC_000083.7 Mus musculus strain B10.RIII chromosome 17
>bNC_000084.7 Mus musculus strain B10.RIII chromosome 18
>bNC_000085.7 Mus musculus strain B10.RIII chromosome 19
>bNC_000086.8 Mus musculus strain B10.RIII chromosome X
>bNC_000087.8 Mus musculus strain B10.RIII chromosome Y
>bNT_165789.3 Mus musculus strain B10.RIII chromosome X unlocalized genomic scaffold
>bNT_187064.1 Mus musculus strain B10.RIII unplaced unlocalized genomic scaffold
>bNW_023337853.1 Mus musculus strain B10.RIII unplaced unlocalized genomic scaffold
>Super-Scaffold_100310 Mus musculus strain B10.RIII unplaced unlocalized genomic scaffold
>Super-Scaffold_100970 Mus musculus strain B10.RIII unplaced unlocalized genomic scaffold
>Super-Scaffold_17 Mus musculus strain B10.RIII unplaced unlocalized genomic scaffold
>Super-Scaffold_187 Mus musculus strain B10.RIII unplaced unlocalized genomic scaffold
>Super-Scaffold_191 Mus musculus strain B10.RIII unplaced unlocalized genomic scaffold
>Super-Scaffold_202 Mus musculus strain B10.RIII unplaced unlocalized genomic scaffold
>Super-Scaffold_295 Mus musculus strain B10.RIII unplaced unlocalized genomic scaffold
>Super-Scaffold_329 Mus musculus strain B10.RIII unplaced unlocalized genomic scaffold
>Super-Scaffold_40 Mus musculus strain B10.RIII unplaced unlocalized genomic scaffold
>Super-Scaffold_45 Mus musculus strain B10.RIII unplaced unlocalized genomic scaffold
>Super-Scaffold_587 Mus musculus strain B10.RIII unplaced unlocalized genomic scaffold
>Super-Scaffold_61 Mus musculus strain B10.RIII unplaced unlocalized genomic scaffold
>Super-Scaffold_71 Mus musculus strain B10.RIII unplaced unlocalized genomic scaffold
>Super-Scaffold_74 Mus musculus strain B10.RIII unplaced unlocalized genomic scaffold
>Super-Scaffold_77 Mus musculus strain B10.RIII unplaced unlocalized genomic scaffold
>Super-Scaffold_78 Mus musculus strain B10.RIII unplaced unlocalized genomic scaffold
>Super-Scaffold_79 Mus musculus strain B10.RIII unplaced unlocalized genomic scaffold


sed -i 's/bNW_023337853.1/RC853 Mus musculus strain B10RIII unplaced unlocalized genomic scaffold/g' /data/CaspiWGSData/b10riii/Results/b10riii.fa
sed -i 's/bNT_187064.1/RC064 Mus musculus strain B10RIII unplaced unlocalized genomic scaffold/g' /data/CaspiWGSData/b10riii/Results/b10riii.fa
sed -i 's/bNT_165789.3/RC789 Mus musculus strain B10RIII chromosome X unlocalized genomic scaffold/g' /data/CaspiWGSData/b10riii/Results/b10riii.fa
sed -i 's/bNC_000087.8/RC087 Mus musculus strain B10RIII chromosome Y/g' /data/CaspiWGSData/b10riii/Results/b10riii.fa
sed -i 's/bNC_000086.8/RC086 Mus musculus strain B10RIII chromosome X/g' /data/CaspiWGSData/b10riii/Results/b10riii.fa
sed -i 's/bNC_000085.7/RC085 Mus musculus strain B10RIII chromosome 19/g' /data/CaspiWGSData/b10riii/Results/b10riii.fa
sed -i 's/bNC_000084.7/RC084 Mus musculus strain B10RIII chromosome 18/g' /data/CaspiWGSData/b10riii/Results/b10riii.fa
sed -i 's/bNC_000083.7/RC083 Mus musculus strain B10RIII chromosome 17/g' /data/CaspiWGSData/b10riii/Results/b10riii.fa
sed -i 's/bNC_000082.7/RC082 Mus musculus strain B10RIII chromosome 16/g' /data/CaspiWGSData/b10riii/Results/b10riii.fa
sed -i 's/bNC_000081.7/RC081 Mus musculus strain B10RIII chromosome 15/g' /data/CaspiWGSData/b10riii/Results/b10riii.fa
sed -i 's/bNC_000080.7/RC080 Mus musculus strain B10RIII chromosome 14/g' /data/CaspiWGSData/b10riii/Results/b10riii.fa
sed -i 's/bNC_000079.7/RC079 Mus musculus strain B10RIII chromosome 13/g' /data/CaspiWGSData/b10riii/Results/b10riii.fa
sed -i 's/bNC_000078.7/RC078 Mus musculus strain B10RIII chromosome 12/g' /data/CaspiWGSData/b10riii/Results/b10riii.fa
sed -i 's/bNC_000077.7/RC077 Mus musculus strain B10RIII chromosome 11/g' /data/CaspiWGSData/b10riii/Results/b10riii.fa
sed -i 's/bNC_000076.7/RC076 Mus musculus strain B10RIII chromosome 10/g' /data/CaspiWGSData/b10riii/Results/b10riii.fa
sed -i 's/bNC_000075.7/RC075 Mus musculus strain B10RIII chromosome 9/g' /data/CaspiWGSData/b10riii/Results/b10riii.fa
sed -i 's/bNC_000074.7/RC074 Mus musculus strain B10RIII chromosome 8/g' /data/CaspiWGSData/b10riii/Results/b10riii.fa
sed -i 's/bNC_000073.7/RC073 Mus musculus strain B10RIII chromosome 7/g' /data/CaspiWGSData/b10riii/Results/b10riii.fa
sed -i 's/bNC_000072.7/RC072 Mus musculus strain B10RIII chromosome 6/g' /data/CaspiWGSData/b10riii/Results/b10riii.fa
sed -i 's/bNC_000071.7/RC071 Mus musculus strain B10RIII chromosome 5/g' /data/CaspiWGSData/b10riii/Results/b10riii.fa
sed -i 's/bNC_000070.7/RC070 Mus musculus strain B10RIII chromosome 4/g' /data/CaspiWGSData/b10riii/Results/b10riii.fa
sed -i 's/bNC_000069.7/RC069 Mus musculus strain B10RIII chromosome 3/g' /data/CaspiWGSData/b10riii/Results/b10riii.fa
sed -i 's/bNC_000068.8/RC068 Mus musculus strain B10RIII chromosome 2/g' /data/CaspiWGSData/b10riii/Results/b10riii.fa
sed -i 's/bNC_000067.7/RC067 Mus musculus strain B10RIII chromosome 1/g' /data/CaspiWGSData/b10riii/Results/b10riii.fa



#### Annotation
cd /data/CaspiWGSData/b10riii/Results/annotation
https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/annotation_releases/current/GCF_000001635.27-RS_2023_04/
https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/annotation_releases/current/GCF_000001635.27-RS_2023_04/GCF_000001635.27_GRCm39_genomic.fna.gz
https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/annotation_releases/current/GCF_000001635.27-RS_2023_04/GCF_000001635.27_GRCm39_genomic.gff.gz
https://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_mammalian/Mus_musculus/latest_assembly_versions/GCA_030265425.1_NEI_Mmus_1.0/GCA_030265425.1_NEI_Mmus_1.0_genomic.fna.gz
https://www.gencodegenes.org/mouse/

feature types prep
cut -f 3 GCF_000001635.27_GRCm39_genomic.gff | sort | uniq > feature_types.txt
cut -f 3 gencode.vM33.chr_patch_hapl_scaff.annotation.gff3 | sort | uniq > feature_types_gencode.txt

The ncbi feature_types.txt has more feature types.

#liftoff -g GCF_000001635.27_GRCm39_genomic.gff -o b10.gff3 -polish -chroms chromosomes.txt -mm2_options="-t 55" -f feature_types.txt -copies -p 55 -m /usr/local/apps/minimap2/2.26/minimap2 GCA_030265425.1_NEI_Mmus_1.0_genomic.fna GCF_000001635.27_GRCm39_genomic.fna

# run as swarm
liftoff -g GCF_000001635.27_GRCm39_genomic.gff -o b10.gff3 -polish -chroms chromosomes.txt -f feature_types.txt -copies -p 55 -m /usr/local/apps/minimap2/2.26/minimap2 GCA_030265425.1_NEI_Mmus_1.0_genomic.fna GCF_000001635.27_GRCm39_genomic.fna

# about 2 days
/data/CaspiWGSData/b10riii/Results/annotation/swarmannotation/liftoff.swarm

# comparison

liftofftools all -r GCF_000001635.27_GRCm39_genomic.fna -t GCA_030265425.1_NEI_Mmus_1.0_genomic.fna -rg GCF_000001635.27_GRCm39_genomic.gff -tg b10_swarm.gff3_db

module load python/3.8
import gffutils
fn = gffutils.example_filename('/vf/users/CaspiWGSData/b10riii/Results/annotation/swarmannotation/b10_swarm.gff3_polished')
#db = gffutils.create_db(fn,'b10_swarm.gff3_db',force=True, keep_order=True,merge_strategy='create_unique', sort_attribute_values=True)
db = gffutils.create_db(fn,'b10_swarm.gff3_polished_db',merge_strategy="warning")

module load agat 0.8.0
agat_sp_statistics.pl --gff b10_swarm.gff3






