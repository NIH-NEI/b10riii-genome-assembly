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

## 6. Quast hicanu
`module load quast`\
#Loading quast, version 5.0.2
`quast.py -o quast/hicanu pacbio/hicanu/hicanudefault3/hicanudefault3.contigs.fasta pacbio/hicanu/hicanudefault4/hicanudefault4.contigs.fasta ../RawData/GCA_000001635.9_GRCm39_genomic.fna.gz`\
`quast.py -o quast/hicanu/vsmm39 -r ../RawData/GCA_000001635.9_GRCm39_genomic.fna.gz -g ../RawData/GCA_000001635.9_GRCm39_genomic.gff.gz pacbio/hicanu/hicanudefault3/hicanudefault3.contigs.fasta pacbio/hicanu/hicanudefault4/hicanudefault4.contigs.fasta`

## 7. pbipa assembly
`ipa local --nthreads 20 --njobs 1 --run-dir /data/CaspiWGSData/b10riii/Results/pacbio/pbipa/default -i /b10riii/RawData/pacbio/pacbio_ccs/F1_1.ccs.fastq -i /b10riii/RawData/pacbio/pacbio_ccs/F1_2.ccs.fastq -i /b10riii/RawData/pacbio/pacbio_ccs/F1_3.ccs.fastq -i /b10riii/RawData/pacbio/pacbio_ccs/F1_4.ccs.fastq`

### 7.1. pbipa assembly filtered
`ipa local --nthreads 20 --njobs 1 --run-dir /b10riii/Results/pacbio/pbipa/defaultfiltered -i /b10riii/Results/pacbio/hifiadapterfilt/F1_1.ccs.filt.fastq.gz -i /b10riii/Results/pacbio/hifiadapterfilt/F1_2.ccs.filt.fastq.gz -i /b10riii/Results/pacbio/hifiadapterfilt/F1_3.ccs.filt.fastq.gz -i /b10riii/Results/pacbio/hifiadapterfilt/F1_4.ccs.filt.fastq.gz`\
#to check snakemake command\
`ipa local -i ../Results/pacbio/pbipa/default/input.fofn --only-print --nthreads 20 --njobs 1`

### 7.2. pbipa assembly quality
`quast.py -o ../Results/quast/pbipa ../Results/pacbio/pbipa/default/assembly-results/final.p_ctg.fasta ../Results/pacbio/pbipa/default/assembly-results/final.a_ctg.fasta ../Results/pacbio/pbipa/default2/assembly-results/final.p_ctg.fasta ../Results/pacbio/pbipa/default2/assembly-results/final.a_ctg.fasta ../RawData/GCA_000001635.9_GRCm39_genomic.fna.gz --threads 60`

## 8. hifiasm
`hifiasm -o ../Results/hifiasm/default.asm -t 600 ../RawData/pacbio/pacbio_ccs/F1_1.ccs.fastq ../RawData/pacbio/pacbio_ccs/F1_2.ccs.fastq ../RawData/pacbio/pacbio_ccs/F1_3.ccs.fastq ../RawData/pacbio/pacbio_ccs/F1_4.ccs.fastq`\
`sbatch --cpus-per-task=60 --mem=1507g --time=10-00:00:00 --partition=largemem /b10riii/Tools/hifiasm60filtered.sbatch`
#convert gfa to fa
`awk '/^S/{print ">"$2;print $3}' default-60-filtered.asm.bp.p_ctg.gfa > default-60-filtered.asm.bp.p_ctg.fa`\
`awk '/^S/{print ">"$2;print $3}' default.asm.bp.p_ctg.gfa > default.asm.bp.p_ctg.fa`

### 8.1 quast with hicanu, pbipa, hifiasm
`quast.py -o /b10riii/Results/quast/hicanu-pbipa-hifiasm /b10riii/Results/pacbio/hicanu/hicanudefault3/hicanudefault3.contigs.fasta /b10riii/Results/pacbio/hicanu/hicanudefault4/hicanudefault4.contigs.fasta /b10riii/Results/pacbio/pbipa/default/assembly-results/final.p_ctg.fasta /b10riii/Results/pacbio/pbipa/default2/assembly-results/final.p_ctg.fasta /b10riii/Results/pacbio/hifiasm/default-60/default.asm.bp.p_ctg.fa /b10riii/Results/pacbio/hifiasm/default-60-filtered/default-60-filtered.asm.bp.p_ctg.fa /b10riii/Results/chromosome/test/GCF_000001635.27_GRCm39_genomic.fna.asmdefaultfiltered.bp.p_ctg.fa.split.reconciled.fa /b10riii/RawData/ncbi_dataset/data/GCF_000001635.27/GCF_000001635.27_GRCm39_genomic.fna --threads 120`

## 9. Polishing using long reads
`source /b10riii/Tools/conda/etc/profile.d/conda.sh`\
`conda activate pbmm2`\
`cd /b10riii/Results/polished/pacbiopolished`\
`echo "/b10riii/RawData/pacbio/pacbio_subreads/F1_1.subreads.bam" > subreadbams.fofn`\
`echo "/b10riii/RawData/pacbio/pacbio_subreads/F1_2.subreads.bam" >> subreadbams.fofn`\
`echo "/data/CaspiWGSData/b10riii/RawData/pacbio/pacbio_subreads/F1_3.subreads.bam" >> subreadbams.fofn`\
`echo "/data/CaspiWGSData/b10riii/RawData/pacbio/pacbio_subreads/F1_4.subreads.bam" >> subreadbams.fofn`\
`pbmm2 align /b10riii/Results/pacbio/assemblies/asmdefault.bp.p_ctg.fa subreadbams.fofn asmaligned.bam --sort -j 80 -J 60 -m 32G --preset SUBREAD
-j, -J number of threads for alignment and sorting`\

#parallelize alignment step\
#below generates the swarm file\
`for i in `cat sequence_heads.txt`; do echo "source /b10riii/Tools/conda/etc/profile.d/conda.sh ; conda activate pbmm2 ; export TMPDIR=/lscratch/\$SLURM_JOB_ID ; cd /b10riii/Results/polished/pacbiopolished/gcpp-parallel/alignments ; pbmm2 align /b10riii/Results/polished/pacbiopolished/gcpp-parallel/"$i".fa /b10riii/Results/polished/pacbiopolished/subreadbams.fofn asmdefaultfiltered-norm-aligned-"$i".bam --sort -j 56 -J 30 -m 8G --preset SUBREAD --log-level DEBUG --log-file asmdefaultfiltered-norm-aligned-"$i".log"; done`\
`swarm -f /b10riii/Tools/pbmm2-alignment-norm-parallel.swarm -g 247 -t 56 --gres=lscratch:800`\
`ls -lh | grep "samtools" | cut -d'.' -f4 | sort | uniq -c`\

### 9.1. GCPP arrow - polishing
`source /b10riii/Tools/conda/etc/profile.d/conda.sh`\
`conda activate gcpp`\
#parallelize gcpp polishing \
`cd /b10riii/Results/polished/pacbiopolished/gcpp-parallel`\
`python ../../../../Tools/Fasta_splitter.py /b10riii/Results/pacbio/assemblies/asmdefaultfiltered.bp.p_ctg.fa > sequence_heads.txt`

#below generates the swarm file
`for i in `cat /b10riii/Results/polished/pacbiopolished/gcpp-parallel/sequence_heads.txt`; do echo "source /b10riii/Tools/conda/etc/profile.d/conda.sh ; TMPDIR="/b10riii/" ; conda activate gcpp ; cd /b10riii/Results/polished/pacbiopolished/gcpp-parallel/polished_seqs ; gcpp -r /b10riii/Results/polished/pacbiopolished/gcpp-parallel/"$i".fa -o asmdefaultfiltered-norm-gcpp.polished."$i".fasta asmfilteredaligned-norm.bam --log-level TRACE --log-file asmdefaultfiltered-norm-gcpp-"$i".log"; done`\
`swarm -f /b10riii/Tools/gcpp-polish-parallel.swarm -g 121 -t 20 --gres=lscratch:300`\

## 10. Hybrid assembly
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
#bionano uses 2.7 in code, but install requires 3.7.7
`ln -s /b10riii/Tools/conda/envs/python37/bin/python /b10riii/Tools/conda/envs/python37/bin/python2.7`\

`perl /b10riii/Tools/perl-bionano/tools/pipeline/1.0/HybridScaffold/1.0/hybridScaffold.pl -n /b10riii/Results/pacbio/assemblies/defaultfiltered.asm.bp.p_ctg.fa -b /b10riii/RawData/bionano/Assembly_data_delivery/output/contigs/exp_refineFinal1/EXP_REFINEFINAL1.cmap -c /b10riii/Tools/perl-bionano/tools/pipeline/1.0/HybridScaffold/1.0/hybridScaffold_DLE1_config.xml -r /b10riii/Tools/perl-bionano/tools/pipeline/1.0/RefAligner/1.0/RefAligner -o /b10riii/Results/hybridscaffold/hifiasm-filtered -f -g -B 2 -N 2`\

cd /data/CaspiWGSData/b10riii/Tools/
singularity shell -B /usr/lib/locale/:/usr/lib/locale/ --bind /data/CaspiWGSData/ singularity-bionano-perl-python-R-time-zip.sif
source /data/CaspiWGSData/b10riii/Tools/conda/etc/profile.d/conda.sh
conda activate python37
cd perl-bionano

## 11. Polishing using short reads
`cd /b10riii/Results/polished/pacbiopolished/gcpp-parallel/polished_seqs`\
`samtools faidx /b10riii/Results/pacbio/assemblies/asmdefaultfiltered.bp.p_ctg.fa`\
`cut -f1 /b10riii/Results/pacbio/assemblies/asmdefaultfiltered.bp.p_ctg.fa.fai > sequence_heads.txt`\
`cd /b10riii/Results/polished/pacbiopolished/gcpp-parallel`\
`for i in `cat sequence_heads.txt`; do samtools faidx "$i".fa ; done`\

#below generates the swarm file\
`for i in `cat sequence_heads.txt`; do echo "source /b10riii/Tools/conda/etc/profile.d/conda.sh ; module load bamtools ; TMPDIR=/b10riii/ ; cd /b10riii/Results/polished/pacbiopolished/gcpp-parallel/polished_seqs ; bamtools merge -in asmfilteredaligned-norm.bam -out "$i".bam -region "$i" ; conda activate pbindex ; pbindex "$i".bam"; done > bam-split-index.swarm`\

`swarm -f /b10riii/Tools/bam-split-index.swarm`\ 

conda activate gcpp
gcpp -r "$i".fasta -o polished_seqs/"$i".polished.fastq -o polished_seqs/"$i".polished.vcf -o polished_seqs/"$i".polished.gff -n 5 --annotateGFF --reportEffectiveCoverage "$i".bam"'
gcpp -r /data/CaspiWGSData/b10riii/Results/polished/pacbiopolished/gcpp-parallel/ptg000111l.fa -o /data/CaspiWGSData/b10riii/Results/polished/pacbiopolished/gcpp-parallel/polished_seqs/ptg000111l.polished.fastq --reportEffectiveCoverage /data/CaspiWGSData/b10riii/Results/polished/pacbiopolished/gcpp-parallel/polished_seqs/ptg000111l.bam

Test polishing using gcpp -w windows option
cut -f1-2 /data/CaspiWGSData/b10riii/Results/pacbio/assemblies/asmdefaultfiltered.bp.p_ctg.fa.fai | awk ' { print $1 ":0-" $2 } ' > contigs.txt
awk "NR > 0 && NR <= 100" contigs.txt | paste -s -d ,
gcpp -w ptg000111l:0-16631 -r /data/CaspiWGSData/b10riii/Results/pacbio/assemblies/asmdefaultfiltered.bp.p_ctg.fa -o /data/CaspiWGSData/b10riii/Results/polished/pacbiopolished/gcpp-parallel/polished_seqs/ptg000111l.polished.fasta /data/CaspiWGSData/b10riii/Results/polished/pacbiopolished/gcpp-parallel/polished_seqs/ptg000111l.bam


