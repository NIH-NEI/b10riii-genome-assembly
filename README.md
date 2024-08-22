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
ipa local --nthreads 20 --njobs 1 --run-dir /data/CaspiWGSData/b10riii/Results/pacbio/pbipa/default -i /data/CaspiWGSData/b10riii/RawData/pacbio/pacbio_ccs/F1_1.ccs.fastq -i /data/CaspiWGSData/b10riii/RawData/pacbio/pacbio_ccs/F1_2.ccs.fastq -i /data/CaspiWGSData/b10riii/RawData/pacbio/pacbio_ccs/F1_3.ccs.fastq -i /data/CaspiWGSData/b10riii/RawData/pacbio/pacbio_ccs/F1_4.ccs.fastq

# filtered
ipa local --nthreads 20 --njobs 1 --run-dir /data/CaspiWGSData/b10riii/Results/pacbio/pbipa/defaultfiltered -i /data/CaspiWGSData/b10riii/Results/pacbio/hifiadapterfilt/F1_1.ccs.filt.fastq.gz -i /data/CaspiWGSData/b10riii/Results/pacbio/hifiadapterfilt/F1_2.ccs.filt.fastq.gz -i /data/CaspiWGSData/b10riii/Results/pacbio/hifiadapterfilt/F1_3.ccs.filt.fastq.gz -i /data/CaspiWGSData/b10riii/Results/pacbio/hifiadapterfilt/F1_4.ccs.filt.fastq.gz

# use --resume --unlock, unlock first and resume again
# Mon May 23 20:26:02 EDT 2022 - need more than 4 days

# to check snakemake command
ipa local -i ../Results/pacbio/pbipa/default/input.fofn --only-print --nthreads 20 --njobs 1

# below is the above commands output
INFO: /data/CaspiWGSData/b10riii/Tools/conda/envs/ipa/bin/ipa local -i ../Results/pacbio/pbipa/default/input.fofn --only-print --nthreads 20 --njobs 1
INFO: ipa.py ipa (wrapper) version=1.8.0 ... Checking dependencies ...
INFO: Dependencies
/data/CaspiWGSData/b10riii/Tools/conda/envs/ipa/bin/python3
/data/CaspiWGSData/b10riii/Tools/conda/envs/ipa/bin/ipa2-task
/data/CaspiWGSData/b10riii/Tools/conda/envs/ipa/bin/falconc
/data/CaspiWGSData/b10riii/Tools/conda/envs/ipa/bin/minimap2
/data/CaspiWGSData/b10riii/Tools/conda/envs/ipa/bin/nighthawk
/data/CaspiWGSData/b10riii/Tools/conda/envs/ipa/bin/pancake
/data/CaspiWGSData/b10riii/Tools/conda/envs/ipa/bin/pblayout
/data/CaspiWGSData/b10riii/Tools/conda/envs/ipa/bin/racon
/data/CaspiWGSData/b10riii/Tools/conda/envs/ipa/bin/samtools
/data/CaspiWGSData/b10riii/Tools/conda/envs/ipa/bin/ipa_purge_dups
/data/CaspiWGSData/b10riii/Tools/conda/envs/ipa/bin/ipa_purge_dups_split_fa
snakemake version=7.7.0
ipa2-task 1.8.0 (commit 89f9124a7aa07d170db0eb490c12a5aeb36bcf64)
 Machine name: 'Linux'
Copyright (C) 2004-2021     Pacific Biosciences of California, Inc.
This program comes with ABSOLUTELY NO WARRANTY; it is intended for
Research Use Only and not for use in diagnostic procedures.

falconc version=1.15.0+git., nim-version=1.2.0
minimap2 version=2.24-r1122
Nighthawk 1.0.0 (commit SL-release-10.2.0-11-g332ccaa)
pancake 1.6.0 (commit SL-release-10.2.0-175-gb750c42)
pblayout 1.4.0 (commit SL-release-10.2.0-64-gd59326a)
racon version=1.5.0
samtools 1.15.1
Using htslib 1.15.1
ipa_purge_dups Version: 1.2.5


To run this yourself:
/data/CaspiWGSData/b10riii/Tools/conda/envs/ipa/bin/python3 -m snakemake -j 1 -d RUN -p -s /data/CaspiWGSData/b10riii/Tools/conda/envs/ipa/etc/ipa.snakefile --configfile RUN/config.yaml --reason

snakemake -s /data/CaspiWGSData/b10riii/Tools/conda/envs/ipa/etc/ipa.snakefile --list-params-changes

# Repeated pbipa with output in default2 - no change in parameters
quast.py -o ../Results/quast/pbipa ../Results/pacbio/pbipa/default/assembly-results/final.p_ctg.fasta ../Results/pacbio/pbipa/default/assembly-results/final.a_ctg.fasta ../Results/pacbio/pbipa/default2/assembly-results/final.p_ctg.fasta ../Results/pacbio/pbipa/default2/assembly-results/final.a_ctg.fasta ../RawData/GCA_000001635.9_GRCm39_genomic.fna.gz --threads 60

#hifiasm (hifiasm-0.16.1) install through conda
#https://hifiasm.readthedocs.io/en/latest/pa-assembly.html
hifiasm -o ../Results/hifiasm/default.asm -t 600 ../RawData/pacbio/pacbio_ccs/F1_1.ccs.fastq ../RawData/pacbio/pacbio_ccs/F1_2.ccs.fastq ../RawData/pacbio/pacbio_ccs/F1_3.ccs.fastq ../RawData/pacbio/pacbio_ccs/F1_4.ccs.fastq

#interactive session didnt work for hifiasm, trying sbatch
#sbatch script in Tools folder
sbatch --cpus-per-task=60 --mem=1507g --time=10-00:00:00 --partition=largemem hifiasm60.sbatch

sbatch --cpus-per-task=60 --mem=1507g --time=10-00:00:00 --partition=largemem /data/CaspiWGSData/b10riii/Tools/hifiasm60filtered.sbatch

#convert gfa to fa
awk '/^S/{print ">"$2;print $3}' default-60-filtered.asm.bp.p_ctg.gfa > default-60-filtered.asm.bp.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' default.asm.bp.p_ctg.gfa > default.asm.bp.p_ctg.fa

# https://hifiasm.readthedocs.io/en/latest/faq.html
# https://github.com/chhylp123/hifiasm

# quast with hicanu, pbipa, hifiasm

quast.py -o ../Results/quast/hicanu-pbipa-hifiasm /data/CaspiWGSData/b10riii/Results/pacbio/hicanu/hicanudefault3/hicanudefault3.contigs.fasta /data/CaspiWGSData/b10riii/Results/pacbio/hicanu/hicanudefault4/hicanudefault4.contigs.fasta ../Results/pacbio/pbipa/default/assembly-results/final.p_ctg.fasta ../Results/pacbio/pbipa/default2/assembly-results/final.p_ctg.fasta /data/CaspiWGSData/b10riii/Results/pacbio/hifiasm/default-60/default.asm.bp.p_ctg.fa /data/CaspiWGSData/b10riii/Results/pacbio/hifiasm/default-60-filtered/default-60-filtered.asm.bp.p_ctg.fa ../RawData/GCA_000001635.9_GRCm39_genomic.fna.gz --threads 60

quast.py -o /data/CaspiWGSData/b10riii/Results/quast/hicanu-pbipa-hifiasm /data/CaspiWGSData/b10riii/Results/pacbio/hicanu/hicanudefault3/hicanudefault3.contigs.fasta /data/CaspiWGSData/b10riii/Results/pacbio/hicanu/hicanudefault4/hicanudefault4.contigs.fasta /data/CaspiWGSData/b10riii/Results/pacbio/pbipa/default/assembly-results/final.p_ctg.fasta /data/CaspiWGSData/b10riii/Results/pacbio/pbipa/default2/assembly-results/final.p_ctg.fasta /data/CaspiWGSData/b10riii/Results/pacbio/hifiasm/default-60/default.asm.bp.p_ctg.fa /data/CaspiWGSData/b10riii/Results/pacbio/hifiasm/default-60-filtered/default-60-filtered.asm.bp.p_ctg.fa /data/CaspiWGSData/b10riii/Results/chromosome/test/GCF_000001635.27_GRCm39_genomic.fna.asmdefaultfiltered.bp.p_ctg.fa.split.reconciled.fa /data/CaspiWGSData/b10riii/RawData/ncbi_dataset/data/GCF_000001635.27/GCF_000001635.27_GRCm39_genomic.fna --threads 120


# all pacbio assemblies in /data/CaspiWGSData/b10riii/Results/pacbio/assemblies

# Polishing using arrow
# align using pbmm2 - conda
# polish using arrow/gcpp - conda

source /data/CaspiWGSData/b10riii/Tools/conda/etc/profile.d/conda.sh
conda activate pbmm2

/data/CaspiWGSData/b10riii/RawData/pacbio/pacbio_subreads/
F1_1.subreads.bam
F1_2.subreads.bam
F1_3.subreads.bam
F1_4.subreads.bam

cd /data/CaspiWGSData/b10riii/Results/polished/pacbiopolished
echo "/data/CaspiWGSData/b10riii/RawData/pacbio/pacbio_subreads/F1_1.subreads.bam" > subreadbams.fofn
echo "/data/CaspiWGSData/b10riii/RawData/pacbio/pacbio_subreads/F1_2.subreads.bam" >> subreadbams.fofn
echo "/data/CaspiWGSData/b10riii/RawData/pacbio/pacbio_subreads/F1_3.subreads.bam" >> subreadbams.fofn
echo "/data/CaspiWGSData/b10riii/RawData/pacbio/pacbio_subreads/F1_4.subreads.bam" >> subreadbams.fofn
cat 

/data/CaspiWGSData/b10riii/Results/pacbio/assemblies/

asmdefault.bp.p_ctg.fa
asmdefaultfiltered.bp.p_ctg.fa
pbipadefault2.final.p_ctg.fasta
pbipadefaultfiltered.final.p_ctg.fasta
hicanudefault3.contigs.fasta
hicanudefault4.contigs.fasta

/data/CaspiWGSData/b10riii/Results/polished/pacbiopolished

pbmm2 align /data/CaspiWGSData/b10riii/Results/pacbio/assemblies/asmdefault.bp.p_ctg.fa subreadbams.fofn asmaligned.bam --sort -j 80 -J 60 -m 32G --preset SUBREAD
-j, -J number of threads for alignment and sorting

SBATCH_PARTITION=largemem swarm -t 72 -g 1507 pbmm2-alignment.swarm

#parallelize alignment step
source /data/CaspiWGSData/b10riii/Tools/conda/etc/profile.d/conda.sh ; TMPDIR="/data/CaspiWGSData/b10riii/" ; conda activate pbmm2 ; cd /data/CaspiWGSData/b10riii/Results/polished/pacbiopolished ; pbmm2 align /data/CaspiWGSData/b10riii/Results/pacbio/assemblies/asmdefaultfiltered.bp.p_ctg.fa subreadbams.fofn asmfilteredaligned-norm.bam --sort -j 56 -J 30 -m 8G --preset SUBREAD --log-level INFO --log-file asmdefaultfiltered-norm.log

#below generates the swarm file
for i in `cat sequence_heads.txt`; do echo "source /data/CaspiWGSData/b10riii/Tools/conda/etc/profile.d/conda.sh ; conda activate pbmm2 ; export TMPDIR=/lscratch/\$SLURM_JOB_ID ; cd /data/CaspiWGSData/b10riii/Results/polished/pacbiopolished/gcpp-parallel/alignments ; pbmm2 align /data/CaspiWGSData/b10riii/Results/polished/pacbiopolished/gcpp-parallel/"$i".fa /data/CaspiWGSData/b10riii/Results/polished/pacbiopolished/subreadbams.fofn asmdefaultfiltered-norm-aligned-"$i".bam --sort -j 56 -J 30 -m 8G --preset SUBREAD --log-level DEBUG --log-file asmdefaultfiltered-norm-aligned-"$i".log"; done


swarm -f /data/CaspiWGSData/b10riii/Tools/pbmm2-alignment-norm-parallel.swarm -g 247 -t 56 --gres=lscratch:800

ls -lh | grep "samtools" | cut -d'.' -f4 | sort | uniq -c


# GCPP arrow - polishing
source /data/CaspiWGSData/b10riii/Tools/conda/etc/profile.d/conda.sh
conda activate gcpp
#https://www.biostars.org/p/273447/

sbatch --cpus-per-task=56 --mem=247g --time=10-00:00:00 --partition=quick /data/CaspiWGSData/b10riii/Tools/gcpp-polish.swarm

### parallelize gcpp polishing - 1 day for 12 contigs
#https://raw.githubusercontent.com/harish0201/General_Scripts/master/Fasta_splitter.py
cd /data/CaspiWGSData/b10riii/Results/polished/pacbiopolished/gcpp-parallel
python ../../../../Tools/Fasta_splitter.py /data/CaspiWGSData/b10riii/Results/pacbio/assemblies/asmdefaultfiltered.bp.p_ctg.fa > sequence_heads.txt

#below generates the swarm file
for i in `cat /data/CaspiWGSData/b10riii/Results/polished/pacbiopolished/gcpp-parallel/sequence_heads.txt`; do echo "source /data/CaspiWGSData/b10riii/Tools/conda/etc/profile.d/conda.sh ; TMPDIR="/data/CaspiWGSData/b10riii/" ; conda activate gcpp ; cd /data/CaspiWGSData/b10riii/Results/polished/pacbiopolished/gcpp-parallel/polished_seqs ; gcpp -r /data/CaspiWGSData/b10riii/Results/polished/pacbiopolished/gcpp-parallel/"$i".fa -o asmdefaultfiltered-norm-gcpp.polished."$i".fasta asmfilteredaligned-norm.bam --log-level TRACE --log-file asmdefaultfiltered-norm-gcpp-"$i".log"; done

swarm -f /data/CaspiWGSData/b10riii/Tools/gcpp-polish-parallel.swarm -g 121 -t 20 --gres=lscratch:300






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
singularity shell -B /usr/lib/locale/:/usr/lib/locale/ --bind /data/CaspiWGSData/ singularity-bionano-perl-python-R-time-zip.sif
source /data/CaspiWGSData/b10riii/Tools/conda/etc/profile.d/conda.sh
conda activate python37
cd perl-bionano

/data/CaspiWGSData/b10riii/Results/hybridscaffold/hifiasm-60-filtered/hybrid_scaffolds/EXP_REFINEFINAL1_bppAdjust_cmap_default-60-filtered_asm_bp_p_ctg_fa_NGScontigs_HYBRID_SCAFFOLD.fasta





--------------------------------
https://www.biostars.org/p/273447/
April 22
arrow - bam sliced
Create index
cd /data/CaspiWGSData/b10riii/Results/polished/pacbiopolished/gcpp-parallel/polished_seqs
samtools faidx /data/CaspiWGSData/b10riii/Results/pacbio/assemblies/asmdefaultfiltered.bp.p_ctg.fa
cut -f1 /data/CaspiWGSData/b10riii/Results/pacbio/assemblies/asmdefaultfiltered.bp.p_ctg.fa.fai > sequence_heads.txt
# also faidx index split fasta files
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
cut -f1-2 /data/CaspiWGSData/b10riii/Results/pacbio/assemblies/asmdefaultfiltered.bp.p_ctg.fa.fai | awk ' { print $1 ":0-" $2 } ' > contigs.txt
awk "NR > 0 && NR <= 100" contigs.txt | paste -s -d ,
gcpp -w ptg000111l:0-16631 -r /data/CaspiWGSData/b10riii/Results/pacbio/assemblies/asmdefaultfiltered.bp.p_ctg.fa -o /data/CaspiWGSData/b10riii/Results/polished/pacbiopolished/gcpp-parallel/polished_seqs/ptg000111l.polished.fasta /data/CaspiWGSData/b10riii/Results/polished/pacbiopolished/gcpp-parallel/polished_seqs/ptg000111l.bam



#1
ccs /data/Sen_Lab_NEI/LRseq/NTLRseq/RawData/DATA_20210315/NS3R189BNTLR/AXI0428_1/SMRTCell_1546/m54313U_210107_193853.subreads.bam /data/Sen_Lab_NEI/LRseq/NTLRseq/Results/NS3R189BNTLR/NS3R189BNTLR_1.ccs.bam --min-rq 0.9 --report-file /data/Sen_Lab_NEI/LRseq/NTLRseq/Results/NS3R189BNTLR/ccs_report.txt
# --min-rq                  FLOAT  Minimum predicted accuracy in [0, 1]. [0.99]
# About 9 hrs to complete

#2
lima /data/Sen_Lab_NEI/LRseq/NTLRseq/Results/NS3R189BNTLR/NS3R189BNTLR_1.ccs.bam /data/Sen_Lab_NEI/LRseq/NTLRseq/Protocols/primer.fasta /data/Sen_Lab_NEI/LRseq/NTLRseq/Results/NS3R189BNTLR/NS3R189BNTLR_1.fl.bam --isoseq --peek-guess --log-file /data/Sen_Lab_NEI/LRseq/NTLRseq/Results/NS3R189BNTLR/lima_log.txt
# About 10 minutes to complete

#3
isoseq3 refine /data/Sen_Lab_NEI/LRseq/NTLRseq/Results/NS3R189BNTLR/NS3R189BNTLR_1.fl.NEB_5p--NEB_Clontech_3p.bam /data/Sen_Lab_NEI/LRseq/NTLRseq/Protocols/primer.fasta /data/Sen_Lab_NEI/LRseq/NTLRseq/Results/NS3R189BNTLR/NS3R189BNTLR_1.flnc.bam --require-polya
# About 5 minutes to complete

#4
isoseq3 cluster  /data/Sen_Lab_NEI/LRseq/NTLRseq/Results/NS3R189BNTLR/NS3R189BNTLR_1.flnc.bam /data/Sen_Lab_NEI/LRseq/NTLRseq/Results/NS3R189BNTLR/NS3R189BNTLR_1.clustered.bam --verbose --split-bam 24
# About 45 minutes to complete

#5
isoseq3 polish /data/Sen_Lab_NEI/LRseq/NTLRseq/Results/NS3R189BNTLR/NS3R189BNTLR_1.clustered.bam /data/Sen_Lab_NEI/LRseq/NTLRseq/RawData/DATA_20210315/NS3R189BNTLR/AXI0428_1/SMRTCell_1546/m54313U_210107_193853.subreads.bam /data/Sen_Lab_NEI/LRseq/NTLRseq/Results/NS3R189BNTLR/NS3R189BNTLR_1_polished.bam
# About 1 hour for one split file
# Run as a swarm job
# https://github.com/graphanalytics/IsoSeq/blob/master/isoseq-clustering.md

#6 Process files for cupcake input
gunzip -c NS3R189BNTLR_1_polished.*.hq.fastq.gz > hq_isoforms.fastq
cat NS3R189BNTLR_1_polished.*.cluster_report.csv > combined.cluster_report.csv

#7 CupCake
ml minimap2
minimap2 -ax splice -t 30 -uf --secondary=no -C5 /data/Sen_Lab_NEI/LRseq/NTLRseq/Protocols/FILES_20220314/hg38.fasta hq_isoforms.fastq > hq_isoforms.sam
# About 2 minutes
sort -k 3,3 -k 4,4n hq_isoforms.sam > hq_isoforms.sorted.sam




collapse_isoforms_by_sam.py --fq --input hq_isoforms.fastq -s hq_isoforms.sorted.sam --dun-merge-5-shorter -o transcripts
get_abundance_post_collapse.py transcripts.collapsed combined.cluster_report.csv
