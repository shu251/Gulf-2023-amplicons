#!/bin/bash
##ENVIRONMENT SETTINGS; CHANGE WITH CAUTION
#SBATCH --export=NONE        #Do not propagate environment
#SBATCH --get-user-env=L     #Replicate login environment

##NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=18S-gom
#SBATCH --time=72:00:00            # 72 hours
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24        # 32 threads
#SBATCH --mem=340GB               # lots of memory, but never the maximum (which is 384G)
#SBATCH --output=18S-qiime2-log.%j

# Load qiime2
module load QIIME2/2024.10-Amplicon

echo "import sequences"

qiime tools import \
 --type 'SampleData[PairedEndSequencesWithQuality]' \
 --input-path /scratch/group/hu-lab/data/tag-seq-data/GoM-2023-18S_2025-05/manifest-gom-2023 \
 --output-path /scratch/group/hu-lab/GoM-amplicon-analysis/GoM-2023-18S_output/pe-gom-18s.qza \
 --input-format PairedEndFastqManifestPhred33V2

qiime demux summarize \
 --i-data /scratch/group/hu-lab/GoM-amplicon-analysis/GoM-2023-18S_output/pe-gom-18s.qza \
 --o-visualization /scratch/group/hu-lab/GoM-amplicon-analysis/GoM-2023-18S_output/pe-gom-18s.qzv

echo "Use cutadapt to remove primer sequences"

qiime cutadapt trim-paired \
 --i-demultiplexed-sequences /scratch/group/hu-lab/GoM-amplicon-analysis/GoM-2023-18S_output/pe-gom-18s.qza \
 --p-cores $SLURM_CPUS_PER_TASK \
 --p-front-f CCAGCASCYGCGGTAATTCC \
 --p-front-r ACTTTCGTTCTTGATYRA \
 --p-error-rate 0.1 \
 --p-overlap 3 \
 --p-match-adapter-wildcards \
 --o-trimmed-sequences /scratch/group/hu-lab/GoM-amplicon-analysis/GoM-2023-18S_output/trimmed-gom-18s.qza

qiime demux summarize \
  --i-data /scratch/group/hu-lab/GoM-amplicon-analysis/GoM-2023-18S_output/trimmed-gom-18s.qza \
  --o-visualization /scratch/group/hu-lab/GoM-amplicon-analysis/GoM-2023-18S_output/trimmed-gom-18s.qzv

echo "Starting DADA2 denoising, error prediction, chimera removal, and ASV determination..."

qiime dada2 denoise-paired \
        --i-demultiplexed-seqs /scratch/group/hu-lab/GoM-amplicon-analysis/GoM-2023-18S_output/trimmed-gom-18s.qza \
        --p-trunc-len-f 260 \
        --p-trunc-len-r 225 \
        --p-max-ee-f 2 \
        --p-max-ee-r 2 \
        --p-min-overlap 10 \
        --p-pooling-method independent \
        --p-n-reads-learn 1000000 \
        --p-n-threads $SLURM_CPUS_PER_TASK \
        --p-chimera-method pooled \
	--o-table /scratch/group/hu-lab/GoM-amplicon-analysis/GoM-2023-18S_output/asv-table-gom-18s.qza \
	--o-representative-sequences /scratch/group/hu-lab/GoM-amplicon-analysis/GoM-2023-18S_output/refseqs-gom-18s.qza \
	--o-denoising-stats /scratch/group/hu-lab/GoM-amplicon-analysis/GoM-2023-18S_output/statsDADA2-gom-18s.qza

echo "dada2 step completed"

qiime tools export \
	--input-path /scratch/group/hu-lab/GoM-amplicon-analysis/GoM-2023-18S_output/asv-table-gom-18s.qza \
	--output-path /scratch/group/hu-lab/GoM-amplicon-analysis/GoM-2023-18S_output/18S-output-2023

biom convert \
	-i /scratch/group/hu-lab/GoM-amplicon-analysis/GoM-2023-18S_output/18S-output-2023/feature-table.biom \
	-o /scratch/group/hu-lab/GoM-amplicon-analysis/GoM-2023-18S_output/gom-2023-18s-asv-table.tsv \
	--to-tsv

echo "asv table is prepped"

qiime metadata tabulate \
	--m-input-file /scratch/group/hu-lab/GoM-amplicon-analysis/GoM-2023-18S_output/statsDADA2-gom-18s.qza \
	--o-visualiation /scratch/group/hu-lab/GoM-amplicon-analysis/GoM-2023-18S_output/statsDADA2-gom-18s.qzv


echo "assign taxonomy"

qiime feature-classifier classify-consensus-vsearch \
	--i-query /scratch/group/hu-lab/GoM-amplicon-analysis/GoM-2023-18S_output/refseqs-gom-18s.qza \
        --i-reference-reads /scratch/group/hu-lab/databases/PR2/pr2_version_5.0.0_seqs.qza \
        --i-reference-taxonomy /scratch/group/hu-lab/databases/PR2/pr2_version_5.0.0_tax.qza \
	--output-dir /scratch/group/hu-lab/GoM-amplicon-analysis/taxonomy-0.9_0.8_PR2-5.0.0 \
        --p-threads $SLURM_CPUS_PER_TASK \
        --p-maxaccepts 10 \
        --p-perc-identity 0.90 \
        --p-min-consensus 0.80

echo "taxonomy assigned, convert output"

qiime tools export \
        --input-path /scratch/group/hu-lab/GoM-amplicon-analysis/taxonomy-0.9_0.8_PR2-5.0.0/classification.qza \
        --output-path /scratch/group/hu-lab/GoM-amplicon-analysis/taxonomy-0.9_0.8_PR2-5.0.0/tax-output
