#!/bin/sh

cd /mnt/d/VERDI/Mathis

#cat *.fa > combined_consensus.fasta

clustalo -i combined_consensus.fasta -o sequences_aligned.fasta 

# Compute the IQTREE phylogenic tree
iqtree -s sequences_aligned.fasta 

# Compute the genetic distance matrix
distmat -sequence sequences_aligned.fasta -outfile distance_matrix.txt -nucmethod 2

# compute a date phylogeny
treetime --tree iqtree/output.treefile --aln clustalo_sequences_aligned.fasta --dates treetime_data.csv 


# ------------- WITH NEXSTRAIN ALIGNMENT --------------------

# https://iqtree.github.io/doc/Dating
iqtree -s 122_nextclade_aligned.fasta --date iqtree_data.txt -pre results/122_nextclade


# compute a date phylogeny
treetime --tree iqtree/output.treefile --aln 122_nextclade_aligned.fasta --dates treetime_data.csv


# ----------------------------------- ROOTED AND WITH BOOTSTRAP ----------------
cd /mnt/d/VERDI/Mathis/122_sample
iqtree -s 122_nextclade_aligned.fasta -o NC_045512.2  -bb 1000

treetime --tree bootstrap_iqtree/122_nextclade_aligned.fasta.contree --aln 122_nextclade_aligned.fasta --dates ../treetime_data.csv


# ---------------------------------- FLOUK ALIGNMENT ------------------------
cd /mnt/d/VERDI/Mathis/flouk_analysis
iqtree -s flouk_align_seq.fa -o Wuhan/Hu-1/2019  -bb 1000

treetime --tree iqtree/flouk_align_seq.fa.contree --aln flouk_align_seq.fa --dates ../treetime_data.csv --reroot Wuhan/Hu-1/2019


cd /mnt/d/VERDI/Mathis/final_seq
clustalo -i verdi_cat.fasta -o final_sequences_aligned.fasta 


cd /mnt/d/VERDI/analysed/verdi_01
cat *.fa > mahidol_seq.fa

cd /mnt/d/VERDI/Mathis/final_seq
mafft --auto final_seq.fasta > aligned_final_seq.fasta


# ------------------------------ MAHIDOL WORKFLOW --------------------
treetime --tree verdiGISAID/verdiGISAID_iqtree/lsd_iqtree.contree --aln data/aligned_verdiGISAIDref_seq.fasta --dates data/treetime_verdiGISAIDref_df.csv --reroot NC_045512.2

seqkit grep -v -f treetime_outliers.txt aligned_verdiGISAIDref_seq.fasta > no_outliers_aligned_verdiGISAIDref_seq.fasta

grep -v -F -f treetime_outliers.txt treetime_verdiGISAIDref_df.csv > no_outliers_treetime_verdiGISAIDref_df.csv


gotree prune \
  -i verdiGISAID/verdiGISAID_iqtree/lsd_iqtree.contree \
  -f data/treetime_outliers.txt \
  -o verdiGISAID/verdiGISAID_iqtree/no_outliers_tree.nwk


treetime --tree verdiGISAID/verdiGISAID_iqtree/no_outliers_tree.nwk --aln data/no_outliers_aligned_verdiGISAIDref_seq.fasta --dates data/no_outliers_treetime_verdiGISAIDref_df.csv --reroot NC_045512.2 --stochastic-resolve



# ------------------------ Try for the Molecular clock of Nexstrain ----------------
treetime --tree verdiGISAID/verdiGISAID_iqtree/no_outliers_tree.nwk --aln data/no_outliers_aligned_verdiGISAIDref_seq.fasta --dates data/no_outliers_treetime_verdiGISAIDref_df.csv --reroot NC_045512.2 --stochastic-resolve --clock-rate 0.001 --clock-std-dev 0.0005

treetime --tree verdiGISAID/verdiGISAID_iqtree/no_outliers_tree.nwk --aln data/no_outliers_aligned_verdiGISAIDref_seq.fasta --dates data/no_outliers_treetime_verdiGISAIDref_df.csv --reroot NC_045512.2 --stochastic-resolve --clock-rate 0.001 --clock-std-dev 0.001


# ------------------------ Keep the polytomies ----------------
treetime --tree verdiGISAID/verdiGISAID_iqtree/no_outliers_tree.nwk --aln data/no_outliers_aligned_verdiGISAIDref_seq.fasta --dates data/no_outliers_treetime_verdiGISAIDref_df.csv --reroot NC_045512.2 --keep-polytomies --clock-rate 0.001 --clock-std-dev 0.0005


# WITH HT 0802 CORRECTED SAMPLE
treetime --tree verdiGISAID/verdiGISAID_iqtree/corrected_no_outliers_tree.nwk --aln data/corrected_no_outliers_aligned_verdiGISAIDref_seq.fasta --dates data/corrected_no_outliers_treetime_verdiGISAIDref_df.csv --reroot NC_045512.2 --stochastic-resolve --clock-rate 0.001 --clock-std-dev 0.0005
