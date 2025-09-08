# Mask Positions based on a guideline (https://virological.org/t/issues-with-sars-cov-2-sequencing-data/473)
library(ape)
library(seqinr)
final_seq <- read.fasta("aligned_verdiref_seq.fasta", as.string = FALSE, seqtype = "DNA")

# Define homoplasic positions to mask
homoplasic_positions <- c(187, 1059, 2094, 3037, 3130, 6990, 8022, 10323, 10741, 11074, 13408, 14786, 19684, 20148, 21137, 24034, 24378, 25563, 26144, 26461, 26681, 28077, 28826, 28854, 29700)

mask_extremities_and_homoplasies <- function(seq,
                                             start1 = 1, end1 = 55,
                                             start2 = 29804, end2 = 29919,
                                             homoplasic_positions = NULL) {
  seq <- as.character(seq)
  
  # Masquer extrémités
  seq[start1:end1] <- "n"
  seq[start2:end2] <- "n"
  
  # Masquer positions homoplasiques si elles existent dans la séquence
  valid_pos <- homoplasic_positions[homoplasic_positions <= length(seq)]
  seq[valid_pos] <- "n"
  
  return(seq)
}


# Mask HOMOPLASIC SITES 


# Apply masking
masked_alignment <- lapply(
  final_seq,
  mask_extremities_and_homoplasies,
  homoplasic_positions = homoplasic_positions
)

write.fasta(
  sequences = masked_alignment,
  names = names(final_seq),
  file.out = "filtered_aligned_verdiref_seq.fasta"
)