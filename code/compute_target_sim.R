require("Biostrings")
require("seqinr")
require("Rcpp")

seqs = read.fasta(file = "../data/kiba_targets.fasta", seqtype = "AA")
seqs = read.fasta(file = "sequence.fasta.txt", seqtype = "AA")

self_aligned_score = rep(NA, length(seqs))

for (i in 1:length(seqs)){
  cat('aligning ',i,'/',length(seqs),' with itself..\n')
  seq = paste(toupper(getSequence(seqs[i][[1]])), sep="", collapse="")
  align_score = pairwiseAlignment(seq, seq, scoreOnly=TRUE, gapExtension=0.5, type="local", substitutionMatrix = "BLOSUM50")
  self_aligned_score[i] = align_score
}

# now compute the alignments for all pairs of targets.. note that I choose alignment parameters type="local" to perform SW alignment
# and gapExtension = 0.5, type="local", substitutionMatrix = "BLOSUM50" to get the same results as the target-target similarities here
# http://staff.cs.utu.fi/~aatapa/data/DrugTarget/

all_combinations = combn(length(seqs),2)

sim_mat = matrix(NA, nrow = length(seqs), ncol = length(seqs))

for (i in 1:ncol(all_combinations)){
  target_a = all_combinations[1,i]
  target_b = all_combinations[2,i]
  cat('computing similarity for ',target_a,', ',target_b,'(',i,'/',ncol(all_combinations),')\n')
  seq_1 = paste(toupper(getSequence(seqs[target_a][[1]])), sep="", collapse="")
  seq_2 = paste(toupper(getSequence(seqs[target_b][[1]])), sep="", collapse="")
  align_score_ab = pairwiseAlignment(seq_1, seq_2, scoreOnly=TRUE, gapExtension = 0.5, type="local", substitutionMatrix = "BLOSUM50")
  similarity = align_score_ab/(sqrt(self_aligned_score[target_a])*sqrt(self_aligned_score[target_b]))
  sim_mat[target_a,target_b] = similarity
  cat('similarity: ',similarity,'\n')
}

for (i in 1:nrow(sim_mat)){
  for (k in 1:ncol(sim_mat)){
    if (!is.na(sim_mat[i,k])){
      sim_mat[k,i] = sim_mat[i,k]
    }
    if (i==k){
      sim_mat[i,k] = 1
    }
  }
}
