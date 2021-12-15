scan_Max_sequence <- function(PSSM, seq){
  motifLength=ncol(PSSM)
  SubScore=c()
  for(i in 1:(nchar(seq)-motifLength+1)){
    SubSeq_char=strsplit(subseq(seq, i, i+motifLength-1), '')[[1]]
    NA_char=setdiff(SubSeq_char,row.names(PSSM))
    SubSeq_char=SubSeq_char[!SubSeq_char %in% NA_char]
    
    weight_vector=rep(1,length(SubSeq_char)) # < new
    SubScore[i]=sum((diag(PSSM[SubSeq_char,])*weight_vector+.Machine$double.eps)) # < new

  }
  list(seq=seq,seq_score=SubScore,max_score=max(SubScore), max_loc=which(SubScore==max(SubScore)))
}
