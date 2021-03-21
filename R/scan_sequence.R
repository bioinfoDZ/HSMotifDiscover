scan_sequence <- function(PSSM, seq, CharData){
  motifLength=ncol(PSSM)
  SubScore=c()
  for(i in 1:(nchar(seq)-motifLength+1)){
    SubSeq_char=strsplit(subseq(seq, i, i+motifLength-1), '')[[1]]
    NA_char=setdiff(SubSeq_char,row.names(PSSM))
    SubSeq_char=SubSeq_char[!SubSeq_char %in% NA_char]
    
    weight_vector=rep(1,length(SubSeq_char)) # < new
    
    u=CharData$UniqueCharGrps
    
    if(length(u)>1){
      for(z in 1:length(u)){
        weight_vector[SubSeq_char %in% CharData$mapped_char[which(CharData$charGrps==u[z])]]= 1/CharData$proportion[as.character(u[z])]# < new
      }
    }
    
    SubScore[i]=sum(log(diag(PSSM[SubSeq_char,])*weight_vector+.Machine$double.eps)) # < new
    
    #SubScore[i]=sum(log(diag(PSSM[SubSeq_char,])+.Machine$double.eps))
  }
  list(seq=seq,seq_score=SubScore,max_score=max(SubScore), max_loc=which(SubScore==max(SubScore)))
}
