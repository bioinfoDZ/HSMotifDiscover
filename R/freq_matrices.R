freq_matrices <- function(remaining_seqs_f,motif_width_f,remaining_mystart, remaining_weight, ND_backGround_prob_mat_f, unique_letters_f, StrandDir, if_RC=FALSE)
{
  
  myend=remaining_mystart+(motif_width_f-1)
  rand_motif_seqs=subseq(remaining_seqs_f, start=remaining_mystart, end=myend)
  
  if(if_RC==TRUE){
      RC_ind=which(StrandDir==1)
      if(length(RC_ind>0)){
        rand_motif_seqs=DNAStringSet(rand_motif_seqs)
        rand_motif_seqs[RC_ind]=reverseComplement(rand_motif_seqs[RC_ind])
      }
  }
  
  
   if(missing(remaining_weight)){
     freq_mat=consensusMatrix(rand_motif_seqs)[unique_letters_f,]
   }else{
    freq_mat=sapply(seq(1:motif_width_f),function(x) (colSums(sweep(letterFrequency(subseq(rand_motif_seqs, start=x, width=1),letters=unique_letters_f), MARGIN=1,remaining_weight,'*')))) # library('Biostrings')
   }
  
  
  freq_mat=freq_mat+.Machine$double.eps
  #prob_mat=t(t(freq_mat)/colSums(freq_mat)) # OR
  prob_mat=sweep(freq_mat,2,colSums(freq_mat),'/')
  
  #readline(prompt="Press [enter] to continue")
  
  
  ## alternate ##
  # freq_mat=matrix(.5, nrow=length(unique_letters_f), ncol=motif_width_f)
  # rownames(freq_mat)=unique_letters_f
  # colnames(freq_mat)=seq(1:motif_width_f)
  # 
  # temp_freq_mat=consensusMatrix(rand_motif_seqs)+0.5
  # freq_mat[rownames(temp_freq_mat),]=temp_freq_mat
  # prob_mat=t(t(freq_mat)/colSums(freq_mat))
  ###
  
  backGround_freq_mat=(colSums(letterFrequency(remaining_seqs_f,letters=unique_letters_f))-rowSums(freq_mat))+.Machine$double.eps
  backGround_prob_mat=backGround_freq_mat/sum(backGround_freq_mat)
  backGround_prob_mat=backGround_prob_mat+.Machine$double.eps  # new line
  
  if(sum(ND_backGround_prob_mat_f)!=0)
  {
    backGround_prob_mat=(backGround_prob_mat+ND_backGround_prob_mat_f)/2
  }
  
  #kl_div=mean(sapply(seq(1,motif_width_f), function(x) sum(freq_mat[,x]*log2(prob_mat[,x]/backGround_prob_mat))))
  kl_div=mean(sapply(seq(1,motif_width_f), function(x)  sum(sapply(rownames(prob_mat), function(b) prob_mat[b,x]*log2(prob_mat[b,x]/backGround_prob_mat[b]))))) # changed
  ##kl_div=mean(sapply(seq(1,motif_width_f), function(x)  sum(sapply(rownames(prob_mat), function(b) freq_mat[b,x]*log2(prob_mat[b,x]/backGround_prob_mat[b])))))
  
  #readline(prompt="Press [enter] to continue")
  
  return(list(motif_PSSM=prob_mat, bkg_mat=backGround_prob_mat, F_score=kl_div))
}


