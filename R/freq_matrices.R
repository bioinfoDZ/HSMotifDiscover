freq_matrices <- function(remaining_seqs,motif_width,remaining_mystart, remaining_weight, ND_backGround_prob_mat, unique_letters)
{
  
  #remaining_mystart=sapply(width(remaining_seqs)-motif_width, sample, 1)
  #remaining_mystart=motif_loc
  myend=remaining_mystart+(motif_width-1)
  rand_motif_seqs=subseq(remaining_seqs, start=remaining_mystart, end=myend)
  
  # if(sum(ND_backGround_prob_mat)!=0)
  # {
  #   remaining_weight=(remaining_weight-min(remaining_weight))/(max(remaining_weight)-min(remaining_weight))# normalise weight between 0 and 1
  # }
  
    #freq_mat=sapply(seq(1:motif_width),function(x) (colSums(letterFrequency(subseq(rand_motif_seqs, start=x, width=1),letters=unique_letters)))) # library('Biostrings')
    #prob_mat=sapply(seq(1:motif_width),function(x) (colSums(letterFrequency(subseq(rand_motif_seqs, start=x, width=1),letters=unique_letters))+0.5)/(length(remaining_seqs)+0.5*length(unique_letters))) # library('Biostrings')
    #freq_mat_=sapply(seq(1:motif_width),function(x) colSums(letterFrequency(subseq(rand_motif_seqs, start=x, width=1),letters=unique_letters)))
  
  freq_mat=sapply(seq(1:motif_width),function(x) (colSums(sweep(letterFrequency(subseq(rand_motif_seqs, start=x, width=1),letters=unique_letters), MARGIN=1,remaining_weight,'*')))) # library('Biostrings')
    #m_=letterFrequency(subseq(rand_motif_seqs, start=1, width=1),letters=unique_letters)
    #m=sweep(letterFrequency(subseq(rand_motif_seqs, start=1, width=1),letters=unique_letters), MARGIN=1,remaining_weight,'*')
    #print(head(m_))
    #print(head(m))
  
  
  freq_mat=freq_mat+.Machine$double.eps
  #prob_mat=t(t(freq_mat)/colSums(freq_mat)) # OR
  prob_mat=sweep(freq_mat,2,colSums(freq_mat),'/')
    #print(freq_mat)
    #print(prob_mat)
  
    #readline(prompt="Press [enter] to continue")
  

  ## alternate ##
  # freq_mat=matrix(.5, nrow=length(unique_letters), ncol=motif_width)
  # rownames(freq_mat)=unique_letters
  # colnames(freq_mat)=seq(1:motif_width)
  # 
  # temp_freq_mat=consensusMatrix(rand_motif_seqs)+0.5
  # freq_mat[rownames(temp_freq_mat),]=temp_freq_mat
  # prob_mat=t(t(freq_mat)/colSums(freq_mat))
  ###
  
  backGround_freq_mat=(colSums(letterFrequency(remaining_seqs,letters=unique_letters))-rowSums(freq_mat))+.Machine$double.eps
  backGround_prob_mat=backGround_freq_mat/sum(backGround_freq_mat)
  backGround_prob_mat=backGround_prob_mat+.Machine$double.eps  # new line
  
  if(sum(ND_backGround_prob_mat)!=0)
  {
    backGround_prob_mat=(backGround_prob_mat+ND_backGround_prob_mat)/2
  }
  
  #kl_div=mean(sapply(seq(1,motif_width), function(x) sum(freq_mat[,x]*log2(prob_mat[,x]/backGround_prob_mat))))
  kl_div=mean(sapply(seq(1,motif_width), function(x)  sum(sapply(rownames(prob_mat), function(b) prob_mat[b,x]*log2(prob_mat[b,x]/backGround_prob_mat[b])))))
  #kl_div=mean(sapply(seq(1,motif_width), function(x)  sum(sapply(rownames(prob_mat), function(b) freq_mat[b,x]*log2(prob_mat[b,x]/backGround_prob_mat[b])))))
  
  #readline(prompt="Press [enter] to continue")
  
  return(list(motif_PSSM=prob_mat, bkg_mat=backGround_prob_mat, F_score=kl_div))
}


