calculate_MaxScoreInd_Of_SelectedSeq <- function(matrices_f, selected_seq_f, CharDataf, motif_width_f, StrandDir, if_RC=FALSE)
{
  #print(' in > calculate_MaxScoreInd_Of_SelectedSeq')

  PSSM=matrices_f$motif_PSSM
  bkg_mat=matrices_f$bkg_mat

  loc_ind_subSeq=seq(1:motif_width_f)
  colnames(PSSM)=loc_ind_subSeq

  score_vec=c()

  if(if_RC==TRUE){
    RC_selected_seq_f=reverseComplement(DNAStringSet(selected_seq_f)) # NN
    RC_score_vec=c() # NN
  }

  lim=width(selected_seq_f)- motif_width_f
  for(i in 1:lim)
  {
      sub_seq=subseq(selected_seq_f, start=i, width=motif_width_f)
      split_sub_seq=strsplit(as.character(sub_seq),"")[[1]]

      if(if_RC==TRUE){
        RC_sub_seq=subseq(RC_selected_seq_f, start=i, width=motif_width_f) # NN
        RC_split_sub_seq=strsplit(as.character(RC_sub_seq),"")[[1]] # NN
      }

      #print(split_sub_seq)
      weight_vector=rep(1,length(split_sub_seq)) # < new

      u=CharDataf$UniqueCharGrps

      if(length(u)>1){
        for(z in 1:length(u)){
          weight_vector[split_sub_seq %in% CharDataf$mapped_char[which(CharDataf$grp==u[z])]]= 1/CharDataf$proportion[as.character(u[z])]# < new
        }
      }

      score_subSeq=prod(diag(PSSM[split_sub_seq,loc_ind_subSeq])*weight_vector) # < new
      score_background_subseq=prod(bkg_mat[split_sub_seq]*weight_vector) # < new
      score_vec[i]=score_subSeq/score_background_subseq

      if(if_RC==TRUE){
        RC_score_subSeq=prod(diag(PSSM[RC_split_sub_seq,loc_ind_subSeq])*weight_vector) #  NN
        RC_score_background_subseq=prod(bkg_mat[RC_split_sub_seq]*weight_vector) #  NN
        RC_score_vec[i]=RC_score_subSeq/RC_score_background_subseq # NN
      }

      #score_subSeq=prod((diag(PSSM[split_sub_seq,loc_ind_subSeq]))) #<
      #score_background_subseq=prod((bkg_mat[split_sub_seq])) #<

      #BString(toString(sub_seq))[1]


  }

  #norm_score = (score_vec-min(score_vec))/(max(score_vec)-min(score_vec))
  norm_score = score_vec/sum(score_vec)
  if(if_RC==TRUE){
    #RC_norm_score = (RC_score_vec-min(RC_score_vec))/(max(RC_score_vec)-min(RC_score_vec))
    RC_norm_score=RC_score_vec/sum(RC_score_vec)
  }


  #max_ind=which(score_vec==max(score_vec))
  #print(norm_score)

  sample_ind = sample(seq(1:lim), size = 1, replace = TRUE, prob = norm_score)

  ret=list()
  if(if_RC==TRUE){
    RC_sample_ind = sample(seq(1:lim), size = 1, replace = TRUE, prob = RC_norm_score)
    normSum=sum(score_vec[sample_ind], RC_score_vec[RC_sample_ind])
    #strandType = sample(c(0,1), size=1, replace=FALSE, prob= c(score_vec[sample_ind]/normSum, RC_score_vec[RC_sample_ind]/normSum))
    strandType = sample(c(0,1), size=1, replace=FALSE, prob= c(score_vec[sample_ind]/normSum, RC_score_vec[RC_sample_ind]/normSum))

    if(strandType==0 & StrandDir==0){
      ret$sample_ind=sample_ind
      ret$RC=0
    }else if(strandType==1 & StrandDir==0){
      ret$sample_ind=RC_sample_ind
      ret$RC=1
    }else if(strandType==1 & StrandDir==1){
      ret$sample_ind=RC_sample_ind
      ret$RC=0
    }else if(strandType==0 & StrandDir==1){
      ret$sample_ind=RC_sample_ind
      ret$RC=1
    }
  }else{
    ret$sample_ind=sample_ind
  }


  return(ret)
}


