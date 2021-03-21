calculate_MaxScoreInd_Of_SelectedSeq <- function(matrices, selected_seq, CharDataf, motif_width)
{
  PSSM=matrices$motif_PSSM
  bkg_mat=matrices$bkg_mat
  loc_ind_subSeq=seq(1:motif_width)
  colnames(PSSM)=loc_ind_subSeq
  # x=c(0,1,1.58,2)
  # y=c(0,2.585,3.16,3.585)
  # sp=spline(x, y)
  # upper_letters_col_ind=which(apply(PSSM[which(rownames(PSSM) %in% LETTERS), ], 2, function(x) sum(x))>0.8)
  # lower_letters_col_ind=which(apply(PSSM[which(rownames(PSSM) %in% letters), ], 2, function(x) sum(x))>0.8)
  
  #approx(x,y,xout=5.019802)
  
  score_vec=c()
  lim=width(selected_seq)- motif_width
  for(i in 1:lim)
      {
        sub_seq=subseq(selected_seq, start=i, width=motif_width)
        split_sub_seq=strsplit(as.character(sub_seq),"")[[1]]
        
        #print(split_sub_seq)
        weight_vector=rep(1,length(split_sub_seq)) # < new
        
        u=CharDataf$UniqueCharGrps
        
        if(length(u)>1){
          for(z in 1:length(u)){
            weight_vector[split_sub_seq %in% CharDataf$mapped_char[which(CharDataf$charGrps==u[z])]]= 1/CharDataf$proportion[as.character(u[z])]# < new
          }
        }
        
        score_subSeq=prod(diag(PSSM[split_sub_seq,loc_ind_subSeq])*weight_vector) # < new
        score_background_subseq=prod(bkg_mat[split_sub_seq]*weight_vector) # < new
        
        #score_subSeq=prod((diag(PSSM[split_sub_seq,loc_ind_subSeq]))) #<
        #score_background_subseq=prod((bkg_mat[split_sub_seq])) #<
        
        #BString(toString(sub_seq))[1]
        score_vec[i]=score_subSeq/score_background_subseq
      }
      
  #max_ind=which(score_vec==max(score_vec))
  sample_ind = sample(seq(1:lim), size = 1, replace = TRUE, prob = score_vec)
  
}