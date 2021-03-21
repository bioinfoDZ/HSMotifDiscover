gibbs_sampling <- function(seq_BString=seq_BString, motif_width=motif_width, des_Alpha_seqs=des_Alpha_seqs, weight_vec2=weight_vec2, unique_letters=unique_letters, CharData=CharData, ND_seq_BString=ND_seq_BString, itr=itr)
{
  
  if(!missing(ND_seq_BString)){
    ND_backGround_freq_mat=colSums(letterFrequency(ND_seq_BString,letters=unique_letters))
    ND_backGround_prob_mat=ND_backGround_freq_mat/sum(ND_backGround_freq_mat)
  }else if(missing(ND_seq_BString)){
    ND_backGround_prob_mat=rep(0,length(unique_letters))
    names(ND_backGround_prob_mat)=unique_letters
  }
    
    ## Gibbs samplig start ##
   
    
    mystart=sapply(width(seq_BString)-motif_width, sample, 1)
    mystart_initial=mystart
    MAX_motifScore=0
    MAX_motifLoc=0
    seq_lengths=width(seq_BString)
    motifScore_rep_count=0
    motifScore_circular_variable=seq(1,100)
    #while(sd(motifScore_circular_variable) > .2)
    for(i in seq(1:itr))
    {
      removed_seq_ind=sample(seq(1:length(seq_BString)),1)
      remaining_seqs=seq_BString[-removed_seq_ind]
      selected_seq=seq_BString[removed_seq_ind]
      matrices=freq_matrices(remaining_seqs, motif_width , mystart[-removed_seq_ind], weight_vec2[-removed_seq_ind,], ND_backGround_prob_mat, unique_letters)
      MaxScoreInd=calculate_MaxScoreInd_Of_SelectedSeq(matrices,selected_seq,CharData, motif_width )
      mystart[removed_seq_ind]=MaxScoreInd
      #print(err)
      motifScore=matrices$F_score
      #print(paste0(as.character(i), ' ',as.character(motifScore_rep_count),' ',as.character(motifScore),' ', as.character(sd(motifScore_circular_variable))))
      motifScore_circular_variable=c(motifScore, motifScore_circular_variable)[1:length(motifScore_circular_variable)]
      if(motifScore>MAX_motifScore){
        MAX_motifScore=motifScore
        MAX_motifLoc=mystart
      }
      ## randomly move discovered motif left or right by the random length after certain interval ##
      motifScore_rep_count=motifScore_rep_count+1
      
      if(motifScore_rep_count>3000)
      {
        print(paste0('====',motifScore_rep_count,'-', motif_width))
        shift_length <- c(2,1,1,2)
        tempStartLst=list()
        F_score=c()
        tempError=c()
        for(k in 1:length(shift_length))
        {
          tempStart=mystart
          if(k>(length(shift_length)/2))  # shift right
          {
            good_ind=which((tempStart+shift_length[k]+ motif_width) < seq_lengths)
            tempStart[good_ind] = tempStart[good_ind] + shift_length[k] ##
            tempMatrices=freq_matrices(remaining_seqs, motif_width , tempStart[-removed_seq_ind], weight_vec2[-removed_seq_ind,], ND_backGround_prob_mat, unique_letters)
            F_score[k]=tempMatrices$F_score
            tempStartLst[[k]]=tempStart
          }else if(k <= (length(shift_length)/2)) # shift left
          {
            good_ind=which(tempStart > shift_length[k])
            tempStart[good_ind] = tempStart[good_ind]-shift_length[k]
            tempMatrices=freq_matrices(remaining_seqs, motif_width , tempStart[-removed_seq_ind], weight_vec2[-removed_seq_ind,], ND_backGround_prob_mat, unique_letters)
            F_score[k]=tempMatrices$F_score
            tempStartLst[[k]]=tempStart
          }
        }
        req_pos_inds=which(F_score>motifScore)
        #shift_FScore_vec=F_score/sum(F_score)
        #random_shift=sample(seq(1:length(shift_FScore_vec)), size = 1, replace = TRUE, prob = shift_FScore_vec)
        #mystart=tempStartLst[[random_shift]]
        
        #req_pos_inds=seq(1:length(shift_length))
        if(length(req_pos_inds)>0)
        {
          sampled_pos_ind=sample(req_pos_inds,1)
          mystart=tempStartLst[[sampled_pos_ind]]
        } else
        {
          mystart=mystart
        }
        motifScore_rep_count=0
      }
    }
    ####
    #final_parameter_data=freq_matrices(seq_BString, motif_width , mystart, unname(unlist(weight_vec2)))
    #final_parameter_data=freq_matrices(seq_BString, motif_width , mystart, unname(unlist(weight_vec2)), ND_backGround_prob_mat, unique_letters)
    final_parameter_data=freq_matrices(seq_BString, motif_width , MAX_motifLoc, unname(unlist(weight_vec2)), ND_backGround_prob_mat, unique_letters)
    

    final_PSSM_matrix=final_parameter_data$motif_PSSM
    final_bkg_vec=final_parameter_data$bkg_mat
    
   
   
    discovered_motifs=as.vector(subseq(seq_BString,start= MAX_motifLoc, width=motif_width))
    
 
     b=lapply(discovered_motifs, function(x) paste0(mapvalues(strsplit(x, '')[[1]],from=CharData$mapped_char, to=CharData$orginal_char, warn_missing = FALSE), collapse=''))
     orgChar_discoveredMotifs=as.character(b)
  
    return(list(orgChar_discoveredMotifs=orgChar_discoveredMotifs,mapChar_discoveredMotifs=discovered_motifs, MAX_motifScore=MAX_motifScore, MAX_motifLoc=MAX_motifLoc, CharDataf=CharData, final_bkg_vec=final_bkg_vec))

}
