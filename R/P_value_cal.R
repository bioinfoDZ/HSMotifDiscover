P_value_cal <- function(final_PSSM_matrix,final_bkg_vec,motif_width, CharDataf )
{
  print(' in > P_value_cal')

    batch_MeanScore_vector=c()
    no_of_batches=6000
    batch_size=40

    final_PSSM_DF=as.data.frame.matrix(final_PSSM_matrix)
    # add missing alphabets in the matrix
    NA_alpha=setdiff(CharDataf$mapped_char,rownames(final_PSSM_DF))
    NA_df=as.data.frame(matrix(0,nrow = length(NA_alpha), ncol = ncol(final_PSSM_DF)))
    rownames(NA_df)=NA_alpha
    colnames(NA_df)=colnames(final_PSSM_DF)
    final_PSSM_matrixM=rbind(final_PSSM_DF,NA_df)
    rm(final_PSSM_matrix)
    rm(final_PSSM_DF)


    NA_alpha2=setdiff(CharDataf$mapped_char,names(final_bkg_vec))
    NA_vec=rep(0,length(NA_alpha2))
    names(NA_vec)=NA_alpha2
    final_bkg_vec=c(final_bkg_vec,NA_vec)
    #print(final_bkg_vec)

   # rand seq ##
    split_sub_seq=rep(NA, motif_width)

    loc_ind_subSeq=colnames(final_PSSM_matrixM)

    for(l in 1:no_of_batches)
    {
      RandMotifSeq_score_vec=c()
      for(m in 1:batch_size)
      {
        #split_sub_seq=sample(c(second_name,first_name), size = motif_width, prob=final_bkg_vec, replace = TRUE)
        u=CharDataf$UniqueCharGrps

        weight_vector=rep(1,length(split_sub_seq)) # < new
        desChars=CharDataf$mapped_char[which(CharDataf$UniqueCharGrps==u[1])]
        split_sub_seq=sample(desChars, size = motif_width, prob=final_bkg_vec[names(final_bkg_vec) %in% desChars], replace = TRUE)

        if(length(u)>1){
          for(w in 2:length(u))
          {
            desChars2=CharDataf$mapped_char[which(CharDataf$UniqueCharGrps==u[w])]
            other_ind=which(colSums(final_PSSM_matrixM[which(rownames(final_PSSM_matrixM) %in% desChars2),])> .90) # <
            split_sub_seq[other_ind]=sample(desChars2, size = length(other_ind), prob=final_bkg_vec[names(final_bkg_vec) %in% desChars2], replace = TRUE)
            weight_vector[other_ind]=1/CharDataf$proportion[as.character(u[w])]
          }
        }

        ##RandMotifSeq_score_vec[m]=prod((diag(as.matrix(final_PSSM_matrixM[split_sub_seq,loc_ind_subSeq])+.Machine$double.eps)))
        #RandMotifSeq_score_vec[m]=sum(log(diag(as.matrix(final_PSSM_matrixM[split_sub_seq,loc_ind_subSeq]))+.Machine$double.eps))
        RandMotifSeq_score_vec[m]=sum(log(diag(as.matrix(final_PSSM_matrixM[split_sub_seq,loc_ind_subSeq])*weight_vector)+.Machine$double.eps))

      }
      batch_MeanScore_vector[l]=mean(RandMotifSeq_score_vec)
    }

    x=batch_MeanScore_vector
    #normalized_x = (x-min(x))/(max(x)-min(x))
    #hist(normalized_x, xlim=c(0,1), breaks=seq(0,1,by=1/10))

    mini=min(x)
    maxi=max(x)
    #hist(x)


    fit <- fitdistr(x, "normal")
    #class(fit)
    para <- fit$estimate
    #hist(x, prob = TRUE)
    #curve(dnorm(x, para[1], para[2]), col = 2, add = TRUE)

    #discovered_motifs=as.vector(subseq(seq_BString,start= mystart, width=motif_width))
    #splitted_discovered_motifs=lapply(discovered_motifs, function(x) strsplit(as.character(x),"")[[1]])
    #motifs_score=sapply(splitted_discovered_motifs, function(x) sum(log(diag(final_PSSM_matrixM[x,loc_ind_subSeq]))))

    #E_value=sapply(motifs_score, function(x) dnorm(x, para[1], para[2]))

    print(' out > P_value_cal')

    return(para)
}
