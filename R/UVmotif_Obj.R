UVmotif_Obj <- function(final_PSSM_matrix,final_bkg_vec,CharDataf )
{
  

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

    motifObject=create_motif(input=as.matrix(final_PSSM_matrixM), alphabet=paste0(CharDataf$mapped_char,collapse = ''), type = "PPM", name = "motif",
                 pseudocount = .Machine$double.eps, bkg=final_bkg_vec)
    
    motifObject@motif=as.matrix(final_PSSM_matrixM[rownames(motifObject@motif),])
    
    return(motifObject)
}
