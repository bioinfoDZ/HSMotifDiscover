#range_motif_data = function(input_HSseq_file,motif_length,out_dir, weight_name, affinity_threshold)
CompileResults = function(input_HSseq_file=input_HSseq_file, motif_length=motif_length, charGrpFile=charGrpFile, seq_weight_file=seq_weight_file, affinity_threshold=affinity_threshold, itr=itr, if_revCompStrand=if_revCompStrand)
{

  print(' in > CompileResults')


  if( missing(seq_weight_file))
  {
    motif_data=mainCode(input_HSseq_file=input_HSseq_file, motif_length=motif_length, characterGrpFile=charGrpFile, itr=itr, if_revCompStrand=if_revCompStrand)
    tableDataMotif=motif_data$seqsPval_TableData

  }else if( !missing(seq_weight_file)) {
    motif_data=mainCode(input_HSseq_file=input_HSseq_file, motif_length=motif_length,  characterGrpFile=charGrpFile, seq_weight_file=seq_weight_file, affinity_threshold=affinity_threshold+.Machine$double.eps,itr=itr, if_revCompStrand=if_revCompStrand)
    tableDataMotif=motif_data$seqsPval_TableData
  }



  if( missing(seq_weight_file))
  {
    max_count=length(motif_data$motif_DesiredSeqsName)
    IC_W=motif_data$IC

    df=data.frame(header=tableDataMotif$header2,seq=tableDataMotif$seqs2,seqLength=tableDataMotif$seqLength, max_score=tableDataMotif$max_score, P_value=tableDataMotif$p_value, MotifLoc=tableDataMotif$MotifLoc, score=tableDataMotif$score_str )

  }else if( !missing(seq_weight_file)) {
    weight_name=basename(file_path_sans_ext(seq_weight_file))

    max_count=length(motif_data$motif_DesiredSeqsName)
    IC_W=motif_data$IC

    df=data.frame(header=tableDataMotif$header2,seq=tableDataMotif$seqs2,seqLength=tableDataMotif$seqLength,seq_weight_g_th=(tableDataMotif$seq_weights2 > tableDataMotif$weight_threshold), max_score=tableDataMotif$max_score, P_value=tableDataMotif$p_value, MotifLoc=tableDataMotif$MotifLoc, seq_weights=tableDataMotif$seq_weights2, score=tableDataMotif$score_str )
    colnames(df)[which(colnames(df)=='seq_weight_g_th')]=paste0('seq_weight_g_th-',tableDataMotif$weight_threshold )
  }


  if(if_revCompStrand==TRUE){
    df=cbind(df,RC_strand=tableDataMotif$RC_strand)
  }

  motif_data$df=df


  # assign(ret_var,list())
  # ret_var[[1]]=motif_data
  # names(ret_var)=paste0('motifData_L',motif_length)
  # return(ret_var)
  motif_data$mapC_PSSM=NULL
  motif_data$mapChar_discoveredMotifs=NULL
  motif_data$mapC_Entropy=NULL
  motif_data$mapC_IC=NULL
  motif_data$mapC_Entropy =NULL
  motif_data$motif_Obj =NULL
  motif_data$seqsPval_TableData =NULL
  motif_data$mapC_MotifLogo=NULL
  motif_data$bkgProb_vec = motif_data$final_bkg_vec
  motif_data$final_bkg_vec= NULL
  motif_data$P_value_para = NULL
  motif_data$mapC_MotifEntropy  = NULL


  motif_data$MAX_motifLoc = NULL
  motif_data$MAX_motifScore = NULL

  motif_data$resultsTable_df = motif_data$df
  motif_data$df = NULL
  motif_data$motif_Pval=mean(motif_data$resultsTable_df$P_value)


  print(' out > CompileResults')


  return(list(motifData=motif_data))



}
