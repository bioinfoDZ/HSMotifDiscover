seqs_TableData <- function(input_HSseq_file, motif_data, seq_weight_file, weight_threshold, if_revCompStrand )
{
  ##### Write outputs to file ###
  print(' in > seqs_TableData')
  HS_seq_data=readtext(input_HSseq_file)
  #HS_seq_data=readtext('HS_seqs_formatted.txt')
  HS_seqs=strsplit(HS_seq_data$text, '\n')[[1]]
  header=HS_seqs[seq(1,length(HS_seqs),2)]
  seqs=HS_seqs[seq(2,length(HS_seqs),2)]
  des_ind=which(nchar(seqs)>= ncol(motif_data$PSSM))
  seqs2=seqs[des_ind]
  header2=header[des_ind]

  temp_PSSM=as.matrix(motif_data$PSSM)
  scanRes=lapply(seqs2, function(x) scan_Max_sequence(temp_PSSM, x))
  names(scanRes)=seq(1,length(seqs2))

  if(if_revCompStrand==TRUE){
    RC_strand=motif_data$RC_strand[des_ind]
    RC_ind=which(RC_strand==1)
    temp_RC_PSSM=reverseComplement(as.matrix(motif_data$PSSM))
    seqsRes_rc=lapply(seqs2[RC_ind], function(x) scan_Max_sequence(temp_RC_PSSM, x))
    names(seqsRes_rc)=RC_ind
    scanRes= modifyList(scanRes, seqsRes_rc[intersect(names(seqsRes_rc), names(scanRes))])
  }


  #names(scanRes)=header2
  score_str=as.vector(do.call(rbind,lapply(seq(1,length(scanRes)), function(x) toString(paste0(scanRes[[x]]$seq_score, sep='')))))
  max_score=sapply(seq(1,length(scanRes)), function(x) scanRes[[x]]$max_score[[1]])
  max_loc=sapply(seq(1,length(scanRes)), function(x) scanRes[[x]]$max_loc[[1]])
  motif_present=sapply(header2, function(x) is.element(x, motif_data$motif_DesiredSeqsName))
  seqLength=nchar(seqs2)
  #p_value=sapply(max_score, function(x) dnorm(x, motif_data$P_value_para[1], motif_data$P_value_para[2]))


  p_value=suppressWarnings({sapply(motif_data$mapChar_discoveredMotifs, function(x) motif_pvalue(motifs=motif_data$motif_Obj, score=score_match(motif_data$motif_Obj, x)) )})
  p_value=p_value[des_ind]



  tableData=c()
  tableData$p_value=p_value
  tableData$header2=header2
  tableData$seqs2=seqs2
  tableData$motif_present=motif_present
  tableData$score_str=score_str
  tableData$MotifLoc=max_loc
  tableData$max_score=max_score
  tableData$seqLength=seqLength

  if(if_revCompStrand==TRUE){
    tableData$RC_strand=motif_data$RC_strand
  }


  if(!missing(weight_threshold) & !missing(seq_weight_file))
  {

    weight_df=read.table(file=seq_weight_file,header = TRUE)
    rownames(weight_df)=weight_df$SeqID
    seq_weights2=weight_df[header2,2]

    tableData$seq_weights2=seq_weights2
    tableData$weight_threshold=weight_threshold
  }

  print(' out > seqs_TableData')
  return(tableData)
}
