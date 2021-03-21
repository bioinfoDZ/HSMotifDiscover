seqs_TableData <- function(input_HSseq_file, motif_data, seq_weight_file, weight_threshold )
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
  #table(nchar(seqs[des_ind]))
  scanRes=lapply(seqs2, function(x) scan_sequence(motif_data$PSSM, x, motif_data$CharData))
  #names(scanRes)=header2
  score_str=as.vector(do.call(rbind,lapply(seq(1,length(scanRes)), function(x) toString(paste0(scanRes[[x]]$seq_score, sep='')))))
  max_score=sapply(seq(1,length(scanRes)), function(x) scanRes[[x]]$max_score[[1]])
  max_loc=sapply(seq(1,length(scanRes)), function(x) scanRes[[x]]$max_loc[[1]])
  motif_present=sapply(header2, function(x) is.element(x, motif_data$motif_DesiredSeqsName))
  seqLength=nchar(seqs2)
  p_value=sapply(max_score, function(x) dnorm(x, motif_data$P_value_para[1], motif_data$P_value_para[2]))


  tableData=c()
  tableData$p_value=p_value
  tableData$header2=header2
  tableData$seqs2=seqs2
  tableData$motif_present=motif_present
  tableData$score_str=score_str
  tableData$MotifLoc=max_loc
  tableData$max_score=max_score
  tableData$seqLength=seqLength

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
