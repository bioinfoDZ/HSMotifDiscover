mainCode <- function(input_HSseq_file=input_HSseq_file, motif_length=motif_length, characterGrpFile=characterGrpFile, seq_weight_file=seq_weight_file, affinity_threshold=affinity_threshold, itr=itr, if_revCompStrand=FALSE){

  print(' in > mainCode')
  if(missing(affinity_threshold) & !missing(seq_weight_file))
  {
    w=(read.table(file=seq_weight_file,header = TRUE))$weight
    affinity_threshold=summary(w)[5]
    print(paste0('aff th:',affinity_threshold))

  }

  if(missing(seq_weight_file) )
  {
    weight=FALSE
  }else if(!missing(seq_weight_file) | (!missing(seq_weight_file) & !missing(affinity_threshold)) ) {
    weight=TRUE
  }



  ### Read HS_seq_data ##

  #input_filenam=gsub(pattern = "\\.txt$", "", input_HSseq_file)
  HS_seq_data=readtext(input_HSseq_file)
  #HS_seq_data=readtext('HS_seqs_formatted.txt')
  HS_seqs=strsplit(HS_seq_data$text, '\n')[[1]]
  header=HS_seqs[seq(1,length(HS_seqs),2)]
  seqs=HS_seqs[seq(2,length(HS_seqs),2)]


  ### generate GREEK letters ##

  GreekLettersVec=greek_vector
  names(GreekLettersVec)=NULL
  GL=GreekLettersVec[seq(1,32)]
  GL=GL[-c(15,25,26,29,30,31)]
  ###### Important variables ###

  CharData=list()
  if(!missing(characterGrpFile))
  {
    #characterGrpFile='mappedChar.txt'
    map_df=read.table(file=characterGrpFile,  sep = '\t', header = TRUE)
    CharData$orginal_char=as.character(map_df$Chars)
    CharData$grp=map_df$CharGroup
    GrpFreq=table(CharData$grp)
    if(any(GrpFreq>26)==TRUE)
    {
      stop("A group can't have more than 26 characters")
    }

   letters_type=c('letters' , 'LETTERS', 'GL')
    mapChar=rep(NA, length(CharData$orginal_char))
    for(x in 1:length(GrpFreq)){
      ind=which(CharData$grp==as.numeric(names(GrpFreq[x])))
      mapChar[ind]= eval(parse(text=letters_type[x]))[1:length(ind)]
      rm(ind)
    }

    CharData$mapped_char=as.character(mapChar)
    CharData$UniqueCharGrps=unique(CharData$grp)
    CharData$alphabets=paste0(CharData$mapped_char,collapse='')
    CharData$proportion=table(CharData$grp)/min(table(CharData$grp))
  } else if(missing(characterGrpFile))
  {
    oneSeq=paste(seqs, collapse = '')
    uniqchars <- function(x) unique(strsplit(x, "")[[1]])
    CharData$orginal_char=sort(uniqchars(oneSeq))
    CharData$mapped_char=CharData$orginal_char
    CharData$grp=rep(1,4)
    CharData$UniqueCharGrps=unique(CharData$grp)
    CharData$alphabets=paste0(CharData$mapped_char,collapse='')
    CharData$proportion=table(CharData$grp)/min(table(CharData$grp))
  }


  ###

  ### map HS_seq data to alphabet and write to fasta formatted file ###

  #mapped_seqs=lapply(seqs, mapToAlphabets, first_filt, second_filt, orginal_char, mapped_char)
  mapped_seqs=lapply(seqs, mapToAlphabets, CharData$orginal_char, CharData$mapped_char)



  #mapped_seqs=lapply(mapped_seqs, function(x) add_pad_seq(x))  # source('add_pad_seq.R')

  intermediate_mapped_file=paste0(format(Sys.time(), "%d%b%Y%H%M%S"),'_ML_',motif_length,'.txt')

  if (file.exists(intermediate_mapped_file))
    file.remove(intermediate_mapped_file)

  sink(intermediate_mapped_file)
  for(i in 1:length(mapped_seqs))
  {
    cat(paste0(header[i],'\n'))
    cat(paste0(mapped_seqs[i],'\n'))
  }
  sink()

  ### Read mapped seq data ##
  HSseqAlpha_data=readtext(intermediate_mapped_file)
  #HSseqAlpha_data=readtext('simulated_seqs.txt')
  HS_Alpha_seqs=strsplit(HSseqAlpha_data$text, '\n')[[1]]
  rm(HSseqAlpha_data)

  header=HS_Alpha_seqs[seq(1,length(HS_Alpha_seqs),2)]
  Alpha_seqs=HS_Alpha_seqs[seq(2,length(HS_Alpha_seqs),2)]
  names(Alpha_seqs)=header

  rm(HS_Alpha_seqs)

  unlink(intermediate_mapped_file)


  all_seq_BString=BStringSet(Alpha_seqs)  # library('Biostrings')
  motif_width=motif_length

  #weight=TRUE
  if(weight==TRUE)
  {
    weight_df=read.table(file=seq_weight_file, header = TRUE)
    rownames(weight_df)=weight_df$SeqID
    weight_vec=weight_df[names(all_seq_BString),'weight']

    rm_ind=which( width(all_seq_BString) <= motif_width )
    if(length(rm_ind)>0)
    {
      weight_vec=weight_df[-rm_ind,'weight']
      longSeq_BString=all_seq_BString[-rm_ind]
    }
    else{
      longSeq_BString=all_seq_BString
    }
    if(floor(affinity_threshold) > 0)
    {
      des_ind=which(weight_vec > affinity_threshold )
      #print(des_ind)
      des_Alpha_seqs=longSeq_BString[des_ind]
      seq_BString=BStringSet(des_Alpha_seqs)  # library('Biostrings')
      weight_vec=data.frame(weight_vec)
      weight_vec2=data.frame(weight_vec[des_ind,]) # modified
      #print(weight_vec2)
      non_des_Alpha_seqs=all_seq_BString[-des_ind]

    } else if(floor(affinity_threshold) == 0){
      #print(paste0('AT:',affinity_threshold))

      des_Alpha_seqs=longSeq_BString
      seq_BString=BStringSet(des_Alpha_seqs)
      weight_vec2=data.frame(weight_vec)
      non_des_Alpha_seqs=all_seq_BString[1]   # artificially put a sequence to get rid of nul background sequence
    }
    colnames(weight_vec2)='V1'
    #print(weight_vec2)


  }else if(weight==FALSE){
    weight_vec=data.frame(rep(1,length(Alpha_seqs)))   # create a weight vector of equal weights
    colnames(weight_vec)='V1'

    des_ind=which( width(all_seq_BString) > motif_width )

    des_Alpha_seqs=all_seq_BString[des_ind]
    seq_BString=BStringSet(des_Alpha_seqs)  # library('Biostrings')

    weight_vec2=data.frame(weight_vec[des_ind,])
    colnames(weight_vec2)='V1'

    non_des_Alpha_seqs=NULL
  }

  rm(Alpha_seqs)
  #print(weight_vec)
  #motifs_seqs=subseq(seq_BString,start= motif_loc, width=7)

  unique_letters=uniqueLetters(seq_BString) # # library('Biostrings')

  if( !is.null(non_des_Alpha_seqs) & weight==TRUE)
  {
    ND_seq_BString=BStringSet(non_des_Alpha_seqs)
    #print(ND_seq_BString)
    #readline(prompt="Press [enter] to continue")
    #weight_vec2=(weight_vec2-min(weight_vec2))/(max(weight_vec2)-min(weight_vec2))# normalise weight between 0 and 1
    weight_vec2=weight_vec2/sum(weight_vec2)# normalise weight between 0 and 1

    print(seq_weight_file)
    print(length(seq_BString))

    motif_data=gibbs_sampling(seq_BString=seq_BString, motif_width=motif_width, des_Alpha_seqs=des_Alpha_seqs, weight_vec2=weight_vec2, unique_letters=unique_letters, CharData=CharData, ND_seq_BString=ND_seq_BString, itr=itr, if_revCompStrand=if_revCompStrand)

    #}else if(length(non_des_Alpha_seqs)==0){
  }else if(weight==FALSE){

    #motif_data=gibbs_sampling(seq_BString=seq_BString, motif_width=motif_width, des_Alpha_seqs=des_Alpha_seqs, weight_vec2=weight_vec2, unique_letters=unique_letters, CharData=CharData, itr=itr, if_revCompStrand=if_revCompStrand)
    motif_data=gibbs_sampling(seq_BString=seq_BString, motif_width=motif_width, des_Alpha_seqs=des_Alpha_seqs, unique_letters=unique_letters, CharData=CharData, itr=itr, if_revCompStrand=if_revCompStrand)

    }

  motif_data=summary_motif_data(motif_data=motif_data, if_revCompStrand=if_revCompStrand)

  if(weight==TRUE)
  {
    motif_data$seqsPval_TableData=seqs_TableData(input_HSseq_file=input_HSseq_file, motif_data=motif_data, seq_weight_file=seq_weight_file, weight_threshold=affinity_threshold, if_revCompStrand=if_revCompStrand)
  }else if( weight==FALSE){
    motif_data$seqsPval_TableData=seqs_TableData(input_HSseq_file=input_HSseq_file, motif_data=motif_data, if_revCompStrand=if_revCompStrand)
  }

  print('out > mainCode')
  return(motif_data)
}
