mapToAlphabets <- function(seq, orginal_char, mapped_char)
{
    #seq='s00g0s00i2s00i2s00i2s00i2s00i2s00i2s00i2s00g0'
    library('stringr')
    seq2=strsplit(seq, split='')[[1]]
    new_str=rep(NA,length(seq2))
    orgCharLengthsVec=unique(nchar(orginal_char))
    
    
    for(orgCharLengths in orgCharLengthsVec)
    {
        xLength_filt=orginal_char[which(nchar(orginal_char)==orgCharLengths)]
        #print(orgCharLengths)
        #print(xLength_filt)
        xLenght_ind_list=str_locate_all(seq, xLength_filt)
        names(xLenght_ind_list)=xLength_filt
        new_str=alternate_string(xLenght_ind_list,new_str)
        
    }
    
    if(length(orgCharLengthsVec)>1)
    {
        new_final_str=rle(new_str)$values
    }else if(orgCharLengthsVec==1)
    {
        new_final_str=new_str
    }
    
    
    library(plyr)
    mapped_str=paste0(mapvalues(new_final_str, from=orginal_char, to=mapped_char, warn_missing = FALSE), collapse = "")
    return(mapped_str)
}
