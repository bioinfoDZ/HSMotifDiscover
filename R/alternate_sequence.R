alternate_sequence <- function(seq)
{

    alternate_string <- function(x_ind_list, new_str)
    {
      for(i in 1:length(x_ind_list))
      {
        name=names(x_ind_list[i])
        ind_data=x_ind_list[[i]]
        if(nrow(ind_data)>0)
        {
          for(j in 1:nrow(ind_data))
          {
            #print(ind_data[j,])
            new_str[ind_data[j,1]: ind_data[j,2]]=name
          }
        }
        
      }
      return(new_str)
    }
    
    
    first_filt=c('I0','G0','I2','G2')
    second_filt=c('A00','A06','A30','A36','H00','H06','H30','H36','S00','S06','S30','S36')
    orginal_char=c(first_filt, second_filt)
    
    first_name=letters[1:4]
    second_name=LETTERS[seq(5,5+11)]
    mapped_char=c(first_name, second_name)
    
    #seq='s00g0s00i2s00i2s00i2s00i2s00i2s00i2s00i2s00g0'
    library('stringr')
    seq2=strsplit(seq, split='')[[1]]
    new_str=rep(NA,length(seq2))
    
    first_ind_list=str_locate_all(seq, first_filt)
    names(first_ind_list)=first_filt
    new_str=alternate_string(first_ind_list,new_str)
    
    second_ind_list=str_locate_all(seq, second_filt)
    names(second_ind_list)=second_filt
    new_str=alternate_string(second_ind_list,new_str)
    new_final_str=rle(new_str)$values
    
    library(plyr)
    mapped_str=mapvalues(new_final_str, from=orginal_char, to=mapped_char, warn_missing = FALSE)
    return(mapped_str)

}
