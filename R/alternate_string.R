alternate_string <- function(x_ind_list, new_str)
{
  for(i in 1:length(x_ind_list))
  {
    nam=names(x_ind_list[i])
    ind_data=x_ind_list[[i]]
    if(nrow(ind_data)>0)
    {
      for(j in 1:nrow(ind_data))
      {
        #print(ind_data[j,])
        new_str[ind_data[j,1]: ind_data[j,2]]=nam
      }
    }
    
  }
  return(new_str)
}