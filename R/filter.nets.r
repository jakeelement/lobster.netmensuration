filter.nets = function(DS, x, probs=c(0.05, 0.95) ){
  
  if(DS == "doorspread.range")  {
    # doorspread sanity check
    i = which( (x > 90) | (x < 0) ) 
    if (length(i) > 0 )  x[i]=NA
    return(x)
  }
  
  if(DS == "wingspread.range") {  
  # wingspread sanity check
    i = which( (x > 25) | (x < 0) ) 
    if (length(i) > 0 ) x[i] = NA
    return(x)
  }
  
  if(DS == "opening.range") {
    # opening sanity check
    i = which( (x > 10) | (x < 0) ) 
    if (length(i) > 0 )  x[i] = NA
    return(x)
  }

  if(DS == "clearance.range") {
    # clearance sanity check
    i = which( (x > 5) | (x < 0) ) 
    if (length(i) > 0 ) x[i] = NA
    return(x)
  }
  
  if(DS == "depth.range") {  
  # depth sanity check
    i = which( (x > 750) | (x < 0) ) 
    if (length(i) > 0 )  x[i] = NA
    return(x)
  }
  
}
