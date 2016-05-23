
bottom.contact.linear = function( sm, O, bcmethods, bcp )  {
  ## Smooth method: smooth data to remove local trends and compute first derivatives 

  res = list( bc0=NA, bc1=NA)

  if(0) {
    load("~/ecomod/groundfish/data/nets/Scanmar/bottom.contact/results/bc.NED2014102.8.rdata")
    sm =data.frame( Z=bc$Z)
    sm$timestamp=bc$timestamp
    sm$ts=bc$ts
    good = bc$good 
    sm$Z[ !good] = NA
    O = bc
  }

  # use only the subset with data for this step
  names( sm) = c("Z", "timestamp", "ts" ) # Z is used to abstract the variable name such that any variable can be passed

 
  # must be careful as timestamp is being converted to text and will lose tzone ... need to reset to correct timezone:
  bcm0 = paste(bcmethods, "0", sep="")
  bcm1 = paste(bcmethods, "1", sep="")

  tt0 =   tt1 =  t( as.data.frame(O[ bcm1 ]) )

  tzone = tz (sm$timestamp[1] )
  # recovert to UTC as time zone is lost with the transpose
  tmp0 = ymd_hms( t( as.data.frame(O[ bcm0 ]) ), tz=tzone )
  tmp1 = ymd_hms( t( as.data.frame(O[ bcm1 ]) ), tz=tzone )

  bottom0.mean =  mean(tmp0, na.rm=TRUE) 
  bottom1.mean =  mean(tmp1, na.rm=TRUE) 
  

  ## at least one solution required to continue  (2 means a valid start and end)
  ## ID best model based upon time .. furthest up a tail is best 
  
  sols = which( sm$timestamp >= bottom0.mean & sm$timestamp <= bottom1.mean )

  left =  min(sols)
  right = max(sols)

  if (0) {
    plot( Z~ts, sm, pch="." )
    abline( v= sm$ts[ c( left, right)], col="red" )  # initial breaks
  } 

  if ( any( !is.finite( c(left, right) ) ) ) return(res)
  
  # remember sm is is the bounds delimited by variance indices so 1=start of variance on left, N is end of variance right
  N = nrow(sm)  
  down = 1:left 
  up =  right:N
  dl = left
  dr = N-right
  # line segments are approximately same length as the up,down segments
  bot_left = c(left:(left+4*dl))
  bot_right = c( (right-4*dr):right)

  # find the row numbers which corresponds to each part of the curve, decent, bottom and ascent
  if(any(down < 1))  down = down [ -which(down < 1) ] 
  if(any(up > N))  up = up [ -which( up > N) ] 
  #  compute linear models for each section and their intersection

  #  right tail
  ri2 = right
  if (( length( which( is.finite( sm$Z[bot_right] ))) > 5) & (length( which( is.finite(sm$Z[up]))) > 5 )) {
    uplm2 = lm( Z ~ ts, sm[up, ], na.action = "na.exclude")
    botlm2r = lm( Z ~ ts, sm[bot_right, ], na.action = "na.exclude")
    cm <- rbind(coef(botlm2r),coef(uplm2)) # Coefficient matrix
    i2=c(-solve(cbind(cm[,2],-1)) %*% cm[,1])
    ii = which(sm$ts >i2[1]) 
    if (length(ii)>0) ri2 = min( ii )
  }

  # left tail
  ri1 =  left
  if (( length( which( is.finite( sm$Z[bot_left] ))) > 5 ) & ( length( which( is.finite( sm$Z[down]))) > 5 )) {
    downlm2 =lm( Z ~ ts, sm[down,], na.action = "na.exclude")
    botlm2l = lm( Z ~ ts, sm[bot_left, ], na.action = "na.exclude")
    # find where the sections intersect
    cm <- rbind(coef(botlm2l),coef(downlm2)) # Coefficient matrix
    i1=c(-solve(cbind(cm[,2],-1)) %*% cm[,1])
    ii = which(sm$ts < i1[1])
    if (length(ii)>0) { 
      ri1 = max(ii)
    } 
  } 
   
  res =  list( bc0=sm$timestamp[ri1], bc1=sm$timestamp[ri2] )

  return(res)
}


