
bottom.contact.maxdepth = function( sm, O, bcmethods, bcp ) {
  # using previously computed estimates of touchdown (bc0) and lift-off (bc1) and associated SD's
  # find maximum depths in the neighbourhoods

  res = list( bc0=NA, bc1=NA)

  names( sm) = c("Z", "timestamp", "ts" ) # Z is used to abstract the variable name such that any variable can be passed
  N = nrow(sm)

  # must be careful as timestamp is being converted to text and will lose tzone ... need to reset to correct timezone:
  bcm0 = paste(bcmethods, "0", sep="")
  bcm1 = paste(bcmethods, "1", sep="")

  tt0 =   tt1 =  t( as.data.frame(O[ bcm1 ]) )

  tzone = lubridate::tz(sm$timestamp[1])
  # reconvert to UTC as time zone is lost with the transpose
  tmp0 = lubridate::ymd_hms( t( as.data.frame(O[ bcm0 ]) ), tz=tzone )
  tmp1 = lubridate::ymd_hms( t( as.data.frame(O[ bcm1 ]) ), tz=tzone )

  bottom0.mean =  mean(tmp0, na.rm=TRUE)
  bottom1.mean =  mean(tmp1, na.rm=TRUE)

  bottom0.sd = max( 10, sd( tmp0, na.rm=TRUE  ), na.rm=TRUE ) # in seconds
  bottom1.sd = max( 10, sd( tmp1, na.rm=TRUE  ), na.rm=TRUE )

  bc0.range = bottom0.mean + bcp$maxdepth.sd.multiplier * bottom0.sd * c(-1, 1)
  bc1.range = bottom1.mean + bcp$maxdepth.sd.multiplier * bottom1.sd * c(-1, 1)

  bc0 = which( sm$timestamp > bc0.range[1] & sm$timestamp < bc0.range[2] )
  bc1 = which( sm$timestamp > bc1.range[1] & sm$timestamp < bc1.range[2] )
  zx0 = abs( sm$Z[bc0] - max( sm$Z[bc0] , na.rm=TRUE ) )
  zx1 = abs( sm$Z[bc1] - max( sm$Z[bc1] , na.rm=TRUE ) )

  zx0i = which( zx0 < bcp$eps.depth )
  zx1i = which( zx1 < bcp$eps.depth )

  if (length( zx0i) > 0 ) {
    z0 = min ( zx0i )
  } else {
    z0 = which.max( sm$Z[bc0] )
  }

  if (length( zx1i) > 0 ) {
    z1 = max ( zx1i )
  } else {
    z1 = which.max( sm$Z[bc1] )
  }

  bc0i = bc0[ z0 ] # choose left-most (earliest)
  bc1i = bc1[ z1 ] # choose right most (latest)

  res =  list( bc0=sm$timestamp[bc0i], bc1=sm$timestamp[bc1i] )

  return(res)
}


