
bottom.contact.incremental = function( sm, bcp=bcp )  {
  ## Algorithm: this is essentially a faster discrete version of the "smooth method"
  ## .. using an equally spaced x-interval, determine locations where change in slope is maximal
 stop( "Not working well enough for full time use")

  if(0) {
    load("~/ecomod/groundfish/data/nets/Scanmar/bottom.contact/results/bc.NED2014102.8.rdata")
    sm =data.frame( Z=bc$Z)
    sm$timestamp=bc$timestamp
    sm$ts=bc$ts
    good = bc$good
    sm$Z[ !good] = NA
  }

  res = c(NA, NA)

  # use only the subset with data for this step
  names( sm) = c("Z", "timestamp", "ts" ) # Z is used to abstract the variable name such that any variable can be passed
  N = nrow(sm)
  if ( N < 50 ) return(res)  # insufficient data

  # interpolate missing data
  sm$Z.cleaned = interpolate.xy.robust( sm[, c("ts", "Z")], method="sequential.linear" , probs=bcp$incremental.quants )

  sm$slope = NA
  sm$rsq = NA
  # n.target = n.target   # nrequired to enter into a linear regression model
  dta = -bcp$incremental.windowsize:bcp$incremental.windowsize
  for ( i in (bcp$incremental.windowsize+1):(N-bcp$incremental.windowsize) ) {
    j = i+dta
    lmres = lm( Z.cleaned ~ ts, data=sm[j,], na.action="na.omit" )
    # sm$rsq[i] = Rsquared(lmres)
    sm$slope[i] =  coefficients(lmres)["ts"]
  }

  sm$slope.cleaned = interpolate.xy.robust( sm[, c("ts", "slope")], method="sequential.linear" , probs=bcp$incremental.quants )

  slope.imax = which.max( sm$slope )  # descent has a positive slope start search from maximal slope
  slope.imin = which.min( sm$slope )  # ascending limb has a negative slope .. start search from min slope
  slope.mid = trunc( (slope.imin + slope.imax) /2 )
  slope.modes = modes( sm$slope )
  # sm$slope.smoothed = interpolate.xy.robust( sm[, c("ts", "slope")], method="inla" )

  if ( !is.finite(slope.modes$sd ) || slope.modes$sd < 1e-3) {
    return(res)
  }

  i0 = slope.imax
  for (i0 in slope.imax:slope.mid) {
    if (sm$slope[i0] < slope.modes$ub2 ) break()
  }

  i1 = slope.imin
  for (i1 in slope.imin:slope.mid) {
    if (sm$slope[i1] > slope.modes$lb2 ) break()
  }

  res =  c( sm$timestamp[i0], sm$timestamp[i1] )

  if (0) {
    plot(Z~ts,sm, pch=".")
    points(Z.cleaned~ts,sm, col="red", pch=20)
    abline( v=sm$ts[ c(i0, i1) ], col="blue" )

    x11()
    plot(slope~ts, sm, pch=20, col="green" )
    abline( h=slope.modes, col="red" )
    abline( v=sm$ts[ c(i0, i1) ], col="blue" )
  }

  sm$slope.smoothed = interpolate.xy.robust( sm[, c("ts", "slope.cleaned")], method="inla" , probs=bcp$incremental.quants, target.r2=0.9 )

  return(res)
}


