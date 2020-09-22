
bottom.contact.gating.variance = function( x, good, bcp ) {

  ## -----------------------------------
  ## A variance based-criterion for gating
  # compute SD in the area of interest and compare with a lagged process to
  # start from centre and move left and continue until sd of residuals begins to deviate sustantially
  # from centre to left

  aoi.range = range( which( good & is.finite(x$depth) ) )
  aoi.min = aoi.range[1]
  aoi.max = aoi.range[2]
  aoi.mid = floor( mean( aoi.range ) ) # approximate midpoint
  aoi = aoi.min:aoi.max

  # determine approximate midpoint of the bottom track
  depths.threshold = mean( x$depth[aoi], na.rm=TRUE ) ## some deep value to initiate search for modal depth
  aoi.depths = x$depth[aoi]
  depths = modes( aoi.depths [ which( aoi.depths > depths.threshold )]  )
  oo = which( x$depth > depths$lb2 & x$depth < depths$ub2 )


  if (length(oo)> 10 ) {
    bot.range = range( oo) ## first pass (fast/rough) estimate of likely range of bottom indices
    aoi.mid = floor( mean( range( oo)) )
  }

  aoi.sd = sd( x$depth, na.rm=TRUE )  ## SD
  buffer = floor(length(aoi)/10) # additional points to add beyond midpoint to seed initial SD estimates
  if(any((aoi.mid-buffer)<1)) buffer = aoi.mid-1 #keep the index positive
  duration = 0

  target.sd = aoi.sd * bcp$noisefilter.sd.multiplier
  Ndata0 = length(aoi)
  Ndata = 0

  # browser()

  while (Ndata < Ndata0 ) {
    for ( j0 in aoi.mid:aoi.min  ) {#  begin from centre to right
      sdtest = sd(( x$depth[ (aoi.mid + buffer):j0]), na.rm=T)
      if ( is.na(sdtest) ) next()
      if ( sdtest  >= target.sd ) break()
    }
    for ( j1 in aoi.mid: aoi.max ) {  #  begin from centre to right
       sdtest =  sd(( x$depth[ (aoi.mid - buffer):j1]), na.rm=T)
      if ( is.na(sdtest) ) next()
      if ( sdtest >= target.sd ) break()
    }
    Ndata = length(unique( c( aoi.mid:j1, j0:aoi.mid) ) )
  }

  return( list( bc0=x$timestamp[ j0 ], bc1=x$timestamp[ j1 ]) )
}


