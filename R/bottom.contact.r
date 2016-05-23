
bottom.contact = function( x, bcp, debugrun=FALSE ) {

  #require(lubridate)
  #require( numDeriv )

  if (debugrun) {
    debug.plot = TRUE
    browser()
  }

  debug.plot = FALSE
  n.min.required = 30
  nx = nrow(x)

  O = list()  # output list
  O$id = bcp$id
  O$error.flag = NA
  O$good = rep(TRUE, nx) # rows that will contain data that passes each step of data quality checks

  if (length (which( (!is.na(x$depth)))) < n.min.required ) return( NULL )

  ##--------------------------------
  # sort in case time is not in sequence
  # timestamps have frequencies higher than 1 sec .. duplciates are created and this can pose a problem
  x = x[order( x$timestamp ) ,]
  x$ts = as.numeric( difftime( x$timestamp, min(x$timestamp), units="secs" ) )

  O$plotdata = x   # save incoming data after creation of time stamps and ordering

  ## PRE-FILTER 1
  if ( exists( "double.depth.sensors", bcp) ) {
    # two depth sensors were used simultaneously but they are not calibrated!
    # remarkably hard to filter this out with any reliability at a later stage ..
    om = modes( x$depth )
    oaoi = range( which( x$depth <= om$ub2 & x$depth >= om$lb2))
    oj = oaoi[1]:oaoi[2]
    oi  = interpolate.xy.robust( x[oj, c("ts", "depth")], method="loess", trim=0.05, probs=c(0.05, 0.95),  target.r2=0.5 )
    oo = which( (x$depth[oj] - oi) < 0)
    x$depth[oj[oo]] = NA
    if (any( is.na( x$depth))) x$depth = interpolate.xy.robust( x[, c("ts", "depth")], method="simple.linear" )
  }

  ## PRE-FILTER 2
  if (all(!is.na(bcp$time.gate))) {
    # simple time-based gating if a time range is provided..
    bcp$trange.max = (bcp$tdif.max+5)
    O$good = bottom.contact.gating.time ( Zt=x$timestamp, good=O$good, bcp=bcp )
    x$depth[ !O$good ] = NA
  }

  ## PRE-FILTER 3
  # simple depth-based gating
  if( !any(x$depth > bcp$depth.min,na.rm=T)) return( NULL )
  O$good = bottom.contact.gating.depth ( Z=x$depth, good=O$good, bcp=bcp)
  x$depth[ !O$good ] = NA


  ## PRE-FILTER 4
  # time and depth-based gating
  mm = modes( x$depth ) # naive first estimate of location of most data depths
  if (mm$lb == mm$ub || is.na(mm$sd)) return(NULL)  # there is no depth variation in the data ... likely a bad data series
  mm.i = which( x$depth > (mm$lb2+bcp$depth.range[1]) & x$depth < (mm$ub2 + bcp$depth.range[2]) )
  O$good[ setdiff(1:nx, mm.i)] = FALSE
  O$good[mm.i] = TRUE

  ## RANGE CHECKS
  x.timerange = range( x$timestamp[O$good], na.rm=TRUE )
  x.dt = difftime( x.timerange[2], x.timerange[1], units="mins" )

  if ( x.dt < (bcp$tdif.min ) ) {
      O$error.flag = "Time track is too short?"
      return(O)
  }

  if ( x.dt > (bcp$tdif.max+10) ) {   # +10 = 5 min on each side
    # data vector is too long ... truncate
    mm.r = range( mm.i )
    mm.dt = difftime( x$timestamp[mm.r[2]], x$timestamp[mm.r[1]], units="mins" )
    if ( mm.dt > (bcp$tdif.max+10) ) {
      O$error.flag = "Time track is too long?"
    }
  }

  ## SANITY CHECKS
  # sometimes multiple tows exist in one track ...
  # over-smooth depth to capture strange tows
  x$dsm = interpolate.xy.robust( x[, c("ts", "depth")], method="sequential.linear", trim=0.05 )
  x$dsm = interpolate.xy.robust( x[, c("ts", "dsm")], method="local.variance", trim=0.05 )
  x$dsm = interpolate.xy.robust( x[, c("ts", "dsm")], method="sequential.linear", trim=0.05 )
  x$dsm = interpolate.xy.robust( x[, c("ts", "dsm")], method="loess", trim=0.05 )

  zz = rep( 0, nx )
  zz[ which(  x$dsm < ( mm$mode / 3 )) ] = 1  # flag shallow areas
  inc.depth = abs( diff( x$dsm ) )

  # also capture strong noise - very obviously wrong data
  rapid.depth.changes = which( inc.depth > bcp$maxdepthchange )
  if ( length( rapid.depth.changes ) > 0 ) {
    zz[ rapid.depth.changes ] = 1
    zz[ rapid.depth.changes-1 ] = 1  # include adjecent points to remove
    zz[ rapid.depth.changes+1 ] = 1
  }
  dzz = diff(zz)
  bnds = c(1, which( dzz != 0 ), nx )
  if ( length (bnds) > 2 ) {
    # i.e. , contaminated by noise or multiple tows
    segs = diff(bnds) # length of each segment
    longest = which.max( segs )
    gg = bnds[longest]:bnds[(longest+1)]
    bad = setdiff( 1:nx, gg)
    O$good[bad] = FALSE
  }

  # fix tails
  mm = modes( x$depth[ O$good ] ) # recompute after gating
  mm.bot = which ( x$depth >= mm$lb2 & x$depth <= mm$ub2 ) # first est of bottom
  mm.med = trunc( median( mm.bot, na.rm=TRUE ) )
  # fix tails when the tails are not a single smooth tail:
  llim = mm$mode + bcp$depth.range[1]
  for ( ll in mm.med:1 ) {
    if ( is.finite( x$depth[ll]  ) ) {
      if ( x$depth[ll] < llim ) break()
    }
  }
  for ( uu in mm.med:nrow(x)) {
    if ( is.finite( x$depth[uu]  ) ) {
      if ( x$depth[uu] < llim ) break()
    }
  }

  todrop = c( 1:ll, uu:nrow(x) )
  x$depth[ todrop ] = NA
  O$good[ todrop] = FALSE


  ## ------------------------------
  ## MAIN NOISE/INTERPOLATION FILTER
  # Some filtering of noise from data and further focus upon area of interest based upon time and depth if possible

  res = NULL
  res = try( bottom.contact.filter.noise ( x=x, good=O$good, bcp=bcp, debug=debugrun ), silent =TRUE )
  if ( "try-error" %in% class(res) ) {
    x$depth = jitter( x$depth )
    res = try( bottom.contact.filter.noise ( x=x, good=O$good, bcp=bcp, debug=debugrun ), silent =TRUE )
  }

  if ( !"try-error" %in% class(res) ) {
    if ( cor( x$depth, res$depth.smoothed, use="pairwise.complete.obs") > 0.999 ) {
      bcp$noisefilter.target.r2 = bcp$noisefilter.target.r2 - 0.1
      x$depth = jitter( x$depth )
      res = try( bottom.contact.filter.noise ( x=x, good=O$good, bcp=bcp, debug=debugrun ), silent =TRUE )
  }}

  x$depth.smoothed = x$depth
  if ( ! "try-error" %in% class( res) )  {
    # retain good/bad information only for the interval where there is data
    # .. this effectively does the filtering in the area of interest only and retains the tails
    # for further influence later
    ig = range( which( res$good ))
    igg = ig[1]:ig[2]
    O$good[igg] = res$good[igg]
    x$depth[ !O$good ] = NA
    x$depth.smoothed= res$depth.smoothed
  }

  if(sum(x$depth-min(x$depth,na.rm=T),na.rm=T)==0) return( NULL )
  if(sum(O$good)==0) return( NULL )

  x.timerange = range( x$timestamp[O$good], na.rm=TRUE )
  x.dt = difftime( x.timerange[2], x.timerange[1], units="mins" )

  if ( x.dt < ( bcp$tdif.min ) ) {
      O$error.flag = "Not enough data?"
      return(O)
  }

  if ( x.dt > (bcp$tdif.max + 10 ) ) {
      O$error.flag = "Too much data?"
      return(O)
  }

  ## FINAL DATA GATING
  # variance gating attempt after the data has been cleaned as much as possible
  O$variance.method0 = NA
  O$variance.method1 = NA
  O$variance.method.indices = NA
  res = NULL
  res = try( bottom.contact.gating.variance ( x, O$good, bcp ), silent =TRUE )

  if ( ! "try-error" %in% class( res) )  {
    if (exists( "bc0", res )) {
    if ( all(is.finite( c(res$bc0, res$bc1 )) ) ) {
      DT = abs( as.numeric( difftime( res$bc0, res$bc1, units="mins" ) ) )
      if ( length(DT) == 1 ) {
        if ( is.finite(DT) &&  DT > bcp$tdif.min & DT < bcp$tdif.max ) {
        O$variance.method0 = res$bc0
        O$variance.method1 = res$bc1
        O$variance.method.indices = which( x$timestamp >= res$bc0 &  x$timestamp <= res$bc1 )
        bad = which( x$timestamp < res$bc0 |  x$timestamp > res$bc1 )
        if (length( bad) > 0) O$good[ bad ] = FALSE
        x$depth[ !O$good ] = NA
      } }
    }}
  }

  if(debug.plot) {
    trange = range( x$ts[O$good], na.rm=TRUE )
    drange = c( quantile( x$depth, c(0.05, 0.975), na.rm=TRUE) , median( x$depth, na.rm=TRUE ) * 1.05 )
    plot(depth~ts, x, ylim=c(drange[2],drange[1]), xlim=c(trange[1],trange[2]), pch=20, cex=0.1, col="gray" )
    mcol = "gray"
    points( depth~ts, x[ O$variance.method.indices, ], pch=20, col=mcol, cex=0.2)
  }


  ## NOTE::: From this point on, O$good is now complete
  ## -- it contains indices of time and depth-based gating as well as the variance based gating


  # finalize selection of area of interest (based upon gating, above)
  aoi.range = range( which( O$good )  )
  aoi.mid = trunc( mean( aoi.range ) ) # approximate midpoint
  aoi.min = aoi.range[1]
  aoi.max = aoi.range[2]
  O$aoi = aoi.min:aoi.max  # stored for use in the linear method where we need to recover the left-most index


  ##--------------------------------
  # Modal method: gating by looking for modal distribution and estimating sd of the modal group in the data
  # first by removing small densities ( 1/(length(i)/nb)  ) and by varying the number of breaks in the histogram
  # until a target number of breaks, nbins with valid data are found
  O$modal.method0 = NA  #### NOTE:: using the 'c' operator on posix strips out the timezone info! this must be retained
  O$modal.method1 = NA  #### NOTE:: using the 'c' operator on posix strips out the timezone info! this must be retained
  O$modal.method.indices = NA

  sm0 = x[ O$aoi, c("depth.smoothed", "timestamp", "ts" ) ]  # send filtered data ... continuity not important .. order is important
  res = NULL
  res = try( bottom.contact.modal( sm=sm0, bcp ), silent=TRUE )
    if ( ! "try-error" %in% class( res) ) {
      if (exists( "bc0", res )) {
      if ( all(is.finite( c(res$bc0, res$bc1 )) ) ) {
        DT =  abs( as.numeric( difftime( res$bc0, res$bc1, units="mins" ) ) )
        if ( length(DT) == 1 ) {
          if ( is.finite(DT) &&  DT > bcp$tdif.min & DT < bcp$tdif.max ) {
          O$modal.method0 = res$bc0 #### NOTE:: using the 'c' operator on posix strips out the timezone info! this must be retained
          O$modal.method1 = res$bc1 #### NOTE:: using the 'c' operator on posix strips out the timezone info! this must be retained
          O$modal.method.indices = which( x$timestamp >= res$bc0  &  x$timestamp <= res$bc1  ) # x correct
        } }
      }}
    }

  if(debug.plot) {
    trange = range( x$ts[O$good], na.rm=TRUE )
    drange = c( quantile( x$depth, c(0.05, 0.975), na.rm=TRUE) , median( x$depth, na.rm=TRUE ) * 1.05 )
    plot(depth~ts, x, ylim=c(drange[2],drange[1]), xlim=c(trange[1],trange[2]), pch=20, cex=0.1, col="gray" )
    mcol = "green"
    points( depth~ts, x[ O$modal.method.indices, ], pch=20, col=mcol, cex=0.2)
  }



  ## ----------------------------
  ## Smooth method: using smoothed data (slopes are too unstable with raw data),
  ## compute first derivatives to determine when the slopes inflect

  O$smooth.method0 = NA #### NOTE:: using the 'c' operator on posix strips out the timezone info! this must be retained
  O$smooth.method1 = NA #### NOTE:: using the 'c' operator on posix strips out the timezone info! this must be retained
  O$smooth.method.indices = NA
  sm0 = x[ O$aoi, c("depth.smoothed", "timestamp", "ts")]  # Send all data within the aoi --- check this .. order is important

  res = NULL
  res = try(
    bottom.contact.smooth( sm=sm0, bcp=bcp ) , silent =TRUE)
    if ( ! "try-error" %in% class( res) ) {
      if (exists( "bc0", res )) {
      if ( all(is.finite( c(res$bc0, res$bc1 )) ) ) {
        DT =  abs( as.numeric( difftime( res$bc0, res$bc1, units="mins" ) ) )
        if ( length(DT) == 1) {
          if ( is.finite(DT) &&  DT > bcp$tdif.min & DT < bcp$tdif.max ) {
          O$smooth.method0 = res$bc0 #### NOTE:: using the 'c' operator on posix strips out the timezone info! this must be retained
          O$smooth.method1 = res$bc1 #### NOTE:: using the 'c' operator on posix strips out the timezone info! this must be retained
          O$smooth.method.indices = which( x$timestamp >= res$bc0 &  x$timestamp <= res$bc1 ) # x correct
        } }
      }}
    }

  if(debug.plot) {
    trange = range( x$ts[O$good], na.rm=TRUE )
    drange = c( quantile( x$depth, c(0.05, 0.975), na.rm=TRUE) , median( x$depth, na.rm=TRUE ) * 1.05 )
    plot(depth~ts, x, ylim=c(drange[2],drange[1]), xlim=c(trange[1],trange[2]), pch=20, cex=0.1, col="gray" )
    mcol = "blue"
    points( depth~ts, x[ O$smooth.method.indices, ], pch=20, col=mcol, cex=0.2)
  }


  ## ----------------------------
  ## maxdepth method: looking for the max depth near the areas of interest (left, right)
  O$maxdepth.method0 =  NA #### NOTE:: using the 'c' operator on posix strips out the timezone info! this must be retained
  O$maxdepth.method1 =  NA #### NOTE:: using the 'c' operator on posix strips out the timezone info! this must be retained
  O$maxdepth.method.indices = NA
  sm0 = x[, c("depth.smoothed", "timestamp", "ts")]  # Send all data within the aoi --- check this .. order is important
  sm0$depth[ !O$good ] = NA
  sm0 = sm0[ O$aoi, ]


  bcmethods=c( "smooth.method", "modal.method" )
  res = NULL
  res = try( bottom.contact.maxdepth( sm=sm0, O=O, bcmethods=bcmethods, bcp=bcp ) , silent=TRUE )

  if ( ! "try-error" %in% class( res) ) {
    if (exists( "bc0", res )) {
    if ( all(is.finite( c(res$bc0, res$bc1 )) ) ) {
      DT =  abs( as.numeric( difftime( res$bc0, res$bc1, units="mins" ) ) )
      if ( length(DT) == 1 ) {
        if ( is.finite(DT) && DT > bcp$tdif.min & DT < bcp$tdif.max ) {
        O$maxdepth.method0 = res$bc0 #### NOTE:: using the 'c' operator on posix strips out the timezone info! this must be retained
        O$maxdepth.method1 = res$bc1 #### NOTE:: using the 'c' operator on posix strips out the timezone info! this must be retained
        O$maxdepth.method.indices = which( x$timestamp >= res$bc0 &  x$timestamp <= res$bc1 ) # x correct
      }}
    }}
  }

  if(debug.plot) {
    trange = range( x$ts[O$good], na.rm=TRUE )
    drange = c( quantile( x$depth, c(0.05, 0.975), na.rm=TRUE) , median( x$depth, na.rm=TRUE ) * 1.05 )
    plot(depth~ts, x, ylim=c(drange[2],drange[1]), xlim=c(trange[1],trange[2]), pch=20, cex=0.1, col="gray" )
    mcol = "yellow"
    points( depth~ts, x[ O$maxdepth.method.indices, ], pch=20, col=mcol, cex=0.2)
  }



  ## ---------------------------
  ## Linear method: looking at the intersection of three lines (up, bot and down)

  O$linear.method0 = NA #### NOTE:: using the 'c' operator on posix strips out the timezone info! this must be retained
  O$linear.method1 = NA #### NOTE:: using the 'c' operator on posix strips out the timezone info! this must be retained
  O$linear.method.indices = NA
  sm0 = x[, c("depth.smoothed", "timestamp", "ts")]  # Send all data within the aoi --- check this .. order is important
  sm0$depth[ !O$good ] = NA
  sm0 = sm0[ O$aoi, ]

  bcmethods=c( "smooth.method", "modal.method", "linear.method" )

  res = NULL
  res = try( bottom.contact.linear( sm=sm0, O=O, bcmethods=bcmethods, bcp=bcp ) , silent=TRUE )

  if ( ! "try-error" %in% class( res) ) {
    if (exists( "bc0", res )) {
    if ( all(is.finite( c(res$bc0, res$bc1 )) ) ) {
      DT =  abs( as.numeric( difftime( res$bc0, res$bc1, units="mins" ) ) )
      if (  length(DT) == 1) {
          if ( is.finite(DT) && DT > bcp$tdif.min & DT < bcp$tdif.max ) {
        O$linear.method0 = res$bc0 #### NOTE:: using the 'c' operator on posix strips out the timezone info! this must be retained
        O$linear.method1 = res$bc1 #### NOTE:: using the 'c' operator on posix strips out the timezone info! this must be retained
        O$linear.method.indices = which( x$timestamp >= res$bc0 &  x$timestamp <= res$bc1 ) # x correct
      } }
    }}
  }


  if(debug.plot) {
    trange = range( x$ts[O$good], na.rm=TRUE )
    drange = c( quantile( x$depth, c(0.05, 0.975), na.rm=TRUE) , median( x$depth, na.rm=TRUE ) * 1.05 )
    plot(depth~ts, x, ylim=c(drange[2],drange[1]), xlim=c(trange[1],trange[2]), pch=20, cex=0.1, col="gray" )
    mcol = "red"
    points( depth~ts, x[ O$linear.method.indices, ], pch=20, col=mcol, cex=0.2)
  }


  ## ---------------------------

  O$manual.method0 = NA #### NOTE:: using the 'c' operator on posix strips out the timezone info! this must be retained
  O$manual.method1 = NA  ### NOTE:: using the 'c' operator on posix strips out the timezone info! this must be retained
  if ( bcp$user.interaction  ) {
    print( "Click with left mouse on start and stop locations and then press ESC or right click to continue.")
    trange = range( x$ts[O$good], na.rm=TRUE )
    drange = c( quantile( x$depth, c(0.05, 0.975), na.rm=TRUE) , median( x$depth, na.rm=TRUE ) * 1.05 )
    plot(depth~ts, x, ylim=c(drange[2],drange[1]), xlim=c(trange[1],trange[2]), pch=20, cex=0.1, col="gray" )
    points( depth~ts, x[O$good,], pch=20, col="orange", cex=0.2)

    useridentified = locator( n=2, type="o", col="cyan")
    u.ts0 = which.min( abs( x$ts-useridentified$x[1] ))
    u.ts1 = which.min( abs( x$ts-useridentified$x[2] ))
    O$manual.method0 = x$timestamp[u.ts0]
    O$manual.method1 = x$timestamp[u.ts1]
    O$manual.method.indices = which( x$timestamp >= O$manual.method0 &  x$timestamp <= O$manual.method1 )
  }


  O$means0 = NA  ### NOTE:: using the 'c' operator on posix strips out the timezone info! this must be retained
  O$means1 = NA  ### NOTE:: using the 'c' operator on posix strips out the timezone info! this must be retained

  bcmethods = c("manual.method", "variance.method", "smooth.method", "modal.method", "maxdepth.method", "linear.method", "means" )
  standard =  which( bcmethods=="manual.method") # gold standard
  direct = which( bcmethods %in%  c("smooth.method", "modal.method", "linear.method", "maxdepth.method" ) )
  imeans = which( bcmethods == "means" )


  # must be careful as timestamp is being converted to text and will lose tzone ... need to reset to correct timezone:
  bcm0 = paste(bcmethods, "0", sep="")
  bcm1 = paste(bcmethods, "1", sep="")

  tzone = tz (x$timestamp[1] )

  # recovert to time zone of incoming data as time zone is lost with the transpose
  tmp0 = ymd_hms( t( as.data.frame(O[ bcm0 ]) ), tz=tzone )
  tmp1 = ymd_hms( t( as.data.frame(O[ bcm1 ]) ), tz=tzone )

  bottom0.mean =  mean(tmp0, na.rm=TRUE)
  bottom1.mean =  mean(tmp1, na.rm=TRUE)

  bottom0.sd = max( 1, sd( tmp0, na.rm=TRUE  ) ) # in seconds
  bottom1.sd = max( 1, sd( tmp1, na.rm=TRUE  ) )

  dflag0 = rowSums(as.matrix(dist(tmp0[direct], upper=TRUE )), na.rm=TRUE)  # which are most extreme
  dflag1 = rowSums(as.matrix(dist(tmp1[direct], upper=TRUE )), na.rm=TRUE)  # which are most extreme

  tooextreme0 = which.max( dflag0 )
  tooextreme1 = which.max( dflag1 )

  trimmed0 = trimmed1 = direct
  if (length( which(dflag0 > 0) ) > 1 ) trimmed0 = trimmed0[-tooextreme0 ]
  if (length( which(dflag1 > 0) ) > 1 ) trimmed1 = trimmed1[-tooextreme1 ]

  tmp0[imeans] = mean( tmp0[ trimmed0 ], na.rm=TRUE )
  tmp1[imeans] = mean( tmp1[ trimmed1 ], na.rm=TRUE )

  O$bottom0 = NA
  O$bottom1 = NA
  O$bottom0.sd = NA
  O$bottom1.sd = NA
  O$bottom0.n = NA
  O$bottom1.n = NA
  O$bottom.diff = NA
  O$bottom.diff.sd = NA

  if ( any (is.na( c( tmp0[ standard ],  tmp1[ standard ]) ) ) ) {
    # no manual standard .. use mean as the standard
    O$bottom0 = tmp0[imeans]
    O$bottom1 = tmp1[imeans]
    O$bottom0.sd = sd(  ( tmp0[ trimmed0]), na.rm=TRUE ) # in secconds
    O$bottom1.sd = sd(  ( tmp1[ trimmed1]), na.rm=TRUE )
    O$bottom0.n = length( which( is.finite( tmp0[ trimmed0] )) )
    O$bottom1.n = length( which( is.finite( tmp1[ trimmed1] )) )
    O$bottom.diff =  difftime( O$bottom1, O$bottom0, units="secs" )
    O$bottom.diff.sd = sqrt(O$bottom0.sd ^2 + O$bottom0.sd^2) # sec
  } else {
    # over-ride all data and use the manually determined results
    O$bottom0 = tmp0[ standard ]
    O$bottom1 = tmp1[ standard ]
    O$bottom0.sd = NA
    O$bottom1.sd = NA
    O$bottom0.n = 1
    O$bottom1.n = 1
    O$bottom.diff = difftime( O$bottom1, O$bottom0, units="secs" )
  }

  tmp = data.frame( start=tmp0)
  tmp$end = tmp1

  tmp$diff = difftime( tmp[,"end"], tmp[,"start"], units="secs" )
  tmp$start.bias =  difftime( tmp[,"start"],  O$bottom0, units="secs" )
  tmp$end.bias   = difftime( tmp[,"end"],  O$bottom1, units="secs" )
  rownames(tmp) = bcmethods

  O$summary = tmp

  # finalised data which have been filtered
  fin.all = which( x$timestamp >= O$bottom0 & x$timestamp <= O$bottom1 )
  if (length( fin.all ) == 0 ) fin.all = min( which( O$good)) : max( which( O$good) )
  fin.good = intersect( fin.all, which( O$good)  )
  fin0 = min( fin.all, na.rm=TRUE)
  fin1 = max( fin.all, na.rm=TRUE)

  O$depth.mean = mean( x$depth[ fin.good ], na.rm=TRUE )
  O$depth.sd = sd( x$depth[ fin.good ], na.rm=TRUE)
  O$depth.n = length( fin.good )
  O$depth.n.bad = length( fin.all) - O$depth.n
  O$depth.smoothed.mean =  mean( x$depth.smoothed[ fin0:fin1 ], na.rm=TRUE )
  O$depth.smoothed.sd = sd( x$depth.smoothed[ fin0:fin1 ], na.rm=TRUE )
  O$good = O$good
  O$depth.filtered = fin.good
  O$depth.smoothed = x$depth.smoothed
  O$ts = x$ts
  O$timestamp = x$timestamp
  O$signal2noise = O$depth.n / length( fin.all )  # not really signal to noise but rather  % informations

  O$bottom.contact.indices = fin.all
  O$bottom.contact = rep( FALSE, nx )
  O$bottom.contact[ fin.all ] = TRUE


  # estimate surface area if possible using wingspread &/or doorspread
  O$surface.area = NA
  if ( exists ( "wingspread", x ) | exists( "doorspread", x ) ) {
    sa = try( surfacearea.estimate( bcp=bcp, O=O ), silent=TRUE )
    if ( ! "try-error" %in% class( sa ) ) O$surface.area = sa
  }


  # for minilog and seabird data .. we have temperature estimates to make ..
  tmean= NA
  tmeansd = NA
  if (exists( "temperature", x )) {
    tmean = mean( x$temperature[fin.all], na.rm=TRUE )
    tmeansd = sd( x$temperature[fin.all], na.rm=TRUE )
  }
  O$res = data.frame( cbind(z=O$depth.mean, t=tmean, zsd=O$depth.sd, tsd=tmeansd,
                            n=O$depth.n, t0=O$bottom0, t1=O$bottom1, dt=O$bottom.diff ) ) # this is really for the snow crab system

  if(debug.plot) {
    trange = range( x$ts[O$good], na.rm=TRUE )
    drange = c( quantile( x$depth, c(0.05, 0.975), na.rm=TRUE) , median( x$depth, na.rm=TRUE ) * 1.05 )
    plot(depth~ts, x, ylim=c(drange[2],drange[1]), xlim=c(trange[1],trange[2]), pch=20, cex=0.1, col="gray" )
    mcol = "yellow"
    points( depth~ts, x[ O$maxdepth.method.indices, ], pch=20, col=mcol, cex=0.2)
  }

  print( O$summary)

  return( O )

}





