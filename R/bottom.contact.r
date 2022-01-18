#' @title  bottom.contact
#' @description  Calculates indicies of touchdown and liftoff
#' @import tcltk geosphere lubridate
#' @return dataframe with various results
#' @export
bottom.contact = function( x, bcp, debugrun=FALSE ) {

require(tcltk)
  # all timestamps must be in Posix/UTC
  # range comparisons seem to fail when they are not
  tzone = "UTC"

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

  if ( sd(x$depth, na.rm=TRUE) < bcp$eps.depth ) {
    O$error.flag = "not enough variability in data?"
    return(O)
  }

  ## inSANITY CHECKS
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
  mm.med = floor( median( mm.bot, na.rm=TRUE ) )
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

  if ( sd(x$depth, na.rm=TRUE) < bcp$eps.depth ) return(NULL)  # not enough variability in data

  ## ------------------------------
  ## MAIN NOISE/INTERPOLATION FILTER
  # Some filtering of noise from data and further focus upon area of interest based upon time and depth if possible
  res = NULL
  res = try( bottom.contact.filter.noise ( x=x, good=O$good, bcp=bcp, debug=debugrun ), silent =TRUE )
  if ( "try-error" %in% class(res) || is.null(res$depth.smoothed) ) {
    x$depth = jitter( x$depth )
    res = try( bottom.contact.filter.noise ( x=x, good=O$good, bcp=bcp, debug=debugrun ), silent =TRUE )
  }

  if ( !"try-error" %in% class(res) ) {
    if (!is.null(res$depth.smoothed)) {
      if ( cor( x$depth, res$depth.smoothed, use="pairwise.complete.obs") > 0.999 || is.na(cor( x$depth, res$depth.smoothed, use="pairwise.complete.obs") )) {
        bcp$noisefilter.target.r2 = bcp$noisefilter.target.r2 - 0.1
        x$depth = jitter( x$depth )
        res = try( bottom.contact.filter.noise ( x=x, good=O$good, bcp=bcp, debug=debugrun ), silent =TRUE )
  }}}

  x$depth.smoothed = x$depth
  if ( ! "try-error" %in% class( res) )  {
    # retain good/bad information only for the interval where there is data
    # .. this effectively does the filtering in the area of interest only and retains the tails
    # for further influence later
    ig = range( which( res$good ))
    igg = ig[1]:ig[2]
    O$good[igg] = res$good[igg]
    x$depth[ !O$good ] = NA
    if (!is.null(res$depth.smoothed)) x$depth.smoothed= res$depth.smoothed
  }

  if(sum(x$depth-min(x$depth,na.rm=T),na.rm=T)==0) return( NULL )
  if(sum(O$good)==0) return( NULL )
  if (sd(x$depth, na.rm=TRUE) < bcp$eps.depth ) return(NULL)  # not enough variability in data
  if (sd(x$depth.smoothed, na.rm=TRUE) < bcp$eps.depth ) return(NULL)  # not enough variability in data

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
    if ("bc0" %in% names(res)) {
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
  aoi.mid = floor( mean( aoi.range ) ) # approximate midpoint
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
      if ("bc0" %in% names(res)) {
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
    #BC - Plots fail in RStudio graphics device, add clause
    dev.new(noRStudioGD = TRUE)
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
      if ("bc0" %in% names(res)) {
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
    #BC - Plots fail in RStudio graphics device, add clause
    dev.new(noRStudioGD = TRUE)
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
    if ("bc0" %in% names(res)) {
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
    #BC - Plots fail in RStudio graphics device, add clause
    dev.new(noRStudioGD = TRUE)
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
     if ("bc0" %in% names(res)) {
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
    #BC - Plots fail in RStudio graphics device, add clause
    dev.new(noRStudioGD = TRUE)
    plot(depth~ts, x, ylim=c(drange[2],drange[1]), xlim=c(trange[1],trange[2]), pch=20, cex=0.1, col="gray" )
    mcol = "red"
    points( depth~ts, x[ O$linear.method.indices, ], pch=20, col=mcol, cex=0.2)
  }


  ## ---------------------------

  O$manual.method0 = NA #### NOTE:: using the 'c' operator on posix strips out the timezone info! this must be retained
  O$manual.method1 = NA  ### NOTE:: using the 'c' operator on posix strips out the timezone info! this must be retained

  if ( bcp$user.interaction  ) {

    satisified = F
    while(!satisified){

    #open plot window, user can modify the numbers to fit their screen
    x11(width = 25, height = 15)


    #Define Margins for multiple y-axis
    par(mar=c(5, 8, 4, 4) + 0.1)

    grange = min(which(O$good)):max(which(O$good))

    if("wingspread" %in% names(O$plotdata)){
    #Plot the door spread
    ylim=c(min(O$plotdata$wingspread[grange], na.rm = TRUE)-1,max(O$plotdata$wingspread[grange], na.rm = TRUE)+1)
    if(is.finite(ylim)){
    plot(O$plotdata$timestamp[grange][which(!is.na(O$plotdata$wingspread[grange]))], O$plotdata$wingspread[grange][which(!is.na(O$plotdata$wingspread[grange]))], axes=F, ylim=ylim, xlab="", ylab="",type="p",col="#0000FF1A", main="",xlim=c(min(O$plotdata$timestamp[grange]), max(O$plotdata$timestamp[grange])))
      if(length(unique(na.omit(O$plotdata$wingspread[grange]))) > 5 ){
       smoothingSpline = smooth.spline(O$plotdata$timestamp[grange][which(!is.na(O$plotdata$wingspread[grange]))], O$plotdata$wingspread[grange][which(!is.na(O$plotdata$wingspread[grange]))], spar=.5)
       lines(smoothingSpline, col="blue")
      }
    abline( h = c(median(O$plotdata$wingspread[grange], na.rm = TRUE)), col = "blue", lty =2 )
    axis(2,col="blue", col.lab = "blue", col.axis = "blue", lwd=1)
    mtext(2,text="Wing Spread",col = "blue",line=2)
    }
    }


    yl=c(min(O$plotdata$opening[grange], na.rm = TRUE)-1,max(O$plotdata$opening[grange], na.rm = TRUE)+1)


    if(!is.numeric(yl[1]) | is.na(yl[1]) ) yl[1] = 0
    if(!is.numeric(yl[2]) | is.na(yl[2]) ) yl[2] = 10

    #Plot the opening
    if("opening" %in% names(O$plotdata)){
      if(length(which(!is.na(O$plotdata$opening[grange]))) > 10){
    par(new=T)
    plot(O$plotdata$timestamp[grange][which(!is.na(O$plotdata$opening[grange]))], O$plotdata$opening[grange][which(!is.na(O$plotdata$opening[grange]))], axes=F, ylim=yl, xlab="", ylab="",type="p",col="#A52A2A1A", main="",xlim=c(min(O$plotdata$timestamp[grange]), max(O$plotdata$timestamp[grange])))
    smoothingSpline = smooth.spline(O$plotdata$timestamp[grange][which(!is.na(O$plotdata$opening[grange]))], O$plotdata$opening[grange][which(!is.na(O$plotdata$opening[grange]))], spar=.5)
    lines(smoothingSpline, col="brown")
    axis(2, ylim=yl,col = "brown",col.lab = "brown", col.axis = "brown",lwd=1,line=3.5)
    mtext(2,text="Opening", col = "brown", line=5.5)
    }
    }

    #Plot the depth (Clickable)
    par(new=T)

    plot(O$plotdata$timestamp[grange][which(!is.na(O$plotdata$depth[grange]))], O$plotdata$depth[grange][which(!is.na(O$plotdata$depth[grange]))], axes=F, ylim=rev(c(min(O$plotdata$depth[grange], na.rm = TRUE)-1,max(O$plotdata$depth[grange], na.rm = TRUE)+1)), xlab="", ylab="", type="l", main=O$id,xlim=c(min(O$plotdata$timestamp[grange]), max(O$plotdata$timestamp[grange])),lwd=1)
    axis(4, ylim=rev(c(min(O$plotdata$depth[grange], na.rm = TRUE)-1,max(O$plotdata$depth[grange], na.rm = TRUE)+1)),lwd=1)
    mtext(4,text="Depth", line = 2)

    #Set up time axis plotting
    abli = seq(O$plotdata$timestamp[grange][1],O$plotdata$timestamp[grange][length(O$plotdata$timestamp[grange])], "1 min")

    #Draw the time axis
    axis.POSIXct(1, at = abli, format = "%H:%M:%S", labels = TRUE)
    ablit = unique(as.character(lubridate::date(abli)))
    mtext(text = ablit, side=1, col="black", line=2)
    abline( v = abli, col = "lightgrey")


    #Visualize neutral haulback
    #Find positions of going backward toward origin. Nuetral haulback
    disfromorigin = rep(0, 1, length(O$plotdata$latitude[grange]))
    for(k in 1:(length(O$plotdata$latitude[grange])-1)){
      p1 = c(O$plotdata$longitude[grange][1],O$plotdata$latitude[grange][1])
      p2 = c(O$plotdata$longitude[grange][k+1],O$plotdata$latitude[grange][k+1])
      if(!is.null(p1) & !is.null(p2)){
      if(!is.na(p1) & !is.na(p2))
           {
             disfromorigin[k+1] = geosphere::distHaversine(p1 , p2, r=6378137)
    }
    }
    }
    neu = rep(0, 1, length(O$plotdata$latitude[grange]))
    for(k in 1:(length(disfromorigin)-1)){
      neu[k] = disfromorigin[k+1]-disfromorigin[k]
      if(k > 1){
        if(neu[k-1]<0 & neu[k]>=0) neu[k-1] = 0
      }
    }

    neuhaul = F
    indn = which(neu < 0)
    if(length(indn) > 10){
      print(paste("Neutral Haulback! ", O$id))
      neuhaul = T
    }

    if(length(indn)>0){
      abline(v = O$plotdata$timestamp[grange][indn], col = "red")
    }

    id = identify(O$plotdata$timestamp[grange], O$plotdata$depth[grange], labels = c(as.character(O$plotdata$timestamp[grange])),n = 1, pos = TRUE)
    start = O$plotdata$timestamp[grange][id$ind]
    points(O$plotdata$timestamp[grange][id$ind], O$plotdata$depth[grange][id$ind], col = "darkgreen", pch = 19)
    #Code that allow identifying clicks
    id2 = identify(O$plotdata$timestamp[grange], O$plotdata$depth[grange], labels = c(as.character(O$plotdata$timestamp[grange])),n = 1, pos = TRUE)
    end = O$plotdata$timestamp[grange][id2$ind]
    points(O$plotdata$timestamp[grange][id2$ind], O$plotdata$depth[grange][id2$ind], col = "darkgreen", pch = 19)

    ##Code taken from bottom.contact.plot.r
    legendtext = NULL
    legendcol = NULL
    legendpch = NULL
    if (all(is.finite(c( O$variance.method0, O$variance.method1) ) ) ) {
      mcol = "pink"
      # points( depth~ts, x[ O$variance.method.indices, ], pch=20, col=mcol, cex=0.2)
      abline (v=O$plotdata$timestamp[min(O$variance.method.indices)], col=mcol, lty="solid", lwd=1.2)
      abline (v=O$plotdata$timestamp[max(O$variance.method.indices)], col=mcol, lty="solid", lwd=1.2)
      DT = as.numeric( difftime( O$variance.method1, O$variance.method0, units="mins" ) )
      legendtext = c( legendtext, paste( "variance gate:   ", round( DT, 2), "" ) )
      legendcol = c( legendcol, mcol)
      legendpch =c( legendpch, 20 )
    }


    if (all(is.finite( c(O$modal.method0, O$modal.method1 ) ) ) ) {
      mcol = "red" # colour for plotting
      # points( depth~ts, x[O$modal.method.indices,], col=mcol, pch=20, cex=0.2)
      abline (v=O$plotdata$timestamp[min(O$modal.method.indices)], col=mcol, lty="dashed")
      abline (v=O$plotdata$timestamp[max(O$modal.method.indices)], col=mcol, lty="dashed")
      DT = as.numeric( difftime( O$modal.method1, O$modal.method0, units="mins" ) )
      legendtext = c( legendtext, paste( "modal:   ", round( DT, 2) ) )
      legendcol = c( legendcol, mcol)
      legendpch =c( legendpch, 20)
    }


    if ( all(is.finite( c( O$smooth.method0,  O$smooth.method1) ) ) ) {
      mcol = "blue"
      # points( depth~ts, x[O$smooth.method.indices,], col=mcol, pch=20, cex=0.2)
      abline (v=O$plotdata$timestamp[min(O$smooth.method.indices)], col=mcol, lty="dotdash", lwd=1.5)
      abline (v=O$plotdata$timestamp[max(O$smooth.method.indices)], col=mcol, lty="dotdash", lwd=1.5)
      DT = as.numeric( difftime( O$smooth.method1, O$smooth.method0, units="mins" ) )
      legendtext = c(legendtext, paste( "smooth:   ", round(DT, 2)) )
      legendcol = c( legendcol, mcol)
      legendpch =c( legendpch, 20)
    }

    if (all(is.finite( c(O$linear.method0, O$linear.method1)  )) ) {
      mcol ="green"
      # points( depth~ts, x[O$linear.method.indices,], col=mcol, pch=20, cex=0.2)
      abline (v=O$plotdata$timestamp[min(O$linear.method.indices)], col=mcol, lty="twodash")
      abline (v=O$plotdata$timestamp[max(O$linear.method.indices)], col=mcol, lty="twodash")
      DT = as.numeric( difftime( O$linear.method1, O$linear.method0, units="mins" ) )
      legendtext = c( legendtext, paste( "linear: ", round( DT, 2) ) )
      legendcol = c( legendcol, mcol)
      legendpch =c( legendpch, 20)
    }


    if (all(is.finite( c( O$maxdepth.method0,  O$maxdepth.method1) )) ) {
      mcol ="orange"
      # points( depth~ts, x[O$linear.method.indices,], col=mcol, pch=20, cex=0.2)
      abline (v=O$plotdata$timestamp[min(O$maxdepth.method.indices)], col=mcol, lty="solid")
      abline (v=O$plotdata$timestamp[max(O$maxdepth.method.indices)], col=mcol, lty="solid")
      DT = as.numeric( difftime( O$maxdepth.method1, O$maxdepth.method0, units="mins" ) )
      legendtext = c( legendtext, paste( "maxdepth: ", round( DT, 2) ) )
      legendcol = c( legendcol, mcol)
      legendpch =c( legendpch, 20)
    }


    if ( !( is.null( legendtext)))  legend( "top", legend=legendtext, col=legendcol, pch=legendpch, bg = "white" )

    #END code taken from bottom.contact


    #User interface

     goodans = FALSE
      res <- tkmessageBox(title = "ManualTD",
                          message = "Are you happy with these values?", icon = "info", type = "yesno")


      if(grepl("yes", res)){
        satisified = TRUE
        goodans = TRUE
        O$manual.method0 = start
        O$manual.method1 = end
        O$manual.method.indices = which( O$plotdata$timestamp >= O$manual.method0 &  O$plotdata$timestamp <= O$manual.method1 )


       # result = rbind(result, row)

        # if(!file.exists(file)){
        #   write.table(row, file = file, sep = ",", row.names = F, quote = F)
        # }else{
        #   write.table(row, file = file, sep = ",", append = T, col.names = F, row.names = F, quote = F )
        # }
      }
      if(grepl("no", res)){

       satisified = FALSE
        goodans = TRUE
      }
      if(!goodans)print("Did not receive valid command, please try again.")

      dev.off()
  }#END plotting loop while not satisfied

  if(is.na(O$manual.method0) || is.na(O$manual.method1)){
    bcp$user.interaction = FALSE
  }
  else{
    if(is.null(bcp$station)) station = unlist(strsplit(bcp$id, "\\."))[4]
    else station = bcp$station
    bcp$id = bcp$trip
    mf = NULL
    if(!is.null(bcp$from.manual.archive)) mf = file.path(bcp$from.manual.archive, "clicktouchdown_all.csv")
    if(!is.null(bcp$from.manual.file)) mf = bcp$from.manual.file
    if(!is.null(mf)){
      manualclick = NULL
      if(file.exists(mf)){
        manualclick = read.csv(mf, as.is=TRUE)

        if(bcp$datasource == "lobster"){
          sta.ind = which(manualclick$station == station & manualclick$trip == bcp$trip)
        }
        else{
          sta.ind = which(manualclick$station == station & manualclick$year == bcp$YR)
        }
        if(length(sta.ind) == 0){
          print(bcp$trip)
           manualclick = rbind(manualclick, data.frame(station = station, start = unlist(as.character(O$manual.method0)), end = unlist(as.character(O$manual.method1)), depth = mean( x$depth, na.rm=TRUE ), year = bcp$YR, trip = bcp$id))
        }
        else{
          manualclick$station[sta.ind] = station
          manualclick$start[sta.ind] = as.character(O$manual.method0)
          manualclick$end[sta.ind] = as.character(O$manual.method1)
          manualclick$depth[sta.ind] = mean( x$depth, na.rm=TRUE )
          manualclick$year[sta.ind] = bcp$YR
          manualclick$trip[sta.ind] = bcp$id
        }
      }
      else{
        manualclick = data.frame(station = station, start = as.character(O$manual.method0), end = as.character(O$manual.method1), depth =  mean( x$depth, na.rm=TRUE ), year = bcp$YR, trip = bcp$id)
      }
      write.csv(manualclick, mf, row.names = FALSE )
    }
  }

  }
  
  
  #BC- Added condition in case user interaction is desired, no need to do this step.
  if (!is.null(bcp$from.manual.archive) && !bcp$user.interaction) {
     print( "Loading values from previously generated .csv")
    if(file.exists(file.path(bcp$from.manual.archive, "clicktouchdown_all.csv"))){
     manualclick = read.csv(file.path(bcp$from.manual.archive, "clicktouchdown_all.csv"), as.is=TRUE)
     station = unlist(strsplit(bcp$id, "\\."))[4]
     sta.ind = which(manualclick$station == station & manualclick$year == bcp$YR)
     if(length(sta.ind == 1)){
       O$manual.method0 = lubridate::ymd_hms(manualclick$start[sta.ind], tz = "UTC")
       O$manual.method1 = lubridate::ymd_hms(manualclick$end[sta.ind], tz = "UTC")

     }
    }
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

  # recovert to time zone of incoming data as time zone is lost with the transpose
  tmp0 = lubridate::ymd_hms( t( as.data.frame(O[ bcm0 ]) ) )  # UTC
  tmp1 = lubridate::ymd_hms( t( as.data.frame(O[ bcm1 ]) ) )  # UTC

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
  if ( "wingspread" %in% names(x)  | "doorspread" %in% names(x) )  {
   sa = try( surfacearea.estimate( bcp=bcp, O=O ), silent=TRUE )
    if ( ! "try-error" %in% class( sa ) ) O$surface.area = sa
  }


  # for minilog and seabird data .. we have temperature estimates to make ..
  tmean= NA
  tmeansd = NA
  if ("temperature" %in% names(x) ) {
    tmean = mean( x$temperature[fin.all], na.rm=TRUE )
    tmeansd = sd( x$temperature[fin.all], na.rm=TRUE )
  }
  O$res = data.frame( cbind(z=O$depth.mean, t=tmean, zsd=O$depth.sd, tsd=tmeansd,
                            n=O$depth.n, t0=O$bottom0, t1=O$bottom1, dt=O$bottom.diff ) ) # this is really for the snow crab system

  # time format / zone gets reset ..
  O$res$t0 = as.POSIXct( O$res$t0, origin=lubridate::origin, tz="UTC" )
  O$res$t1 = as.POSIXct( O$res$t1, origin=lubridate::origin, tz="UTC" )
  O$res$dt = O$res$t1 - O$res$t0

  if(debug.plot) {
    trange = range( x$ts[O$good], na.rm=TRUE )
    drange = c( quantile( x$depth, c(0.05, 0.975), na.rm=TRUE) , median( x$depth, na.rm=TRUE ) * 1.05 )
    #BC - Plots fail in RStudio graphics device, add clause
    dev.new(noRStudioGD = TRUE)
    plot(depth~ts, x, ylim=c(drange[2],drange[1]), xlim=c(trange[1],trange[2]), pch=20, cex=0.1, col="gray" )
    mcol = "yellow"
    points( depth~ts, x[ O$maxdepth.method.indices, ], pch=20, col=mcol, cex=0.2)
  }

  print( O$summary)

  return( O )

}





