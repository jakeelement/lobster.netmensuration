

surfacearea.estimate = function( bcp, O ) {
  
  nms = O$plotdata[ O$bottom.contact,]

  if (! ( exists( "wingspread", nms) || exists( "doorspread", nms) ) ) {
    warning( "Wingspread and doorspread were not found" )
    return(NA) 
  }

  nms$sm.wingspread = NA
  nms$sm.doorspread = NA 


  nd = nrow(nms)
  if (0) {
    plot( doorspread ~ timestamp, nms )
    plot( wingspread ~ timestamp, nms )
    plot( wingspread ~ doorspread, nms, type="b" )
    plot( wingspread ~ depth, nms, type="l" )
    plot( doorspread ~ depth, nms, type="l" )
  }

  # estimate/predict missing data where possible:
  jj = which( is.na( nms$doorspread ) & !is.na( nms$wingspread) )
  if (length( jj ) > 0 ) {
    kd = lm( doorspread ~ wingspread, nms, na.action="na.omit" )
    if ( ! "try-error" %in% class( kd) ) nms$doorspread[ jj ] = predict( kd, nms[jj,] ) 
  }

  ii = which( is.na( nms$wingspread ) & !is.na( nms$doorspread) )
  if (length( ii) > 0 ) {
    kw = try (lm( wingspread ~ doorspread, nms, na.action="na.omit" ), silent=TRUE) 
    if ( ! "try-error" %in% class( kw) ) nms$wingspread[ ii ] = predict( kw, nms[ii,] ) 
  }


  ## smooth and filter wingspread
  ii = which( is.finite( nms$wingspread ) ) 
  if (length(ii) > 10 ){

    nms$sm = interpolate.xy.robust( nms[, c("ts", "wingspread")], method="sequential.linear",  
        trim=bcp$noisefilter.trim, probs=bcp$noisefilter.quants )
    nms$sm = interpolate.xy.robust( nms[, c("ts", "sm")], method="moving.window", 
        trim=bcp$noisefilter.trim, target.r2=bcp$noisefilter.target.r2, mv.win=bcp$noisefilter.var.window )
      
    kk = nms$wingspread - nms$sm
    i = which.quantile ( kk, probs=bcp$noisefilter.quants, inside=FALSE ) 
    if (length(i) > 0) nms$wingspread [ i ] = NA
    
    # redo interpolation with filtered data removed
    test = NULL
    test = interpolate.xy.robust( nms[, c("ts", "wingspread")], method="sequential.linear",
      probs=bcp$noisefilter.quants, target.r2=bcp$noisefilter.target.r2, trim=bcp$noisefilter.trim )
    if ( "try-error" %in% class(test) ) {
      test = interpolate.xy.robust( nms[, c("ts", "wingspread")], method="moving.window",  
        trim=bcp$noisefilter.trim, probs=bcp$noisefilter.quants )
    }
    nms$sm.wingspread = test 
  }

    ### --- doorspread

  ii = which( is.finite( nms$doorspread ) ) 
  if (length(ii) > 10 ) {

    nms$sm = interpolate.xy.robust( nms[, c("ts", "doorspread")], method="sequential.linear" , 
        trim=bcp$noisefilter.trim, probs=bcp$noisefilter.quants )
    nms$sm = interpolate.xy.robust( nms[, c("ts", "sm")], method="moving.window",
        trim=bcp$noisefilter.trim, target.r2=bcp$noisefilter.target.r2, mv.win=bcp$noisefilter.var.window )
      
    kk = nms$doorspread - nms$sm
    i = which.quantile ( kk, probs=bcp$noisefilter.quants, inside=FALSE )  
    if (length(i) > 0) nms$doorspread [ i ] = NA
    
    # redo interpolation with filtered data removed
    test =NULL
    test  = interpolate.xy.robust( nms[, c("ts", "doorspread")], method="sequential.linear", 
        probs=bcp$noisefilter.quants, target.r2=bcp$noisefilter.target.r2, trim=bcp$noisefilter.trim )
    if ( "try-error" %in% class(test) ) {
      test = interpolate.xy.robust( nms[, c("ts", "doorspread")], method="moving.window" , 
        trim=bcp$noisefilter.trim, probs=bcp$noisefilter.quants )
    }

    nms$sm.doorspread = test
  }

    ## depths -- smooth and filter one last time

    # dh = NULL  # incremental distance horizonatal
    # dv = NULL  # incremental distance vertical 
    # smooth again as there is occasionally v. high freq noise still in the depth data
    depth.smoothed = O$depth.smoothed [ O$bottom.contact ]
#    depth.smoothed = interpolate.xy.robust( cbind( nms$ts, depth.smoothed), method="moving.window",  
#        target.r2=bcp$noisefilter.target.r2, trim=bcp$noisefilter.trim, mv.win=bcp$noisefilter.var.window )  

    dv = c(0, abs( diff( depth.smoothed )) / 1000 ) # km  .. incremental difference in vertical distance
    dv = interpolate.xy.robust( cbind( nms$ts, dv), method="sequential.linear" , 
        trim=bcp$noisefilter.trim, probs=bcp$noisefilter.quants )
    dv = interpolate.xy.robust( cbind( nms$ts, dv), method="moving.window",  
        target.r2=bcp$noisefilter.target.r2, trim=bcp$noisefilter.trim, mv.win=bcp$noisefilter.var.window )  


    # determine the algorithm to use for use to determine tow distance
    coords =  c( "longitude", "latitude" )
    dh = rep(NA, nd)  # km.. incremental distance of  horizontal component  
    npos = length( unique( paste( nms$longitude, nms$latitude) ) )
    if (npos > 30) { 
      # good gps positioning in modern series 
      # otherwise if poor spatial resolution: resort to using ship speeds and to estimate horizontal distance .. 
      # this needs to be done at a project level so, here we set the SA estimates to NA and provide best dist vertical
      # and net spread data
   
      # sometimes GPS loses contact and GPS position stays static though still moving .. 
      # elimiate zero/small values and interpolate where required 
      H = geosphere::distGeo ( nms[ 1:(nd-1), coords ], nms[ 2:nd, coords ] ) / 1000 ## in meters .. convert to km

      oo = which( H < bcp$gps.distance.range.valid.km[1] | H > bcp$gps.distance.range.valid.km[2] )
      if (length( oo) > 1 ) {
        nms$longitude[oo] = NA
        nms$latitude [oo] = NA
        nms$longitude = interpolate.xy.robust( nms[, c("ts", "longitude") ], method="sequential.linear" , 
          trim=bcp$noisefilter.trim, probs=bcp$noisefilter.quants )
        nms$latitude = interpolate.xy.robust( nms[, c("ts", "latitude") ], method="sequential.linear" , 
          trim=bcp$noisefilter.trim, probs=bcp$noisefilter.quants )
        H = geosphere::distGeo ( nms[ 1:(nd-1), coords ], nms[ 2:nd, coords ] ) / 1000 ## in meters .. convert to km
      }

      # recheck
      pp = which( H < bcp$gps.distance.range.valid.km[1] | H > bcp$gps.distance.range.valid.km[2] )
      if (length( pp ) > 1 ) H[pp] = NA
      H = interpolate.xy.robust( cbind( nms$ts[1:(nd-1)], H), method="sequential.linear" , 
          trim=bcp$noisefilter.trim, probs=bcp$noisefilter.quants )
      H = interpolate.xy.robust( cbind( nms$ts[1:(nd-1)], H), method="moving.window",  
          target.r2=bcp$noisefilter.target.r2, trim=bcp$noisefilter.trim, mv.win=bcp$noisefilter.var.window )  
      dh = c( 0, H ) 
    }

    ## integration step
    dd = sqrt( dh^2 + dv^2)  ## total incremental distance, km (path length)
    sa.incremental.door = dd * ( nms$sm.doorspread / 1000 )  # km^2 (path length * path width) 
    sa.ds = sum( abs( sa.incremental.door ) ) 

    sa.incremental.wing = dd * ( nms$sm.wingspread / 1000 )  # km^2
    sa.ws = sum( abs( sa.incremental.wing ) ) 

    out = list( 
      distances.total = cumsum(dd) , 
      distances.vertical = cumsum(dv) ,
      distances.horizontal = cumsum(dh) ,
      door.sa.cum = cumsum( sa.incremental.door ), 
      wing.sa.cum = cumsum( sa.incremental.wing ),
      door.sa=sa.ds, 
      wing.sa=sa.ws, 
      door.mean=mean(nms$sm.doorspread), 
      wing.mean=mean(nms$sm.wingspread),  
      door.sd=sd(nms$sm.doorspread), 
      wing.sd=sd(nms$sm.wingspread),
      wing.smoothed=nms$sm.wingspread, 
      door.smoothed=nms$sm.doorspread )  
    
    return (out)

  }

