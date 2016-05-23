
bottom.contact.filter.noise = function( x, good, bcp, debug=FALSE ) {

  # First, a simple high pass filter to remove noisy data (very large extreme fluctuations)
  # Then, truncation of data based upon variance of central region to permit use of more fine-scale methods
  # to determine anomlaous data

  # do not operate upon depth (keep as reference)... operate upon Z and Z.smoothed
  x$Z = x$depth  # x$Z will be returned ...
  nx = nrow(x)

  # 1. check range of the data and recheck using modal filter of the bottom part of the "good" area
  fr = range( which(good) )
  aoi = fr[1]:fr[2]

  Zmin = min( x$depth, na.rm=TRUE )
  tooshallow = which( x$Z < Zmin)
  if (length( tooshallow)>0) {
    good[ tooshallow ] = FALSE
    x$Z[ tooshallow ] = NA
  }

  mm = modes( x$Z[aoi] )
  mm.aoi =  which ( x$Z >= mm$lb2*0.9 & x$Z <= mm$ub2*1.1  ) # fisrt estimate of bottom

  Zmax = max (x$depth[mm.aoi], na.rm=TRUE ) + 40  # sometimes very steep tows ...
  toodeep = which( x$Z > Zmax )
  if (length( toodeep )>0) {
    good[ toodeep ] = FALSE
    x$Z[ toodeep ] = NA
  }

  aoi = min( mm.aoi ) : max( mm.aoi )
  if (length( aoi) < 30 ) return (NULL)

  # 2. focus upon aoi: identify data that are likely noise and mark them
  #  filter localized (adjacent) deviations
  tm0 = x[aoi, c("ts", "Z")]

  tm = tm0
  tm$Z[ !good[aoi]] = NA
  test = NULL
  test = try( interpolate.xy.robust( tm, method="sequential.linear",
      trim=bcp$noisefilter.trim, probs=bcp$noisefilter.quants ), silent =TRUE  )
  if ( ! ( class( test ) %in% "try-error" | ( length( which(is.finite(test))) < 30)  ) ) {
    kk = x$Z[aoi] - test
    # ll =  which( abs(kk) < bcp$eps.depth )
    # i = intersect( ll, which.quantile ( kk, probs=bcp$noisefilter.quants, inside=FALSE )  )
    i = which.quantile ( kk, probs=bcp$noisefilter.quants, inside=FALSE )
    j = setdiff(1:nx , i )
    if (length(i) > 0) good[aoi[i]] = FALSE
    # if (length(j) > 0) good[aoi[j]] = TRUE
    x$Z[ !good] = NA
  }

  # filter out based upon local-mean estimates
  tm = tm0
  tm$Z[ !good[aoi] ] = NA
  test = NULL
  test = try( interpolate.xy.robust( tm, method="moving.window",
      target.r2=bcp$noisefilter.target.r2, trim=bcp$noisefilter.trim, mv.win=bcp$noisefilter.var.window ), silent =TRUE  )
  if ( ! ( class( test ) %in% "try-error" | ( length( which(is.finite(test))) < 30)  ) ) {
    kk = x$Z[aoi] - test
    # ll =  which( abs(kk) < bcp$eps.depth )
    # i = intersect( ll, which.quantile ( kk, probs=bcp$noisefilter.quants, inside=FALSE )  )
    i = which.quantile ( kk, probs=bcp$noisefilter.quants, inside=FALSE )
    j = setdiff(1:nx, i )
    if (length(i) > 0) good[aoi[i]] = FALSE
    # if (length(j) > 0) good[aoi[j]] = TRUE
    x$Z[ !good] = NA
  }

  # filter out based upon local-variance estimates
  tm = tm0
  tm$Z[ !good[aoi] ] = NA
  test = try( interpolate.xy.robust( tm, method="local.variance",
      target.r2=bcp$noisefilter.target.r2, trim=bcp$noisefilter.trim, mv.win=bcp$noisefilter.var.window ) )
  if ( ! ( class( test ) %in% "try-error" | ( length( which(is.finite(test))) < 30)  ) ) {
    kk = x$Z[aoi] - test
    # ll =  which( abs(kk) < bcp$eps.depth )
    # i = intersect( ll, which.quantile ( kk, probs=bcp$noisefilter.quants, inside=FALSE )  )
    i = which.quantile ( kk, probs=bcp$noisefilter.quants, inside=FALSE )
    j = setdiff(1:nx , i )
    if (length(i) > 0) good[aoi[i]] = FALSE
    # if (length(j) > 0) good[aoi[j]] = TRUE
    x$Z[ !good] = NA
  }

   # filter out based upon local-mean estimates
  tm = tm0
  tm$Z[ !good[aoi] ] = NA
  test = NULL
  test = try( interpolate.xy.robust( tm, method="linear",
      target.r2=bcp$noisefilter.target.r2, trim=bcp$noisefilter.trim, mv.win=bcp$noisefilter.var.window ), silent =TRUE  )
  if ( ! ( class( test ) %in% "try-error" | ( length( which(is.finite(test))) < 30)  ) ) {
    kk = x$Z[aoi] - test
    # ll =  which( abs(kk) < bcp$eps.depth )
    # i = intersect( ll, which.quantile ( kk, probs=bcp$noisefilter.quants, inside=FALSE )  )
    i = which.quantile ( kk, probs=bcp$noisefilter.quants, inside=FALSE )
    j = setdiff(1:nx, i )
    if (length(i) > 0) good[aoi[i]] = FALSE
    # if (length(j) > 0) good[aoi[j]] = TRUE
    x$Z[ !good] = NA
  }

  # finalize update of aoi
  aoi.finite = which( is.finite(x$Z) )
  aoi.range  = range( aoi.finite )
  aoi = aoi.range[1] : aoi.range[2]

  tm2 = x[aoi, c("ts", "Z")]
  test = NULL
  test = try( interpolate.xy.robust( tm2, method=bcp$noisefilter.smoother,
    target.r2=bcp$noisefilter.target.r2, trim=bcp$noisefilter.trim ), silent =TRUE  )
    if ( ! ( class( test ) %in% "try-error" | ( length( which(is.finite(test))) < 30)  ) ) {
      kk = x$depth[aoi] - test
      i = which.quantile ( kk, probs=bcp$noisefilter.quants, inside=TRUE )
      j = which.quantile ( kk, probs=bcp$noisefilter.quants, inside=FALSE )
      # reset "good"
      if (length(i) > 0) good[aoi[i]] = TRUE
      if (length(j) > 0) good[aoi[j]] = FALSE
    }
    rm(tm2); gc()

# 3. finalize id of good data and create a smooth based upon the good data
#   by this point, the Z are the highest probability of good clean data ..
#   operate upon finding a smooth return them to the data for futher processing
  x$Z.smoothed = x$depth  # begin with raw data
  x$Z.smoothed[ aoi.finite ] = x$Z[aoi.finite]
  x$Z.smoothed[!good[aoi]] = NA

# smooth the tails
  buff = 3 ## ~ sec
  leftt = fr[1] : (min(aoi+buff) )  # add a buffer
  rightt = (max(aoi)-buff  ) :  fr[2]   # add a buffer
  tails = c( fr[1]:min(aoi), max(aoi):fr[2] )

  if (length( leftt) > buff ) {
    # force tails to surface .. helps tails to go in right direction when not enough data retained
    x$Z.smoothed[ fr[1] ] = Zmin
    test = NULL
    test = try( interpolate.xy.robust( x[leftt, c("ts", "Z.smoothed")],  method="sequential.linear",
        trim=bcp$noisefilter.trim, probs=bcp$noisefilter.quants ), silent =TRUE  )
    if ( ! ( class( test ) %in% "try-error" | ( length( which(is.finite(test))) < 3)  ) ) x$Z.smoothed[leftt] = test

    x$Z.smoothed[ fr[1] ] = Zmin
    test = NULL
    test = try( interpolate.xy.robust( x[leftt, c("ts", "Z.smoothed")],  method="linear",
        trim=bcp$noisefilter.trim, probs=bcp$noisefilter.quants ), silent =TRUE  )
    if ( ! ( class( test ) %in% "try-error" | ( length( which(is.finite(test))) < 3)  ) ) x$Z.smoothed[leftt] = test

    x$Z.smoothed[ fr[1] ] = Zmin
    test = NULL
    test = try( interpolate.xy.robust( x[leftt, c("ts", "Z.smoothed")],  method="moving.window",
        trim=bcp$noisefilter.trim, probs=bcp$noisefilter.quants ), silent =TRUE  )
    if ( ! ( class( test ) %in% "try-error" | ( length( which(is.finite(test))) < 3)  ) ) x$Z.smoothed[leftt] = test
  }

  if (length( rightt ) > buff ) {
    x$Z.smoothed[ fr[2] ] = Zmin
    test = NULL
    test = try( interpolate.xy.robust( x[rightt, c("ts", "Z.smoothed")],  method="sequential.linear",
        trim=bcp$noisefilter.trim, probs=bcp$noisefilter.quants ), silent =TRUE )
    if ( ! ( class( test ) %in% "try-error" | ( length( which(is.finite(test))) < 3)  ) ) x$Z.smoothed[rightt] = test

    x$Z.smoothed[ fr[2] ] = Zmin
    test = NULL
    test = try( interpolate.xy.robust( x[rightt, c("ts", "Z.smoothed")],  method="linear",
        trim=bcp$noisefilter.trim, probs=bcp$noisefilter.quants ), silent =TRUE )
    if ( ! ( class( test ) %in% "try-error" | ( length( which(is.finite(test))) < 3)  ) ) x$Z.smoothed[rightt] = test

    x$Z.smoothed[ fr[2] ] = Zmin
    test = NULL
    test = try( interpolate.xy.robust( x[rightt, c("ts", "Z.smoothed")],  method="moving.window",
        trim=bcp$noisefilter.trim, probs=bcp$noisefilter.quants ), silent =TRUE )
    if ( ! ( class( test ) %in% "try-error" | ( length( which(is.finite(test))) < 3)  ) ) x$Z.smoothed[rightt] = test

  }

  if (length( tails ) > 0 ) {
    kk = x$depth[tails] - x$Z.smoothed[tails]
    i = which.quantile ( kk, probs=bcp$noisefilter.quants, inside=FALSE )
    if (length(i) > 0 ) {
      good [ tails[i] ]  = FALSE
    }
  }


    # last smooth: with all good data back together ... force tails to be very shallow
  sm = x[, c("ts", "Z")]
  sm$Z[ range(aoi) ] = min( x$depth[aoi], na.rm=TRUE )   # force terminal points to be the most shallow

    # test = NULL
    # test = try( interpolate.xy.robust( sm[aoi,],  method="sequential.linear",
    #     trim=bcp$noisefilter.trim, probs=bcp$noisefilter.quants ), silent =TRUE )
    # if ( ! ( class( test ) %in% "try-error" | ( length( which(is.finite(test))) < 3)  ) ) {
    #     kk = x$depth[aoi] - test
    #     j = which.quantile ( kk, probs=bcp$noisefilter.quants, inside=FALSE )
    #     if (length(j) > 0) {
    #       good[aoi[j]] = FALSE
    #       sm$Z[ !good ] = NA
    #     }
    # }

    # test = NULL
    # test = try( interpolate.xy.robust( sm[aoi,],  method="moving.window",
    #     trim=bcp$noisefilter.trim, probs=bcp$noisefilter.quants ), silent =TRUE )
    # if ( ! ( class( test ) %in% "try-error" | ( length( which(is.finite(test))) < 3)  ) ) {
    # kk = x$depth[aoi] - test
    #     j = which.quantile ( kk, probs=bcp$noisefilter.quants, inside=FALSE )
    #     if (length(j) > 0) {
    #       good[aoi[j]] = FALSE
    #       sm$Z[ !good ] = NA
    #     }
    # }

    # last smooth ... wrap up
    test = NULL
    test = try( interpolate.xy.robust( sm[aoi,], method=bcp$noisefilter.smoother,
      target.r2=bcp$noisefilter.target.r2, trim=bcp$noisefilter.trim ), silent =TRUE  )
    if ( ( class( test ) %in% "try-error" | ( length( which(is.finite(test))) < 30)  ) ) {
      test = try( interpolate.xy.robust( sm[aoi,], method=bcp$noisefilter.smoother,
        target.r2=(bcp$noisefilter.target.r2-0.1), trim=bcp$noisefilter.trim ), silent =TRUE  )
    }
    if ( ! ( class( test ) %in% "try-error" | ( length( which(is.finite(test))) < 30)  ) ) {
      kk = x$depth[aoi] - test
      i = which.quantile ( kk, probs=bcp$noisefilter.quants, inside=TRUE )
      if (length(i) > 0) good[aoi[i]] = TRUE
      x$Z = x$depth
      x$Z[!good] = NA

      sm$Z[!good] = NA
      test = try( interpolate.xy.robust( cbind( sm[aoi,"ts"], test), method=bcp$noisefilter.smoother, # yes smooth again on the smoothed
         target.r2=(bcp$noisefilter.target.r2), trim=bcp$noisefilter.trim ), silent =TRUE  )
      if ( ( class( test ) %in% "try-error" | ( length( which(is.finite(test))) < 30)  ) ) {
        test = try( interpolate.xy.robust( sm[aoi,], method=bcp$noisefilter.smoother,
          target.r2=(bcp$noisefilter.target.r2-0.1), trim=bcp$noisefilter.trim ), silent =TRUE  )
      }
      if ( ! ( class( test ) %in% "try-error" | ( length( which(is.finite(test))) < 30)  ) ) {
        kk = x$depth[aoi] - test
        i = which.quantile ( kk, probs=bcp$noisefilter.quants, inside=TRUE )
        j = which.quantile ( kk, probs=bcp$noisefilter.quants, inside=FALSE )
        if (length(i) > 0) good[aoi[i]] = TRUE
        if (length(j) > 0) good[aoi[j]] = FALSE
        x$Z = x$depth
        x$Z[!good] = NA
        sm$Z =NA
        sm$Z[aoi] = test
      }
    }
    # last resort: loess
    if ( ( class( test ) %in% "try-error" | ( length( which(is.finite(test))) < 30)  ) ) {
      # reset sm ... and do a simple loess
      sm = x[, c("ts", "Z")]
      sm$Z[ range(aoi) ] = min( x$depth[aoi], na.rm=TRUE )   # force terminal points to be the most shallow
      test = try( interpolate.xy.robust( cbind( sm[aoi,"ts"], test), method="loess", # yes smooth again on the smoothed
         target.r2=(bcp$noisefilter.target.r2), trim=bcp$noisefilter.trim ), silent =TRUE  )
      if ( ! ( class( test ) %in% "try-error" | ( length( which(is.finite(test))) < 30)  ) ) {
          kk = x$depth[aoi] - test
          i = which.quantile ( kk, probs=bcp$noisefilter.quants, inside=TRUE )
          j = which.quantile ( kk, probs=bcp$noisefilter.quants, inside=FALSE )
          if (length(i) > 0) good[aoi[i]] = TRUE
          if (length(j) > 0) good[aoi[j]] = FALSE
          x$Z = x$depth
          x$Z[!good] = NA
          sm$Z =NA
          sm$Z[aoi] = test
      }
    }

    if(debug) {
      plot( depth~ ts, x, pch=20, col="black", cex= 0.5, ylim=rev(range(x$depth, na.rm=T)) )
      points( Z ~ ts, x, pch=20, col="red", cex=0.8 )
      points ( Z ~ ts, sm, col="green")
      lines ( Z ~ ts, sm[aoi,], col="orange", lwd=3)
      browser()
    }

  return( list( depth.smoothed=sm$Z, good=good, depth.filtered=x$Z ) )
}

