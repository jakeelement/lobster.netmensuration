
bottom.contact.smooth = function( sm, bcp )  {
  ## Smooth method: smooth data to remove local trends and compute first derivatives

  res = list( bc0=NA, bc1=NA)

  names( sm) = c("Z", "timestamp", "ts" ) # Z is used to abstract the variable name such that any variable can be passed
  N = nrow(sm)

  # expand tails a bit more:
  Z0 = modes( sm$Z )
  aoi.i = which( sm$Z > (Z0$lb2 + bcp$depth.range[1]/2) & sm$Z < (Z0$ub2 + bcp$depth.range[2]/2 ) )

  if (length( which( is.finite( aoi.i))) < 10  ) return( res) # need some data

  sm$dZ = grad( approxfun( sm$ts, sm$Z, rule=2 ), sm$ts, method="simple" )

  dZ.modes = modes( sm$dZ[aoi.i] )
  ii = which( abs(sm$dZ) < dZ.modes$sd * bcp$smooth.zeros ) # remove discretization induced zero-values
  sm$dZ[ii] = NA

  ub =  dZ.modes$ub2
  lb =  dZ.modes$lb2

  # remove high freq variation
  sm$dZ.smoothed = interpolate.xy.robust( sm[, c("ts", "dZ")],
        trim=bcp$noisefilter.trim, probs=bcp$smooth.filter.quants, method="sequential.linear" )

  sm$dZ.smoothed = interpolate.xy.robust( sm[, c("ts", "dZ.smoothed")],
        trim=bcp$noisefilter.trim, probs=bcp$smooth.filter.quants, method="moving.window" )

  ml0 = which.max( sm$dZ.smoothed )
  mr0 = which.min( sm$dZ.smoothed )

  dZ.sd = sd( sm$dZ.smoothed[aoi.i], na.rm=TRUE )
  ub =  dZ.modes$mode + dZ.sd * bcp$smooth.sd.multiplier
  lb =  dZ.modes$mode - dZ.sd * bcp$smooth.sd.multiplier

  # now using (smoothed) first derivative determine inflection points (crossing of the (near)-zero line)
  iib = intersect( which( sm$dZ.smoothed <= ub & sm$dZ.smoothed >= lb ), c(ml0:mr0))
  ib = min(iib): max(iib)

  ng = trunc( length(ib) / 10 ) # small chunk close to candidate area of interest
  ml = ib[ ng ]
  mr = ib[ (length(ib) - ng) ]

  smooth.threshold = 0.5
  mwin = - bcp$smooth.windowsize : bcp$smooth.windowsize

  # left side
  if ( abs( min(sm$Z[ml0:ml], na.rm=TRUE) - Z0$mode ) < Z0$sd ) {
    # assuming it is truncated without the left tail ..
    k0 = 1
  } else {
    for( k0 in ml0:ml) {
      ri = k0+mwin
      ri = ri[ ri >= 1 & ri <=N ]
      nri = length(ri)
      if ( nri >= 1.5*(bcp$smooth.windowsize) ) {
        nx = length( which( sm$dZ.smoothed[ri] < ub ) )
        prx = nx / nri
        if (prx > smooth.threshold) break()
      }
    }
  }

  if ( abs( min(sm$Z[mr0:mr], na.rm=TRUE) - Z0$mode ) < Z0$sd ) {
    # assuming it is truncated without the right tail ..
    k1 = N
  } else {
    for( k1 in mr0:mr) {
      ri = k1+mwin
      ri = ri[ ri >= 1 & ri <=N ]
      nri = length(ri)
      if ( nri >= 1.5*(bcp$smooth.windowsize) ) {
        nx = length( which( sm$dZ.smoothed[ri] > lb ) )
        prx = nx / nri
        if (prx > smooth.threshold) break()
      }
    }
  }

  res =  list( bc0=sm$timestamp[k0], bc1=sm$timestamp[k1] )

  return(res)
}


