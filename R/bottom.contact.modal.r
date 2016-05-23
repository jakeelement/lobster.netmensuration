
bottom.contact.modal = function( sm, bcp ) {
  # Modal method: gating by looking for modal distribution and estimating sd of the modal group
  # until a target number of breaks, nbins with valid data are found
  # see density and bw.SJ for more options

  res = list( bc0=NA, bc1=NA)

  names(sm) = c("Z", "timestamp", "ts" ) # Z is used to abstract the variable name such that any variable can be passed
  N = nrow(sm)
  if ( N < 50 ) return(res)  # insufficient data

  sm$Z = sm$Z + runif( N, min=-0.5, max=0.5 ) * bcp$eps.depth   # add some noise as data are too "flat"

  # first pass to crudely determine location of mode
  Z0 = modes( sm$Z )

  r0 =1
  r1 =N
  # catch situations where there is no variance in data due to only bottom-track being recorded and truncated:
  # no discrimination possible
  if ( !is.finite(Z0$sd) || Z0$sd < 1e-3) {
    oo = which ( abs( sm$Z - Z0$mode ) < bcp$eps.depth )  # anything within eps m of the "mode"
    if (length(oo) > 50) {
      rr = range(oo)
      res = c( sm$timestamp[rr[1]], sm$timestamp[rr[2]] )
    }
    return(res)
  }

  # trim tails a bit more:
  fa = range( which( sm$Z < (Z0$mode + bcp$depth.range[2]/2) & sm$Z > (Z0$mode + bcp$depth.range[1]/2) ) )
  if ( abs( diff( fa)) > 30 ) {
    Z0 = modes( sm$Z[ fa[1]: fa[2] ] )
  }

  aoi.i = which( sm$Z > Z0$lb2 &  sm$Z < Z0$ub2 )
  aoi = min(aoi.i) : max(aoi.i)
  aoi.n = length(aoi)
  aoi.trunc = trunc( aoi.n / 5 )
  aoi = aoi[ aoi.trunc] : aoi[ aoi.n-aoi.trunc ]  # trim off 1/5 off each tail

  aoi.inner = trunc( length(aoi) / 3 )

  left = min(aoi.i) + c( 0 : aoi.inner )
  right = max(aoi.i) + c( - aoi.inner : 0 )

  modal.threshold = 0.5
  mwin = -bcp$modal.windowsize: bcp$modal.windowsize

  #final pass using the initial estimates as a filter .. operating on residuals
  # create a variable with any linear trend in the depths removed as this can increase the precision of the method
  if (length( aoi ) > 0 ) {

    botlm = lm( Z ~ ts , data=sm[ left,], weights=Z^2, na.action="na.omit")
    sm$Zl = sm$Z - predict( botlm, newdata=sm )
    Zrmodes.left = modes( sm$Zl[ left ] )
    i.left = which( sm$Zl[left] > Zrmodes.left$lb2 & sm$Zl[left] < Zrmodes.left$ub2    )
    i.mod = left[i.left]
    ill = 5:trunc(median(i.mod) )

    if ( abs( min(sm$Zl[ill], na.rm=TRUE) - Zrmodes.left$mode ) <   bcp$modal.sd.multiplier* Zrmodes.left$sd ) {
      r0 = 1 # nothing to do .. assuming it is truncated without the left tail ..
    } else {
      for ( r0 in ill){
        ri = r0 + mwin
        ri = ri[ ri >= 1 & ri <=N ]
        nri = length(ri)
        if ( nri > 3 ) {
          nx = length( which( sm$Zl[ri] > Zrmodes.left$lb2 & sm$Zl[ri] < Zrmodes.left$ub2 )  )
          prx = nx / nri
          if (prx > modal.threshold) break()
        }
      }
    }

    botlm = lm( Z ~ ts , data=sm[ right,], weights=Z^2, na.action="na.omit")
    sm$Zr = sm$Z - predict( botlm, newdata=sm )
    Zrmodes.right = modes( sm$Zr[ right ] )
    i.right = which( sm$Zr[right] > Zrmodes.right$lb2 & sm$Zr[right] < Zrmodes.right$ub2 )
    i.mod = right[i.right]
    irr = (N-5):trunc(median(i.mod))

    if ( abs( min(sm$Zr[irr], na.rm=TRUE) - Zrmodes.right$mode ) <  bcp$modal.sd.multiplier * Zrmodes.right$sd ) {
      r1 = N # nothing to do ... assuming it is truncated without the right tail ..
    } else {

      for ( r1 in irr ) {
        ri = r1+mwin
        ri = ri[ ri >= 1 & ri <=N ]
        nri = length(ri)
        if ( nri > 3 ) {
          nx = length( which( sm$Zr[ri] > Zrmodes.right$lb2 & sm$Zr[ri] < Zrmodes.right$ub2 )  )
          prx = nx / nri
          if (prx > modal.threshold) break()
        }
      }
    }
  }

  res = list( bc0=sm$timestamp[r0], bc1=sm$timestamp[r1] )

  return(res)
}

