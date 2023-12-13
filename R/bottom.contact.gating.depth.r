
bottom.contact.gating.depth = function( Z, good, bcp ) {
  
  ## Preliminary gating: simple range limits (gating) of depths to remove real extremes
  # eliminate records with NA for depth  
  i = which(!is.finite(Z))
  if (length(i) > 0) good[i] = FALSE


  ##--------------------------------
  # eliminiate shallow records due to operation at sea level
  i = which(Z< bcp$depth.min )
  if (length(i) > 0) good[i] = FALSE



  zm = modes( Z[good] )
  zsd = zm$sd
  ir = which( Z < (zm$lb2 +bcp$depth.range[1]) | Z > (zm$ub2 + bcp$depth.range[2])  )  # a bit more inclusive as these are crude gates at an early stage
  if (length(ir) > 0) good[ir] =FALSE
  depth.bottom = zm$mode


  ##--------------------------------
  # eliminiate records that are shallower than a given percentage of the bottom depth
  i = which(  Z< ( bcp$depthproportion * depth.bottom ) )
  if (length(i) > 0) good[i] = FALSE
  # points( Z[good] ~ iz[good] , pch=20, col="green" )

  ##--------------------------------
  ## record filter
  depth.diff =  Z - depth.bottom 
  i = which( depth.diff < bcp$depth.range[1] | depth.diff> bcp$depth.range[2])
  if (length(i) > 0) good[i] = FALSE
  # points( Z[good] ~ iz[good] , pch=20, col="cyan" )

     
  return(good)

}


