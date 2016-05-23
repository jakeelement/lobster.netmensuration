

look.for.multiple.sets = function( x,  tdiff.min=20 ) {
  
  # some logs contain multiple sets in one file
  # break them down into separate groups

  if (0) {
    tzone="UTC"
    yr=2009
    fn="/home/jae/ecomod/groundfish/data/nets/Scanmar/datalogs/2009/2009-Aug20-024348.SET.LOG"
    x = load.scanmar.rawdata( fn, yr=yr)
     threshold.depth=20  
    tdiff.min=20 
  }
  
  ii = which( x$depth > 20 )
  if (length(ii) > 30) {
    zm = modes( x$depth[ii ] ) #  
    threshold.depth = zm$mode / 3  # retain bottom 2/3 of data 
  } else {
    threshold.depth = 20 
  }

  shallow = which( x$depth < threshold.depth )
  if (length(shallow)>0) x$depth[shallow] = NA
  x = x[order( x$timestamp ) ,]
  fn0 = unique(x$nm_id)
  set = list()
  startdata = which( is.finite(x$depth) )
  if (length(startdata) < 30) return(x) # nothing to do
  
  ir = range( startdata)  # start and ending index
  nr = nrow(x)

  active = rep(0, nr)
  active[ which( is.finite(x$depth) ) ] = 1
  
  da = rep( 0, nr )
  da[1:(nr-1)] = diff( active)

  i0 = which( da ==  1 )  # starting 
  i1 = which( da == -1 ) # ending

  tdiff0 = diff( x$timestamp[ i0 ] ) / 60  # minutes
  tdiff1 = diff( x$timestamp[ i1 ] ) / 60  # minutes
  tdiffs = difftime( x$timestamp[ i1 ], x$timestamp[ i0 ], units="mins" )

  i = which( tdiffs > tdiff.min )
  if (length(i) %in% c( 0,1) ) {
    # nothing to do if there is no data or if a single tow 
    # plot( depth~timestamp, x)
    return(x)  # do nothing
  }
  
  print( "Multiple tows in one file, breaking into separate data streams and adding a numeric suffix to file name ...")

  tn = 0
  for ( ti in i) {
    tn = tn + 1
    j = which( x$timestamp >= x$timestamp[i0[ti]] &  x$timestamp <= x$timestamp[i1[ti]] )
    fnj =  paste(fn0, tn, sep=".")
    x$nm_id[j] = fnj
    print( fnj )
  }
 
  todrop = which(  x$nm_id == fn0 )
  if (length(todrop) > 0 ) x=x[-todrop,]

  return(x)
}



