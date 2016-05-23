
parse.marport.config = function( cfgfn ) {
  
  cfg = readLines( cfgfn )
  cfg = cfg[ -grep("^;.+$", cfg ) ]
  cfg = cfg[ grep("=", cfg ) ]

  res = as.data.frame( matrix( unlist(strsplit( cfg, split="=")), ncol=2, byrow=TRUE), stringsAsFactors=FALSE )
  names(res) = c( "sensorid", "variable" )

  res$units =NA
  res$units[ grep("ROL", res$variable) ] ="Roll angle: degrees"
  res$units[ grep("PIT", res$variable) ] ="Pitch angle: degrees"
  res$units[ grep("TMP", res$variable) ] ="Temperature: degrees Celcius"
  res$units[ grep("DPT", res$variable) ] ="Depth: meters"
  res$units[ grep("DST", res$variable) ] ="Distance: meters"
  res$units[ grep("HGT", res$variable) ] ="Clearance: meters"
  
  return(res)

}


