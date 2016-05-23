

load.marport.rawdata = function( fnroot, fncfg, tzone="America/Halifax" ) {

  require(lubridate)
  if (FALSE) {
    fncfg=sensorconfig
    fnroot=fl
  }

# fn="/home/jae/Downloads/net/Alfred\ Needler-Ned_2013_002-045"
  # cfgfn = "/home/jae/Downloads/net/Config_Old.txt"

  sensorcodes = as.data.frame( matrix(
    c( "01PIT", "PitchDoorPort",
       "02PIT", "PitchDoorStar",
       "03PIT", "PitchWingPort",
       "01ROL", "RollDoorPort",
       "02ROL", "RollDoorStar",
       "03ROL", "RollWingPort",
       "02TMP", "TempDoorStar",
       "03TMP", "TempWingPort",
       "04TMP", "TempWingStar",
       "01DPT", "DepthDoorPort",
       "02DPT", "DepthDoorStar",
       "03DPT", "DepthWingPort",
       "04DPT", "DepthWingStar",
       "05DPT", "DepthScanmar",
       "01DST", "doorspread",
       "03DST", "wingspread",
       "03HGT", "opening",
       "06HGT", "opening.scanmar"
    ), ncol=2,  byrow=TRUE),  stringsAsFactors =FALSE )
  colnames(sensorcodes) = c("sensorcode", "variablename")

  cfg = parse.marport.config( fncfg )

  header = cfg[ which( is.na( cfg$units)) , ]
  bottom = cfg[ which( !is.na( cfg$units)), ]

  lookuptable = merge( sensorcodes, bottom, by.x="sensorcode", by.y="variable", all.x=TRUE, all.y=TRUE)
  names( header) = c( "sensorid", "sensorcode", "units")
  header$variablename=NA

  cfg = rbind( header, lookuptable[,names(header)] )

  fn.gps = paste(fnroot, "gps", sep=".")
  fn.sgp = paste(fnroot, "sgp", sep=".")
  fn.dlog =  paste(fnroot, "log", sep=".")

  if (!( file.exists(fn.gps) & file.exists(fn.sgp) & file.exists(fn.dlog) ) ) return(NULL)

  gps = read.table( fn.gps, sep=",", as.is=TRUE, header=FALSE )
  names(gps) = c("Vessel", "Cruise", "set", "timestamp", "latitude", "longitude" )

  sgp = read.table( fn.sgp, sep=",", as.is=TRUE, header=FALSE )
  names(sgp) = c("Vessel", "Cruise", "set",  "timestamp", "sensor", "value" )

  dlog = read.table( fn.dlog, sep=",", as.is=TRUE, header=FALSE )
  names(dlog) =  c("Vessel", "Cruise", "set", "timestamp", "event" )

  sgp = sgp[ , c( "timestamp", "sensor", "value" ) ]
  sgp$sensor = as.character( sgp$sensor)

  marport = gps
  vnames1 = names(gps)
  newvnames = cfg$variablename[ which(!is.na(cfg$variablename))]
  outputvnames = c(vnames1, newvnames )

  for (i in 1:nrow(cfg)) {
    if ( is.na( cfg$variablename[i])) next()
    if ( is.na( cfg$sensorid[i])) {
      marport[[cfg$variablename[i]]] = NA
      marport[[cfg$variablename[i]]] = as.numeric( marport[[cfg$variablename[i]]])
      next()
    }
    matchingsensordata = which( sgp$sensor == cfg$sensorid[i] )
    if (length(matchingsensordata)==0) {
      marport[[cfg$variablename[i]]] = NA
      marport[[cfg$variablename[i]]] = as.numeric( marport[[cfg$variablename[i]]])
      next()
    }
    sdata = NULL
    sdata = sgp[ matchingsensordata ,]
    newvariablename = cfg$variablename[i]
    names(sdata) = c( "timestamp", "sensor", newvariablename )
    marport = merge( marport, sdata[,c("timestamp",newvariablename)], by="timestamp", all.x=TRUE, all.y=FALSE )
  }

  marport = marport[, outputvnames]
  marportId = paste( marport$Vessel, marport$Cruise, marport$set )

  marport$timestamp = mdy_hms( marport$timestamp ) ## need to check if mdy or dmy ...
  dlog$timestamp = mdy_hms(dlog$timestamp)
  tz(marport$timestamp) = tzone
  tz(dlog$timestamp) = tzone

    threshold.seconds = 1*60*60 # 1hr in seconds
    nmids = unique( marportId )
    for (ii in nmids) {
      jj = which( marportId == ii)
      tstamp = marport$timestamp[jj]
      r = range(tstamp, na.rm=TRUE)
      y = as.numeric( difftime(r[2], r[1]), units="secs")  # in seconds

      if ( y > threshold.seconds ) { # as duration is in seconds
        # if there is a timestamp problem, the problematic records are those with hour values that are soon after midnight
        # .. assume any values from midnight to 2 AM need to be recoded to the next day's value
        hrs = hour( tstamp )
        z = which( hrs < 2 )  # 2AM is the cutoff
        if ( length(z) > 0 ) {
          day( tstamp[z]) = day( tstamp[z])+1
        }
        marport$timestamp[jj] = tstamp
      }
    }

  return(marport)
}

