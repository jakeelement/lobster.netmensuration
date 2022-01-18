convert.marport.sds_csv2netmind = function( fnroot, yr, redo.marport_conversion = F ) {
  fn_marport_proc = file.path(project.datadirectory(),"bio.snowcrab", "data", "marport", paste("marport_proc_", yr, ".RDATA", sep = ""))
  if(redo.marport_conversion){
  
    fnroot = file.path("S:", "Survey", "Annual Files by Year", paste("ENS Snow Crab ",yr, "Survey", sep = ""), "completed", paste("Snow Crab Survey ",yr, sep = ""), "SCFINAL.csv")
    con = file(fnroot, "r")
    line = readLines(con, n = -1)
    ind = grep("GPRMC", line)
    inda = grep("192.168", line) 
    indb = grep("RAW", line)  
    indc = intersect(inda, indb)
    indall = sort(c(indc, ind))
    lines = line[indall]
    close(con)
    depth = ""
    wing = ""
    height = ""
    pitch = ""
    roll = ""
    temp = ""
    
    marport = data.table::data.table()
    for(i in 1:length(lines)) {
      print(i)
      lin = lines[i]
      lin = gsub(";", ",", lin)
      
      
      if(grepl("GPRMC", lin)){
        gps = unlist(strsplit(lin, ","))[c(1, 3, 5, 7, 9, 11)]
        marport = rbind(marport, cbind(gps[1], gps[6], gps[2], gps[3], gps[4], gps[5], height, NA, wing, depth, temp))
        depth = NA
        wing = NA
        height = NA
        pitch = NA
        roll = NA
        temp = NA
      }
      if(grepl("DEPTH", lin)){
        depth = unlist(strsplit(lin, ","))[5]
      }
      if(grepl("DISTANCE", lin)){
        wing = unlist(strsplit(lin, ","))[5]
      }
      if(grepl("ROLL", lin)){
        roll = unlist(strsplit(lin, ","))[5]
      }
      if(grepl("PITCH", lin)){
        pitch = unlist(strsplit(lin, ","))[5]
      }
      if(grepl("TEMPERATURE", lin)){
        temp = unlist(strsplit(lin, ","))[5]
      }
      if(grepl("SENSORDTB", lin)){
        height = unlist(strsplit(lin, ","))[5]
      }
      
    }
    
    names(marport) = c("localtime", "Date",	"Time",	"Latitude",	"Longitude",	"Speed",	"Primary",	"Secondary",	"DoorSpread",	"Depth",	"Temperature")
    
    marport$localtime = lubridate::as_datetime(as.numeric(marport$localtime)/86400000 + lubridate::as_date("1970/01/01"), tz = "UTC")
    marport$localtime = lubridate::with_tz(marport$localtime, "America/Halifax")
    marport = tidyr::separate(data = marport, col = "Latitude", into = c('latdeg', 'latmin'), sep = 2, remove = F)
    marport = tidyr::separate(data = marport, col = "Longitude", into = c('londeg', 'lonmin'), sep = 3, remove = F)
    marport$lat = as.numeric(marport$latdeg) + as.numeric(marport$latmin)/60
    marport$lon = (as.numeric(marport$londeg) + as.numeric(marport$lonmin)/60)*-1
    
   
    save(marport, file =  fn_marport_proc)
  }
  load(fn_marport_proc)
  
  fn_netmind_arc = file.path(project.datadirectory(),"bio.snowcrab", "data", "netmind", "archive", yr)
  if(!dir.exists(fn_netmind_arc)){ 
    dir.create(fn_netmind_arc)
  }
  fl = list.files(fn_netmind_arc)
  set = snowcrab.db( DS="setInitial" ) 
  set = set[which(set$yr == yr),]
  
  set$station[which(nchar(set$station) == 1)] = paste("00", set$station[which(nchar(set$station) == 1)], sep = "")
  set$station[which(nchar(set$station) == 2)] = paste("0", set$station[which(nchar(set$station) == 2)], sep = "")
  missing.set = set[which(!paste("ep", set$station, ".txt", sep = "") %in% fl),]
  missing.set = missing.set[order(missing.set$timestamp),]
  
  
  for(i in 1:nrow(missing.set)){
    curset = missing.set[i,]
     if(curset$trip == "S06102021" && curset$station == "130") curset$station = "305"
    sind = which(geosphere::distm(c(curset$lon, curset$lat),cbind(marport$lon, marport$lat), fun = distHaversine) < 2000)
    if(length(sind) == 0){
      warning(paste("Could Not find Marport data for station: ", curset$station, sep = ""))
      next()
    }
    start = min(sind)
    end = max(which(geosphere::distm(c(curset$lon1, curset$lat1),cbind(marport$lon, marport$lat), fun = distHaversine) < 2000))
    marsub = marport[start:end,]
    
    
    range = which((abs(dmy_hms(paste(marsub$Date, marsub$Time, sep = " ")) - curset$timestamp)) < minutes(15))
    marsub = marsub[range,]
    
    
    if(nrow(marsub) < 10){
      warning(paste("Could Not find Marport data for station: ", curset$station, sep = ""))
      next()
    }

    header =  paste("FileName: ",yr, "SnowCrabSurvey ep", curset$station, ".txt \nLocal Time: ", format(marsub$localtime[1], format = "%a %b %d  %H:%M:%S %Y"),"\nShip:  Trip: ", curset$trip, "  Tow:  ", curset$set, "\nComments: ", sep = "")
    
    marsub$Latitude = paste(marsub$latdeg, marsub$latmin, "N", sep = " ")
    marsub$Longitude = paste(marsub$londeg, marsub$lonmin, "W", sep = " ")
    marsub$Date = format(dmy(marsub$Date), format = "%y%m%d")
    
    marsub$localtime = NULL
    marsub$latdeg = NULL 
    marsub$latmin = NULL
    marsub$londeg = NULL      
    marsub$lonmin = NULL
    marsub$lat = NULL      
    marsub$lon = NULL
    
    file = file.path(project.datadirectory(),"bio.snowcrab", "data", "netmind", "archive", yr, paste("ep", curset$station, ".txt", sep = ""))
    
    cat(header, '\n',  file = file)
    suppressWarnings(write.table(marsub, file, append = T, quote = F, row.names = F, col.names = T, sep = "  "))
        
  }

  

}

load.marport.rawdata = function( fnroot, fncfg ) {

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

  tzone="America/Halifax"
  marport$timestamp = mdy_hms( marport$timestamp, tz=tzone ) ## need to check if mdy or dmy ...
  dlog$timestamp = mdy_hms(dlog$timestamp, tz=tzone )

  # convert to internal TZ
  marport$timestamp = with_tz(marport$timestamp, "UTC")
  dlog$timestamp = with_tz(dlog$timestamp, "UTC" )

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

