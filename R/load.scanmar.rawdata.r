
load.scanmar.rawdata = function( fn, yr=NULL ) {
  # Scanmar is always UTC!!!

  scanmar=NULL
  header = readLines(fn, n=10, encoding="UTF-8", skipNul=TRUE)
  datestring = basename(fn)

  ds = unlist( strsplit( datestring, "-") )

  header = iconv(enc2utf8(header),sub="byte")  # fix mutlibyte characters

  if (nchar(ds[1]) == 4 ) {
    yr0 = as.numeric( as.character( ds[1] ) )
    if ( yr0 >= 1990 & yr0 <= 2014 ) {
      # most likely year .. from data prior to 2015
      # for 1990 to 2014 .. naming method is consistent
      mon = substring( ds[2], 1, 3 )
      day = gsub( mon, "", ds[2] )
    } else {
      print( fn)
      stop( "Error reading data: unexpected filename encountered" )
    }
  } else {
    # more info added to name in 2015 ... altered parsing
    yr0 = as.numeric( unlist( strsplit( ds[1], "\\." ) )[3])
    if ( yr0 >= 2014 ) {
      mon = substring( ds[2], 1, 3 )
      day = gsub( mon, "", ds[2] )

    } else {
      print( fn)
      print(header)
      stop( "Error reading data: unexpected filename/data encountered" )
    }
  }

  if ( !is.null(yr) ) {
    if ( yr != yr0 ) {
      print (fn)
      stop( "Error in parsing year" )
    }
  }

  # extract timestamp of "Start Set" --- these are not always reliable .. ignore
  # line.start = grep("Start Set", header, ignore.case=T  )

  # timestring = header[ line.start]
  # time = substring(timestring, 1, 8)

  if (length(header) < 8 ) return( NULL )

  tmpfile = file.path (getwd(), "tmp.scanmar")
  # remove "*" -- hard to do internally without fileswapping out and reloading
  tmp = readLines(fn, encoding="UTF-8", skipNul=TRUE)
  tmp = tmp[-c(1:5)] # headers sometimes long
  tmp = iconv(enc2utf8(tmp),sub="byte")  # fix mutlibyte characters

  u = nchar(tmp)
  v = which(u>40)
  if (length(v)<100) {return(NULL)}

  tmp = tmp[v]

  write( tmp, file=tmpfile )

  scanmar = read.table( file=tmpfile, sep="", as.is=T, colClasses="character", header=F)
  file.remove(tmpfile)
  nr0 = nrow(scanmar)
  if (nr0 < 30 ) return(NULL)

  # nc0 = ncol(scanmar)
  # if ( nc0 < 12 ) return ( out ) # no net metrics stored ( only position and speed )  -- ignored
  # if ( nc0 > 14 ) stop( fn)  # should not be the case unless new data streams start
  #if ( nc0 == 12) {
    # older files do not have depths, add a dummy column
    # more modern data have 13 columns 2000 +
   # scanmar$depth= NA
  #}

  colnames(scanmar) = c("time", "depth", "doorspread", "wingspread", "opening", "clearance",
                        "ltspeed", "gyro", "latitude", "longitude")

  numbers = c("depth", "doorspread", "wingspread", "opening", "clearance",
              "ltspeed", "gyro", "latitude", "longitude")

  scanmar = factor2number(scanmar, numbers)

  scanmar$doorspread = filter.nets("doorspread.range", scanmar$doorspread)
  scanmar$wingspread = filter.nets("wingspread.range", scanmar$wingspread)
  scanmar$clearance = filter.nets("clearance.range", scanmar$clearance)
  scanmar$opening = filter.nets("opening.range", scanmar$opening)
  scanmar$depth = filter.nets("depth.range", scanmar$depth)

  scanmar$timestamp= paste(yr0, mon, day, scanmar$time, sep="-" )
  scanmar$timestamp=gsub(":","-",scanmar$timestamp)
  scanmar$timestamp = lubridate::ymd_hms(scanmar$timestamp, tz="UTC" ) # Scanmar is always UTC!!! .. which is the default of ymd_hms

  scanmar$nm_id = basename(fn)

    threshold.seconds = 1*60*60 # 1hr in seconds
    nmids = unique( scanmar$nm_id )
    for (ii in nmids) { # in case there are multiple sets
      jj = which( scanmar$nm_id == ii)
      tstamp = scanmar$timestamp[jj]
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
        scanmar$timestamp[jj] = tstamp
      }
    }

  scanmar = look.for.multiple.sets ( scanmar )
  scanmar$id = NA  # filled in below if from 2015 and onwards, otherwise this will be filled in later once the position/time is matched to gsinf (or, in perley .. Perley.mission.set )

  if (yr0 >= 2015) {
    hdr = gsub( "<", "", header[1] )  # when multibyte chars exist they are surrounded by "<>" .. remove
    hdr = gsub( ">", "", hdr )
    h1 = unlist( strsplit( hdr, "[[:punct:]]" ))
    m1 = nchar(h1[1])
    mission = substring( h1[1], m1-9, m1 )

    set =  as.character( h1[2] )
    header.id.fn = paste( mission, set, sep=".")
    header.id = paste( mission, as.character(as.numeric(set)), sep=".")

    f0 = unlist( strsplit( basename( fn), "[[:punct:]]" ))
    filename.id = paste( f0[1], f0[2], sep=".")

    scanmar$id = header.id

    if (header.id.fn != filename.id){
      print( paste( header.id.fn, " -- ", filename.id) )
      print( "Header and filename do not indicate the same mission & set" )
      print( "Using header id over filename ... please check the file." )
    }
  }

  return(scanmar)

}



