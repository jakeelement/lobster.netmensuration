
  serialbus.connect = function( params ) {
    #\\ connect to a serial bus for reading
    
    if (!exists( "OS", params)) params$OS="windows" 
    if (!exists( "baudrate", params)) params$baudrate=384000
    if (!exists( "mode", params)) params$mode="r"  # read only by default, see ?file Modes
   
    if (params$OS == "linux") {
      # see man stty for more option
      if (!exists( "port", params)) params$port = "/dev/ttyS0"
      if (!exists( "stopbit", params)) params$stopbit = "-cstopb" #-cstop = 1, cstop=2
      if (!exists( "databits", params)) params$databits = "cs7" # csN , where N is number of bits
      if (!exists( "parity", params)) params$parity = "parodd" # parodd is odd, -parodd is even
      connection.command = paste( "stty", "-F", params$port, params$baudrate, 
        "parenb", params$parity,  params$databits, params$stopbit )
    }

    if (params$OS == "windows" ) {
      if (!exists( "port", params)) params$port = "COM1"
      if (!exists( "stopbit", params)) params$stopbit = "1"
      if (!exists( "databits", params)) params$databits = "7" # csN , where N is number of bits
      if (!exists( "parity", params)) params$parity = "E" # parodd is odd, -parodd is even
      connect0 = paste( params$baudrate, params$parity, params$databits ,params$stopbit, sep=",")
      connection.command = paste( "mode", " ",  params$port, ":", connect0, sep="") 
    }
    system( connection.command )  # initiate connection
    con = file( params$port, open=params$mode )  
    return(con)
  }


