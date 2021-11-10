

# example of how to connect to a serial port to stream data from sensors:

con = serialbus.connect( list( port="COM1", baudrate=9600, parity="E", databits=8, stopbit=1 ) )

# to work with MSWindows:
# system( "mode COM1:4800,E,7,1" ) # look up the options for "mode"

con = file( "COM1", open="r" )
outfile = "test.txt"

if (isOpen(con)) {
  cat( paste( "Starting logging of serial bus::", date()), file=outfile )
  while (1) {
    res = scan( con, n=1, quiet=TRUE, what="character")
    cat( res, file=outfile, append=TRUE )
  }
}


