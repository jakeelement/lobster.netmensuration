
#// template for using net mensuration methods

### -----------------------------------------------------------------------
### -- example usage directly upon a single file, perhaps on board/in situ
### -----------------------------------------------------------------------

# file list can be created and then a gui can be
# created to use mouse to choose file and run etc ... be my guest .. :)

# set up libs and bio functions
p = aegis.survey::groundfish_parameters()

datadir = file.path( p$scanmar.dir, "datalogs", "2015", "NED2015002" )  # storage location
fn = "NED2015002.028.2015-Mar21-162742.SET.LOG" # filename of data file to examine
fl = file.path( datadir, fn)
rawdata = load.scanmar.rawdata( fl )
bcp = list( id=fn, nr=nrow(rawdata), tdif.min=9, tdif.max=45, user.interaction=TRUE )  ### yes some are as short as 9 min .. turn user interaction off if you would like automatic solutions only
bcp = bottom.contact.parameters( bcp ) # add other default parameters
bc =  bottom.contact(rawdata, bcp )
bottom.contact.plot( bc, netspread=TRUE )
str(bc)


#### -----------------------------------------------------------------------
### -- example usage directly upon a a directory of scanmar data
### -----------------------------------------------------------------------


p = aegis.survey::groundfish_parameters()
datadir = file.path( p$scanmar.dir, "datalogs", "2014" )  # storage location for raw data logs

fn.base = "basedata.rdata"  # storage of outputs
fn.meta = "metadata.rdata"

# names of variables in logs
varnames = c( "id", "nm_id", "ltspeed", "ctspeed", "wingspread", "doorspread", "clearance", "opening",
                  "latitude", "longitude", "depth", "gyro", "timestamp")
filelist = list.files(path=datadir, pattern="\\.log$", full.names=T, recursive=TRUE, ignore.case=TRUE)
basedata = NULL
metadata = NULL
for ( fl in filelist ) {
  print(fl)
  rawdata = load.scanmar.rawdata( fl )  # read file one at a time
  if (!is.null(rawdata)) {
    rawdata$ctspeed = NA  # missing in modern data so add to permit merging with historical data
    rawdata = rawdata[ , varnames ]
    bcp = list( id=fn, nr=nrow(rawdata), tdif.min=9, tdif.max=45, user.interaction=TRUE )
    bcp = bottom.contact.parameters( bcp ) # add other default parameters
    bc =  bottom.contact(rawdata, bcp )
    bottom.contact.plot( bc, netspread=TRUE )
    str(bc)
    stats = as.data.frame( cbind(
      id=bc$id, bc$res,
      t0=bc$summary[ "means" ,"start"] ,  # choose "manual" instead of "means" if you want manually id'ed start/end times
      t1=bc$summary[ "means" ,"end"]
    ))
    basedata = rbind( basedata, rawdata)
    metadata = rbind( metadata, stats )
  }
}
save(metadata, file=fn.meta, compress=TRUE)
save(basedata, file=fn.base, compress= TRUE)

