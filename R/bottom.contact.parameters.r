
  bottom.contact.parameters = function( bcp=list() ) {

    # basic defaults for any that have not yet been defined
    if ( !exists("id", bcp)) bcp$id="noid"
    if ( !exists("nr", bcp)) bcp$nr = NA

    if ( !exists("tdif.min", bcp)) bcp$tdif.min=8 # min time difference (minutes)
    if ( !exists("tdif.max", bcp)) bcp$tdif.max=52  # max time difference (minutes) .. including tails
    if ( !exists("depthproportion", bcp)) bcp$depthproportion=0.6  # depthproportion controls primary (coarse)gating
    if ( !exists("depth.min", bcp)) bcp$depth.min= 15
    if ( !exists("depth.range", bcp)) bcp$depth.range=c(-60,60)
    if ( !exists("time.gate", bcp)) bcp$time.gate=NA
    if ( !exists("eps.depth", bcp)) bcp$eps.depth = 2 # m

    if ( !exists("maxdepthchange", bcp)) bcp$maxdepthchange = 15 # max fluctuation in depth (m) between sensor pings

    # expected distance between GPS pings .. i.e. 10 cm to 0.1 km  --- might need to modify depending upon ping rate
    if ( !exists("gps.distance.range.valid.km", bcp)) bcp$gps.distance.range.valid.km =  c( 1e-4, 1e-1 )

    if ( !exists("noisefilter.trim", bcp)) bcp$noisefilter.trim = 0.1  # proportion of data to remove/trim based upon adjacent differences that are too extreme and local variance windows
    if ( !exists("noisefilter.var.window", bcp)) bcp$noisefilter.var.window = 9  # 1/2 of moving window used to compute local variance and moving window mean
    if ( !exists("noisefilter.inla.h", bcp)) bcp$noisefilter.inla.h = 0.05
    if ( !exists("noisefilter.inla.diagonal", bcp)) bcp$noisefilter.inla.diagonal = 0.05
    if ( !exists("noisefilter.inla.ngroups", bcp)) bcp$noisefilter.inla.ngroups = 500

    if ( !exists("noisefilter.quants", bcp)) bcp$noisefilter.quants = c(0.05, 0.95)
    if ( !exists("noisefilter.target.r2", bcp)) bcp$noisefilter.target.r2 = 0.9 # for noise filtering  .. ignore variations less than this threshold
    if ( !exists("noisefilter.sd.multiplier", bcp)) bcp$noisefilter.sd.multiplier = 5 # for noise filtering  .. ignore variations less than this threshold
    if ( !exists("noisefilter.postfilter.method", bcp)) bcp$noisefilter.smoother = "inla" # for noise filtering  .. ignore variations less than this threshold

    if ( !exists("smooth.filter.quants", bcp)) bcp$smooth.filter.quants=c(0.05, 0.95) # for dZ
    if ( !exists("smooth.windowsize",  bcp))  bcp$smooth.windowsize = 5  # number of data points to assess for consistent passage into the non-modal (slopes) area
    if ( !exists("smooth.sd.multiplier", bcp)) bcp$smooth.sd.multiplier = 1  # multiplier to the SD to id slopes that are low vs high
    if ( !exists("smooth.zeros", bcp)) bcp$smooth.zeros = 0.1  # multiplier to the SD to id slopes that should be ignored before smoothing

    if ( !exists("modal.sd.multiplier", bcp)) bcp$modal.sd.multiplier=2 # to detect if end point has been prefiltered/truncated
    if ( !exists("modal.trim", bcp)) bcp$modal.trim = 0.001  # keep low as it is really only for interpolation of NAs
    if ( !exists("modal.filter.quants", bcp)) bcp$modal.filter.quants = c(0.05, 0.95) # as above
    if ( !exists("modal.windowsize",  bcp))  bcp$modal.windowsize = 5  # number of data points to assess for consistent passage into the non-modal (depths) area

    if ( !exists("maxdepth.sd.multiplier", bcp)) bcp$maxdepth.sd.multiplier = 2 #

    if ( !exists("linear.sd.multiplier", bcp)) bcp$linear.sd.multiplier=2 # to detect if end point has been prefiltered/truncated
    if ( !exists("linear.trim", bcp)) bcp$linear.trim = 0.001  # keep low as it is really only for interpolation of NAs
    if ( !exists("linear.filter.quants", bcp)) bcp$linear.filter.quants = c(0.05, 0.95) # as above

    if ( !exists("user.interaction", bcp)) bcp$user.interaction=FALSE  # if you want to try to manually determine end points too
    if ( !exists("from.manual.archive", bcp)) bcp$from.manual.archive = project.datadirectory("bio.snowcrab", "data", "touchdown", "results")

    return(bcp)

  }


