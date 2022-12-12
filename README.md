# netmensuration

Utilities to interoperate with net mensuration data streams and determine touchdown/lift off and associated statistics.

Installation:

```
install.packages( "devtools", ask=F, dependencies=TRUE )
require( "devtools")
install_github( "jakeelement/netmensuration" )
```

External dependency: INLA (http://www.r-inla.org)  OR use this method: 
remotes::install_version("INLA", version="20.06.29",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)  ##change version number to get version of INLA built for whatever R version you're running (version list here: https://groups.google.com/g/r-inla-discussion-group/c/c0__KycX1j4)


This version of netmensuration works for lobster ILTS data. For original package written by Jae Choi, see jae0/netmensuration
