# netmensuration 

Utilities to interoperate with net mensuration data streams and determine touchdown/lift off and associated statistics. Requires ecomodUtils and ...


```
# to enable inter-operability with github
require( devtools ) # or install.packages( "devtools", dependencies=TRUE )

# this is to bootstrap the ecomod suite of tools
install_github( "jae0/ecomodUtils" ) 

# to use some of the functionality:
require( ecomodUtils ) # this should ideally be placed into your .Rprofile

# once the above is loaded, you can load standard R libraries with:
RLibrary( "mgcv", "sp", "nlme" ) 

# to load other ecomod-related packages from github: 
ecomodLibrary ( "ecomodUtils", "netmensuration" )   

# loadfunctions still operates properly for alternate locations (here the orgininal ecomod was moved to ecomod0) 
loadfunctions("netmensuration", alternate.directory="~/ecomod0")  # "/home/jae/ecomod0/netmensuration/src/_Rfunctions/"
loadfunctions("netmesuration")  # "/home/jae/ecomod/netmensuration/R"


# list of currently available ecomod packages on github:
ecomodLibraryList()

# The above list is hard coded into the function. If you have a project to add to ecomod, it would have to be updated in this file and below.

```

The currently available list of packages that make up ecomod include:

  * ecomodTemplate <https::/github.com/jae0/ecomodTemplate> (to use as base file to start a project quickly) 
  * ecomodUtils <https::/github.com/jae0/ecomodUtils> (formerly _ecomodSetup) 
  * snowcrab <https::/github.com/jae0/snowcrab>
  * groundfish <https::/github.com/jae0/groundfish> 
  * netmensuration <https::/github.com/jae0/netmensuration> 


This project preserves some of the original functionality of ecomod to source structured directories, except the default is now to source the directory: */R/ rather than */src/_Rfunctions/ . 


To set up the environment, modify your Rprofile to include:

```
  # initialize the ecomod environment
	ecomod.workdirectory = file.path( homedir, "tmp" )		 ### replace with correct path
	ecomod.directory = file.path( homedir, "ecomod" )   ### replace with correct path
	ecomod.datadirectory = file.path( homedir, "ecomod_data" )   ### replace with correct path
  pkgsInstalled = .packages(all.available = TRUE)
  if (!"ecomodUtils" %in% pkgsInstalled ) {
    if (!"devtools" %in% pkgsInstalled ) install.packages("devtools", dependencies=TRUE, ask=FALSE)
    require( devtools)
    install_github( "jae0/ecomodUtils")
  }
  require( ecomodUtils )

```


#### Useful links:

To make your own package, look at the structure of the ecomodUtils package. It is the mininimal set required for creating a package.

Conventions: naming of your package -- library name exists in the same namespace as regular R libraries and so you need to be careful about name conflicts. I suggest using ecomodXXX just to be consistent. 

#### Details on expected directory structure to interoperate with devtools::install_github() 

  http://r-pkgs.had.co.nz/description.html 

#### Other ways of installing devtools:

  https://github.com/hadley/devtools




