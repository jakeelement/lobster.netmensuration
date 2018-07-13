
bottom.contact.plot = function ( O, netspread=FALSE ) {
  
  if ( exists("error.flag", O) && !is.na( O$error.flag ) )  {
    x = O$plotdata
    trange = range( x$ts, na.rm=TRUE )
    drange = range( x$depth, na.rm=TRUE) 
    #BC - Plots fail in RStudio graphics device, add clause
    dev.new(noRStudioGD = TRUE)
    plot(depth~ts, x, ylim=c(drange[2] +2.5, drange[1] -5), xlim=c(trange[1]-20,trange[2]+20), type="p", pch=".", col="black", cex=2 )
    title( sub=paste(O$id, O$error.flag) ) 
    return() 
  }

    x = O$plotdata

  # all data within range of gated limits
    dtrg = O$variance.method.indices
    if (length( which (is.finite( dtrg))) == 0 ) dtrg = O$good
    trange = range( x$ts[dtrg], na.rm=TRUE )
    drange = range( O$depth.smoothed[dtrg], na.rm=TRUE) 
    #BC - Plots fail in RStudio graphics device, add clause
    dev.new(noRStudioGD = TRUE)
    plot(depth~ts, x, ylim=c(drange[2] + 5, drange[1] -10 ), xlim=c(trange[1]-60,trange[2]+60), type="n", xlab="" )
    legendtext = NULL
    legendcol = NULL
    legendpch = NULL
    points( depth~ts, x[O$good,], pch=21, col="slategray", cex=0.5)
    points( depth~ts, x[!O$good,], pch=21, col="red", cex=0.4)


    if (all(is.finite(c( O$variance.method0, O$variance.method1) ) ) ) {
        mcol = "gray"
        # points( depth~ts, x[ O$variance.method.indices, ], pch=20, col=mcol, cex=0.2)
        abline (v=x$ts[min(O$variance.method.indices)], col=mcol, lty="solid", lwd=1.2)
        abline (v=x$ts[max(O$variance.method.indices)], col=mcol, lty="solid", lwd=1.2)
        DT = as.numeric( difftime( O$variance.method1, O$variance.method0, units="mins" ) )
        legendtext = c( legendtext, paste( "variance gate:   ", round( DT, 2), "" ) )
        legendcol = c( legendcol, mcol)
        legendpch =c( legendpch, 20 ) 
    }
    

    if (all(is.finite( c(O$modal.method0, O$modal.method1 ) ) ) ) {
        mcol = "red" # colour for plotting
        # points( depth~ts, x[O$modal.method.indices,], col=mcol, pch=20, cex=0.2)       
        abline (v=x$ts[min(O$modal.method.indices)], col=mcol, lty="dashed")
        abline (v=x$ts[max(O$modal.method.indices)], col=mcol, lty="dashed")
        DT = as.numeric( difftime( O$modal.method1, O$modal.method0, units="mins" ) )
        legendtext = c( legendtext, paste( "modal:   ", round( DT, 2) ) )
        legendcol = c( legendcol, mcol)
        legendpch =c( legendpch, 20) 
    }

    
    if ( all(is.finite( c( O$smooth.method0,  O$smooth.method1) ) ) ) {
          mcol = "blue"
          # points( depth~ts, x[O$smooth.method.indices,], col=mcol, pch=20, cex=0.2)   
          abline (v=x$ts[min(O$smooth.method.indices)], col=mcol, lty="dotdash", lwd=1.5)
          abline (v=x$ts[max(O$smooth.method.indices)], col=mcol, lty="dotdash", lwd=1.5)
          DT = as.numeric( difftime( O$smooth.method1, O$smooth.method0, units="mins" ) )
          legendtext = c(legendtext, paste( "smooth:   ", round(DT, 2)) )
          legendcol = c( legendcol, mcol)
          legendpch =c( legendpch, 20) 
    }
  
    if (all(is.finite( c(O$linear.method0, O$linear.method1)  )) ) {
          mcol ="green"
          # points( depth~ts, x[O$linear.method.indices,], col=mcol, pch=20, cex=0.2)      
          abline (v=x$ts[min(O$linear.method.indices)], col=mcol, lty="twodash")
          abline (v=x$ts[max(O$linear.method.indices)], col=mcol, lty="twodash")
          DT = as.numeric( difftime( O$linear.method1, O$linear.method0, units="mins" ) )
          legendtext = c( legendtext, paste( "linear: ", round( DT, 2) ) )
          legendcol = c( legendcol, mcol)
          legendpch =c( legendpch, 20) 
     }


    if (all(is.finite( c( O$maxdepth.method0,  O$maxdepth.method1) )) ) {
      mcol ="orange"
          # points( depth~ts, x[O$linear.method.indices,], col=mcol, pch=20, cex=0.2)      
          abline (v=x$ts[min(O$maxdepth.method.indices)], col=mcol, lty="solid")
          abline (v=x$ts[max(O$maxdepth.method.indices)], col=mcol, lty="solid")
          DT = as.numeric( difftime( O$maxdepth.method1, O$maxdepth.method0, units="mins" ) )
          legendtext = c( legendtext, paste( "maxdepth: ", round( DT, 2) ) )
          legendcol = c( legendcol, mcol)
          legendpch =c( legendpch, 20) 
     }

    if (all(is.finite( c(O$bottom0, O$bottom1) )) ) {
      mcol ="purple"
        td0 = abs( difftime( O$bottom0, x$timestamp)) 
        td1 = abs( difftime( O$bottom1, x$timestamp) )
        ii0  =which.min( td0)
        ii1  =which.min( td1)
          abline (v=x$ts[ ii0], col=mcol, lty="solid", lwd=1.6)
          abline (v=x$ts[ ii1], col=mcol, lty="solid", lwd=1.6)
          DT = as.numeric( difftime( O$bottom1, O$bottom0, units="mins" ) )
          legendtext = c( legendtext, paste( "Trimmed mean: ", round( DT, 2) ) )
          legendcol = c( legendcol, mcol)
          legendpch =c( legendpch, 20) 
     }


     if (all(is.finite( c(O$manual.method0, O$manual.method1 ) )) ) {
        mcol ="cyan"
        # points( depth~ts, x[O$manual.method.indices,], col=mcol, pch=20, cex=0.2)      
        abline (v=x$ts[min(O$manual.method.indices)], col=mcol, lty="dotdash")
        abline (v=x$ts[max(O$manual.method.indices)], col=mcol, lty="dotdash")
        tdif = abs( as.numeric(difftime(O$manual.method1, O$manual.method0 , units="mins")) )
        tdiffman = round( tdif, 2)
        legendtext = c( legendtext, paste( "manual: ", tdiffman  ) ) 
        legendcol = c( legendcol, mcol)
        legendpch =c( legendpch, 20) 
     }
   
    if ( !( is.null( legendtext)))  legend( "top", legend=legendtext, col=legendcol, pch=legendpch )
     
    #plot.new(); plot( slopes ~ ts, x2 )
    lines( O$depth.smoothed ~ x$ts, col="brown" )
    # points( depth0~ts, x[!O$good,], col="red", cex=1 )   ## points dropped from filters
    
    title( sub=paste( O$id,  
      " :: time mean=", signif( as.numeric(O$bottom.diff)/60, 3), 
      "; time sd=", signif( as.numeric(O$bottom.diff.sd)/60, 3) , 
      "; n=", O$depth.n,  
      ";", "\n", " depth mean=", signif( O$depth.mean, 3), 
      "; depth sd=", signif( O$depth.sd, 3), 
      "; signal=", signif(O$signal2noise, 3) , sep="")  )

    if (netspread & !is.na(O$surface.area) ) {
      #BC - Plots fail in RStudio graphics device, add clause
      dev.new(noRStudioGD = TRUE)
      par(mfrow=c(2,2))
      bts = O$ts[ O$bottom.contact ] 

      dr = range( c( O$plotdata$wingspread[ O$bottom.contact ] , O$surface.area$wing.smoothed, O$surface.area$wing.mean - 2*O$surface.area$wing.sd,  O$surface.area$wing.mean + 2*O$surface.area$wing.sd ) , na.rm=TRUE )
      plot(  O$plotdata$wingspread[ O$bottom.contact ] ~ bts, pch=20, col="lightgray", cex=0.5, ylim=dr, xlab="", ylab="wingspread, m" ) 
      abline( h=O$surface.area$wing.mean, lwd=2, col="gray") 
      abline( h=O$surface.area$wing.mean - 2*O$surface.area$wing.sd, lty="dotted", col="lightgray" ) 
      abline( h=O$surface.area$wing.mean + 2*O$surface.area$wing.sd, lty="dotted", col="lightgray" ) 
      lines( O$surface.area$wing.smoothed ~ bts, col="red" )
      legend( "bottom", legend=c( paste("Mean:", round(O$surface.area$wing.mean,2), "\n", 
                                        "SD:", round(O$surface.area$wing.sd,2) ))) 

      dr = range( c( O$plotdata$doorspread[ O$bottom.contact ] , O$surface.area$door.smoothed, O$surface.area$door.mean - 2*O$surface.area$door.sd,  O$surface.area$door.mean + 2*O$surface.area$door.sd ) , na.rm=TRUE )
      plot(  O$plotdata$doorspread[ O$bottom.contact ] ~ bts, pch=20, col="lightgray", cex=0.5, ylim=dr, xlab="", ylab="doorspread, m" ) 
      abline( h=O$surface.area$door.mean, lwd=2, col="gray") 
      abline( h=O$surface.area$door.mean - 2*O$surface.area$door.sd, lty="dotted", col="lightgray" ) 
      abline( h=O$surface.area$door.mean + 2*O$surface.area$door.sd, lty="dotted", col="lightgray" ) 
      lines( O$surface.area$door.smoothed ~ bts, col="blue" )
      legend( "bottom", legend=c( paste("Mean:", round(O$surface.area$door.mean,2), "\n", 
                                        "SD:", round(O$surface.area$door.sd,2) ))) 

   
      plot( O$surface.area$distances.total ~ bts, type="l", col="red", xlab="", ylab="Cummulative distance (km)" ) 
      lines( O$surface.area$distances.vertical ~ bts,  col="orange" ) 
      lines( O$surface.area$distances.horizontal ~ bts,  col="cyan" ) 
      legend( "topleft", legend=c( paste("Total:", max( round(O$surface.area$distances.total,2) ), "\n", 
                                         "Horizontal:",  max( round(O$surface.area$distances.horizontal,2)), "\n",
                                         "Vertical:",  max( round(O$surface.area$distances.vertical,2) )))) 


      plot( c( O$surface.area$door.sa.cum, O$surface.area$wing.sa.cum ) ~ c( bts, bts) , type="n", xlab="", ylab="SA (km^2)" ) 
      lines( O$surface.area$door.sa.cum ~ bts,  col="blue" ) 
      lines( O$surface.area$wing.sa.cum ~ bts,  col="red" ) 
      legend( "topleft", legend=c( paste("Door:", max( round(O$surface.area$door.sa.cum, 4) ), "\n", 
                                         "Wing:",  max( round(O$surface.area$wing.sa.cum,4) )))) 
    
    }
}


