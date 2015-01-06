#' Plot overlapping longitudinal hormone data 
#' 
#' @param x hormLong object (produced from hormBaseline) [required]
#' @param hormone_var name for the hormone variable.  It must be listed in by_var during hormBaseline [required]
#' @param colors list of colors for the fills. It needs to have as many colors as hormone types  [default='red,blue']
#' @param add_fill should the area under the curve be filled with color (TRUE or FALSE) [default=TRUE]
#' @param two_axes if TRUE then a second y-axis will be added to the right side of plot.  It can only be used with two hormones.
#' If there are more than two hormones, function will give an error. [default=FALSE]
#' @param date_format the format of the date variable on x-axis. Default is 01-Jan format. See Appendix 1 in help manual 
#' for examples of other formats [default = '\%d-\%b']
#' @param xscale  determines if x-axis is free ('free') to change for each panel or remain the same ('fixed') for all panels  [default = 'free']
#' @param yscale  determines if y-axis is free ('free') to change for each panel or remain the same ('fixed') for all panels [default = 'free']
#' @param plot_per_page the number of plot panels per page, by row. [default = 4]
#' @param plot_height  the height of individual plot panels (in inches).  Pdf page height is determined by both plot_per_page and plot_height. [default = 2]
#' @param plot_width  the width of the pdf page. [default = 6]
#' @param save_plot indicates whether to save plot as a file [default = TRUE]
#' @param ...   generic plotting options [optional]  
#' @return nothing  Produces a pdf file saved at current working directory
#' @export
#' @examples
#' 
#' result <- hormBaseline(data=hormLynx, criteria=2, by_var='AnimalID, Hormone', time_var='datetime', 
#'            conc_var='Conc', event_var='Events' )
#' hormPlotOverlap( result, hormone_var='Hormone', colors='red, dark green, blue' )
#' 
#'# no color fill...
#' hormPlotOverlap( result, hormone_var='Hormone', colors='red, dark green, blue', add_fill=FALSE )
#'
#'# Two-axes option is only vaalid if there are only two hormones 
#' hormPlotOverlap( result, hormone_var='Hormone', two_axes=T) # this will produce an error
#'
#'# Let's try again wtih hormElphant that only has two hormones
#' result <- hormBaseline(data=hormElephant, criteria=2, by_var='Ele, Hormone', time_var='Date', 
#'              conc_var='Conc_ng_ml', event_var='Event'  )
#' hormPlotOverlap( result, hormone_var='Hormone', colors='red, dark green', two_axes=TRUE )
#' 

hormPlotOverlap <- function(x, hormone_var, colors='red, blue', date_format='%d-%b', 
                     add_fill=TRUE, two_axes=FALSE,
                     xscale='free', yscale='free', 
                     plot_per_page=4, plot_height=2, plot_width=6, save_plot=TRUE,...){
  
#--- main check ---#
  graphics.off() # just to make sure no devices are open

  checkClass(x, 'hormLong')
  
  checkPlotOpts(plot_per_page, plot_width, plot_height, save_plot, xscale, yscale, date_format)
  
  checkPlotOverlap(hormone_var, x$by_var)

#-- set-up ---#
  colors <- cleanByvar(colors) # make a vector   
  by_var_v <- cleanByvar(x$by_var) 
  by_var_v <- by_var_v[by_var_v!=hormone_var]
  time_var <- x$time_var
  conc_var <- x$conc_var
  data <- x$data
  data <- data[ do.call(order, data[c(by_var_v,time_var,hormone_var)]), ]
  
  data$plot_title <- getPlotTitle(data, by_var=by_var_v)

  #--- check the two hormone scenario ---#
  if( two_axes ){ 
    if( length(unique(data[,hormone_var]))!=2 ){
      stop(paste0('Number of hormones must be two. Currently, there are ',length(unique(data[,hormone_var])) ) )
    }
  }
    
#--- create plots ---#
  if( save_plot ){
    pdf('hormPlotOverlap.pdf', height=plot_per_page * plot_height, width = plot_width )
  }

  if( two_axes ){par(mfrow=c(plot_per_page,1), mar=c(2,4,2,4),oma=c(2,2,2,2))
    }else{ par(mfrow=c(plot_per_page,1), mar=c(2,4,2,0.5),oma=c(2,2,2,0.5))}
  for( i in unique(data$plot_title) ){
    ds_sub <- data[data$plot_title==i, ]

    events <- getEventInfo(ds_sub, x$event_var, time_var)
    
    ds_sub <- ds_sub[!is.na(ds_sub[,conc_var]),]

    #--- set up scales ---#
      y_lim <- getPlotlim(d_s=ds_sub, d_f=data, var=conc_var, scale=yscale )
      x_lim <- getPlotlim(d_s=ds_sub, d_f=data, var=time_var, scale=xscale )
 
      if( two_axes & yscale=='free' ){
        y_lim[1] <- 0
        y_lim[2] <- max(ds_sub[ ds_sub[hormone_var]==sort(unique(ds_sub[,hormone_var]))[1],conc_var],na.rm=T)*1.2
        ymax2 <-  max(ds_sub[ ds_sub[hormone_var]==sort(unique(ds_sub[,hormone_var]))[2],conc_var],na.rm=T)*1.2
      }  
      if( two_axes & yscale=='fixed' ){
        y_lim[1] <- 0
        y_lim[2] <- max(data[ data[hormone_var]==sort(unique(data[,hormone_var]))[1],conc_var],na.rm=T)*1.2
        ymax2 <-  max(data[ data[hormone_var]==sort(unique(data[,hormone_var]))[2],conc_var],na.rm=T)*1.2
      }  
    
    #--- main plot
      require(lubridate)
      loop <- 0
      for(h in sort(unique(ds_sub[,hormone_var]))){
        loop <- loop + 1
        ds_sub1 <- ds_sub[ds_sub[,hormone_var]==h,]
        if( loop==1 ){
          y_lab = conc_var
          if( two_axes ){ y_lab <- h }
           plot(ds_sub1[,conc_var] ~ ds_sub1[,time_var], type='l',xlim=x_lim, ylim=y_lim, 
                        xlab=NA, ylab=y_lab, xaxt='n', col=colors[loop] )
        }
        if(loop>1 & two_axes==FALSE){lines(ds_sub1[,time_var],ds_sub1[,conc_var],col=colors[loop] )}
        if(loop>1 & two_axes){
          par( new=T )
          plot(ds_sub1[,conc_var] ~ ds_sub1[,time_var], type='l',xlim=x_lim, ylim=c(y_lim[1],ymax2), 
                    xlab=NA, ylab=NA, xaxt='n', yaxt='n', col=colors[loop] )
            axis(4)
            mtext(h,4,line=2.5)
        }
        t_order <- c(ds_sub1[,time_var], rev(ds_sub1[,time_var]) )
        c_order <- c(rep(0,length(ds_sub1[,conc_var])),rev(ds_sub1[,conc_var])) 
        if( add_fill ) polygon(t_order,c_order,col=adjustcolor(colors[loop], alpha=0.25),border=colors[loop]) 
        
        if( nrow(events) > 0 ){ 
          events1 <- events[events[,hormone_var]==h,  ] 
          plotEventInfo(events1,t=time_var, e=x$event_var) 
        }

      legend('topleft',legend=sort(unique(ds_sub[,hormone_var])),fill=adjustcolor(colors, alpha=0.25),
              bty='n',cex=0.9,bg=NA, pt.cex=0.6)

      plotAxes(ds_sub, time_var, x_lim, d_f=date_format)

    } # end sort(unique())
    mtext(unique(ds_sub$plot_title),side=3,line=0.25)
 
  } # end i unique(data$plot_title)
  if( save_plot ){  dev.off()
      cat( paste0('\n *********\nNote: plots are saved at: \n', getwd(),'/hormPlotOverlap.pdf \n***** \n\n')  )
  }
} # end of function
