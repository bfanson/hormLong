#' Plot overlapping longitudinal hormone data 
#' 
#' @param x hormLong object (produced from hormBaseline) [required]
#' @param hormone_var name for the hormone variable.  It must be listed in by_var during hormBaseline [required]
#' @param colors list of colors for the fills. It needs to have as many colors as hormone types  [default='red,blue']
#' @param add_fill should the area under the curve be filled with color (TRUE or FALSE) [default=TRUE]
#' @param two_axes if TRUE then second axis will added to right side of plot.  Can only with two hormones
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
#' result <- hormBaseline(data=hormLynx, criteria=2, by_var='AnimalID, Hormone', time_var='datetime', conc_var='Conc' )
#' hormPlotOverlap( result, hormone_var='Hormone', colors='red, dark green, blue' )

#'# no color fill...
#' hormPlotOverlap( result, hormone_var='Hormone', colors='red, dark green, blue', add_fill=FALSE )
#'
#'# Two-axes option is only vaalid if there are only two hormones 
#' hormPlotOverlap( result, hormone_var='Hormone', two_axes=T ) # this will produce an error
#'
#'# Let's try again wtih hormElphant that only has two hormones
#' result <- hormBaseline(data=hormElephant, criteria=2, by_var='Ele, Hormone', time_var='Date', 
#'              conc_var='Conc_ng_ml' )
#' hormPlotOverlap( result, hormone_var='Hormone', colors='red, dark green', two_axes=TRUE )
#' 

hormPlotOverlap <- function(x, hormone_var='horm_type', colors='red, blue', date_format='%d-%b', 
                     add_fill=TRUE, two_axes=FALSE,
                     xscale='free',yscale='free', 
                     plot_per_page=4, plot_height=2, plot_width=6, save_plot=TRUE,...){
  
#--- main check ---#
  graphics.off() # just to make sure no devices are open

  if( class(x)!='hormLong'){
      stop('Object needs to be hormLong.  Run hormBaseline() first')
  }
  if( !is.numeric(plot_per_page) | plot_per_page<1 ){
    stop('plot_per_page needs to be numeric and greater than 0')
  }
  plot_per_page <- as.integer(plot_per_page ) # make sure whole number
    
  if( !is.logical(save_plot)  ){
    stop( paste('save_plot needs to be either TRUE or FALSE, not:',save_plot) )
  }

  if( !(yscale %in% c('free','fixed') ) ){
    stop( paste('yscale must be either "free" or "fixed", not:',yscale) )
  }
  if( !(xscale %in% c('free','fixed') ) ){
    stop( paste('xscale must be either "free" or "fixed", not:',xscale) )
  }
  colors <- cleanByvar(colors) # make a vector   
  
#-- set-up ---#
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
    if(!is.null(x$event_var)){
      events <- unique( ds_sub[ !is.na(ds_sub[,x$event_var]) & ds_sub[,x$event_var]!='',c(x$event_var,time_var)])
      }else{events <- data.frame()}
    ds_sub <- ds_sub[!is.na(ds_sub[,conc_var]),]

    #--- set up scales
      if( yscale=='free'){
        ymin <- min(ds_sub[,conc_var])
        ymax <- max(0, max(ds_sub[,conc_var]) )*1.1
      }else{
        ymin <- min(data[,conc_var],na.rm=T)
        ymax <- max(data[,conc_var],na.rm=T)*1.1
      }
      if( xscale=='free'){
        xmin <- min( ds_sub[,time_var],na.rm=T )
        xmax <- max( ds_sub[,time_var],na.rm=T )
      }else{
        xmin <- min( data[,time_var],na.rm=T)
        xmax <- max( data[,time_var],na.rm=T)
      }    
      if( two_axes & yscale=='free' ){
        ymin <- 0
        ymax  <-  max(ds_sub[ ds_sub[hormone_var]==sort(unique(ds_sub[,hormone_var]))[1],conc_var],na.rm=T)*1.1
        ymax2 <-  max(ds_sub[ ds_sub[hormone_var]==sort(unique(ds_sub[,hormone_var]))[2],conc_var],na.rm=T)*1.1
      }  
      if( two_axes & yscale=='fixed' ){
        ymin <- 0
        ymax  <-  max(data[ data[hormone_var]==sort(unique(data[,hormone_var]))[1],conc_var],na.rm=T)*1.1
        ymax2 <-  max(data[ data[hormone_var]==sort(unique(data[,hormone_var]))[2],conc_var],na.rm=T)*1.1
      }  
    
    #--- main plot
      require(lubridate)
      loop <- 0
      for(h in sort(unique(ds_sub[,hormone_var]))){
        loop <- loop+1
        ds_sub1 <- ds_sub[ds_sub[,hormone_var]==h,]
        if(loop==1){
          y_lab=conc_var
          if( two_axes ){ y_lab <- h}
           plot(ds_sub1[,conc_var] ~ ds_sub1[,time_var], type='l',xlim=c(xmin,xmax), ylim=c(ymin,ymax), 
            xlab=NA, ylab=y_lab, xaxt='n', col=colors[loop] )
        }
        if(loop>1 & two_axes==FALSE){lines(ds_sub1[,time_var],ds_sub1[,conc_var],col=colors[loop] )}
        if(loop>1 & two_axes){
          par( new=T )
          plot(ds_sub1[,conc_var] ~ ds_sub1[,time_var], type='l',xlim=c(xmin,xmax), ylim=c(ymin,ymax2), 
            xlab=NA, ylab=NA, xaxt='n', yaxt='n', col=colors[loop] )
          axis(4)
          mtext(h,4,line=2.5)
        }
        t_order <- c(ds_sub1[,time_var], rev(ds_sub1[,time_var]) )
        c_order <- c(rep(0,length(ds_sub1[,conc_var])),rev(ds_sub1[,conc_var])) 
        if( add_fill) polygon(t_order,c_order,col=adjustcolor(colors[loop], alpha=0.25),border=colors[loop]) 
      }
      legend('topleft',legend=sort(unique(ds_sub[,hormone_var])),fill=adjustcolor(colors, alpha=0.25),
              bty='n',cex=0.9,bg=NA, pt.cex=0.6)
      if(is.numeric(ds_sub[,time_var])){ axis(1)
      }else if( is.Date(ds_sub[,time_var]) ){
          ats <- seq( xmin, xmax, length.out = 5)
          axis.Date(1,at=ats, format=date_format)
      }else if( is.POSIXct(ds_sub[,time_var]) ){
          ats <- seq( xmin, xmax, length.out = 5)
          axis.POSIXct(1,at=ats,format=date_format)
      } 
    
      mtext(unique(ds_sub$plot_title),side=3,line=0.25)
       if( nrow(events)>0 ){
        for(l in 1:nrow(events)){
          arrows(x0=events[l,time_var],x1 =events[l,time_var],y0=ymax*0.95,y1=ymax*0.8, length = 0.1)
          text(x=events[l,time_var],y=ymax*0.99, events[l,x$event_var])
        }
      }
  }
  if( save_plot ){  dev.off()
      cat( paste0('\n *********\nNote: plots are saved at: \n', getwd(),'/hormPlotOverlap.pdf \n***** \n\n')  )
  }
}
