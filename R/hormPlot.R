#' Plot longitudinal hormone data with baseline information
#' 
#' @param x hormLong object (produced from hormBaseline) [required]
#' @param date_format the format of the date variable on x-axis. See help for examples of other formats [default = '%d-%b']
#' @param log_scale determines if y-axis is log10-scale or not. log-scale='y' makes log scale [default='n']  
#' @param plot_per_page the number of plot panels per page, by row. [default = 4]
#' @param save_plot indicates whether to save plot as a file [default = TRUE]
#' @param plot_height  the height of individual plot panels (in inches).  Pdf page height is determined by both plot_per_page and plot_height. [default = 2]
#' @param plot_width  the width of the pdf page. [default = 6]
#' @param yscale  determines if y-axis should be free ('free') to change for each panel or remain the same ('fixed') for all panels [default = 'free']
#' @param xscale  determines if x-axis should be free ('free') to change for each panel or remain the same ('fixed') for all panels  [default = 'free']
#' @param ...   generic plotting options [optional]  
#' @return nothing  Produces a pdf file saved at current working directory
#' @export
#' @examples
#' 
#' result <- hormBaseline(data=hormone, criteria=1, by_var='sp, sex, id', time_var='date', conc_var='conc' )
#' hormPlot(result, yscale='fixed',xscale='fixed' )
#' 

hormPlot <- function(x, date_format='%d-%b', log_scale='n', plot_per_page=4, save_plot=TRUE, 
                     plot_height=2, plot_width=6, yscale='free', xscale='free',...){
  
#--- main check ---#
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
  
  
#-- set-up ---#
  by_var_v <- cleanByvar(x$by_var) 
  time_var <- x$time_var
  conc_var <- x$conc_var
  data <- x$data
  data <- data[ do.call(order, data[c(by_var_v,time_var)]), ]
  
  data$plot_title <- getPlotTitle(data, by_var=by_var_v)

    
#--- create plots ---#
  if( save_plot ){
    pdf('hormPlot.pdf', height=plot_per_page * plot_height, width = plot_width )
  }

  par(mfrow=c(plot_per_page,1), mar=c(2,4,2,0.5),oma=c(2,2,2,1))
  for( i in unique(data$plot_title)){
    ds_sub <- data[data$plot_title==i, ]
    if(!is.null(x$criteria)){
       baseline <- getCutoff( ds_sub[ds_sub$conc_type=='base',conc_var], criteria=x$criteria )
      }else{ baseline <- 1 }
    if(!is.null(x$event_var)){
      events <- ds_sub[ !is.na(ds_sub[,x$event_var]) & ds_sub[,x$event_var]!='',c(x$event_var,time_var)]
      }else{events <- data.frame()}
    ds_sub <- ds_sub[!is.na(ds_sub[,conc_var]),]

    #--- set up scales
      if( yscale=='free'){
        ymin <- min(ds_sub[,conc_var])*0.95
        ymax <- max(baseline, max(ds_sub[,conc_var]) )*1.1
      }else{
        ymin <- min(data[,conc_var],na.rm=T)*0.95
        ymax <- max(data[,conc_var],na.rm=T)*1.1
      }
      if( xscale=='free'){
        xmin <- min( ds_sub[,time_var],na.rm=T )
        xmax <- max( ds_sub[,time_var],na.rm=T )
      }else{
        xmin <- min( data[,time_var],na.rm=T)
        xmax <- max( data[,time_var],na.rm=T)
      }    
    #--- main plot
    if( is.null(x$y_lab) ){x$y_lab <- conc_var}
    if( log_scale=='y' ){
        plot(ds_sub[,conc_var] ~ ds_sub[,time_var], type='l',xlim=c(xmin,xmax), ylim=c(ymin,ymax), 
            xlab=NA, ylab=x$y_lab, xaxt='n', log=log_scale)
    }else{
        plot(ds_sub[,conc_var] ~ ds_sub[,time_var], type='l',xlim=c(xmin,xmax), ylim=c(ymin,ymax), 
            xlab=NA, ylab=x$y_lab, xaxt='n')
    }
      if(is.numeric(ds_sub[,time_var])){ axis(1)
      }else if( is.Date(ds_sub[,time_var]) ){
          axis.Date(1,ds_sub[,time_var], format=date_format)
      }else if( is.POSIXct(ds_sub[,time_var]) ){
          axis.POSIXct(1,ds_sub[,time_var],format=date_format)
      } 
    
      points(ds_sub[,time_var], ds_sub[,conc_var],pch=19)
      mtext(unique(ds_sub$plot_title),side=3,line=0.25)
      abline(h = baseline, lty=2)
      if( nrow(events)>0 ){
        for(l in 1:nrow(events)){
          arrows(x0=events[l,time_var],x1 =events[l,time_var],y0=ymax*0.95,y1=ymax*0.8, length = 0.1)
          text(x=events[l,time_var],y=ymax*0.99, events[l,x$event_var])
        }
      }
  }
  if( save_plot ){  dev.off()
      cat( paste0('\n *********\nNote: plots are saved at: \n', getwd(),'/hormPlot.pdf \n***** \n\n')  )
  }
}
