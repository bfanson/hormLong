#' Plot longitudinal hormone data with baseline information
#' 
#' @param x hormLong object (produced from hormBaseline)
#' @param plot_per_page the number of plot panels per page, by row.
#' @param save_plot indicates whether 
#' @param plot_height  the height of individual plot panels (in inches).  Pdf page height is determined
#' by both plot_per_page and plot_height.
#' @param plot_width  the width of the pdf page.
#' @param yscale  determines if y-axis should be free ('free') to change for each panel or
#' remain the same ('fixed') for all panels 
#' @param xscale  determines if x-axis should be free ('free') to change for each panel or
#' remain the same ('fixed') for all panels 
#' @return NA
#' @export
#' @examples
#' 
#' result <- hormBaseline(data=hormone, criteria=1, by_var='sp, sex, id', time_var='date', conc_var='conc' )
#' hormPlot(result, yscale='fixed',xscale='fixed' )

hormPlot <- function(x, plot_per_page=4, save_plot=TRUE, plot_height=2, plot_width=6,
                     yscale='free', xscale='free',...){
  
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
  by_var_v <- by_var_v <- cleanByvar(x$by_var) 
  time_var <- x$time_var
  conc_var <- x$conc_var
  data <- x$data
  data <- data[ do.call(order, data[c(by_var_v,time_var)]), ]
  
  data$plot_title <- getPlotTitle(data, by_var=by_var_v)

    
#--- create plots ---#
  if( save_plot ){
    pdf('hormPlot.pdf', height=plot_per_page * plot_height, width = plot_width )
  }

  par(mfrow=c(plot_per_page,1), mar=c(2,4,2,0.5),oma=c(2,2,2,0))
  for( i in unique(data$plot_title)){
    ds_sub <- data[data$plot_title==i, ]
    baseline <- hormCutoff( ds_sub[ds_sub$conc_type=='base',conc_var], criteria=x$criteria )
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
    plot(ds_sub[,conc_var] ~ ds_sub[,time_var], type='l',xlim=c(xmin,xmax), ylim=c(ymin,ymax), 
          xlab=NA, ylab=conc_var)
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



#' Boxplot of individual concentrations 
#' 
#' @param data dataset 
#' @param id_var name of the individual variable
#' @return NA
#' @export
#' @examples
#' 
#' hormBoxplot(data=hormone, conc_var='conc',id_var='id')
#' 

hormBoxplot <- function(data, conc_var, id_var, plot_height=4, plot_width=6, save_plot=T){
  if( save_plot ){
    pdf('hormBoxplot.pdf', height=plot_height, width = plot_width )
    cat( paste0('\n *********\nNote: plots are saved at: \n', getwd(),'/hormBoxplot.pdf \n***** \n\n')  )
  }
  par(mar=c(3,4,1,1))
  boxplot( data[,conc_var] ~ data[,id_var],ylab=conc_var )  
    points(as.numeric(as.factor(data[,id_var])), data[,conc_var] )
  if( save_plot ){  dev.off() }
}


