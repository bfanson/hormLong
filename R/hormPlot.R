#' Plot longitudinal hormone data with baseline information
#' 
#' @param x hormLong object (produced from hormBaseline) [required]
#' @param date_format the format of the date variable on x-axis.  Default is 01-Jan format. See Appendix 1 in help manual 
#' for examples of other formats [default = '\%d-\%b']
#' @param xscale  determines if x-axis is free ('free') to change for each panel or remain the same ('fixed') for all panels  [default = 'free']
#' @param yscale  determines if y-axis is free ('free') to change for each panel or remain the same ('fixed') for all panels [default = 'free']
#' @param plot_per_page the number of plot panels per page, by row. [default = 4]
#' @param plot_height  the height of individual plot panels (in inches).  
#' Pdf page height is determined by both plot_per_page and plot_height. [default = 2]
#' @param plot_width  the width of the pdf page. [default = 6]
#' @param filename  the filename of the pdf file [default = 'hormPlot']
#' @param save_plot indicates whether to save plot as a file [default = TRUE]
#' @return nothing  Produces a pdf file saved at current working directory
#' @export
#' @examples
#' 
#' result <- hormBaseline(data=hormLynx, criteria=2, by_var='AnimalID, Hormone', time_var='datetime', 
#'                        conc_var='Conc', event_var='Events' )
#' hormPlot( result )
#'
#'# Sometimes you may want to fix y-axis to be same across all plots so heights are comparable
#' hormPlot( result, yscale='fixed' ) 
#' 

hormPlot <- function(x, date_format='%d-%b',
                     xscale='free', yscale='free',
                     plot_per_page=4, plot_height=2, plot_width=6, 
                     filename='hormPlot',  save_plot=TRUE){
  
#--- main check ---#
  graphics.off() # just to make sure no devices are open
  
  checkClass(x, 'hormLong')
  checkPlotOpts(plot_per_page, plot_width, plot_height, save_plot, xscale, yscale, date_format)

  
#-- set-up ---#
  by_var_v <- cleanByvar(x$by_var) 
  time_var <- x$time_var
  conc_var <- x$conc_var
  data <- x$data
  data <- data[ do.call(order, data[c(by_var_v,time_var)]), ]
  
  data$plot_title <- getPlotTitle(data, by_var=by_var_v)

#--- create plots ---#
  require(lubridate)
  if( save_plot ){
    pdf(paste0(filename,'.pdf'), height=plot_per_page * plot_height, width = plot_width )
  }

  par(mfrow=c(plot_per_page,1), mar=c(2,4,2,0.5),oma=c(2,2,2,1))
  for( i in unique(data$plot_title) ){
    ds_sub   <- data[data$plot_title==i, ]
    baseline <- getBaseline(ds_sub, x$criteria, conc_var)

    events <- getEventInfo(ds_sub, x$event_var, time_var)

    ds_sub <- ds_sub[!is.na(ds_sub[,conc_var]),]

    #--- set up scales ---#
      y_lim <- getPlotlim(d_s=ds_sub, d_f=data, var=conc_var, scale=yscale, base=baseline)
      x_lim <- getPlotlim(d_s=ds_sub, d_f=data, var=time_var, scale=xscale)
    
    #--- main plot
    if( is.null(x$y_lab) ){x$y_lab <- conc_var}
      plot(ds_sub[,conc_var] ~ ds_sub[,time_var], type='l',xlim=x_lim, ylim=y_lim, 
            xlab=NA, ylab=x$y_lab, xaxt='n')
  
      plotAxes(ds_sub, time_var, x_lim, date_format)

      points(ds_sub[,time_var], ds_sub[,conc_var],pch=19)
      mtext(unique(ds_sub$plot_title),side=3,line=0.25)
      abline(h = baseline, lty=2)

      plotEventInfo(events, t=time_var, e=x$event_var) 
  }
  if( save_plot ){  dev.off()
      cat( paste0('\n *********\nNote: plots are saved at: \n', getwd(),'/',filename,'.pdf \n***** \n\n')  )
  }
}
