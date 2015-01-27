#' Plot longitudinal hormone data with baseline information
#' 
#' @param x hormLong object (produced from hormBaseline) [required]
#' @param date_format the format of the date variable on x-axis.  Default is 01-Jan format. See Appendix 1 in help manual 
#' for examples of other formats [default = '\%d-\%b']
#' @param break_cutoff the maximum number of days between consecutive points. 
#' Above this cutoff value, a break is created.  Default is  [default = Inf]
#' @param color colour of the line and points [default='black']
#' @param symbol number to indicate point symbol. e.g. 1=open circle, 2=open triangle, 15=closed square, 19=closed circle  [default=19]
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
#' hormPlot( result, yscale='fixed', line_color='red', symbol=1 ) 
#'
#'# You may want to just include both date and time on the x-axis instead of date
#' hormPlot( result, date_format='%d/%m/%y %H:%m' )   # does not makes sense here but the code shows how 
#'
#' 
#'# You may want to not join points that are temporally separated by large gaps
#'result <- hormBaseline(data=hormElephant, criteria=2, by_var='Ele, Hormone', time_var='Date', 
#'              conc_var='Conc_ng_ml' , event_var='Event')
#' hormPlot( result )  # everything connected
#' hormPlot( result, yscale='fixed', break_cutoff=30 ) # do not join any points more than 20 days
#' 

hormPlot <- function(x, date_format='%d-%b', break_cutoff=Inf,
                     xscale='free', yscale='free',
                     color='black', symbol=19,
                     plot_per_page=4, plot_height=2, plot_width=6, 
                     filename='hormPlot', save_plot=TRUE){
  
#--- main check ---#
  graphics.off() # just to make sure no devices are open
  
  checkClass(x, 'hormLong')
  checkPlotOpts(plot_per_page, plot_width, plot_height, save_plot, xscale, yscale, date_format)

  checkPlotBreaks(break_cutoff, 10) 

#-- set-up ---#
  by_var_v <- cleanByvar(x$by_var) 
  time_var <- x$time_var
  conc_var <- x$conc_var
  data <- x$data
  data <- data[ do.call(order, data[c(by_var_v,time_var)]), ]
  
  data$plot_title <- getPlotTitle(data, by_var=by_var_v)

#-- check for missing data in by_var_v, time_var, hormone_var --#
   data <- checkPlotMissing(data, var_list=c(by_var_v,time_var) )

#--- get break information ---#
  require(lubridate)
  if( is.Date(data[,time_var]) ){ data[,time_var] <- as.numeric( data[,time_var] ) * 3600*24 
    }else{  data[,time_var] <- as.numeric( data[,time_var] ) } 
  data$row_id <- 1:nrow(data)
    tmp <- do.call(rbind,by(data,data$plot_title,FUN=function(x, time=time_var) cbind(row_id=x$row_id,diff=c(-99, diff(x[,time]))  ) ))
    data <- merge(data,tmp,all.x=T)
    data$brk <- 0
    data <- data[ do.call('order', data[c('plot_title',time_var)]), ] #need to resort the data
    for(i in 1:nrow(data)){
      if(data$diff[i] == -99){brk_num=1}
      if(data$diff[i]/(3600*24) > break_cutoff){brk_num=brk_num+1}
      data$brk[i] <- brk_num
    }
   data[,time_var] <- as.POSIXct( data[,time_var], origin='1970-01-01 00:00:00', tz='UTC' )

#--- create plots ---#
  require(lubridate)
  if( save_plot ){
    pdf(paste0(filename,'.pdf'), height=plot_per_page * plot_height, width = plot_width )
  }

  par(mfrow=c(plot_per_page,1), mar=c(2,4,2,0.5),oma=c(2,2,2,1))
  for( i in unique(data$plot_title) ){
    ds_sub   <- data[data$plot_title==i, ]
    baseline <- getBaseline(ds_sub, x$criteria, conc_var)

    events <- getEventInfo(ds_sub , x$event_var, time_var)

    ds_sub <- ds_sub[!is.na(ds_sub[,conc_var]),]

    #--- set up scales ---#
      y_lim <- getPlotlim(d_s=ds_sub, d_f=data, var=conc_var, scale=yscale, base=baseline)
      x_lim <- getPlotlim(d_s=ds_sub, d_f=data, var=time_var, scale=xscale)
    
    #--- main plot
    if( is.null(x$y_lab) ){x$y_lab <- conc_var}
      plot(ds_sub[,conc_var] ~ ds_sub[,time_var], type='n',xlim=x_lim, ylim=y_lim, 
            xlab=NA, ylab=x$y_lab, xaxt='n')
  
      plotLines(ds_sub, conc_var, time_var, color)

      plotAxes(ds_sub, time_var, x_lim, date_format)

      points(ds_sub[,time_var], ds_sub[,conc_var], pch=symbol, col=color)
      mtext(unique(ds_sub$plot_title),side=3,line=0.25)
      abline(h = baseline, lty=2)

      plotEventInfo(events, t=time_var, e=x$event_var) 
  }
  if( save_plot ){  dev.off()
      cat( paste0('\n *********\nNote: plots are saved at: \n', getwd(),'/',filename,'.pdf \n***** \n\n')  )
  }
}
