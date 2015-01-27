#' Plot longitudinal hormone ratios  
#' 
#' @param x hormLong object (produced from hormBaseline) [required]
#' @param hormone_var name of the hormone variable.  This must have been included 
#' as part of by_var in hormBaseline  [required]
#' @param hormone_num name of the hormone that is the numerator in the ratio. [required]
#' @param hormone_denom name of the hormone that is the denominator in the ratio. [required]
#' @param ... arguments for hormPlot [optional]  
#' @return nothing  Produces a pdf file saved at current working directory
#' @export
#' @examples
#' 
#' result <- hormBaseline(data=hormElephant, criteria=2, by_var='Ele, Hormone', time_var='Date', 
#'                        conc_var='Conc_ng_ml', event_var='Event' )
#' hormPlotRatio( x=result, hormone_var='Hormone', hormone_num='Cortisol', hormone_denom='Progesterone' )
#' 

hormPlotRatio <- function(x, hormone_var, hormone_num, hormone_denom,...){
#                       date_format='%d-%b',xscale='free', yscale='free',
#                      plot_per_page=4, plot_height=2, plot_width=6, 
#                      filename='hormPlot',  save_plot=TRUE){
  
#--- main check ---#
  graphics.off() # just to make sure not devices are open
  
  checkClass(x, 'hormLong')

  checkPlotRatio(hormone_var, hormone_num, hormone_denom, x$by_var)
  
  
#-- set-up a new hormLong object ---#
  by_var_v <- cleanByvar(x$by_var) 
  new_by_var <- by_var_v[ by_var_v!=hormone_var ]

  ds_new <- x$data

  ds_new <- checkPlotMissing(ds_new, var_list=c(by_var_v,time_var) ) #get rid of missing data

  ds_new <- ds_new[ ds_new[,hormone_var] %in% c(hormone_num,hormone_denom), c(by_var_v,x$time_var,x$conc_var) ]
  ds_new_n <- ds_new[ ds_new[,hormone_var]==hormone_num,   ]
    ds_new_n[,hormone_var] <- NULL
  ds_new_d <- ds_new[ ds_new[,hormone_var]==hormone_denom, ] 
    ds_new_d[,hormone_var] <- NULL
  names(ds_new_n)[ which( names(ds_new_n)==x$conc_var) ] <- 'num_conc'
  names(ds_new_d)[ which( names(ds_new_d)==x$conc_var) ] <- 'denom_conc'
  ds_new <- merge( ds_new_n, ds_new_d) # inner merge
  ds_new$ratio <- ds_new$num_conc/ds_new$denom_conc
  ds_new <- ds_new[ !is.na(ds_new$ratio), ]
    ds_new$num_conc <- NULL 
    ds_new$denom_conc <- NULL 

  #--- remove any undefined ratios ---#
    if( max(ds_new$ratio)==Inf ){
      cat('WARNING: ratio contains Infinity which is due to denominator being zero.  All Infinity ratios
            will be removed. \n\n')
      ds_new <- subset( ds_new, ratio < Inf )
    }

#--- add in events ---#
  if(!is.null(x$event_var) ){
    ds_event <- x$data 
    ds_event <- ds_event[ !is.na(ds_event[,x$event_var]) , c(new_by_var,x$time_var,x$event_var)]
    ds_new[,x$event_var] <- NA
    ds_event[,'ratio']   <- NA
    ds_new <- rbind(ds_new, ds_event)        
  }

#-- set-up a new hormLong object ---#
  new_x <- x
  new_x$by_var   <-  paste0(new_by_var, collapse=',')
  new_x$criteria <- NULL 
  new_x$conc_var <- 'ratio'
  new_x$data <- ds_new
  new_x$y_lab <- paste0( hormone_num, '/', hormone_denom )


#--- call hormPlot with new hormLong object ---#
  hormPlot(x=new_x, filename='hormRatio',...)

}
