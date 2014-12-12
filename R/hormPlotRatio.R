#' Plot longitudinal hormone ratios  
#' 
#' @param x hormLong object (produced from hormBaseline) [required]
#' @param horm_var name of the hormone variable.  This must have been included 
#' as part of by_var in hormBaseline  [required]
#' @param horm_num name of the hormone that is the numerator in the ratio. [required]
#' @param horm_denom name of the hormone that is the denominator in the ratio. [required]
#' @param ... arguments for hormPlot [optional]  
#' @return nothing  Produces a pdf file saved at current working directory
#' @export
#' @examples
#' 
#' ds <- hormDate(hormone2,date_var = 'date', date_order = 'ymd')
#' res <- hormBaseline(data=ds,criteria=1,by_var='sp,horm_type,id',conc_var = 'conc',time_var='date',
#'                      event_var='event')
#' hormPlotRatio( x=res, horm_var='horm_type', horm_num='E', horm_denom='prog' )
#' 

hormPlotRatio <- function(x, horm_var, horm_num, horm_denom, ...){
#--- main check ---#
  if( class(x)!='hormLong'){
      stop('Object needs to be hormLong.  Run hormBaseline() first')
  }
  if( missing(horm_var) | missing(horm_num) | missing(horm_denom)){
      stop('horm_var and horm_num and horm_denom must be specified')
  }
  
  by_var_v <- cleanByvar(x$by_var) 
  if( !(horm_var %in% by_var_v) ){
    stop('horm_var must have been specified in by_var when running hormBaseline')
  }
  if( !horm_num %in% unique(x$data[,horm_var])){
    stop('horm_num must be a level in horm_var.  check spelling and capitalization')
  }  
  if( !horm_denom %in% unique(x$data[,horm_var]) ){
    stop('horm_denom must be a level in horm_var.  check spelling and capitalization')
  }  

#-- set-up a new hormLong object ---#
  new_by_var <- by_var_v[ by_var_v!=horm_var ]
  ds_new <- x$data
  ds_new <- ds_new[ ds_new[,horm_var] %in% c(horm_num,horm_denom), c(by_var_v,x$time_var,x$conc_var) ]
  ds_new_n <- ds_new[ ds_new[,horm_var]==horm_num,   ]
    ds_new_n[,horm_var] <- NULL
  ds_new_d <- ds_new[ ds_new[,horm_var]==horm_denom, ] 
    ds_new_d[,horm_var] <- NULL
  names(ds_new_n)[ which( names(ds_new_n)==x$conc_var) ] <- 'num_conc'
  names(ds_new_d)[ which( names(ds_new_d)==x$conc_var) ] <- 'denom_conc'
  ds_new <- merge( ds_new_n, ds_new_d) # inner merge
  ds_new$ratio <- ds_new$num_conc/ds_new$denom_conc
  ds_new <- ds_new[ !is.na(ds_new$ratio), ]
    ds_new$num_conc <- NULL 
    ds_new$denom_conc <- NULL 

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
  new_x$y_lab <- paste0( horm_num, '/', horm_denom )


#--- call hormPlot with new hormLong object ---#
  hormPlot(x=new_x, ...)

}
