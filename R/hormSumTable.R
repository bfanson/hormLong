#' Summarise results from baseline calculation.
#' 
#' @param x hormLong object (produced from hormBaseline)
#' @param file_type choose either rtf or csv
#' @param num_decimals number of decimal places
#' @return hormLong object.
#' @export
#' @examples
#' 
#' result <- hormBaseline(data=hormone, by_var='sp, sex, id', time_var='date', conc_var='conc' )
#' hormSumTable(result) 

hormSumTable <- function(x,file_type='rtf', num_decimals=2){
#-- checks --#  
  if( class(x)!='hormLong'){
      stop('Object needs to be hormLong.  Run hormBaseline() first')
    }
  if( !(file_type %in% c('rtf','csv') ) ){
      stop( paste('xscale must be either "rtf" or "csv", not:',file_type) )
    }
#--- set-up info ---#
  ds <- x$data
  conc_var <- x$conc_var
  by_var_v <- cleanByvar(x$by_var) 

#--- make table ---#
  ds_b <- aggregate(ds[,conc_var], by = ds[c(by_var_v,'conc_type')], FUN = mean, na.rm=T )
    ds_b <- merge( subset(ds_b,conc_type=='base',-conc_type), subset(ds_b,conc_type=='peak',-conc_type),
                   by=by_var_v,all=T )
    names(ds_b)[(ncol(ds_b)-1):ncol(ds_b)] <- c('base_mean','peak_mean')

  ds2 <- getSumStat(data=ds,name='mean', func=function(x)mean(x,na.rm=T), add_ds=ds_b  )
  ds2 <- getSumStat(data=ds,name='sd', func=function(x)sd(x,na.rm=T), add_ds=ds2 )
  ds2 <- getSumStat(data=ds,name='percent_cv', func= function(x) sd(x,na.rm=T)/mean(x,na.rm=T)*100, add_ds=ds2  )
  ds2 <- getSumStat(data=ds,name='cutoff', func= function(y) hormCutoff(y, criteria=x$criteria ), add_ds=ds2 )
  ds2 <- getSumStat(data=ds,name='min', func= function(x) min(x,na.rm=T), add_ds=ds2 )
  ds2 <- getSumStat(data=ds,name='max', func= function(x) max(x,na.rm=T), add_ds=ds2 )
  ds2 <- getSumStat(data=ds,name='median', func= function(x) median(x,na.rm=T), add_ds=ds2 )

  ds_out <- cbind( ds2[,by_var_v], sapply(ds2[,(length(by_var_v)+1):ncol(ds2)],round,num_decimals) )

#--- output table ---#
  if( file_type=='csv'){
    write.csv(ds_out,file='hormSumTable.csv',quote=F, row.names=F)
    cat( paste0('\n *********\nNote: table saved at: \n', getwd(),'/hormSumTable.csv \n***** \n\n')  )
  }
  if( file_type=='rtf'){
    write.rtf(ds_out,file='hormSumTable.rtf')
    cat( paste0('\n *********\nNote: table saved at: \n', getwd(),'/hormSumTable.rtf \n***** \n\n')  )
  }
}