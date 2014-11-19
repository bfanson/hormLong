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
  ds_m <- aggregate(ds[,conc_var], by = ds[c(by_var_v)], FUN = mean, na.rm=T )
    names(ds_m)[ncol(ds_m)] <- 'mean'
  ds_s <- aggregate(ds[,conc_var], by = ds[c(by_var_v)], FUN = sd, na.rm=T )
    names(ds_s)[ncol(ds_s)] <- 'sd'
    ds_out <- merge(ds_m, ds_s,all=T)
  ds_p <- aggregate(ds[,conc_var], by = ds[c(by_var_v)], FUN = function(x) sd(x,na.rm=T)/mean(x,na.rm=T)*100 )
    names(ds_p)[ncol(ds_p)] <- 'percent_cv'
    ds_out <- merge(ds_out, ds_p,all=T)
  ds_b <- aggregate(ds[,conc_var], by = ds[c(by_var_v,'conc_type')], FUN = mean, na.rm=T )
    ds_b <- merge( subset(ds_b,conc_type=='base',-conc_type), subset(ds_b,conc_type=='peak',-conc_type),
                   by=by_var_v,all=T )
    names(ds_b)[(ncol(ds_b)-1):ncol(ds_b)] <- c('base_mean','peak_mean')
    ds_out <- merge(ds_out,ds_b, all=T) 
  ds_c <- aggregate(ds[,conc_var], by = ds[c(by_var_v)], FUN = hormCutoff, criteria=x$criteria )
    names(ds_c)[ncol(ds_c)] <- 'cuttoff'
    ds_out <- merge(ds_out, ds_c, all=T)

  ds_out <- cbind( ds_out[,by_var_v], sapply(ds_out[,(length(by_var_v)+1):ncol(ds_out)],round,num_decimals) )
#--- output table ---#
  if( file_type=='csv'){
    write.csv(ds_out,file='hormSumTable.csv',quote=F, row.names=F)
    cat( paste0('\n *********\nNote: table saved at: \n', getwd(),'/hormSumTable.csv \n***** \n\n')  )
  }
  if( file_type=='csv'){
    write.rtf(ds_out,file='hormSumTable.rtf')
    cat( paste0('\n *********\nNote: table saved at: \n', getwd(),'/hormSumTable.rtf \n***** \n\n')  )
  }
}