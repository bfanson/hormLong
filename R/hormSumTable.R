#' Summarise results from baseline calculation.
#' 
#' @param x hormLong object (produced from hormBaseline) [required]
#' @param num_decimals number of decimal places [default = 2]
#' @return nothing  Produces a file saved to current working directory
#' @export
#' @examples
#' 
#' result <- hormBaseline(data=hormLynx, criteria=2, by_var='AnimalID, Hormone', time_var='datetime', conc_var='Conc' )
#' hormSumTable(result) 

hormSumTable <- function(x, num_decimals=2){
#-- checks --#  
  if( class(x)!='hormLong'){
      stop('Object needs to be hormLong.  Run hormBaseline() first')
    }

#--- set-up info ---#
  ds <- ridFactor(x$data)
  conc_var <- x$conc_var
  by_var_v <- cleanByvar(x$by_var) 

#--- make table ---#
  ds_b <- aggregate(ds[,conc_var], by = ds[c(by_var_v,'conc_type')], FUN = mean, na.rm=T )
    ds_b <- merge( subset(ds_b,conc_type=='base',-conc_type), subset(ds_b,conc_type=='peak',-conc_type),
                   by=by_var_v,all=T )
    names(ds_b)[(ncol(ds_b)-1):ncol(ds_b)] <- c('base_mean','peak_mean')
    ds_b$peak_base <- ds_b$peak_mean/ds_b$base_mean

  ds2 <- getSumStat(data=ds,name='n', func= function(x) length( x[!is.na(x)] ), add_ds=ds_b,c_var=conc_var,by_var = by_var_v )
  ds2 <- getSumStat(data=ds,name='mean', func=function(x)mean(x,na.rm=T), add_ds=ds2,c_var=conc_var,by_var = by_var_v  )
  ds2 <- getSumStat(data=ds,name='sd', func=function(x)sd(x,na.rm=T), add_ds=ds2,c_var=conc_var,by_var = by_var_v )
  ds2$sem <- ds2$sd/sqrt(ds2$n)  
  ds2 <- getSumStat(data=ds,name='percent_cv', func= function(x) sd(x,na.rm=T)/mean(x,na.rm=T)*100, add_ds=ds2,c_var=conc_var,by_var = by_var_v  )
  ds2 <- getSumStat(data=ds[ds$conc_type=='base',],name='cutoff', func= function(y) getCutoff(y, criteria=x$criteria), add_ds=ds2,c_var=conc_var,by_var = by_var_v )
  ds2 <- getSumStat(data=ds,name='min', func= function(x) min(x,na.rm=T), add_ds=ds2,c_var=conc_var,by_var = by_var_v )
  ds2 <- getSumStat(data=ds,name='max', func= function(x) max(x,na.rm=T), add_ds=ds2,c_var=conc_var,by_var = by_var_v )
  ds2 <- getSumStat(data=ds,name='median', func= function(x) median(x,na.rm=T), add_ds=ds2,c_var=conc_var,by_var = by_var_v )
    
  ds2 <- merge(ds2,ds_b[,c(by_var_v, 'peak_base')],all=T)

  ds2[,(length(by_var_v)+1):ncol(ds2)] <- sapply(ds2[,(length(by_var_v)+1):ncol(ds2)],round,num_decimals)  
  ds_out <- ds2[, c(by_var_v, 'n', 'mean','median','sd','sem','percent_cv','min','max','cutoff',
                              'base_mean','peak_mean','peak_base')]

#--- output table ---#
    write.csv(ds_out, file='hormSumTable.csv', quote=F, row.names=F,na='')
    cat( paste0('\n *********\nNote: table saved at: \n', getwd(),'/hormSumTable.csv \n***** \n\n')  )
}