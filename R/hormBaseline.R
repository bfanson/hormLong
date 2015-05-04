#' Run iterative process to calculate baseline.
#' 
#' @param data dataset containing the hormone data. [required]
#' @param by_var names of the grouping variables.  Column names should be separated by a comma:
#'  e.g. 'name1' for one column, 'name1, name2' for two columns, 'name1, name2, name3' for three columns,
#'  etc.  Remember, capitalization matters. [optional]
#' @param time_var name of the time variable (e.g. 'date', 'datetime', 'days') [required].
#' @param conc_var name of the concentration variable (response variable) [required].
#' @param event_var name of event variable  [optional]
#' @param criteria baseline criteria (mean + criteria * SD) [default=2]
#' @param save_data determines if the output dataset should be saved to a csv file [default = T ]
#' @return hormLong object.
#' @export
#' @examples
#' 
#' result <- hormBaseline(data=hormLynx, criteria=2, by_var='AnimalID, Hormone', time_var='datetime', conc_var='Conc' )
 

hormBaseline <- function(data, by_var, conc_var, time_var, criteria=2, event_var=NULL, save_data=T ) {
#--- initial checks of arguments ---#
  checkCapitalization(data, var_list = c(by_var,conc_var,time_var,event_var) )

  if(missing(data)){
      stop('data must be specified')
    }
  if(!is.data.frame(data)){
      stop('data must be a data.frame object')
    }
  data <- ridFactor(data) # get rid of all factors
  
  if(missing(by_var)){
    message('Warning: No by_var included ... baseline value is for whole dataset')
    data$hold_id <- 1    # set-up a new by_var
    by_var <- 'hold_id'  
    }
  if(!is.character(by_var)){
    stop('by_var must be a character string')
    }
  if( is.null(criteria) | is.na(criteria) ){
    stop('criteria must be specified')
    }
  if(!is.numeric(criteria)){
    stop('criteria must be numeric')
    }

  if(missing(conc_var)){
      stop('conc_var (e.g. the response variable containing the concentration) must be specified')
    }
  if(!is.numeric(data[,conc_var])){
    cat('WARNING: conc_var should be numeric. Program will convert automatically.  Non-numeric
         values will be become missing.  Any non_numeric values are listed below:\n\n')
    print('------')
    suppressWarnings(print( data[ is.na(as.numeric( data[,conc_var])),] ))
    cat('\n')
    suppressWarnings(data[,conc_var] <- as.numeric(data[,conc_var])  )
    }

  if(missing(time_var)){
      stop('time_var must be specified')
    }
  if( !( class(data[,time_var])[1] %in% c('Date','numeric','POSIXct') ) ){
    time_class <- class(class(data[,time_var]))
    stop(paste0('time_var must be numeric or Date (your time_var is "',time_class,'" variable') )
    }

  by_var_v <- cleanByvar(by_var)
  
#--- checks of data.frame ---#
  data <- data[ do.call('order', data[c(by_var_v,time_var)]), ] 
  data$row_id <- 1:nrow(data) 
  data1 <- data[,c('row_id',by_var_v,conc_var)] # keep only columns needed

  if( any( is.na(data1[,by_var_v]) ) ){
    cat("ERROR: The following rows have missing by_var or hormone_var data and will be removed from plotting function:\n")
    print( data1[ apply(data1[,by_var_v], 1, function(x){ any( is.na(x) ) }  ), ])
    stop("Remove/fix the above rows that have missing by_var information")
  }


  #--- remove na --#
    data1 <- na.omit(data1)

  #--- check for n=1 issues ---#
    hold <- aggregate(data1[,conc_var], by = data1[c(by_var_v)], FUN = function(y) length(y) )
    if( any( hold$x==1 ) ){
      row_id <- which( hold$x==1 )
      message('Warning: only one sample per byvar: sd does not exist. The following group \n have been skipped:')
      out <- merge(data1, hold[row_id, by_var_v])
      print( out )
      data1 <- data1[!(data1$row_id %in% out$row_id),] 
    }
    rm(hold)

#--- iterative algorithm ---#
  hold <- data1
  hold$exclude <- 0
  hold$x <- 0
  total <- 99
  loop <- 0
  cat('\n')
  print('*---  Iteration History   ----*')
  while(total>0){
    hold <- subset( hold, exclude==0, -x ) 
    hold_add <- aggregate(hold[,conc_var], by = hold[c(by_var_v)], FUN = function(x) getCutoff(x,criteria=criteria) )
    hold <- merge(hold, hold_add )
    hold$exclude <- ifelse( hold[,conc_var] > hold$x, 1,0 ) 
    total <- sum(hold$exclude)
    loop <- loop + 1
    print(paste('Iteration = ', loop,':  total removed is ',total) )
  }

#--- add in base/peak information and clean up dataset ---#
  data1 <- merge(data1, hold[,c('row_id','exclude')], by='row_id', all.x=T) 
  data1$exclude[is.na(data1$exclude)] <- 1
  data1$conc_type <- ifelse(data1$exclude==1,'peak','base') 
 
  data2 <- merge(data, data1[,c('row_id','conc_type')],by='row_id',all.x=T)
  data_class <- list( data=data2, by_var=by_var, time_var=time_var,
                      conc_var=conc_var, event_var=event_var, criteria=criteria) 
  class(data_class) <- 'hormLong'

    #--- write dataset ---#
  if( save_data ){
      write.csv(data2,file='hormBaseData.csv',quote=F, row.names=F)
      cat( paste0('\n *********\nNote: table saved at: \n', getwd(),'/hormBaseData.csv \n***** \n\n')  )
  }

  return(data_class)

}

