#' Run iterative process to calculate baseline.
#' 
#' @param data dataset.
#' @param by_var grouping variable.
#' @param time_var grouping variable.
#' @param conc_var grouping variable.
#' @param criteria grouping variable.
#' @return hormLong object.
#' @export
#' @examples
#' 
#' result <- hormBaseline(data=hormone, by_var='sp, sex, id', time_var='date', conc_var='conc', event_var='Event' )
 

hormBaseline <- function(data, by_var,conc_var, time_var,criteria=2, event_var ) {
#--- initial checks of arguments ---#
  if(missing(data)){
      stop('data must be specified')
    }
  if(!is.data.frame(data)){
      stop('data must be a data.frame object')
    }
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
    stop('conc_var must be numeric')
    }

  if(missing(time_var)){
      stop('time_var must be specified')
    }
  if( !( class(data[,time_var]) %in% c('Date','numeric') ) ){
    time_class <- class(class(data[,time_var]))
    stop(paste0('time_var must be numeric or Date (your time_var is "',time_class,'" variable') )
    }
  
  if(missing(event_var)){
    event_var <- NULL    
  }

  
  by_var_v <- cleanByvar(by_var) # make by_var a vector
  if( sum(!( c(by_var_v,conc_var,time_var)  %in% names(data) ))>0 ){
      stop('not all variables are present in data.set.  check your column names')
    }


#--- checks of data.frame ---#
  data$row_id <- 1:nrow(data) 
  data1 <- data[,c('row_id',by_var_v,conc_var,time_var)] # keep only columns needed

  #--- remove na --#
    data1 <- na.omit(data1)

  #--- check for n=1 issues ---#
    hold <- aggregate(data1[,conc_var], by = data1[c(by_var_v)], FUN = function(y) length(y) )
    if( any( hold$x==1 ) ){
      row_id <- which( hold$x==1 )
      message('Error: only one sample per byvar: sd does not exist')
      print(data1[row_id,])
      stop('program stopped. Remove these groups')
    }
    rm(hold)

#--- iterative algorithm ---#
  hold <- data1
  hold$exclude <- 0
  total <- 99
  loop <- 0
  while(total>0){
    hold <- hold[hold$exclude==0,] 
    hold_add <- aggregate(data1[,conc_var], by = data1[c(by_var_v)], FUN = function(x)getCutoff(x,criteria=criteria) )
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
  return(data_class)
}

