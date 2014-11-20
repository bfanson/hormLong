#' Write data from a hormLong object to a file
#' 
#' @param x hormLong object (produced from hormBaseline)
#' @param filename name of file to be saved. Unless path is specified, file will be
#' saved in your current working directory 
#' @param file_type  determines type of file (e.g. csv, tab-delimited) 
#' @param ...  options for write.csv 
#' @return nothing
#' @export
#' @examples
#' 
#' result <- hormBaseline(data=hormone, by_var='sp, sex, id', time_var='date', conc_var='conc' )
#' hormWrite(result, filename='data_out.csv', file_type='csv' )

hormWrite <- function(x, filename, file_type='csv',... ){
#--- main check ---#
  if( missing(filename) ){
      stop('please specify a filename')
  }
  if( file_type!='csv' ){
      stop('only csv file type has currently been implemented')
  }
  
  write.csv( x$data, file=filename, quote=FALSE, row.names=FALSE,...) 
  cat(paste0('*****\nNote: If no file location was specified, then file is at:\n',getwd(),'\n***') )
}


#' Helper function to calculate baseline cutoff
#' 
#' @param x hormone concentration
#' @param criteria the number of standard deviations above baseline 
#' @return nothing
#' @export
#' @examples
#' 
#' conc <- rnorm(100)
#' hormCutoff(conc, 2)

hormCutoff <- function(x, criteria){ 
  mean(x, na.rm=T)+sd(x, na.rm=T)*criteria 
}


#' Helper function to work with dates
#' 
#' @param data data set with the variables to convert
#' @param date_var variable name for date column in format '2014-jan-01'
#' @param time_var variable name for the time column in format '10:10 PM'
#' @return new datetime variable
#' @export
#' @examples
#' 
#' date <- c('2014-jan-01','2014-jan-02')
#' time <- c('10:10 Am', '10:30 pm')
#' hormDate(date,time)

hormDate <- function(data,date_var,time_var){
  if( missing(date_var) ){stop('please provide date_var. Time_var is optional ')}
  if( !missing(time_var) & any( grepl('-',data[,date_var])==F ) ){
    stop('check format of date_var.  It should be in from: 01-jan-2014')
  }
  if( !missing(time_var)){
    if(any( (grepl('pm',tolower(data[,time_var]))==F |  
                              grepl('am',tolower(time_var))==F)==F ) ){
    stop('check format of date_var.  It should be in from: 10:00 AM')
  }}
  if(!missing(date_var) & !missing(time_var)){
    datetime <- as.POSIXct(paste(data[,date_var],data[,time_var]), format="%Y-%m-%d %I:%M %p", tz='GMT' )
  }
  if(!missing(date_var) & missing(time_var)){
    datetime <- as.Date(data[,date_var])
  }
  data$datetime<-datetime  
  return(data)
}
  

#' Read hormone data from a csv file
#' 
#' @param none  no arguments
#' @return imported dataset
#' @export
#' @examples
#' 
#' ds <- hormRead()

hormRead <- function( ){
  ds<- read.csv(file.choose() )
  cat('Below shows the first 6 lines of the imported dataset.  check that it imported correctly.
      If not, check help file.\n\n')
  print( head(ds) )
  return(ds)
}


#' Automated function that runs baseline calculation, make plots, outputs dataset as csv
#' 
#' @param file_method  if 'choose' then user is prompted. if 'write' then user writes the name 
#' @return nothing in R.  Output hormPlot.pdf and horm_data_out.csv
#' @export
#' @examples
#' 
#' ds <- hormAuto(by_var='sp, sex, id', time_var='date', conc_var='conc')

hormAuto<- function(file_method='write',by_var, time_var, conc_var, filename='test' ){
  ds  <- hormRead()
  ds <- hormone
  res <- hormBaseline(data=ds, by_var=by_var,time_var=time_var,conc_var=conc_var)
  hormPlot(res)
  hormWrite(res,'horm_data_out.csv')
#   if( file_method=='write'){ 
#     if(missing(filename){ }
}


#' Run iterative process to calculate baseline.
#' 
#' @param data dataset.
#' @param by_var grouping variable.
#' @param time_var grouping variable.
#' @param conc_var grouping variable.
#' 
#' @return hormLong object.
#' @export
#' @examples
#' 
#' result <- hormACTH(x=result)
 
hormACTH <- function(x, start=0, end=10) {
  if( class(x)!='hormLong'){
      stop('Object needs to be hormLong.  Run hormBaseline() first')
  }
  ds <- x$data
  by_var_v <- cleanByvar(x$by_var)
  time_var <- x$time_var
  conc_var <- x$conc_var
  ds <- ds[ start <= ds[,time_var] & ds[,time_var] <= start, ]
  ds_base <- getSumStat( ds[ ds[,'conc_type']=='base',], name='baseline', func=function(x) mean(x,na.rm=T) )
  ds <- merge(ds, ds_base, all.x=T)
  ds$conc_baseline <- ds[,conc_var]/ x$baseline
  
}



#' Splits and trims by_var

cleanByvar <- function(x){
  by_var_v <- gsub("^\\s+|\\s+$", "", unlist(strsplit(x,',')) )
  return(by_var_v)
}

#' Write results to rtf file
#' 

write.rtf <- function(x,file){
  require(rtf)
  rtf<-RTF(file, width=7, height=11, font.size=9, omi=c(3,1,1,1))
      addTable.RTF(rtf, x)
  done(rtf)
}

#' Helper to get aggregate statistics
#' 

  getSumStat <- function(data=ds,name='mean', func=function(x)mean(x,na.rm=T), add_ds ){
    ds1 <- aggregate(data[,conc_var], by = data[c(by_var_v)], FUN = func )
      names(ds1)[ncol(ds1)] <- name
      if(!missing(add_ds)){ ds1 <- merge(add_ds,ds1,all=T)}
      return(ds1)
  }    
