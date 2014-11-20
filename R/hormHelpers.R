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
#' @param date date column in format '2014-jan-01'
#' @param time time column in format '10:10 PM'
#' @return new datetime variable
#' @export
#' @examples
#' 
#' date <- c('2014-jan-01','2014-jan-02')
#' time <- c('10:10 Am', '10:30 pm')
#' hormDate(date,time)

hormDate <- function(date,time){
  if( any( grepl('-',date)==F ) ){
    stop('check format of date.  It should be in from: 2014-jan-01')
  }
  if( any( (grepl('pm',tolower(time))==F |  grepl('am',tolower(time))==F)==F ) ){
    stop('check format of date.  It should be in from: 10:00 AM')
  }
  datetime <- as.POSIXct(paste(date,time), format="%Y-%b-%d %I:%M %p", tz='GMT' )
  return(datetime)
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