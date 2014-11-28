#' Write data from a hormLong object to a file
#' 
#' @param x hormLong object (produced from hormBaseline)
#' @param filename name of file to be saved. Unless path is specified, file will be
#' saved in your current working directory 
#' @param file_type  determines type of file (only csv currently)  
#' @return nothing
#' @export
#' @examples
#' 
#' result <- hormBaseline(data=hormone, by_var='sp, sex, id', time_var='date', conc_var='conc' )
#' head(result)

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
#' @return cutoff value
#' @export
#' @examples
#' 
#' conc <- rnorm(100)
#' getCutoff(conc, 2)

getCutoff <- function(x, criteria){
  mean(x, na.rm=T)+sd(x, na.rm=T)*criteria 
}


#' Helper function to work with dates
#' 
#' @param data data set with the variables to convert
#' @param date_var variable name for date column in format '2014-01-01'
#' @param time_var variable name for the time column in format '10:10 PM'
#' @param name name of the new date variable created by the function
#' @return new datetime variable
#' @export
#' @examples
#' 
#' date <- c('2014-01-01','2014-01-01')
#' time <- c('10:10 AM', '10:30 PM')
#' ds <- data.frame(date=date,time=time)
#' hormDate(ds,date_var='date',time_var='time',name='datetime' )

hormDate <- function(data,date_var,time_var,name='datetime'){
  if( missing(date_var) ){stop('please provide date_var. Time_var is optional ')}
  if( !missing(time_var) & any( grepl('-',data[,date_var])==F ) ){
    stop('check format of date_var.  It should be in from: 2014-01-01')
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
  data$rename_this<-datetime  
  names(data)[ncol(data)] <- name
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

#' Splits and trims by_var
#' 
#' @param x by_var to be cleans (e.g. 'species, hormone' )
#' @return x split and trimmed of whitespace
#' @export
#' @examples
#' 
#' ds <- hormRead()

cleanByvar <- function(x){
  by_var_v <- gsub("^\\s+|\\s+$", "", unlist(strsplit(x,',')) )
  return(by_var_v)
}

#' Write results to rtf file
#'
#' @param x dataframe to be writen to a file
#' @param file filename 
#' @return x split and trimmed of whitespace
#' @export
#' @examples
#' 
#' write.rtf(cars,'cars.rtf')

write.rtf <- function(x,filename){
  require(rtf)
  rtf<-RTF(filename, width=7, height=11, font.size=9, omi=c(3,1,1,1))
      addTable.RTF(rtf, x)
  done(rtf)
}


#' Combined the columns listed in the by_var_v
#'
#' @param data dataframe containing the variables in by_var_v
#' @param by_var the by_var_v (vector of names) 
#' @return the combined columns as a single vector, separated by a semicolon
#' @export
#' @examples
#' 
#' (plot_title <- getPlotTitle(hormone,by_var=c('sp','sex')))


getPlotTitle <- function(data, by_var=by_var_v){
  #-- create titles--#
    if( length(by_var)>1){
        plot_title <- do.call(paste1, data[,c(by_var)] )
     }else{plot_title <- data[,c(by_var)] }
}
paste1 <- function(...) paste(...,sep='; ')


#' Get the area under a curve using trapezoid method
#'
#' @param t time vector
#' @param c response vector (e.g. hormone concentration) 
#' @return calculates the area under the curve for the curve
#' @export
#' @examples
#' 
#' time <- 1:10
#' conc <- c( rep(0,4),2,2,rep(0,4))
#' plot(conc~time,type='l')
#' getAUC(time,conc)

getAUC <- function(t,c){
    require(zoo)
    (AUC <- sum(diff(t)*rollmean(c,2)))
  }  

#' Assigns peak number and get area under curve for each peak
#'
#' @param k dataframe (broken first by by_var)
#' @param date time variable
#' @param conc response variable (e.g. concentration)
#' @return return dataset with peaks identified and AUC calculated
#' @export
#' @examples
#'   
#' head(hormone) 


getPeakInfo <- function(k, date=time_var,conc=conc_var){
  if( !any(k$conc_type=='peak')){
    return(data.frame(peak_num=0,t=0, c=0, cutoff=k$cutoff[1], AUC=0))
  }else{
    cutoff <- k$cutoff[1] 
    t<-as.numeric(k[,date])
    c<-k[,conc] - cutoff
    c_t <- k$conc_type

  #-- assign peak number ---#
    num <- 0
    peak_num <- rep(0,length(t) )
    if(c_t[1]=='peak'){ num <- num+1; peak_num[1] <- num} 
    for(i in 2:length(t) ){
        if( c_t[i]=='peak' & c_t[i-1]=='base' ){ num <- num+1}
        if( c_t[i]=='peak'){ peak_num[i] <- num} 
      }
  
  #-- extract peak AUC info ---#
   if(max(peak_num)>0){
      peaks <- list( )
      for( p in 1:max(peak_num) ){
        st <-  max(1, min( which(peak_num==p) )-1 )
        end <- min(length(t), max(which(peak_num==p))+1 )
        st_ind <- ifelse( min( which(peak_num==p))==1,1,0)
        
      #-- calculate point that slope crosses cutoff --#
        m <- diff(c[st:end])/diff(t[st:end])
        b <- c[st:(end-1)] - m*t[st:(end-1)]
        t0 <- -b[unique(c(1,length(b)))]/m[unique(c(1,length(m)))]
        if( length(m)==1 & st==1){ t0 <-c(t[st],t0)}   
        if( length(m)==1 & st>1){ t0  <-c(t0,t[end])}   
        if( st_ind==1){ t0[1] <-c(t[st])  }
      
        t_new <- t[c(st,st:end,end)]
        t_new[c(1:2,(length(t_new)-1):length(t_new))] <- rep(t0,each=2) 
          
        c_new <- c[c(st,st:end,end)]
        c_new[c(1,length(c_new))] <- 0
        c_new[c(2,(length(c_new)-1))] <- pmax(c_new[c(2,(length(c_new)-1))],0)
        
        AUC <- getAUC(t_new,c_new)
        peaks[p] <- list( data.frame(peak_num=p,t=t_new, c=c_new, cutoff=cutoff, AUC=AUC) )
        #plot(c_new~t_new)
      }
    }
  return( do.call(rbind, peaks) )
  }
}  

#' Helper function to get aggregate statistics
#'
#' @param data dataframe
#' @param name new name for variable being produced
#' @param func function to be used in aggregate()
#' @param add_ds optional dataframe - if included, the new variable will be added to the supplied dataframe 
#' @param c_var name of the conc_var (e.g. concentration)
#' @param by_var vector of by_var (usually use by_var_v )
#' @return return either new variable by itself or in the dataframe listed in add_ds 
#' @export
#' @examples
#' 
#' ds <- getSumStat(hormone, name='mean_var', func= function(x) mean(x,na.rm=T), by_var=c('sp','sex'), c_var='conc' )
#' ds

getSumStat <- function(data=ds,name='mean', func=function(x)mean(x,na.rm=T), add_ds, c_var=conc_var,by_var=by_var_v ){
    ds1 <- aggregate(data[,c_var], by = data[c(by_var)], FUN = func )
      names(ds1)[ncol(ds1)] <- name
      if(!missing(add_ds)){ ds1 <- merge(add_ds,ds1,all=T)}
      return(ds1)
  }    

#' Helper function that converts all factors to string in a dataset
#'
#' @param data dataframe
#' @return dataframe with all factors converted to strings 
#' @export
#' @examples
#'
#' sapply(iris,class)
#' ds <- ridFactor(iris)
#' sapply(ds,class)
#' 

ridFactor <- function(data){ 
  i <- sapply(data, is.factor)
  data[i] <- lapply(data[i], as.character)
  return(data)
}
