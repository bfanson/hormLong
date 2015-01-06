#' Helper function to calculate baseline cutoff
#' 
#' @param  x hormone concentration
#' @param  criteria the number of standard deviations above baseline 
#' @return cutoff value
#' @export
#' @examples
#' 
#' conc <- rnorm(100)
#' getCutoff(conc, 2)

getCutoff <- function(x, criteria){
  mean(x, na.rm=T)+sd(x, na.rm=T)*criteria 
}


#' Helper function to work with dates.  
#'  
#' Helper function to work with dates.  User needs to specify at least date_var and date_order. 
#' If time_var is included, then it will output a date_time variable
#' 
#' @param data dataset with the variables to convert [required]
#' @param date_var variable name for date column.  Form should be a standard date form 
#' (e.g. '12AUG2014','01/01/2014') [required]
#' @param time_var variable name for the time column in 24-hour format (e.g. '20:10:01', '01:04:01') [optional]
#' @param name name of the new date variable created by the function [default = 'datetime']
#' @param date_order specifies the day, month, year order in abbreviated form (e.g. 'ymd','mdy','dmy'). 
#' [default='dmy']
#' @return dataset with the new date/time variable
#' @export
#' @examples
#'
#'# Combining date and time  
#' date <- c('2014-01-01','2014-Mar-01')
#' time <- c('10:10:01', '20:30:23')
#' ds <- data.frame(date=date,time=time)
#' hormDate(ds,date_var='date',time_var='time', name='datetime', date_order='ymd' )
#' 
#'# Function is robust to date formatting as long as day, month, year order are the same  
#' date <- c('01-01-2014','01-Mar-2014','01/April/2014')
#' ds <- data.frame(date=date)
#' hormDate(ds,date_var='date',name='datetime', date_order='dmy' )


hormDate <- function(data, date_var, time_var, name='datetime', date_order='dmy') {
  require(lubridate)
  date_order <- tolower(date_order)
  if( is.numeric(data[,date_var]) ) data[,date_var] <- as.character(data[,date_var]) 
  if( missing(date_var)   ){stop('please provide date_var. time_var is optional ')}
  if( missing(date_order) ){stop('please provide date_order...e.g. ymd, dmy, mdy')}
  if( (date_order %in% c('ymd','ydm','mdy','dmy' ))==F ){stop('date_order must be ymd, ydm, mdy or dmy ')}

  if( !missing(time_var)){
    if(any( (grepl('pm',tolower(data[,time_var]))==T |  
                                grepl('am',tolower(data[,time_var]))==T)==T ) ){
      stop('check format of time_var.  It should be in 24-hour format with hours, minutes and seconds (e.g. 20:10:01)')
    }
  }
  
  #--- calculate new date --#
   if(!missing(date_var) & missing(time_var)){
      datetime <- do.call( date_order, list(data[,date_var]))
    }else if(!missing(date_var) & !missing(time_var) ){
      datetime <- as.POSIXct(paste( do.call( date_order, list(data[,date_var])),data[,time_var]), format="%Y-%m-%d %H:%M:%S", tz='UTC' )
    }
  
   data$rename_this <- datetime  
   names(data)[ncol(data)] <- name
   return(data)
}
  


#' Read hormone data from a csv file
#' 
#' @param none no arguments
#' @return imported dataset
#' @export
#' @examples
#' 
#' ds <- hormRead()

hormRead <- function( ){
  ds<- read.csv(choose.files(caption='Select hormone data file') ) 
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
        end_ind <- ifelse( k$conc_type[end]=='peak' ,1,0)
        
      #-- calculate point that slope crosses cutoff --#
        m <- diff(c[st:end])/diff(t[st:end])
        b <- c[st:(end-1)] - m*t[st:(end-1)]
        t0 <- -b[unique(c(1,length(b)))]/m[unique(c(1,length(m)))]
        if( length(m)==1 & st==1){ t0 <- c(t[st],t0)}   
        if( length(m)==1 & st>1){ t0  <- c(t0,t[end])}   
        if( st_ind==1){ t0[1] <- c(t[st])  }
        if( end_ind ){ t0[length(t0)] <- as.numeric(k[end,date])}
      
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


#' Helper function that plots Event info
#'
#' @param d_e events dataframe
#' @param t time_var name
#' @param e event_var name
#' @return nothing plots event text on graphs using par('user)[4] coordinate 
#' @export
#'

plotEventInfo <- function( d_e=events, t=x$time_var,  e=x$event_var )
  if( nrow(d_e)>0 ){
    ymax_e <- par('usr')[4]
    for(l in 1:nrow(d_e)){
      arrows(x0=d_e[l,t],x1 =d_e[l,t],y0=ymax_e*0.92,y1=ymax_e*0.76, length = 0.1)
      text(x=d_e[l,t],y=ymax_e*0.96, d_e[l,e])
    }
}

#' Helper function for get axes limits
#'
#' @param d_s ds_sub dataframe
#' @param d_f data dataframe
#' @param var response variable
#' @param scale xscale or yscale
#' @param base baseline measure which is needed if higher than max(data[,conc_var])
#' @return limits a vector of min and max   
#' @export
#'

getPlotlim <- function(d_s, d_f, var, scale, base=NA){
    if( scale=='free' ){
        a_min <- min( d_s[,var] )
        if( is.na(base) ){ a_max <- max( d_s[,var] ) 
        }else{ a_max <- max(base, max(d_s[,var]) ) }
        
    }else if(scale=='fixed') {
        a_min <- min(d_f[,var],na.rm=T)
        a_max <- max(d_f[,var],na.rm=T)
      }
  
  lim <- c(a_min, a_max)
  return( lim )
}

#' Helper function for get axes limits
#'
#' @param d_s ds_sub dataframe
#' @param d_f data dataframe
#' @param var response variable
#' @param scale xscale or yscale
#' @param base baseline measure which is needed if higher than max(data[,conc_var])
#' @return limits a vector of min and max   
#' @export
#'

plotAxes <- function(d_s, t, x_s, d_f  ){
  require(lubridate)
      if( is.numeric(d_s[,t]) ){ 
        axis(1)
      }else if( is.Date(d_s[,t]) ){
          ats <- seq( x_s[1], x_s[2], length.out = 5)
          axis.Date(1,at=ats, format=d_f)
      }else if( is.POSIXct(d_s[,t]) ){
          ats <- seq( x_s[1], x_s[2], length.out = 5)
          axis.POSIXct(1,at=ats,format=d_f)
      } 
}  
  

#' Helper function for checking class
#'
#' @param var ds_sub dataframe
#' @param class baseline measure which is needed if higher than max(data[,conc_var])
#' @return none
#' @export
#'

checkClass <- function( var, v_class) {
  if( class(var) != v_class ){
    if( v_class == 'hormLong' ){
      stop( 'Object needs to be hormLong.  Run hormBaseline() first' )
    }else{
      stop( paste0( 'Wrong class for ', var, '. It needs to be class ',v_class )   )
    }
  }
}



#' Helper function for checking plot options
#'
#' @param p_p plot_per_page
#' @param w   plot_width
#' @param h   plot_height
#' @param s   save_plot
#' @param x   xscale
#' @param y   yscale
#' @param d   date_format
#' @param l_b lower_bound
#' @return none
#' @export
#'

checkPlotOpts <- function(p_p_p, w, h, s, x=NA, y=NA,d ){
  if( !is.numeric(p_p_p) | p_p_p<1 ){
    stop( paste0('plot_per_page needs to be an numeric and greater than 0.  Currently it is ',p_p_p) )
  }
  if( !is.numeric(w) | w<0 ){
    stop( paste0('plot_width needs to be numeric and greater than 0.  Currently it is ',w) )
  }
  if( !is.numeric(h) | w<0 ){
    stop( paste0('plot_height needs to be numeric and greater than 0.  Currently it is ',h) )
  }
  if( !is.logical(s) ){
    stop( paste0('save_plot needs to be logical (T or F).  Currently it is ',save_plot) )
  }
  if(!is.na(y)){
    if( !(y %in% c('free','fixed') ) ){
      stop( paste0('yscale must be either "free" or "fixed", not "',y,'"') )
    }
  }
  if(!is.na(y)){
    if( !(x %in% c('free','fixed') ) ){
      stop( paste0('xscale must be either "free" or "fixed", not "',x,'"') )
    }
  }
  if( !is.character(d)  ){
    stop( paste0('date_format must be a character field (e.g. "%d/%m"), not "',d,'"') )
  }
  if( !grepl('\\%',d)  ){
    stop( paste0('check your date_format. No % sign detected. It should be like "%d/%m". It is currently "',d,'"') )
  }
}


#' Helper function for checking plot options for hormArea
#'
#' @param d   date_format
#' @param l_b lower_bound
#' @return none
#' @export
#'

checkPlotArea <- function(m, l_b){
  if( m != 'trapezoid' ){ 
      stop('no other method currently implemented. Please write method="trapezoid" ')
  }
  if( !(l_b %in% c('origin','baseline','peak') ) ){ 
      stop(paste0("lower_bound is incorrect.  It must be 'origin', 'baseline', 'peak': you wrote '", l_b,"'") )
    }
}

#' Helper function for checking plot options for hormPlotBreaks
#'
#' @param b_c  break_cutoff
#' @param b_b  break_buffer
#' @return none
#' @export
#'
checkPlotBreaks <- function(b_c, b_b){
  if( !is.numeric(b_c) | b_c < 0 ){
    stop( paste0('break_cutoff must be numeric and nonnegative.  It is currently ',b_c) )
  }
  if( !is.numeric(b_b) | b_b < 0 ){
    stop( paste0('break_buffer must be numeric and nonnegative.  It is currently ',b_b) )
  }
}

#' Helper function for checking plot options for hormPlotOverlap
#'
#' @param h_v  hormone_var
#' @param b_v  by_var
#' @return none
#' @export
#'

checkPlotOverlap <- function(h_v, b_v){
  b_v_v <- cleanByvar(b_v) 
  if( missing(h_v) ){
      stop(' hormone_var must be specified ')
  }
  if( !is.character(h_v) ){
      stop(' hormone_var must be character ')
  }
  if( !(h_v %in% b_v_v) ){
      stop(' hormone_var must be specified in the by_var ')
  }
}


#' Helper function for checking plot options for hormPlotOverlap
#'
#' @param h_v  hormone_var
#' @param b_v  by_var
#' @return none
#' @export
#'

checkPlotRatio <- function(h_v, h_n, h_d, b_v){

  b_v_v <- cleanByvar(b_v) 
  
  if( missing(h_v) ){
      stop(' hormone_var must be specified ')
  }
  if( !is.character(h_v) ){
      stop(' hormone_var must be character ')
  }
  if( !(h_v %in% b_v_v) ){
      stop(' hormone_var must be specified in the by_var ')
  }
  
  if(  missing(h_n) | missing(h_d)){
      stop('hormone_num and hormone_denom must be specified')
  }
  
  if( !(h_v %in% b_v_v) ){
    stop('hormone_var must have been specified in by_var when running hormBaseline')
  }
  if( !(h_n %in% unique(x$data[,h_v])) ){
    stop('hormone_num must be a level in hormone_var.  check spelling and capitalization')
  }  
  if( !h_d %in% unique(x$data[,h_v]) ){
    stop('hormone_denom must be a level in hormone_var.  check spelling and capitalization')
  }  
  
}


#' Helper function for getting Baseline for plotting
#'
#' @param d_s  ds_sub
#' @param crit x$criteria
#' @param conc conc_var
#' @return baseline 
#' @export
#' 
getBaseline <- function(d_s, crit, conc){
      if(!is.null(crit)){
         baseline <- getCutoff( d_s[d_s$conc_type=='base',conc], criteria=crit )
      }else{ baseline <- 1 }
      return(baseline)
    }
