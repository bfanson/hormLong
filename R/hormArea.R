#' Get area under the curve 
#' 
#' @param x hormLong object (produced from hormBaseline)
#' @param plot_per_page the number of plot panels per page, by row.
#' @param save_plot indicates whether 
#' @param plot_height  the height of individual plot panels (in inches).  Pdf page height is determined
#' by both plot_per_page and plot_height.
#' @param plot_width  the width of the pdf page.
#' @param yscale  determines if y-axis should be free ('free') to change for each panel or
#' remain the same ('fixed') for all panels 
#' @param xscale  determines if x-axis should be free ('free') to change for each panel or
#' remain the same ('fixed') for all panels 
#' @return NA
#' @export
#' @examples
#' 
#' Citation: Variation within and between Birds in Corticosterone Responses of Great Tits (Parus major)
#' General and Comparative Endocrinology, Volume 125, Issue 2, February 2002, Pages 197–206
#' J.F. Cockrem, B. Silverin
#' 
#' result <- hormBaseline(data=hormone, criteria=1, by_var='sp, sex, id', time_var='date', conc_var='conc' )
#' hormArea(result)

hormArea <- function(x, method='trapezoid'){
  stop('function under development')
  if( method!='trapezoid'){ 
    stop('no other method currently implemented ')
  }
  by_var_v <- by_var_v <- cleanByvar(x$by_var) 
  time_var <- x$time_var
  conc_var <- x$conc_var
  ds <- x$data
  ds <- ds[!is.na(ds[,conc_var]),]
  i <- sapply(ds, is.factor)
  ds[i] <- lapply(ds[i], as.character)

#--- AUC funciton ---#
  shiftByCutoff <- function(t, c,cutoff){ 
      c <- c - cutoff
      m <- diff(c)/diff(t)
      b <- c[1:(length(c)-1)] - m*t[1:(length(c)-1)]
      t0 <- -b/m
      t_new <- as.character(t)
      c_new <- pmax(c,0)

      if(c[1]<0 & c[1+1]>0 ) t_new[1]   <- as.character(t0[1]) # increasing 
      if(c[1]>0 & c[1+1]<0 ) t_new[1+1] <- as.character(t0[1]) # decreasing

      for(i in 2:(length(t)-1) ){
        if(c[i]<0 & c[i+1]>0 ) t_new[i]   <- as.character(t0[i]) # increasing 
        if(c[i]>0 & c[i+1]<0 ) t_new[i+1] <- as.character(t0[i]) # decreasing
        if( c[i-1]>0 & c[i]<0 & c[i+1]>0 ) t_new[i] <- 
              paste(as.character(t0[i]),as.character(t0[i+1]),sep='-')  # both decresing and increasing
      }
    return( data.frame(t_adj=t_new,c_adj=c_new) )
  }

#--- peak analysis ---#
  #-- shift curve by cutoff ---#
    ds <- ds[ do.call(order, ds[,c(by_var_v,time_var)]),]
    tmp <- do.call("rbind", as.list(by(ds, ds[,by_var_v], 
                            function(k) cbind( shiftByCutoff(t=as.numeric(k[,time_var]),
                                                      c=k[,conc_var],10), row_id=k$row_id ) ) ) 
            )
    ds1 <- merge(ds, tmp, all.x=T )


  #-- number peaks ---#
    calcPeaks <- function(t,c){
        peak_num <- rep(0,length(t))
        num <- 0
        if(c[1] > 0 | (c[1]==0 & c[2]>0) ){num <- 1}
        peak_num[1] <- num
        
        for(i in 2:(length(t)-1) ){
          if( (c[i]==0 & c[i+1]>0) ){ num <- num + 1}
          peak_num[i] <- num
          if( c[i]==0 & c[i-1]==0 & c[i+1]==0  ) peak_num[i] <- NA 
        }
        if( c[length(t)]>0 & c[length(t)-1]==0 )  num <- num + 1 
        peak_num[length(t)] <- num
        return( peak_num )
      }
  tmp <- do.call("rbind", as.list(by(ds1, ds1[,by_var_v], 
                            function(k) cbind( peak_num= calcPeaks(t=k$t_adj,c=k$c_adj), row_id=k$row_id ) ) ) )
  ds1 <- merge(ds1, tmp, all.x=T)

  #--- get AUC ---#
   hold <- ds1[!is.na(ds1$peak_num),]
   tmp <- do.call("rbind", as.list(by(hold, hold[,c(by_var_v,'peak_num')], 
                            function(k) cbind( AUC=getAUC(t=k$t_adj, c=k$c_adj), row_id=k$row_id ) ) ) )
   ds1 <- merge(ds1, tmp, all.x=T)

   paste1 <- function(...) paste(...,sep='; ')
   ds1$plot_title <- do.call(paste1, ds1[,c(by_var_v)] )

  #-- plot out results ---#
    ggplot(ds1) + geom_line(aes_string(time_var,conc_var)) + geom_point(aes(as.Date(t_adj),c_adj+10)) +
      geom_ribbon(aes(x=as.Date(t_adj),ymin=10,ymax=c_adj+10),alpha=0.3) + 
      geom_hline(yintercept=10)+ 
      geom_text(aes(x=as.Date(t_adj),y=c_adj*1.1,label=peak_num),vjust=0.3) +
      facet_grid(plot_title~.)

}

getAUC <- function(t,c){
    require(zoo)
    (AUC <- sum(diff(t)*rollmean(c,2)))
  }  

getPeakInfo <- function(k){
  cutoff=10
  t<-as.numeric(k$date)
  c<-k$conc - cutoff
  c_t <- k$conc_type
  row_id <- 1:nrow(k)
 
  #-- assign peak number ---#
    num <- 0
    peak_num <- rep(0,length(t) )
    if(c_t[1]=='peak'){ num <- num+1; peak_num[1] <- num} 
    for(i in 2:length(t) ){
      if( c_t[i]=='peak' & c_t[i-1]=='base' ){ num <- num+1; peak_num[i] <- num}
    }
  
  #-- grab each peak and calculate stats --#
  if(max(peak_num)>0){
    peaks <- list( )
    for( p in 1:max(k$peak_num) ){
      st <-  max(1, min( which(peak_num==p) )-1 )
      end <-  min(length(t), which(peak_num==p)+1 )  
    #-- cap t and c in case a peak starts or ends the sampling period --#
        
    
    #-- calculate point that slope crosses cutoff --#
      m <- diff(c[st:end])/diff(t[st:end])
      b <- c[st:(end-1)] - m*t[st:(end-1)]
      t0 <- -b[unique(c(1,length(b)))]/m[unique(c(1,length(m)))]
      if( length(m)==1 & st==1){ t0 <-c(t[st],t0)}   
      if( length(m)==1 & st>1){ t0  <-c(t0,t[end])}   
    
      t_new <- rep(t[st:end],each=2)
      t_new[c(1:2,(length(t_new)-1):length(t_new))] <- rep(t0,each=2) 
        
      c_new <- rep(c[st:end],each=2)
      c_new[c(1,length(c_new))] <- 0
      c_new[c(2,(length(c_new)-1))] <- pmax(c_new[c(2,(length(c_new)-1))],0)
      
      AUC <- getAUC(t_new,c_new)
      peaks[p] <- list( data.frame(peak_num=p,t=t_new, c=c_new, cutoff=cutoff, AUC=AUC) )
    }
  }
  return( do.call(rbind, peaks))
}  

