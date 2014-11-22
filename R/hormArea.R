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
#' result <- hormBaseline(data=hormone, criteria=0, by_var='sp, sex, id', time_var='date', conc_var='conc' )
#' hormArea(result)

hormArea <- function(x, method='trapezoid',plot_per_page=4){
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
  ds$plot_title <- getPlotTitle(ds)
  ds <- getSumStat(data=ds,name='cutoff', func= function(y) hormCutoff(y, criteria=x$criteria ), add_ds=ds )


#--- peak analysis ---#
  #-- shift curve by cutoff ---#
    ds <- ds[ do.call(order, ds[,c(by_var_v,time_var)]),]

#   k1 <- ds[ds$plot_title==unique(ds$plot_title)[1], ]
#   getPeakInfo(k1)

    ds_pk <- do.call("rbind", as.list(by(ds, ds[,by_var_v], 
                            function(k) data.frame(getPeakInfo(k),plot_title=unique(k$plot_title )) ) ) ) 
    ds_pk <- ds_pk[ds_pk$peak_num>0,]  # get rid individuals with no peaks       

  data<-ds
  #-- plot out results ---#
    for(i in unique(ds$plot_title)){
      ds_h <- ds[ds$plot_title==i,]
      ds_pk1 <- ds_pk[ds_pk$plot_title==i,]
      ggplot(ds_h) + geom_line(aes_string(time_var,conc_var)) + 
        geom_point(data=ds_pk1, aes(as.Date(t),c+cutoff)) +
        geom_ribbon(data=ds_pk1, aes(as.Date(t),ymin=cutoff,ymax=c+cutoff),alpha=0.3) + 
        geom_hline(aes(yintercept=cutoff))  
  #      geom_text(data=ds_pk,aes(x=t,y=c,label=peak_num),vjust=0.3) 
    }
    
  if( save_plot ){
    pdf('hormArea.pdf', height=plot_per_page * plot_height, width = plot_width )
  }

  par(mfrow=c(plot_per_page,1), mar=c(2,4,2,0.5),oma=c(2,2,2,0))
  for( i in unique(data$plot_title)){
    ds_sub <- data[data$plot_title==i, ]
    pk <- ds_pk[ds_pk$plot_title==i, ]
    baseline <- ds_sub$cutoff[1]
    ds_sub <- ds_sub[!is.na(ds_sub[,conc_var]),]
    
    #--- set up scales
      if( yscale=='free'){
        ymin <- min(ds_sub[,conc_var])*0.95
        ymax <- max(baseline, max(ds_sub[,conc_var]) )*1.1
      }else{
        ymin <- min(data[,conc_var],na.rm=T)*0.95
        ymax <- max(data[,conc_var],na.rm=T)*1.1
      }
      if( xscale=='free'){
        xmin <- min( ds_sub[,time_var],na.rm=T )
        xmax <- max( ds_sub[,time_var],na.rm=T )
      }else{
        xmin <- min( data[,time_var],na.rm=T)
        xmax <- max( data[,time_var],na.rm=T)
      }    
    #--- main plot
    plot(ds_sub[,conc_var] ~ ds_sub[,time_var], type='l',xlim=c(xmin,xmax), ylim=c(ymin,ymax), 
          xlab=NA, ylab=conc_var)
      points(ds_sub[,time_var], ds_sub[,conc_var],pch=19)
      mtext(unique(ds_sub$plot_title),side=3,line=0.25)
      abline(h = baseline, lty=2)
      if( max(pk$peak_num)>0 ){
        for(p in 1:max(pk$peak_num) ){
          pk1 <- pk[ pk$peak_num==p, ]
          t_order <- c(pk1$t,rev(pk1$t))
          c_order <- c(rep(0,length(pk1$t)),rev(pk1$c)) + pk1$cutoff
          with(pk1, polygon(t_order,c_order,col='grey') )
          text( mean(pk1$t), pk1$cutoff, labels=p, adj = c(0,1))
        }
      }
  }
  if( save_plot ){  dev.off()
      cat( paste0('\n *********\nNote: plots are saved at: \n', getwd(),'/hormPlot.pdf \n***** \n\n')  )
  }


  #--- produce table
    


}



getPeakInfo <- function(k){
  if( !any(k$conc_type=='peak')){
    return(data.frame(peak_num=0,t=0, c=0, cutoff=k$cutoff[1], AUC=0))
  }else{
    cutoff <- k$cutoff[1] 
    t<-as.numeric(k$date)
    c<-k$conc - cutoff
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
        end <-  min(length(t), max(which(peak_num==p))+1 )  
        
      #-- calculate point that slope crosses cutoff --#
        m <- diff(c[st:end])/diff(t[st:end])
        b <- c[st:(end-1)] - m*t[st:(end-1)]
        t0 <- -b[unique(c(1,length(b)))]/m[unique(c(1,length(m)))]
        if( length(m)==1 & st==1){ t0 <-c(t[st],t0)}   
        if( length(m)==1 & st>1){ t0  <-c(t0,t[end])}   
        if( length(m)>1 & st==1){ t0[1] <-c(t[st])  }
        if( length(m)>1 & end==length(t)){ t0[length(t0)] <-c(t[end])}   
      
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
  return( do.call(rbind, peaks))
  }
}  

