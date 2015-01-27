#' Longitudinal graph with breaks 
#' 
#' @param x hormLong object (produced from hormBaseline) [required]
#' @param break_cutoff the maximum number of days between consecutive points. 
#' Above this cutoff value, a break is created. [default = 40]
#' @param break_buffer size of the gap (in number of days) between data groups.  Larger values will
#' create larger spaces between data groups. [default = 60]
#' @param date_format the format of the date variable on x-axis. Default is 01-Jan format. See Appendix 1 in help manual 
#' for examples of other formats [default = '\%d-\%b']
#' @param plot_per_page the number of plot panels per page, by row. [default = 4]
#' @param plot_height  the height of individual plot panels (in inches).  Pdf page height is determined
#' by both plot_per_page and plot_height. [default = 2]
#' @param plot_width  the width of the pdf page.[default = 6]
#' @param save_plot indicates whether to save plot as a file [default = TRUE]
#' @return nothing  Produces a pdf file saved at current working directory
#' @export
#' 
#' @examples
#' 
#' 
#' result <- hormBaseline(data=hormElephant, criteria=2, by_var='Ele, Hormone', time_var='Date', 
#'              conc_var='Conc_ng_ml' , event_var='Event')
#' hormPlotBreaks( result ) 
#'# compare to regular hormPlot, especially Ele1; Cortisol
#' hormPlot( result ) 



hormPlotBreaks <- function(x, break_cutoff=40, break_buffer=60, date_format='%d-%b',
                           plot_per_page=4, plot_height=2, plot_width=6, save_plot=TRUE){

#--- checks for missing values ---#
  graphics.off() # just to make sure no devices are open
  
  checkPlotOpts(plot_per_page, plot_width, plot_height, save_plot, d=date_format)

  checkPlotBreaks(break_cutoff, break_buffer)
  
#--- set-up data ---#
  require(lubridate)
  by_var_v <- cleanByvar(x$by_var) 
  time_var <- x$time_var
  conc_var <- x$conc_var
  data <- x$data
  data <-ridFactor(data)
  data$plot_title <- getPlotTitle(data,by_var=by_var_v)

#-- check for missing data in by_var_v, time_var, hormone_var --#
   data <- checkPlotMissing(data, var_list=c(by_var_v,time_var) )

#-- prepare data by geting cutoff and making date seconds 
  data <- getSumStat(data=data,name='cutoff', func= function(y) getCutoff(y, criteria=x$criteria ), add_ds=data, by_var=by_var_v, c_var=conc_var )
  if( is.Date(data[,time_var]) ){ data[,time_var] <- as.numeric( data[,time_var] ) * 3600*24 
    }else{  data[,time_var] <- as.numeric( data[,time_var] ) } 


#--- prep data ---#
  #--- calculate breaks using break_cutoff ---#
    data$row_id <- 1:nrow(data)
    tmp <- do.call(rbind,by(data,data$plot_title,FUN=function(x, time=time_var) cbind(row_id=x$row_id,diff=c(-99, diff(x[,time]))  ) ))
    data <- merge(data,tmp,all.x=T)
    data$brk <- 0
    for(i in 1:nrow(data)){
      if(data$diff[i] == -99){brk_num=1}
      if(data$diff[i]/(3600*24) > break_cutoff){brk_num=brk_num+1}
      data$brk[i] <- brk_num
    }

  #--- calculate adjusted time values --#
    getBreaks <-function(ds,x,y){ 
      min_brk <- aggregate(ds[x],ds['brk'],function(x) min(x) )
        names(min_brk)[2] <- 'x_min'
      rg_brk <- aggregate(ds[x],ds['brk'],function(x) diff( range(x) )) 
        names(rg_brk)[2]<-'x_rg'
      if( nrow(rg_brk) > 1 ){ rg_brk$x_acc <- cumsum(c(0,rg_brk$x_rg[1:(nrow(rg_brk)-1)]) )
        }else{rg_brk$x_acc <- 0} 
       adj_brk <- merge(min_brk,rg_brk)
       adj_brk$x_diff <- adj_brk$x_min - (adj_brk$x_acc + (adj_brk$brk-1)*break_buffer*(3600*24))
       ds <- merge(ds,adj_brk[,c('brk','x_diff')],by='brk', all.x=T )
       ds$x_adj <- ds[,x] - ds$x_diff
      return(ds)
    }
    ds1 <- do.call(rbind,by(data,data$plot_title,FUN=function(k,t=time_var,c=conc_var) getBreaks(k,t,c) ) ) 

#--- plot out the data
  #--- create plots ---#
  if( save_plot ){
    pdf('hormPlotBreak.pdf', height=plot_per_page * plot_height, width = plot_width )
  }

  par(mfrow=c(plot_per_page,1), mar=c(2,4,2,0.5),oma=c(2,2,2,0))
  for( i in unique(data$plot_title) ){
    ds_sub <- ds1[data$plot_title==i, ]
    baseline <- getCutoff( ds_sub[ds_sub$conc_type=='base',conc_var], criteria=x$criteria )

   plot(ds_sub[,conc_var]~ds_sub$x_adj, type='n', xaxt='n',ylab=conc_var, xlab=NA) 
   mtext(unique(ds_sub$plot_title),side=3,line=0.25)
     abline(h = baseline, lty=2)

   for(b in 1:max(ds_sub$brk)){
      ds_sub1 <- ds_sub[ds_sub$brk==b,] 
      points(ds_sub1$x_adj,ds_sub1[,conc_var], pch=19)
      lines(ds_sub1$x_adj,ds_sub1[,conc_var])
    #  text(ds_sub1$x_adj,-2+rnorm(length(ds_sub1$x_adj)),labels=ds_sub1$x,cex=0.6)
        p_at <-pretty( c(min(ds_sub1$x_adj), max(ds_sub1$x_adj), n=ceiling(12*nrow(ds_sub1)/nrow(ds_sub)) ) )
        p_at <- p_at[ p_at >= min(ds_sub1$x_adj)]
        p_at <- p_at[ p_at <= max(ds_sub1$x_adj)]
        if( length(p_at) < 2) p_at <- c( min(ds_sub1$x_adj), max(ds_sub1$x_adj ) )
        axis(1,at=p_at, labels = 
               format( as.POSIXct( p_at+ds_sub1$x_diff[1] ,origin='1970-01-01'), date_format) )

    events <- getEventInfo(ds_sub1, x$event_var, 'x_adj')
    
    plotEventInfo(events, t='x_adj', e=x$event_var) 

    } # then brk loop
  } # end plot_title  

  if( save_plot ){  
      dev.off()
      cat( paste0('\n *********\nNote: plots are saved at: \n', getwd(),'/hormPlotBreak.pdf \n***** \n\n')  )
  }
} # end function


