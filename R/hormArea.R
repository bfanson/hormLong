#' Plot the area under the curve (AUC)  
#' 
#' @param x hormLong object (produced from hormBaseline) [required]
#' @param lower_bound the lower limit to calculate the area under the curve. This can be
#' 'origin', 'baseline', or 'peak'.  Origin is all values above 0.  Baseline is all values above 
#' baseline mean.  Peak is all values above peak cutoff. [default = 'origin']
#' @param method the AUC method to use. Only Trapezoid and spline method has been implemented. Units
#' for AUC is conc units * day. [default = 'trapezoid']
#' @param date_format the format of the date variable on x-axis. Default is 01-Jan format. See Appendix 1 in help manual 
#' for examples of other formats [default = '\%d-\%b']
#' @param xscale  determines if x-axis is free ('free') to change for each panel or
#' remain the same ('fixed') for all panels  [default = 'free']
#' @param yscale  determines if y-axis is free ('free') to change for each panel or
#' remain the same ('fixed') for all panels  [default = 'free']
#' @param plot_per_page the number of plot panels per page, by row. [default = 4]
#' @param save_plot indicates whether [default = TRUE]
#' @param plot_height  the height of individual plot panels (in inches).  Pdf page height is determined
#' by both plot_per_page and plot_height. [default = 2]
#' @param plot_width  the width of the pdf page. [default = 6]
#' @return nothing is returned.  Pdf file is produced and a summary table 
#' @export
#' 
#' @details 
#' Citation for the AUC used in this example:\cr\cr
#'   Cockrem JF and Silverin B. 2002. Variation within and between Birds in Corticosterone 
#'   Responses of Great Tits (Parus major).  General and Comparative Endocrinology 125:197.
#' 
#' @examples 
#' 
#' result <- hormBaseline(data=hormLynx, criteria=2, by_var='AnimalID, Hormone', 
#'             time_var='datetime', conc_var='Conc', event_var='Events' )
#' 
#'# AUC only for values above cutoff threshold 
#'#  defined in hormBaseline
#' hormArea(result,lower_bound='peak')    
#' 
#'# AUC only for values above baseline mean 
#' hormArea(result,lower_bound='baseline') 
#' 
#'# AUC only for values above 0 
#' hormArea(result,lower_bound='origin')   

hormArea <- function(x, lower_bound = 'origin', method='trapezoid', date_format='%d-%b', 
                     xscale='free', yscale='free',
                     plot_per_page=4, plot_height=2, plot_width=6, save_plot=T){
#-- initial checking ---# 
  graphics.off() # just to make sure not devices are open
  
  checkClass(x, 'hormLong')

  checkPlotOpts(plot_per_page, plot_width, plot_height, save_plot, xscale, yscale, date_format)
  checkPlotArea(method, lower_bound)


#--- get hormLong object and break up ---#  
  by_var_v <- cleanByvar(x$by_var) 
  time_var <- x$time_var
  conc_var <- x$conc_var
  data <- x$data
  data <- data[!is.na(data[,conc_var]),]
  data <-ridFactor(data)
  data$plot_title <- getPlotTitle(data,by_var=by_var_v)
  if( lower_bound =='origin'){
    data$cutoff <- 0
    data$conc_type[ data[,conc_var] >=0 ] <- 'peak'
  }else if(lower_bound=='baseline'){
    getBaseLine <- function(c){
      return( mean( c, na.rm=T ) )
    }
    data <- getSumStat(data=data[data$conc_type=='base',], name='cutoff', func= function(y) getBaseLine(y), 
                                                by_var=by_var_v, c_var=conc_var,add_ds=data  )
    data$conc_type[ data[,conc_var] >= data$cutoff ] <- 'peak'
 
  }else if(lower_bound=='peak'){
    data <- getSumStat(data=data[data$conc_type=='base',], name='cutoff', func= function(y) getCutoff(y, criteria=x$criteria ), 
                                                by_var=by_var_v, c_var=conc_var,add_ds=data  )
  }

#--- deal with duplicate time points ---#
  if( any( duplicated( data[,c(by_var_v,time_var)] )==T ) ){
    cat('Warning: There are duplicated date/time, so only the lowest concentration is kept.
          The duplicated rows are listed below: \n\n')
    print( data[ duplicated( data[,c(by_var_v,time_var)] )==T,c(by_var_v,time_var,conc_var)] )
    data <- data[ duplicated( data[,c(by_var_v,time_var)] )==F,]
  }

#--- peak analysis ---#
  #-- shift curve by cutoff ---#
    # k <- data[data$plot_title==unique(data$plot_title)[1],]
    ds_pk <- do.call("rbind", as.list(by(data, data[,by_var_v], 
                          function(k, date=time_var,conc=conc_var,m=method) data.frame(getPeakInfo(k,date,conc,m),
                                        plot_title=unique(k$plot_title ) ) ) ) ) 
    ds_pk <- ds_pk[ds_pk$peak_num>0,]  # get rid individuals with no peaks       

#-- plot out results ---#
  
  if( save_plot ){
    pdf('hormArea.pdf', height=plot_per_page * plot_height, width = plot_width )
  }

  par(mfrow=c(plot_per_page,1), mar=c(2,4,2,0.5),oma=c(2,2,2,0))
  for( pt in unique(data$plot_title) ){
    ds_sub <- data[data$plot_title==pt, ]
    pk <- ds_pk[ds_pk$plot_title==pt, ]
    baseline <- ds_sub$cutoff[1]
    
    events <- getEventInfo(ds_sub, x$event_var, time_var)
    
    ds_sub <- ds_sub[!is.na(ds_sub[,conc_var]),]
    
    #--- set up scales
      ylim_expand=1.1
      if( method=='spline'){ ylim_expand = 1.3}
      y_lim <- getPlotlim(d_s=ds_sub, d_f=data, var=conc_var, scale=yscale, max_expand=ylim_expand, base=baseline)
      x_lim <- getPlotlim(d_s=ds_sub, d_f=data, var=time_var, scale=xscale)
      if( lower_bound=='origin' ) y_lim[1] <- 0

    #--- main plot
    plot(ds_sub[,conc_var] ~ ds_sub[,time_var], type='l',xlim=x_lim, ylim=y_lim, 
          xlab=NA, ylab=conc_var, xaxt='n')
      points(ds_sub[,time_var], ds_sub[,conc_var],pch=19)
      mtext(unique(ds_sub$plot_title),side=3,line=0.25)
      abline(h = baseline, lty=2)
    
      plotAxes(ds_sub, time_var, x_lim,d_f=date_format)
    
      if( length( unique(pk$peak_num) ) > 0 ){
        for(p in 1:max(pk$peak_num) ){
          pk1 <- pk[ pk$peak_num==p, ]
          if( method=='trapezoid'){
            t_order <- c(pk1$t,rev(pk1$t))
            c_order <- c(rep(0,length(pk1$t)),rev(pk1$c)) + pk1$cutoff
          }
          if( method=='spline'){
            hold <- data.frame( x = pk1$t, y = pk1$c) # get rid of double points at the same date-time
            hold1 <-  do.call(rbind, list(by(hold,hold[1],function(x) max(x['y']))) )
            hold <- data.frame( x=as.numeric(colnames(hold1)),y=hold1[1,] )
            new_val <- spline(hold$x, hold$y, method='natural' )
            t_order <- c( new_val$x, rev(new_val$x) )
            c_order <- c( rep(0,length(new_val$x)), rev(new_val$y))  + unique(pk1$cutoff)
          }

          with(pk1, polygon(t_order,c_order,col=adjustcolor('grey', alpha=0.7), border='dark grey') )
          text( mean(pk1$t), max(c_order), labels=p, adj = c(0,-0.5))
        }
      }
      plotEventInfo(events, t=time_var, e=x$event_var) 
  } # end p loop
  if( save_plot ){  
      dev.off()
      cat( paste0('\n *********\nNote: plots are saved at: \n', getwd(),'/hormArea.pdf \n***** \n\n')  )
  }


  #--- produce table ---#
    require(lubridate)
    ds_tbl <- unique( ds_pk[,c('plot_title','peak_num','AUC')])    
    ds_tbl <- merge(ds_tbl, unique(data[,c(by_var_v,'plot_title')]), all.x=T)
    if(  is.Date(data[,time_var]) | is.numeric(data[,time_var])  ){
      ds_tbl$AUC <- round(ds_tbl$AUC,3) # get conc * day
    }else{
      ds_tbl$AUC <- round(ds_tbl$AUC,3) / (3600 * 24)  # get conc * day
    }
    ds_out <- ds_tbl[,c(by_var_v,'peak_num','AUC')]
    names(ds_out)[which(names(ds_out)=='AUC')] <-'peak AUC'    

    #-- merge in total AUC ---#
#     do.call("rbind", as.list(by(data, data[,by_var_v], 
#                   function(k) data.frame(totalAUC=getAUC(k[,time_var],k[,conc_var]),plot_title=unique(k$plot_title ) ) )))

    #--- write dataset ---#
      write.csv(ds_out,file='hormAUCtable.csv',quote=F, row.names=F)
      cat( paste0('\n *********\nNote: table saved at: \n', getwd(),'/hormAUCtable.csv \n***** \n\n')  )

}





