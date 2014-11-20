#' Create longitudinal plot with smoothers
#' 
#' @param data dataset.
#' @param by_var by variable.
#' @param time_var time variable.
#' @param conc_var response variable.
#' @param group_var grouping variable.
#' @return hormLong object.
#' @export
#' @examples
#' 
#' result <- hormSmooth(data=hormone2, by_var='sp, id', group_var='horm_type', time_var='date', conc_var='conc')
 

hormSmooth <- function(data,by_var,group_var,time_var,conc_var, plot_per_page=4,
                       plot_height=2, plot_width=6, yscale='free', xscale='free',
                       smoothness=c(0.7,0.7), shape=c(19,20), colour=c('red','blue'), line_type=c(1,2),
                       line_width=c(1,1),   ...) {
  
#--- add in checks ---#
  
  
#-- set-up ---#
  by_var_v <- by_var_v <- cleanByvar(by_var) 
  time_var <- time_var
  conc_var <- conc_var
  data <- data[ do.call(order, data[c(by_var_v,group_var,time_var)]), ]
  
  #-- create titles--#
    paste1 <- function(...) paste(...,sep='; ')
    data$plot_title <- do.call(paste1, data[,c(by_var_v)] )
    
#--- create plots ---#
  if( save_plot ){
    pdf('hormSmooth.pdf', height=plot_per_page * plot_height, width = plot_width )
    cat( paste0('\n *********\nNote: plots are saved at: \n', getwd(),'/hormSmooth.pdf \n***** \n\n')  )
  }

  par(mfrow=c(plot_per_page,1), mar=c(2,4,2,0.5),oma=c(2,2,2,0))
  for( i in unique(data$plot_title)){
    ds_sub <- data[data$plot_title==i, ]
    events <- ds_sub[ !is.na(ds_sub$event) & ds_sub$event!='',c('event',time_var)]
    ds_sub <- ds_sub[!is.na(ds_sub[,conc_var]),]
    
    #--- set up scales
      if( yscale=='free'){
        ymin <- min(ds_sub[,conc_var])*0.95
        ymax <- max(max(ds_sub[,conc_var]) )*1.2
      }else{
        ymin <- min(data[,conc_var],na.rm=T)*0.95
        ymax <- max(data[,conc_var],na.rm=T)*1.2
      }
      if( xscale=='free'){
        xmin <- min( ds_sub[,time_var],na.rm=T )
        xmax <- max( ds_sub[,time_var],na.rm=T )
      }else{
        xmin <- min( data[,time_var],na.rm=T)
        xmax <- max( data[,time_var],na.rm=T)
      }    
    
    #--- main plot ---#
    plot(ds_sub[,conc_var] ~ ds_sub[,time_var], type='n',xlim=c(xmin,xmax), ylim=c(ymin,ymax), 
          xlab=NA, ylab=conc_var)
      mtext(unique(ds_sub$plot_title),side=3,line=0.25)
      loop <- 1
      for(g in levels(ds_sub[,group_var]) ){
        ds_sub1 <- ds_sub[ds_sub$horm_type==g,]
        points(ds_sub1[,time_var], ds_sub1[,conc_var],pch=shape[loop],col=colour[loop] )
        lines(loess.smooth(ds_sub1[,time_var],ds_sub1[,conc_var], span=smoothness[loop]),col=colour[loop]  ) 
        loop <- loop + 1
      }
      legend('topleft',legend=levels(ds_sub[,group_var]),col=colour,pch=shape, lty=line_type,bty='n')
  }
  if( save_plot ){  dev.off() }
}
