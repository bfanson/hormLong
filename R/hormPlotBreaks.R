#' Longitudinal graph with breaks 
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
#' 
#' 

hormPlotBreaks <- function(x, break_cutoff=40, break_buffer=60,
                           plot_per_page=4, save_plot=TRUE, plot_height=2, plot_width=6){

  #stop('function under development')

  by_var_v <- cleanByvar(x$by_var) 
  time_var <- x$time_var
  conc_var <- x$conc_var
  data <- x$data
  data <- data[!is.na(data[,conc_var]),]
  data <-ridFactor(data)
  data$plot_title <- getPlotTitle(data,by_var=by_var_v)
  data <- getSumStat(data=data,name='cutoff', func= function(y) getCutoff(y, criteria=x$criteria ), add_ds=data, by_var=by_var_v, c_var=conc_var )

#--- prep data ---#
  #--- calculate breaks using break_cutoff ---#
    data$row_id <- 1:nrow(data)
    tmp <- do.call(rbind,by(data,data$plot_title,FUN=function(x, time=time_var) cbind(row_id=x$row_id,diff=c(-99, diff(x[,time]))  ) ))
    data <- merge(data,tmp,all.x=T)
    data$brk<-0
    for(i in 1:nrow(data)){
      if(data$diff[i] == -99){brk_num=1}
      if(data$diff[i] > break_cutoff){brk_num=brk_num+1}
      data$brk[i] <- brk_num
    }
  #--- calculate adjusted time values --#
    getBreaks <-function(ds,x,y){ 
      min_brk <- aggregate(ds[x],ds['brk'],min)
        names(min_brk)[2]<-'x_min'
      rg_brk <- aggregate(ds[x],ds['brk'],function(x) diff( range(x))) 
        names(rg_brk)[2]<-'x_rg'
      if( nrow(rg_brk)>1){rg_brk$x_acc <- cumsum(c(0,rg_brk$x_rg[1:(nrow(rg_brk)-1)]) )
        }else{rg_brk$x_acc <- 0} 
       adj_brk <- merge(min_brk,rg_brk)
       adj_brk$x_diff <- adj_brk$x_min-(adj_brk$x_acc + (adj_brk$brk-1)*break_buffer)
       ds <- merge(ds,adj_brk[,c('brk','x_diff')],by='brk', all.x=T )
       ds$x_adj <- ds[,x] - ds$x_diff
      return(ds)
    }
    ds1 <- do.call(rbind,by(data,data$plot_title,FUN=function(k,t=time_var,c=conc_var) getBreaks(k,t,c) ) ) 

#--- plot out the data#--- create plots ---#
  if( save_plot ){
    pdf('hormPlotBreak.pdf', height=plot_per_page * plot_height, width = plot_width )
  }

  par(mfrow=c(plot_per_page,1), mar=c(2,4,2,0.5),oma=c(2,2,2,0))
  for( i in unique(data$plot_title)){
    ds_sub <- ds1[data$plot_title==i, ]
    baseline <- getCutoff( ds_sub[ds_sub$conc_type=='base',conc_var], criteria=x$criteria )
    ds_sub <- ds_sub[!is.na(ds_sub[,conc_var]),]

   plot(ds_sub[,conc_var]~ds_sub$x_adj, type='n', xaxt='n',ylab=conc_var) 
     mtext(unique(ds_sub$plot_title),side=3,line=0.25)
     abline(h = baseline, lty=2)

   for(b in 1:max(ds_sub$brk)){
      ds_sub1 <- ds_sub[ds_sub$brk==b,] 
      points(ds_sub1$x_adj,ds_sub1[,conc_var])
      lines(ds_sub1$x_adj,ds_sub1[,conc_var])
    #  text(ds_sub1$x_adj,-2+rnorm(length(ds_sub1$x_adj)),labels=ds_sub1$x,cex=0.6)
        p_at <-pretty( c(min(ds_sub1$x_adj), max(ds_sub1$x_adj), n=12*nrow(ds_sub1)/nrow(ds_sub) ) )
        axis(1,at=p_at, labels = c(p_at[1:(length(p_at)-1)]+ds_sub1$x_diff[1],NA) )

#     if(!is.null(x$event_var)){
#       events <- ds_sub1[ !is.na(ds_sub1[,x$event_var]) & ds_sub1[,x$event_var]!='',c(x$event_var,time_var)]
#       }else{events <- data.frame()}
    
    } # then brk loop
  } # end plot_title  

  if( save_plot ){  
      dev.off()
      cat( paste0('\n *********\nNote: plots are saved at: \n', getwd(),'/hormPlotBreak.pdf \n***** \n\n')  )
  }
} # end function


