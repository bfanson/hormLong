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
#' result <- hormBaseline(data=hormone, criteria=2, by_var='sp, sex, id', time_var='date', conc_var='conc' )
#' hormArea(result)

hormArea <- function(x, method='trapezoid',plot_per_page=4,plot_height=2, plot_width=6, save_plot=T, 
                     yscale='free',xscale='free'){
  #stop('function under development')
  if( method!='trapezoid'){ 
    stop('no other method currently implemented ')
  }
  by_var_v <- cleanByvar(x$by_var) 
  time_var <- x$time_var
  conc_var <- x$conc_var
  data <- x$data
  data <- data[!is.na(data[,conc_var]),]
  data <-ridFactor(data)
  data$plot_title <- getPlotTitle(data,by_var=by_var_v)
  data <- getSumStat(data=data,name='cutoff', func= function(y) getCutoff(y, criteria=x$criteria ), add_ds=data, by_var=by_var_v, c_var=conc_var )

#--- peak analysis ---#
  #-- shift curve by cutoff ---#
    # k <- data[data$plot_title==unique(data$plot_title)[1],]
    ds_pk <- do.call("rbind", as.list(by(data, data[,by_var_v], 
                            function(k, date=time_var,conc=conc_var) data.frame(getPeakInfo(k,date,conc),plot_title=unique(k$plot_title ) ) ) ) ) 
    ds_pk <- ds_pk[ds_pk$peak_num>0,]  # get rid individuals with no peaks       

#-- plot out results ---#
  
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
          text( mean(pk1$t), max(c_order), labels=p, adj = c(0,-0.5))
        }
      }
  }
  if( save_plot ){  
      dev.off()
      cat( paste0('\n *********\nNote: plots are saved at: \n', getwd(),'/hormArea.pdf \n***** \n\n')  )
  }


  #--- produce table ---#
    ds_tbl <- unique( ds_pk[,c('plot_title','peak_num','AUC')])    
    ds_tbl <- merge(ds_tbl, unique(data[,c(by_var_v,'plot_title')]), all.x=T)
    ds_tbl$AUC <- round(ds_tbl$AUC,3)
    ds_out <- ds_tbl[,c(by_var_v,'peak_num','AUC')]
    names(ds_out)[which(names(ds_out)=='AUC')] <-'peak AUC'    

    #-- merge in total AUC ---#
#     do.call("rbind", as.list(by(data, data[,by_var_v], 
#                   function(k) data.frame(totalAUC=getAUC(k[,time_var],k[,conc_var]),plot_title=unique(k$plot_title ) ) )))

    #--- write dataset ---#
      write.csv(ds_out,file='hormAUCtable.csv',quote=F, row.names=F)
      cat( paste0('\n *********\nNote: table saved at: \n', getwd(),'/hormAUCtable.csv \n***** \n\n')  )

}





