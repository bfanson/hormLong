#' Plot the area under the curve for peak values 
#' 
#' @param x hormLong object (produced from hormBaseline) [required]
#' @param lower_bound the lower limit to caculate the area under the curve. This can be
#' 'origin', 'baseline', or 'peak' [default = 'origin']
#' @param method the AUC method to use.  Currently,only trapezoid method has been implemented [default = 'trapezoid']
#' @param xscale  determines if x-axis should be free ('free') to change for each panel or
#' remain the same ('fixed') for all panels  [default = 'free']
#' @param yscale  determines if y-axis should be free ('free') to change for each panel or
#' remain the same ('fixed') for all panels  [default = 'free']
#' @param plot_per_page the number of plot panels per page, by row. [default = 4]
#' @param save_plot indicates whether [default = TRUE]
#' @param plot_height  the height of individual plot panels (in inches).  Pdf page height is determined
#' by both plot_per_page and plot_height. [default = 2]
#' @param plot_width  the width of the pdf page. [default = 6]
#' @return nothing is returned.  Pdf file is produced and a summary table 
#' @export
#' @examples
#' 
#' Citation: Cockrem JF and Silverin B. 2002. Variation within and between Birds in Corticosterone 
#' Responses of Great Tits (Parus major).  General and Comparative Endocrinology 125\:197–206.
#' 
#' 
#' result <- hormBaseline(data=hormLynx,  criteria=3, by_var='AnimalID, Hormone', time_var='datetime', conc_var='Conc' )
#' hormArea(result,lower_bound='peak')

hormArea <- function(x, lower_bound = 'origin' , method='trapezoid', xscale='free',yscale='free',
                     plot_per_page=4, plot_height=2, plot_width=6, save_plot=T){
#-- initial checking ---#  
  if( method!='trapezoid'){ 
    stop('no other method currently implemented ')
  }
  if( !(lower_bound %in% c('origin','baseline','peak') ) ){ 
    stop(paste0("lower_bound is incorrect.  It must be 'origin', 'baseline', 'peak': you wrote '", lower_bound,"'") )
  }
  graphics.off() # just to make sure not devices are open

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
    getBaseline <- function(c){
      return( mean( c, na.rm=T ) )
    }
    data <- getSumStat(data=data[data$conc_type=='base',], name='cutoff', func= function(y) getBaseline(y), 
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
                            function(k, date=time_var,conc=conc_var) data.frame(getPeakInfo(k,date,conc),
                                        plot_title=unique(k$plot_title ) ) ) ) ) 
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
      if( length(pk$peak_num)>0 ){
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





