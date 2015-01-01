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
#' result <- hormBaseline(data=hormLynx, criteria=2, by_var='AnimalID, Hormone', time_var='Date', conc_var='Conc' )
 

hormSmooth <- function(data,by_var,group_var,time_var,conc_var, plot_per_page=4,
                       plot_height=2, plot_width=6, yscale='free', xscale='free',
                       smoothness=0.7, shape=1, colour='black', line_type=1,
                       line_width=1, add_points=T,...) {
stop('function under development')
  graphics.off() # just to make sure not devices are open

#--- add in checks ---#
   if(missing(data)){
      stop('data must be specified')
    }
  if(!is.data.frame(data)){
      stop('data must be a data.frame object')
    }
  if(missing(by_var)){
    message('Warning: No by_var included ... baseline value is for whole dataset')
    data$hold_id <- 1    # set-up a new by_var
    by_var <- 'hold_id'  
    }
  if(!is.character(by_var)){
    stop('by_var must be a character string')
    }
  if(missing(group_var)){
    message('Warning: No group_var included ... only one line per plot')
    data$hold_grp <- 1    # set-up a new by_var
    group_var <- 'hold_grp'  
    }
  if(!is.character(group_var)){
    stop('group_var must be a character string')
    }

  by_var_v <- cleanByvar(by_var) # make by_var a vector
  if( sum(!( c(by_var_v,conc_var,time_var)  %in% names(data) ))>0 ){
      stop('not all variables are present in data.set.  check your column names')
    }
  if(missing(conc_var)){
      stop('conc_var (e.g. the response variable containing the concentration) must be specified')
    }
  if(!is.numeric(data[,conc_var])){
    stop('conc_var must be numeric')
    }
  if(missing(time_var)){
      stop('time_var must be specified')
    }
  if( !( class(data[,time_var]) %in% c('Date','numeric') ) ){
    time_class <- class(class(data[,time_var]))
    stop(paste0('time_var must be numeric or Date (your time_var is "',time_class,'" variable') )
    }

  
#-- set-up ---#
  time_var <- time_var
  conc_var <- conc_var
  data <- data[ do.call(order, data[c(by_var_v,group_var,time_var)]), ]
  
  #-- create titles--#
    paste1 <- function(...) paste(...,sep='; ')
    if( length(by_var_v)>1){
        data$plot_title <- do.call(paste1, data[,c(by_var_v)] )
     }else{data$plot_title <- data[,c(by_var_v)] }
  #-- create page numbers ---#
    hold <- unique( data[,by_var_v])
    hold$pg_num <- ((1:nrow(hold))-1)%/%plot_per_page +1
    data<-merge(data,hold, all.x=T)

#--- create plots ---#
  if( save_plot ){
    pdf('hormSmooth.pdf', height=plot_per_page * plot_height, width = plot_width )
    cat( paste0('\n *********\nNote: plots are saved at: \n', getwd(),'/hormSmooth.pdf \n***** \n\n')  )
  }
  for( pg in 1:max(pg_num)){
    ds_sub <- data[data$pg_num==pg,]
    #--- main plot ---#
    for(i in unique(ds_sub$plot_title)){
      ds_sub1 <- ds_sub[ds_sub$plot_title==i,]
      fig1 <- ggplot(data=ds_sub1,aes_string(x=time_var,y=conc_var, group=group_var)) 
          if(add_points) fig1 <- fig1 + geom_point() 
      fig1 <- fig1 + geom_smooth(span=smoothness,se=F) + 
                labs(title=i)
    }
  
    }
  if( save_plot ){  dev.off() }
}



#   par(mfrow=c(plot_per_page,1), mar=c(2,4,2,0.5),oma=c(2,2,2,0))
#   for( i in unique(data$plot_title)){
#     ds_sub <- data[data$plot_title==i, ]
#     events <- ds_sub[ !is.na(ds_sub$event) & ds_sub$event!='',c('event',time_var)]
#     ds_sub <- ds_sub[!is.na(ds_sub[,conc_var]),]
#     
#     #--- set up scales
#       if( yscale=='free'){
#         ymin <- min(ds_sub[,conc_var])*0.95
#         ymax <- max(max(ds_sub[,conc_var]) )*1.2
#       }else{
#         ymin <- min(data[,conc_var],na.rm=T)*0.95
#         ymax <- max(data[,conc_var],na.rm=T)*1.2
#       }
#       if( xscale=='free'){
#         xmin <- min( ds_sub[,time_var],na.rm=T )
#         xmax <- max( ds_sub[,time_var],na.rm=T )
#       }else{
#         xmin <- min( data[,time_var],na.rm=T)
#         xmax <- max( data[,time_var],na.rm=T)
#       }    
#     
#     #--- main plot ---#
#     plot(ds_sub[,conc_var] ~ ds_sub[,time_var], type='n',xlim=c(xmin,xmax), ylim=c(ymin,ymax), 
#           xlab=NA, ylab=conc_var)
#       mtext(unique(ds_sub$plot_title),side=3,line=0.25)
#       loop <- 1
#       for(g in unique(ds_sub[,group_var]) ){
#         ds_sub1 <- ds_sub[ds_sub[,group_var]==g,]
#         if(points)points(ds_sub1[,time_var], ds_sub1[,conc_var],pch=shape[loop],col=colour[loop] )
#         lines(loess.smooth(ds_sub1[,time_var],ds_sub1[,conc_var], col=colour[loop], span=smoothness[loop]))
#         #loop <- loop + 1
#       }
#       legend('topleft',legend=levels(ds_sub[,group_var]),col=colour,pch=shape, lty=line_type,bty='n')

  
