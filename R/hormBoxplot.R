#' Boxplot of individual concentrations 
#' 
#' @param data dataset containing the hormone data. [required]
#' @param conc_var name of the concentration variable [required]
#' @param id_var name of the individual variable [required]
#' @param by_var names of the grouping variables.  Column names should be separated by a comma:
#'  e.g. 'name1' for one column, 'name1, name2' for two columns, 'name1, name2, name3' for three columns,
#'  etc.  Remember, capitalization matters. [optional]
#' @param log_scale determines if y-axis is log10-scale or not. log-scale='y' makes 
#' log scale [default='n']  
#' @param plot_per_page the number of plot panels per page, by row. [default = 1]
#' @param plot_height  the height of individual plot panels (in inches).  Pdf page height is 
#' determined by both plot_per_page and plot_height. [default = 5]
#' @param plot_width  the width of the pdf page. [default = 6]
#' @return nothing Produces a pdf file with the graph in current working directory
#' 
#' @export
#' @examples
#' 
#' hormBoxplot(data=hormLynx, conc_var='Conc', id_var='AnimalID', by_var='Hormone')
#' 
#'# the following will fail because there are zero concentrations.
#' hormBoxplot(data=hormLynx, conc_var='Conc', id_var='AnimalID', by_var='Hormone', log_scale='y')
#'
#'# One option to deal with zeros 
#' ds <- hormLynx            # make a new object
#' ds$Conc1 <- ds$Conc + 1   # add 1 to all values concenration and save as a new column
#' hormBoxplot(data=ds, conc_var='Conc1', id_var='AnimalID', by_var='Hormone', log_scale='y')
#' 

hormBoxplot <- function(data, conc_var, id_var, by_var, log_scale='n', 
                        plot_per_page=1, plot_height=5, plot_width=6, save_plot=T){
  
  graphics.off() # just to make sure no devices are open
  if(missing(conc_var)){
      stop('conc_var (e.g. the response variable containing the concentration) must be specified')
    }
  if(!is.numeric(data[,conc_var])){
    cat('WARNING: conc_var should be numeric. Program will convert automatically.  Non-numeric
         values will be become missing.  Any non_numeric values are listed below:\n\n')
    print('------')
    if(any(is.na(as.numeric( data[,conc_var]))==T)){
      print( data[ is.na(as.numeric( data[,conc_var])),] )
    }
    cat('\n')
    data[,conc_var] <- as.numeric(data[,conc_var])
    }
  
  if(missing(id_var)){
      stop('id_var (e.g. animal name) must be specified')
    }
  if(missing(by_var)){
      by_var <- 'by'
      data$by <- ' '   # create a by_variable so following code works
    }

  if( save_plot ){
    pdf('hormBoxplot.pdf', height=plot_height, width = plot_width )
    cat( paste0('\n *********\nNote: plots are saved at: \n', getwd(),'/hormBoxplot.pdf \n***** \n\n')  )
  }
  
  by_var_v <- cleanByvar(by_var ) 
  data$plot_title <- getPlotTitle(data, by_var=by_var_v) 
  data <- data[ !is.na(data[conc_var]),]
  
  #--- check for zeros ---#
  if( log_scale=='y' & nrow( data[data[,conc_var] == 0, ]) ){
    stop('There are zero concentrations. log of 0 is undefined. You may want to add +1 to all concentrations and rerun program')
  }
  
  
  par(mfrow=c(plot_per_page,1), mar=c(4,4,2,0.5),oma=c(2,2,2,2))
  for(p in unique( data$plot_title ) ){
    data1 <- subset(data, plot_title == p )
    data1 <- data1[ order(data1[,id_var]),]
    if(log_scale=='y'){
      h<-boxplot( data1[,conc_var] ~ data1[,id_var],ylab=conc_var, xaxt='n',log = 'y'  )
    }else{
      h<-boxplot( data1[,conc_var] ~ data1[,id_var],ylab=conc_var, xaxt='n'  )
    }
      mtext(p, side = 3, cex=1.1)
      points(as.numeric(as.factor(data1[,id_var])), data1[,conc_var] )
      axis(1,at=1:length(h$names),labels=NA)
      if( log_scale == 'y'){
      text(1:length(h$names), 10^par("usr")[3] , srt = 315, adj = c(0,2),  
          labels = h$names, xpd = TRUE)
      }else{
      text(1:length(h$names), par("usr")[3] - (par("usr")[4]-par("usr")[3])*0.1 , srt = 315, adj = 0,  
          labels = h$names, xpd = TRUE)
      }
    
  }  
  if( save_plot ){  dev.off() }
}

