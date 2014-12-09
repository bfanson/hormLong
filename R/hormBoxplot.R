#' Boxplot of individual concentrations 
#' 
#' @param data dataset [required]
#' @param id_var name of the individual variable [required]
#' @param conc_car name of the concentration variable [required]
#' @param log_scale determines if y-axis is log10-scale or not. log-scale='y' makes log scale [default='n']  
#' @return nothing...Produces a pdf file with the graph in current working directory
#' @export
#' @examples
#' 
#' hormBoxplot(data=hormone, conc_var='conc',id_var='id')
#' 

hormBoxplot <- function(data, conc_var, id_var, plot_height=4, plot_width=6, save_plot=T, log_scale='n'){
  
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

  if( save_plot ){
    pdf('hormBoxplot.pdf', height=plot_height, width = plot_width )
    cat( paste0('\n *********\nNote: plots are saved at: \n', getwd(),'/hormBoxplot.pdf \n***** \n\n')  )
  }
  
  
  par(mar=c(5,4,1,1))
  data <- data[ order(data[,id_var]),]
  if(log_scale=='y'){
    h<-boxplot( data[,conc_var] ~ data[,id_var],ylab=conc_var, xaxt='n',log = 'y'  )
  }else{
    h<-boxplot( data[,conc_var] ~ data[,id_var],ylab=conc_var, xaxt='n'  )
  }
    points(as.numeric(as.factor(data[,id_var])), data[,conc_var] )
    axis(1,at=1:length(h$names),labels=NA)
    text(1:length(h$names), par("usr")[3] - (par("usr")[4]-par("usr")[3])*0.1 , srt = 315, adj = 0,  
        labels = h$names, xpd = TRUE)
  if( save_plot ){  dev.off() }
}

