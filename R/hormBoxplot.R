#' Boxplot of individual concentrations 
#' 
#' @param data dataset 
#' @param id_var name of the individual variable
#' @return NA
#' @export
#' @examples
#' 
#' hormBoxplot(data=hormone, conc_var='conc',id_var='id')
#' 

hormBoxplot <- function(data, conc_var, id_var, plot_height=4, plot_width=6, save_plot=T){
  
  if(missing(conc_var)){
      stop('conc_var (e.g. the response variable containing the concentration) must be specified')
    }
  if(!is.numeric(data[,conc_var])){
    cat('WARNING: conc_var should be numeric. Program will convert automatically.  Non-numeric
         values will be become missing.  Any non_numeric values are listed below:\n\n')
    print('------')
    print( data[ is.na(as.numeric( data[,conc_var])),] )
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
  
  
  
  par(mar=c(3,4,1,1))
  boxplot( data[,conc_var] ~ data[,id_var],ylab=conc_var )  
    points(as.numeric(as.factor(data[,id_var])), data[,conc_var] )
  if( save_plot ){  dev.off() }
}
