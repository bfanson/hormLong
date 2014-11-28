#' Automated function that runs baseline calculation, make plots, outputs dataset as csv
#' 
#' @param none  
#' @return nothing in R.  Output hormPlot.pdf and horm_data_out.csv
#' @export
#' @examples
#' 
#' ds <- hormAuto(by_var='sp, sex, id', time_var='date', conc_var='conc')

hormAuto<- function(by_var, time_var, conc_var, filename='test' ){
  stop('function under development')
  ds  <- hormRead()
  ds <- hormone
  res <- hormBaseline(data=ds, by_var=by_var,time_var=time_var,conc_var=conc_var)
  hormPlot(res)
  hormWrite(res,'horm_data_out.csv')

}
