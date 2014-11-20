#' Run iterative process to calculate baseline.
#' 
#' @param data dataset.
#' @param by_var grouping variable.
#' @param time_var grouping variable.
#' @param conc_var grouping variable.
#' 
#' @return hormLong object.
#' @export
#' @examples
#' 
#' result <- hormBaseline(data=hormone, by_var='sp, sex, id', time_var='date', conc_var='conc' )
 

hormACTH <- function(data) {
