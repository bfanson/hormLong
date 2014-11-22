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
#' @examples
#' 

hormPlotBreaks <- function(){

x<-c(1:10,101:110,200:240)
brk<-c(rep(1,10),rep(2,10),rep(3,41))
y<-rnorm(length(x))
breaks<-3
ds<-data.frame(x=x,brk=brk,y=y)

par(mfrow=c(1,3),oma=c(3,4,1,1))
for(i in 1:breaks){
  par(mar=c(3,0.2,1,0.2))
  with(ds[ds$brk==i,], plot(y~x, ylim=c(min(ds$y),max(ds$y)), xlim=c(min(x)*0.96,max(x)*1.04),
                            bty='n',yaxt='n',type='l',col=i,xaxs='i' ) )
  if(i==1) axis(2)
  box()
  if( getRidLines==T){
    if(i!=breaks) abline(v = max(ds[ds$brk==i,'x'])*1.04, col='white')
    if(i!=1)abline(v = min(ds[ds$brk==i,'x'])*0.96, col='white')
  }
  axis(1,labels=F)
  #axis.break(1,max(ds[ds$brk==i,'x'])*1.02,style="gap",) 

}

# option 2 
ds<-data.frame(x=x,brk=brk,y=y)
buffer=8
if( breaks>1 ){
  min_brk <- aggregate(ds['x'],ds['brk'],min)
    names(min_brk)[2]<-'x_min'
  rg_brk <- aggregate(ds['x'],ds['brk'],function(x) diff( range(x))) 
    names(rg_brk)[2]<-'x_rg'
  rg_brk$x_acc <- cumsum(c(0,rg_brk$x_rg[1:(nrow(rg_brk)-1)]) ) 
  adj_brk <- merge(min_brk,rg_brk)
  adj_brk$x_diff <- adj_brk$x_min-(adj_brk$x_acc + (adj_brk$brk-1)*buffer)
  ds <- merge(ds,adj_brk[,c('brk','x_diff')],by='brk', all.x=T )
  ds$x_adj <- ds$x - ds$x_diff
}

#par(mfrow=c(1,3),oma=c(3,4,1,1))
with(ds, plot(y~x_adj, type='n', xaxt='n') )
for(i in 1:breaks){
  ds1 <- ds[ds$brk==i,] 
  points(ds1$x_adj,ds1$y,col=i)
  lines(ds1$x_adj,ds1$y,col=i)
#  text(ds1$x_adj,-2+rnorm(length(ds1$x_adj)),labels=ds1$x,cex=0.6)
  p_at <-pretty(c(min(ds1$x_adj),max(ds1$x_adj)))
  axis(1,at=p_at, labels = p_at+ds1$x_diff[1] )
}
}