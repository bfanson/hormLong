
ds <- hormRead()
ds <- hormDate(data=ds, date_var='Date')

res <- hormBaseline(data=ds, by_var='id,hormone',conc_var='preg.ng.ml',time_var='datetime')
hormPlot(res)

hormSmooth(data=ds, by_var='id',group_var='hormone',conc_var='preg.ng.ml',time_var='datetime',
           smoothness=0.2,points=F)
hormBoxplot(data=ds,conc_var='preg.ng.ml',id_var = 'id')

hormSumTable(res)

hormWrite(res, filename='example_out.csv')

