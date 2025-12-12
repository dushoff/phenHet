library(shellpipes)
manageConflicts()

library(ggplot2); theme_set(theme_bw())

pal <- okabe_ito <- c("#E69F00", "#56B4E9", "#009E73") 

loadEnvironments()
Rf <- rdsRead()

iplot <- (ggplot(Rf)
	+ aes(t, obs, color=name)
	+ geom_line()
	+ coord_cartesian(xlim=c(0, 50))
	+ scale_color_manual(values = pal)
)
print(iplot)

cplot <- (ggplot(Rf)
	+ aes(t, value, color=name)
	+ geom_line()
	+ geom_point(aes(size=obs))
	+ coord_cartesian(xlim=c(0, 50))
	+ scale_color_manual(values = pal)
	+ scale_size_area()
	+ geom_hline(yintercept=Ri0, color=pal[[3]])
	+ geom_hline(yintercept=Rc0, color=pal[[2]])
)
print(cplot)
print(cplot + scale_y_log10())

