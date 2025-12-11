library(shellpipes)
manageConflicts()

library(ggplot2); theme_set(theme_bw())

okabe_ito <- c("#E69F00", "#56B4E9", "#009E73") 

loadEnvironments()
Rf <- rdsRead()

cplot <- (ggplot(Rf)
	+ aes(t, value, color=name)
	+ geom_line()
	+ coord_cartesian(xlim=c(0, 50))
	+ scale_color_manual(values = okabe_ito)
)
print(cplot)
print(cplot + scale_y_log10())

