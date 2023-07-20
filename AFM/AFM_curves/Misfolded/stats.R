setwd('C:/Users/jhmva/Documents/AFM_selected_curves/Short')

d.AFM.04 <- read.table('stats_test04.txt', header=TRUE, sep='\t', quote='"')
h.LC.04 <- hist(d.AFM.04$LC, breaks=50)
plot(d.AFM.04$LC, d.AFM.04$Break_force)

d.AFM.07 <- read.table('stats_test07.txt', header=TRUE, sep='\t', quote='"')
h.LC.07 <- hist(d.AFM.07$LC, breaks=50)

setwd('C:/Users/jhmva/Documents/AFM_selected_curves/Long')

d.AFM.04 <- read.table('stats_test04.txt', header=TRUE, sep='\t', quote='"')
h.LC.04 <- hist(d.AFM.04$LC, breaks=80)
plot(d.AFM.04$LC, d.AFM.04$Break_force)

d.AFM.04 <- read.table('last_events04.txt', header=TRUE, sep='\t', quote='"')
h.LC.04 <- hist(d.AFM.04$LC, breaks=80)
plot(d.AFM.04$LC, d.AFM.04$Break_force)
plot(d.AFM.04$LC, d.AFM.04$Fitting_error)

# Read data
setwd('C:/Users/jhmva/Documents/AFM_selected_curves')
d.long <- read.table('Long/stats_final2.txt', header=TRUE, sep='\t', quote="")
d.mis <- read.table('Misfolded/stats_final2.txt', header=TRUE, sep='\t', quote="")

#plot(d.long$Break_extension, d.long$Fitting_error, type='p', pch=16, xlab='Extension length (nm)', ylab='Standard deviation (nm)', ylim=c(0, max(c(d.long$Fitting_error, d.mis$Fitting_error))))
#points(d.mis$Break_extension, d.mis$Fitting_error, col='red', pch=16)

#plot(d.long$Break_extension, d.long$Fitting_error, type='p')
#points(d.mis$Break_extension, d.mis$Fitting_error, col='red')
#boxplot(d.long$Fitting_error, d.mis$Fitting_error, xlab='State', ylab='Standard deviation (nm)', main='Fitting error of the WLC fits of the\nfully extended native and misfolded states')
#axis(1, labels=c('Native', 'Misfolded'), at=c(1,2), las=1)
#t.test(d.long$Fitting_error, d.mis$Fitting_error)

#plot(d.long$Break_extension, type='p')
#points(d.mis$Break_extension, col='red')
#boxplot(d.long$Fitting_error, d.mis$Fitting_error)
# Student's T-test for differences in values
t.test(d.long$Fitting_error, d.mis$Fitting_error)

# Plot distributions and density plots
hist(d.long$Fitting_error, col=rgb(1,0,0,0.5), freq=FALSE, main=NULL, xlab='Standard Deviation (pN)', xlim=c(0,400)) #, main='Fitting Errors of Proteins in Native and Misfolded States'
hist(d.mis$Fitting_error, col=rgb(0,0,1,0.5), add=TRUE, freq=FALSE, breaks=200)
lines(density(d.long$Fitting_error), lwd=2)
lines(density(d.mis$Fitting_error, adjust=0.3), lwd=2)
legend("topright", c("Native", "Misfolded"), fill=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)))

#########
# Read data for Mean errors and number of events per curve
setwd('C:/Users/jhmva/Documents/AFM_selected_curves')
d.long.mean <- read.table('Long/mean_errs.txt', header=TRUE, sep='\t', quote="")
d.mis.mean <- read.table('Misfolded/mean_errs.txt', header=TRUE, sep='\t', quote="")

# Student's T-test for differences in values
t.test(d.long.mean$Mean_error, d.mis.mean$Mean_error)

# Plot distributions and density plots
hist(d.long.mean$Mean_error, col=rgb(1,0,0,0.5), freq=FALSE, main='Fitting Errors of Proteins in Native and Misfolded States', xlab='Standard Deviation (pN)', xlim=c(0,200))
hist(d.mis.mean$Mean_error, col=rgb(0,0,1,0.5), add=TRUE, freq=FALSE, breaks=100)
lines(density(d.long.mean$Mean_error), lwd=2)
lines(density(d.mis.mean$Mean_error), lwd=2)
legend("topright", c("Native", "Misfolded"), fill=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)))

hist(d.long.mean$N_events, col=rgb(1,0,0,0.5), freq=FALSE, main='Number of Events per Curve', xlab='Number of Events', xlim=c(0,10))
hist(d.mis.mean$N_events, col=rgb(0,0,1,0.5), add=TRUE, freq=FALSE)
