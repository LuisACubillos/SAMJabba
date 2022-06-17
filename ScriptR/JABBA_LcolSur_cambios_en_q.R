# Installation
# install.packages(devtools)
# devtools:: install_github("jabbamodel/JABBA")
rm(list = ls())
library(JABBA)

#><>><>><>><>><>><>><>><>><>><>><>
# Langostino colorado Sur con cambios en q (CPUE)
#><>><>><>><>><>><>><>><>><>><>><>
df <- read.csv("Data/lac_sur_1968-20212q.csv",sep=";")
head(df)
plot(df$Year,df$CPUE1,ty="b")
points(df$Year,df$CPUE2,pch=19)
plot(df$Year,df$Bv,ty="b")
File = "~/Rwork/SAMJabba/SAM"
assessment = "LCOL_Sur2021_2q"
output.dir = file.path(File,assessment)
dir.create(output.dir,showWarnings = F)
setwd(output.dir)

# Prepara los datos  ----------------------------------------------
cpue = df[,c(1,3,4,5)]
colnames(cpue) = c("Year","cpue1","cpue2","Bv")
#cpue = df[,c(1,4)]
#colnames(cpue) = c("Year","Bv")
se = df[,c(1,7,8,9)]
colnames(se) = c("Year","cpue1","cpue2","Bv")
#se = df[,c(1,7)]
#colnames(se) = c("Year","Bv")

catch = df[,c(1,2)]
lcolsur <- list()
lcolsur$cpue <- cpue
lcolsur$se <- se
lcolsur$catch <- catch


#------------------------------------------------------
# Ajuste simple de JABBA con Biomasa y SE
#-------------------------------------------------------
# Compile JABBA JAGS model and input object

## PELLA com BmsyK=0.4
jblcolsur_M2 = build_jabba(catch=lcolsur$catch,
                               cpue=lcolsur$cpue,
                               se=lcolsur$se,
                               assessment=assessment,
                               scenario = "M2",
                               model.type = "Pella",
                               BmsyK=0.4,
                               sigma.est = FALSE,
                               r.prior=c(0.5,0.5),
                               K.prior = c(70000,0.3),
                               igamma = c(0.01,0.001),
                               Plim=0.2,
                               psi.prior = c(2,0.5))
# Fit JABBA (here mostly default value - careful)
lcolsurM2 = fit_jabba(jblcolsur_M2,save.jabba=TRUE,output.dir=output.dir,quickmcmc = TRUE) # quick run

# Make individual plots
jbplot_catch(lcolsurM2)
jbplot_catcherror(lcolsurM2)
jbplot_ppdist(lcolsurM2)
jbplot_mcmc(lcolsurM2)
jbplot_residuals(lcolsurM2)
jbplot_cpuefits(lcolsurM2)
jbplot_runstest(lcolsurM2)
jbplot_logfits(lcolsurM2)
jbplot_procdev(lcolsurM2)
jbplot_trj(lcolsurM2,type="BBmsy",add=T)
jbplot_trj(lcolsurM2,type="FFmsy",add=T)

# Plot all
jabba_plots(jabba=lcolsurM2,output.dir = output.dir)

#lamsur1$pars_posterior[,1]
lcolsurM2$estimates
lcolsurM2$pars

## PELLA (m estimado libre)
jblcolsur_M3 = build_jabba(catch=lcolsur$catch,
                           cpue=lcolsur$cpue,
                           se=lcolsur$se,
                           assessment=assessment,
                           scenario = "M3",
                           model.type = "Pella_m",
                           sigma.est = FALSE,
                           r.prior=c(0.5,0.5),
                           K.prior = c(70000,0.3),
                           igamma = c(0.01,0.001),
                           Plim=0.2,
                           psi.prior = c(2,0.5))
# Fit JABBA (here mostly default value - careful)
lcolsurM3 = fit_jabba(jblcolsur_M3,save.jabba=TRUE,output.dir=output.dir,quickmcmc = TRUE) # quick run

# Make individual plots
jbplot_catch(lcolsurM3)
jbplot_catcherror(lcolsurM3)
jbplot_ppdist(lcolsurM3)
jbplot_mcmc(lcolsurM3)
jbplot_residuals(lcolsurM3)
jbplot_cpuefits(lcolsurM3)
jbplot_runstest(lcolsurM3)
jbplot_logfits(lcolsurM3)
jbplot_procdev(lcolsurM3)
jbplot_trj(lcolsurM3,type="BBmsy",add=T)
abline(h=0.5,lty=2)
jbplot_trj(lcolsurM3,type="FFmsy",add=T)

# Plot all
jabba_plots(jabba=lcolsurM3,output.dir = output.dir)


#lamsur1$pars_posterior[,1]
round(lcolsurM3$estimates,3)
round(lcolsurM3$pars,5)


# Organize folders by creating a "retro" subfolder
retro.dir = file.path(output.dir,"retro")
dir.create(retro.dir,showWarnings = F)

# Run hindcasts
hclamsurM2 = jabba_hindcast(jblcolsur_M2,save.hc=F,plotall=F,output.dir = retro.dir,peels = 0:5)

# Retro Analysis Summary plot
jbplot_retro(hclamsurM2,as.png = F,single.plots = F,output.dir = retro.dir)
# Zoom-in
mohnsrho = jbplot_retro(hclamsurM2,as.png = F,single.plots = F,output.dir = retro.dir,Xlim=c(2010,2021))

# Save plot and note Mohn's rho statistic
mohnsrho = jbplot_retro(hclamsurM2,as.png = T,single.plots = F,output.dir = retro.dir)
# eval mohnsrho
mohnsrho
mohnsrho["rho.mu",]



Bt <- lamsurM3$timeseries[,,3]
Ft <- lamsurM3$timeseries[,,4]
plot(Bt[,1],Ft[,1],type = "n",las=1,ylim=c(0,3),xlim=c(0,3),ylab="F/Fmsy",xlab="B/Bmsy",axes = FALSE)
polygon(c(0,0,0.5,0.5),c(0,3,3,0),col = "grey90", border = "grey90")
polygon(c(1.5,1.5,3,3),c(0,0.9,0.9,0),col = "lightgreen", border = "lightgreen")
polygon(c(0.5,0.5,0.8,0.8),c(0,1,1,0),col = "yellow", border = "yellow")
polygon(c(0.5,0.5,3,3),c(1,3,3,1),col = "grey90", border = "grey90")
abline(h=1)
lines(c(0.8,0.8),c(0,1))
lines(c(0.5,0.5),c(0,1))
points(Bt[,1],Ft[,1],pch=19,cex=0.3)
lines(Bt[,1],Ft[,1])
lines(c(0.871,2.358),c(0.570,0.570),lwd=4)
lines(c(1.4892878,1.4892878),c(0.312,1.101),lwd=4)
points(1.4892878,0.570,pch=19,cex=1.4)
points(1.4892878,0.570,pch=19,cex=1,col="white")
axis(side=2,las=1)
axis(side=1)




#################### HASTA AQUI 

# Compila JABBA JAGS y objetos de entrada
jbinput5 = build_jabba(catch=lcolsur$catch,cpue=lcolsur$cpue,se=lcolsur$se,assessment=assessment,
                        scenario = "S3final",
                        model.type = "Pella",BmsyK=0.4,
                        r.prior=c(0.3,0.5),
                        psi.prior = c(2,0.2),
                        add.catch.CV = TRUE,
                        proc.dev.all = FALSE, 
                        igamma=c(0.01,0.001),
)

lcolsur5 = fit_jabba(jbinput5,save.jabba=TRUE,output.dir=output.dir)


# Make individual plots
jbplot_catch(lamsur5)
jbplot_catcherror(lamsur5)
jbplot_ppdist(lamsur5)
jbplot_mcmc(lamsur5)
jbplot_residuals(lamsur5)
jbplot_cpuefits(lamsur5)
jbplot_runstest(lamsur5)
jbplot_logfits(lamsur5)
jbplot_procdev(lamsur5)

jbplot_trj(lamsur5,type="BBmsy",add=F)
jbplot_trj(lamsur5,type="FFmsy",add=T)


jbplot_kobe(lcolsur5,add=F)
#lamsur1$pars_posterior[,1]
round(lcolsur5$estimates,3)
round(lcolsur5$pars,3)
