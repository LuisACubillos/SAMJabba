# Installation
# install.packages(devtools)
# devtools:: install_github("jabbamodel/JABBA")
library(JABBA)

#><>><>><>><>><>><>><>><>><>><>><>
# Langostino amarillo Sur
#><>><>><>><>><>><>><>><>><>><>><>
df <- read.csv("Data/lam_sur_1979-2021.csv",sep=";")
head(df)
plot(df$Year,df$CPUE,ty="b")
plot(df$Year,df$Bv,ty="b")
File = "~/Rwork/JabbaSA/SAM"
assessment = "LAM_Sur2021"
output.dir = file.path(File,assessment)
dir.create(output.dir,showWarnings = F)
setwd(output.dir)

# Prepara los datos  ----------------------------------------------
cpue = df[,c(1,3,4)]
colnames(cpue) = c("Year","cpue","Bv")
se = df[,c(1,6,7)]
colnames(se) = c("Year","cpue","Bv")
catch = df[,c(1,2)]
lamsur <- list()
lamsur$cpue <- cpue
lamsur$se <- se
lamsur$catch <- catch

#############################################
# PELLA-TOMLINSON
# Modelo 1
#
# Compile JABBA JAGS model and input object
jbinputM1 = build_jabba(catch=lamsur$catch,cpue=lamsur$cpue,se=lamsur$se,assessment=assessment,scenario = "M1",model.type = "Pella",BmsyK=0.4,sigma.est = FALSE,igamma = c(0.001,0.001))
# Fit JABBA (here mostly default value - careful)
lamsurM1 = fit_jabba(jbinputM1,save.jabba=TRUE,output.dir=output.dir,quickmcmc = TRUE) # quick run

jbplot_catch(lamsurM1)
jbplot_catcherror(lamsurM1)
jbplot_ppdist(lamsurM1)
jbplot_mcmc(lamsurM1)
jbplot_residuals(lamsurM1)
jbplot_cpuefits(lamsurM1)
jbplot_runstest(lamsurM1)
jbplot_logfits(lamsurM1)
jbplot_procdev(lamsurM1)
jbplot_trj(lamsurM1,type="BBmsy",add=T)
jbplot_trj(lamsurM1,type="FFmsy",add=T)
lamsurM1$pars_posterior[,1]
round(lamsurM1$estimates,3)
round(lamsurM1$pars,3)
jbplot_kobe(lamsurM1,add=F)

# Compile JABBA JAGS model and input object
jbinputM2 = build_jabba(catch=lamsur$catch,cpue=lamsur$cpue,se=lamsur$se,assessment=assessment,scenario = "M2",model.type = "Pella_m",sigma.est = FALSE,igamma = c(0.001,0.001))
# Fit JABBA (here mostly default value - careful)
lamsurM2 = fit_jabba(jbinputM2,save.jabba=TRUE,output.dir=output.dir,quickmcmc = TRUE) # quick run

jbplot_catch(lamsurM2)
jbplot_catcherror(lamsurM2)
jbplot_ppdist(lamsurM2)
jbplot_mcmc(lamsurM2)
jbplot_residuals(lamsurM2)
jbplot_cpuefits(lamsurM2)
jbplot_runstest(lamsurM2)
jbplot_logfits(lamsurM2)
jbplot_procdev(lamsurM2)
jbplot_trj(lamsurM2,type="BBmsy",add=T)
jbplot_trj(lamsurM2,type="FFmsy",add=T)
lamsurM2$pars_posterior[,1]
round(lamsurM2$estimates,3)
round(lamsurM2$pars,3)

# Compila JABBA JAGS y objetos de entrada
jbinputM3 = build_jabba(catch=lamsur$catch,cpue=lamsur$cpue,se=lamsur$se,assessment=assessment,
                      scenario = "M3",
                      model.type = "Pella",
                      BmsyK=0.4,
                      r.prior=c(0.2,0.5),
                      psi.prior = c(0.95,0.1),
                      add.catch.CV = TRUE,
                      proc.dev.all = FALSE, 
                      igamma=c(0.001,0.001),
                      )

lamsurM3 = fit_jabba(jbinputM3,save.jabba=TRUE,output.dir=output.dir)


# Make individual plots
jbplot_catch(lamsurM3)
jbplot_catcherror(lamsurM3)
jbplot_ppdist(lamsurM3)
jbplot_mcmc(lamsurM3)
jbplot_residuals(lamsurM3)
jbplot_cpuefits(lamsurM3)
jbplot_runstest(lamsurM3)
jbplot_logfits(lamsurM3)
jbplot_procdev(lamsurM3)
jbplot_trj(lamsurM3,type="BBmsy",add=T)
jbplot_trj(lamsurM3,type="FFmsy",add=T)
round(lamsurM3$estimates,3)
round(lamsurM3$pars,3)

# Plot all
jabba_plots(jabba=lamsurM3,output.dir = output.dir)

# Organize folders by creating a "retro" subfolder
retro.dir = file.path(output.dir,"retro")
dir.create(retro.dir,showWarnings = F)

# Run hindcasts
hclamsurM3 = jabba_hindcast(jbinputM3,save.hc=F,plotall=F,output.dir = retro.dir,peels = 0:7)

# Retro Analysis Summary plot
jbplot_retro(hclamsurM3,as.png = F,single.plots = F,output.dir = retro.dir)
# Zoom-in
mohnsrho = jbplot_retro(hclamsurM3,as.png = F,single.plots = F,output.dir = retro.dir,Xlim=c(2010,2021))

# Save plot and note Mohn's rho statistic
mohnsrho = jbplot_retro(hclamsurM3,as.png = T,single.plots = F,output.dir = retro.dir)
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



#------------------------------------------------------
# Ajuste simple de JABBA  Solamente Catch
#-------------------------------------------------------
#jbinput = build_jabba(catch=lamsur$catch,model.type = "Fox",assessment=assessment,scenario =  "CatchOnly" ,b.prior=c(0.7,0.2,2009,"bbmsy"),psi.prior = c(1,0.1))
# Ajuste JABBA
#lamsur = fit_jabba(jbinput,save.jabba=TRUE,output.dir=output.dir)
# Make individual plots
#jbplot_catch(lamsur1)
#jbplot_catcherror(lamsur1)
#jbplot_trj(lamsur1,type="BBmsy",add=T)
#jbplot_trj(lamsur1,type="FFmsy",add=T)

#------------------------------------------------------
# Ajuste simple de JABBA con Biomasa y SE
#-------------------------------------------------------
# Compile JABBA JAGS model and input object
jbinput0 = build_jabba(catch=lamsur$catch,cpue=lamsur$cpue,se=lamsur$se,assessment=assessment,scenario = "S1",model.type = "Fox",sigma.est = FALSE,igamma = c(0.001,0.001))
# Fit JABBA (here mostly default value - careful)
lamsur0 = fit_jabba(jbinput1,save.jabba=TRUE,output.dir=output.dir,quickmcmc = TRUE) # quick run

# Make individual plots
jbplot_catch(lamsur0)
jbplot_catcherror(lamsur0)
jbplot_ppdist(lamsur0)
jbplot_mcmc(lamsur0)
jbplot_residuals(lamsur0)
jbplot_cpuefits(lamsur0)
jbplot_runstest(lamsur0)
jbplot_logfits(lamsur0)
jbplot_procdev(lamsur0)
jbplot_trj(lamsur0,type="BBmsy",add=T)
jbplot_trj(lamsur0,type="FFmsy",add=T)
lamsur0$pars_posterior[,1]
round(lamsur0$estimates,3)
round(lamsur0$pars,3)

# SCHAEFER
# Compile JABBA JAGS model and input object
jbinput_schaefer = build_jabba(catch=lamsur$catch,cpue=lamsur$cpue,se=lamsur$se,assessment=assessment,scenario = "S2",model.type = "Schaefer",sigma.est = FALSE,igamma = c(0.001,0.001))
# Fit JABBA (here mostly default value - careful)
lamsur1 = fit_jabba(jbinput_schaefer,save.jabba=TRUE,output.dir=output.dir,quickmcmc = TRUE) # quick run

jbplot_catch(lamsur1)
jbplot_catcherror(lamsur1)
jbplot_ppdist(lamsur1)
jbplot_mcmc(lamsur1)
jbplot_residuals(lamsur1)
jbplot_cpuefits(lamsur1)
jbplot_runstest(lamsur1)
jbplot_logfits(lamsur1)
jbplot_procdev(lamsur1)
jbplot_trj(lamsur1,type="BBmsy",add=T)
jbplot_trj(lamsur1,type="FFmsy",add=T)
lamsur1$pars_posterior[,1]
round(lamsur1$estimates,3)
round(lamsur1$pars,3)
