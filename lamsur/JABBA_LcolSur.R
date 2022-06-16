# Installation
# install.packages(devtools)
# devtools:: install_github("jabbamodel/JABBA")
library(JABBA)

#><>><>><>><>><>><>><>><>><>><>><>
# Langostino colorado Sur
#><>><>><>><>><>><>><>><>><>><>><>
df <- read.csv("Data/lac_sur_1968-2021.csv",sep=";")
head(df)
plot(df$Year,df$CPUE,ty="b")
plot(df$Year,df$Bv,ty="b")
File = "~/Rwork/SAMJabba/SAM"
assessment = "LCOL_Sur2021"
output.dir = file.path(File,assessment)
dir.create(output.dir,showWarnings = F)
setwd(output.dir)

# Prepara los datos  ----------------------------------------------
cpue = df[,c(1,3,4)]
colnames(cpue) = c("Year","cpue","Bv")
#cpue = df[,c(1,4)]
#colnames(cpue) = c("Year","Bv")
se = df[,c(1,6,7)]
colnames(se) = c("Year","cpue","Bv")
#se = df[,c(1,7)]
#colnames(se) = c("Year","Bv")

catch = df[,c(1,2)]
lcolsur <- list()
lcolsur$cpue <- cpue
lcolsur$se <- se
lcolsur$catch <- catch

## PELLA Bmsy/k=0.4
jblcolsur_M1 = build_jabba(catch=lcolsur$catch,
                                cpue=lcolsur$cpue,
                                se=lcolsur$se,
                                assessment=assessment,
                                scenario = "M3",
                                model.type = "Pella",
                                BmsyK=0.4,
                                sigma.est = FALSE,
                                r.prior=c(0.5,0.5),
                                K.prior = c(70000,0.5),
                                igamma = c(0.01,0.001),
                                Plim=0.2,
                                psi.prior = c(2,0.1))
# Fit JABBA (here mostly default value - careful)
lcolsurM1 = fit_jabba(jblcolsur_M1,save.jabba=TRUE,output.dir=output.dir,quickmcmc = TRUE) # quick run

# Make individual plots
jbplot_catch(lcolsurM1)
jbplot_catcherror(lcolsurM1)
jbplot_ppdist(lcolsurM1)
jbplot_mcmc(lcolsurM1)
jbplot_residuals(lcolsurM1)
jbplot_cpuefits(lcolsurM1)
jbplot_runstest(lcolsurM1)
jbplot_logfits(lcolsurM1)
jbplot_procdev(lcolsurM1)

# Plot all
jabba_plots(jabba=lcolsurM1,output.dir = output.dir)

#lamsur1$pars_posterior[,1]
lcolsurM1$estimates
lcolsurM1$pars

jbplot_kobe(lcolsurM3,add=F)

## PELLA (m estimado libre)
jblcolsur_M2 = build_jabba(catch=lcolsur$catch,cpue=lcolsur$cpue,se=lcolsur$se,
                               assessment=assessment,
                               scenario = "M2",
                               model.type = "Pella_m",
                               sigma.est = FALSE,
                               r.prior=c(0.5,0.5),
                               K.prior = c(70000,0.3),
                               igamma = c(0.001,0.001),
                               Plim=0.2,
                               psi.prior = c(2,0.1))


# Fit JABBA (here mostly default value - careful)
lcolsurM2 = fit_jabba(jblcolsur_M2,save.jabba=TRUE,output.dir=output.dir,quickmcmc = TRUE) # quick run

head(lcolsur4$kobe)

# Make individual plots
jbplot_catch(lcolsur4)
jbplot_catcherror(lcolsur4)
jbplot_ppdist(lcolsur4)
jbplot_mcmc(lcolsur4)
jbplot_residuals(lcolsur4)
jbplot_cpuefits(lcolsur4)
jbplot_runstest(lcolsur4)
jbplot_logfits(lcolsur4)
jbplot_procdev(lcolsur4)

jbplot_trj(lcolsur4,type="BBmsy",add=T)
jbplot_trj(lcolsur4,type="FFmsy",add=T)

#lamsur1$pars_posterior[,1]
lcolsur4$estimates
lcolsur4$pars




#------------------------------------------------------
# Ajuste simple de JABBA con Biomasa y SE
#-------------------------------------------------------
# Compile JABBA JAGS model and input object
jblcolsur_input1 = build_jabba(catch=lcolsur$catch,cpue=lcolsur$cpue,se=lcolsur$se,assessment=assessment,scenario = "S1",
                               model.type = "Fox", sigma.est = FALSE, r.prior=c(0.5,0.5),K.prior = c(70000,0.3),
                               igamma = c(0.001,0.001),Plim=0.2,psi.prior = c(2,0.1))
# Fit JABBA (here mostly default value - careful)
lcolsur1 = fit_jabba(jblcolsur_input1,save.jabba=TRUE,output.dir=output.dir,quickmcmc = FALSE) # quick run

head(lcolsur1$kobe)

# Make individual plots
jbplot_catch(lcolsur1)
jbplot_catcherror(lcolsur1)
jbplot_ppdist(lcolsur1)
jbplot_mcmc(lcolsur1)
jbplot_residuals(lcolsur1)
jbplot_cpuefits(lcolsur1)
jbplot_runstest(lcolsur1)
jbplot_logfits(lcolsur1)
jbplot_procdev(lcolsur1)

jbplot_trj(lcolsur1,type="BBmsy",add=T)
jbplot_trj(lcolsur1,type="FFmsy",add=T)

#lamsur1$pars_posterior[,1]
lcolsur1$estimates
lcolsur1$pars





#### SCHAEFER
jblcolsur_input2 = build_jabba(catch=lcolsur$catch,cpue=lcolsur$cpue,se=lcolsur$se,assessment=assessment,scenario = "S2",
                               model.type = "Schaefer",sigma.est = FALSE,r.prior=c(0.5,0.5),K.prior = c(70000,0.3),
                               igamma = c(0.001,0.001),Plim=0.2,psi.prior = c(2,0.1))
# Fit JABBA (here mostly default value - careful)
lcolsur2 = fit_jabba(jblcolsur_input2,save.jabba=TRUE,output.dir=output.dir,quickmcmc = TRUE) # quick run

head(lcolsur2$kobe)

# Make individual plots
jbplot_catch(lcolsur2)
jbplot_catcherror(lcolsur2)
jbplot_ppdist(lcolsur2)
jbplot_mcmc(lcolsur2)
jbplot_residuals(lcolsur2)
jbplot_cpuefits(lcolsur2)
jbplot_runstest(lcolsur2)
jbplot_logfits(lcolsur2)
jbplot_procdev(lcolsur2)

jbplot_trj(lcolsur2,type="BBmsy",add=T)
jbplot_trj(lcolsur2,type="FFmsy",add=T)

#lamsur1$pars_posterior[,1]
lcolsur2$estimates
lcolsur2$pars
