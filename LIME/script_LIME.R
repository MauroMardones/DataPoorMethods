##############################################################
########### Run length based models LBSPR and LIME ###########
########## mpons@uw.edu - September 2018 #####################
##############################################################
rm(list=ls())
########
R.version$os # to check how lucky you are ...
Sys.getenv('PATH')
session_info()

getRversion()


install.packages("LBSPR")
install.packages("devtools", repos='http://cran.us.r-project.org')
install.packages("TMB")
devtools::install_github("kaskr/TMB_contrib_R/TMBhelper", force=T)
devtools::install_github("merrillrudd/LIME")
#devtools::install_github("brodieG/fansi")



find_rtools() # should be TRUE, assuming you have Rtools 3.5
#8920811a4a029297c12eb16745054a29cdbe8b7e

#fb21b779f95a16bf77037218522135694af61e90

library(pkgbuild)
library(devtools)
library("LBSPR")
library("TMB")
library("LIME")
library(dplyr)
library(ggplot2)


### read length data 
Main.dir<-"C:/Users/mauricio.mardones/Documents/IFOP/Cursos/UW_Seattle/LBSPR_LIME" # your directory 
setwd(Main.dir)
Length.comps<-read.csv("length_struc_urchin.csv",header=T) 
head(Length.comps) #First column is the year, the following are the numbers of individuals in each lenght bin
tail(Length.comps)
#Biology
MaxAge=12
L50=43.2 #Jaramillo
L95=65 # verrificar bibliografia
M=0.25
h=0.8
wla=0.0005
wlb=2.97973
K=0.139
t0=(-0.45)
Linf=136
LenCV=0.1  
SigmaR=0.4 
SigmaF=0.2
SigmaC=0.2
SigmaI=0.2
R0=1 
qcoef=1e-5 
start_ages=0 
rho=0 
nseasons=1
binwidth=2
S50=65
S95=70
####################################################################
##### Length based methods #################
####################################################################
minL<-34
maxL<-136
Bins<- seq(from=minL, to= maxL, 
           by = binwidth)

#####################################################################
########## plot parametres #################
####################################################################

x11()
par(mfrow=c(2,2), mar=c(4,4,3,1))
plot(lh$L_a, type="l", lwd=4, xlab="Age", ylab="Length (cm)")
plot(lh$W_a, type="l", lwd=4, xlab="Age", ylab="Weight (g)")
plot(lh$Mat_l, type="l", lwd=4, xlab="Length (cm)", ylab="Proportion mature")
plot(lh$S_l, type="l", lwd=4, xlab="Length (cm)", ylab="Proportion selected to gear")

plot(lh$S_fl[1,], type="l", lwd=4, xlab="Length (cm)", ylab="Proportion selected to gear")
plba <- with(lh, age_length(highs, lows, L_a, CVlen))


tallas <- seq(40,138,2)

x11()
plot(tallas,plba[1,], type ="n", ylab="Probabilidad", xlab="Tallas(mm)")
for(i in 1:12){
  lines(tallas,plba[i,], col=i)
}



#####################################################################
##################### simulate data ################################
######################################################################

true <- generate_data(modpath=NULL,
                      itervec=1,
                      lh=lh,
                      Fdynamics="Ramp",
                      Rdynamics="AR",
                      Nyears=17,
                      Nyears_comp=17,
                      comp_sample=200,
                      init_depl=0.8,
                      seed=1)



## years with length data -- rename with your own years with length data
length_years <- rownames(true$LF)
## length bins -- rename with your upper length bins
length_bins <- colnames(true$LF)
###########################################################################
###### list of parameters to use for LBSPR and LIME (use create_lh_list)
lh <- create_lh_list(vbk=K, linf=Linf, t0=t0,
                     lwa=wla, lwb=wlb, 
                     S50=S50, S95=S95, selex_input="length", 
                     M50=L50, M95=L95, selex_type=c("logistic"),
                     maturity_input="length", M= M, 
                     SigmaR=SigmaR, SigmaF=SigmaF, SigmaC=SigmaC,
                     SigmaI=SigmaI,CVlen=LenCV, 
                     h=h, R0=R0, qcoef=qcoef,
                     start_ages=start_ages, rho=rho,
                     nseasons=nseasons, binwidth=binwidth, 
                     AgeMax= MaxAge,
                     Frate=0.1,
                     Fequil=0.25,
                     nfleets=1)

lfdata<-as.matrix(Length.comps[,-1])
colnames(lfdata)<-Bins

years <- Length.comps[,1]
rownames(lfdata)<-years
lf <- lfdata;lf

## input data
data_list <- list("years"=as.numeric(rownames(lf)), "LF"=lf)

## create input list -- adds some zeros on the end as well to make sure there is room for larger fish
inputs <- create_inputs(lh=lh, input_data=data_list)
mids<- seq(from=minL+binwidth/2, 
           to= maxL+binwidth/2, 
           by = binwidth)
inputs$mids<-mids
## run LIME
res <- run_LIME(modpath=NULL,
                input=inputs,
                data_avail="LC")

## check TMB inputs
Inputs <- res$Inputs
## Report file
Report <- res$Report
## Standard error report
Sdreport <- res$Sdreport

## check convergence
hessian <- Sdreport$pdHess
gradient <- res$opt$max_gradient <= 0.001
hessian == TRUE & gradient == TRUE

#Example: Length data only
lc_only <- run_LIME(modpath = NULL, input = inputs_all, data_avail = "LC")

###################################
## Run LBSPR
###################################

LB_pars <- new("LB_pars")
LB_pars@Species <- ""
LB_pars@MK <- lh$M/lh$vbk 
LB_pars@M <- lh$M
LB_pars@Linf <- lh$linf
LB_pars@CVLinf <- lh$CVlen
LB_pars@L50 <- lh$ML50 
LB_pars@L95 <- lh$ML95
LB_pars@Walpha <- lh$lwa
LB_pars@Wbeta <- lh$lwb
LB_pars@SL50 <- S50 
LB_pars@SL95 <- S95
LB_pars@BinWidth <- binwidth
LB_pars@Steepness <- 0.8
LB_pars@R0 <- 1
LB_pars@L_units<-"mm"
LB_lengths <- new("LB_lengths")
LB_lengths@LMids <- mids
LB_lengths@LData <- t(inputs$LF[,,1])
LB_lengths@Years <- data_list$years
LB_lengths@NYears <- length(data_list$years)

#RUN de model
lbspr_res <- LBSPRfit(LB_pars=LB_pars, LB_lengths=LB_lengths, yrs=NA, Control=list(modtype="GTG"))

plot_LCfits(Inputs=Inputs,
            Report=Report,
            LBSPR=lbspr_res)

plot_output(Inputs=Inputs,
            Report=Report,
            Sdreport=Sdreport,
            lh=lh,
            LBSPR=lbspr_res,
            plot=c("Fish","Rec","SPR","Selex"),
            set_ylim=list("SPR" = c(0,1)),
            true_years=inputs$years)

plot_LCfits(Inputs=list("LF"=true$LF))


plot_LCfits(LF_df = LF_df, Inputs = lc_only$Inputs, Report = lc_only$Report)
