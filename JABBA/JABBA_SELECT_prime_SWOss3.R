##><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
## Prime File for SWOss3 JABBA-SELECT example
## written by Henning Winker
##><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Delete all objects
rm(list=ls())
gc()
# Install required packages if missing
list.of.packages <- c("gplots", "coda","rjags","R2jags","fitdistrplus","reshape","mvtnorm","scales","reshape2","r4ss")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
# Load Packages
library(gplots);library(coda);library(rjags);library(R2jags);library("fitdistrplus");library(reshape);
library(mvtnorm);library(scales);library(reshape2);library(r4ss)


#---------------------------------------------------------------------
# Set Working directory file where to store the results
File = "C:/Work/Research/GitHub/JABBA-SELECTbeta"
# Set working directory for JABBA R source code
JABBA.file = "C:/Work/Research/GitHub/JABBA-SELECTbeta"
# Set Assessment
assessment = "SWOss3"
# Version
version = "v1.2beta"

#-------------------------------------------
# Load ss2jabba rdata for JABBA-Select
#-------------------------------------------

load(paste0(File,"/",assessment,"/ss4js_",assessment,".rdata"),verbose=T)


#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Graphic, Output, Saving (.RData) settings 
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
KOBE.plot = TRUE # Produces JABBA Kobe plot 
KOBE.type = c("ICCAT","IOTC")[2] # ICCAT uses 3 colors; IOTC 4 (incl. orange) 
Biplot= TRUE # Produces a "post-modern" biplot with buffer and target zones (Quinn & Collie 2005)
SP.plot = c("standard","phase")[2] # Produces standard or 'Kobe phase' SP plot  
save.trajectories =TRUE # saves posteriors of P=B/K, B/Bmsy and H/Hmsy as .RData object 
harvest.label = c("Hmsy","Fmsy")[2] # choose label preference H/Hmsy versus Fmsy
CPUE.plot= TRUE # Runs state-tool to produce "alligned" multi-CPUE plot  
meanCPUE = FALSE # Uses averaged CPUE from state-space tool instead of individual indices  
Projection = TRUE # Use Projections: requires to define TACs vectors 
save.projections = TRUE # saves projection posteriors as .RData object 
catch.metric = "(t)" # Define catch input metric e.g. (tons) "000 t" etc 
Reproduce.seed = FALSE # If FALSE a random seed assigned to each run, if TRUE set.seed(123)
runASEM=TRUE  
jabba2FRL = TRUE
# Save entire posterior as .RData object
save.all = FALSE # (if TRUE, a very large R object of entire posterior is saved)  
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>


# Choose Sceanrio name for creating a seperate folder
# Scenario 1: WG SEs
# Scenario 2: SEs ca
# Scenario 3: -JPN
# Scenario 4: -CAN_Late

Scenarios = c("base_wg")
s=1

for(s in 1:length(Scenarios)){
#runASEM=ifelse(s==1,TRUE,FALSE) # run ASEM prior MCMC on first run  
  
Scenario = Scenarios[s] 
Mod.names = c("JS") 

  
 #--------------------------------------------------
 # Read csv files
 #--------------------------------------------------

 # get SELECT info
 select = ss4js$select

 
 # Load assessment data
 catch = read.csv(paste0(File,"/",assessment,"/catch",assessment,".csv"))
 cpue = read.csv(paste0(File,"/",assessment,"/cpue",assessment,".csv"))#

 # Use SEs from csv file for abudance indices (TRUE/FALSE)
 SE.I = TRUE
 
 if(SE.I ==TRUE){
  se = read.csv(paste0(File,"/",assessment,"/se",assessment,".csv"))
 }
 
 # Read select csv
 select = read.csv(paste0(File,"/",assessment,"/select",assessment,".csv"))
 
 # Read selex (selectivity) csv
 selex = read.csv(paste0(File,"/",assessment,"/selex",assessment,".csv"))
 
 
 
 names(cpue)
 names(catch)

  
  #---------------------------------------
  # option to exclude CPUE time series 
  #---------------------------------------
  
  #------------------------------------------------------
  # Option use mean CPUE from state-space cpue averaging
  #-----------------------------------------------------
  # Produce CPUE plot average plot
  CPUE.plot = TRUE
  
  #------------------------------------------------------
  # mean and CV and sd for unfished spawning biomass SB0
  #------------------------------------------------------
  mu.SB0 = 200000; CV.SB0 = 2; sd.SB0=sqrt(log(CV.SB0^2+1)) 
  
  #-----------------------------------------------------------
  # mean and CV and sd for Initial depletion level P1= SB/SB0
  #-----------------------------------------------------------
  # Set the initial depletion prior SB1/SB0 
  # To be converted into a lognormal prior (with upper bound at 1.1)
  # psi.prior = "lnorm"
  # or to be converted into a Beta prior
  # psi.prior = "beta"
  
  psi.prior= "lnorm"
  # specify as mean and CV 
  mu.psi = 1
  CV.psi = 0.05
  
  P_bound = c(0.03,1.2)
  #--------------------------------------------------------------
  # Determine estimation for catchability q 
  #--------------------------------------------------------------
  # Assign q to abundance
  sets.q = 1:(ncol(cpue)-1)  # here 1: South Early+Recent, 2: South-East Early+Recent  
  
  
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
  # Observation Error
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
  
  #To Estimate additional observation variance set sigma.est = TRUE
  sigma.est = c(TRUE,FALSE)[2]
  
  # Assign common variance estimates to abundance indices
  # Here it assumed that same flagged fleets have the same additional observation variance
  sets.var = c(1,2,2,3,3,3,4,5,6)
  
  # As option for data-weighing
  # minimum fixed observation error for each variance set (optional choose 1 value for both)
  fixed.obsE = c(0.01) # Important if SE.I is not availble
  
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
  # Process Error
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
  #Estimate set sigma.proc == True
  sigma.proc = TRUE
  proc.dev.all=1979 # start process error in 1979 as in ss3 assessment
  #------------------------------------------
  if(sigma.proc == TRUE){
    proc.type = c("igamma","lnorm")[1] # choose 1: inverse-gamma or 2: lognormal
    
    if(proc.type=="lnorm"){
    pr.proc = c(log(0.08),0.2) # Option for lognormal process error
    }
    if(proc.type=="igamma"){
    #pr.proc = c(0.001,0.001) # Option for inverse-gamma prior
      #pr.proc = c(10,0.1)
      pr.proc = c(0.001,0.001)
      gamma.check = 1/rgamma(1000,pr.proc[1],pr.proc[2]) # Process error check
    # check mean process error + CV
    mu.proc = sqrt(mean(gamma.check)); CV.proc = sd(sqrt(gamma.check))/mean(sqrt(gamma.check))
    # check CV
    round(c(mu.proc,CV.proc),3)
    quantile(sqrt(gamma.check),c(0.1,0.9))
    }  
    }else{
    sigma.proc = 0.05 #IF Fixed (sigma.est = FALSE): typicallly 0.05-0.15 (see Ono et al. 2012)
  
    }
  #--------------------------------------------
  
  
  
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
  # Prior specification for Model 5: JABBA-SELECT
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
  # The following section is ignored if: Model < 5
    
    # Option to set SBmsy reference point (e.g. SBmsy/SB0 = 0.4)
    # If SBmsy_SB0 = NULL the reference SBmsy_SB0 that produces MSY will be used 
    SBmsy_SB0 = NULL   # Standard Reference setting for SA Linefish assessment reference points
    
  #---------------------------------------------------------------
  # STOCK PARAMETERS for prior generation Hmsy as a fuction of r
  #---------------------------------------------------------------
    
    minage <- 0  																						
    maxage <- ss4js$stock.pars$Amax
    
    # Number of sexes order Female and Males (this model is sex-structured)
    nsexes = 2
    # VBGF parameters
    Linf <- c(ss4js$stock.pars$vbgf_F[1],ss4js$stock.pars$vbgf_M[1])
    kappa <- c(ss4js$stock.pars$vbgf_F[2],ss4js$stock.pars$vbgf_M[2])
    t0 <- c(ss4js$stock.pars$vbgf_F[3],ss4js$stock.pars$vbgf_M[3])
    
    #Length-weight
    aW <- exp(c(ss4js$stock.pars$LW_F[1],ss4js$stock.pars$LW_F[1])) 																						
    bW <- c(ss4js$stock.pars$LW_F[2],ss4js$stock.pars$LW_F[2])
    
    # Maturity 
    #maturity = 3 # a single value is taken as age (knife-edge) 
    maturity =ss4js$stock.pars$mat  #c(Lm50,Lm95) # two values are taken as length-based: Lm50 and Lm95 
    
    # Natural mortality estimate (affects Hmsy)
    M = 0.2
    CV.M = 0.001 # ss3 approach M "fixed"
    
    # steepness B&H SSR (effects Hmsy and determines SBmsy/SB0)
    h = 0.8
    CV.h = 0.06
    
 #------------------------------------------------
 # Selectivity will determine changes in r (Fmsy)
 #------------------------------------------------
    
    # Selectivity SL50 must be sufficiently different (+-5%) between "fleets" to seperate r 

    # only unique SL50 values are permitted (no replicates) 
    SL50 <- as.numeric(selex[1,-1]) 
    SL95 <- as.numeric(selex[2,-1])  # If unknown set to 0.05*SL50 ~ knife-edge
    
    # Define point where descening limb starts (set Linf for logistic)
    SL.desc <-as.numeric(selex[3,-1])  # mean of half-normal 
    # Define rate of decreasing selectivity
    CV.desc <- as.numeric(selex[4,-1]) # CV of half-normal 
    # Define minimum descending limp between 0 and 1
    min.desc =as.numeric(selex[5,-1]) 
    
    # number of different Hmsy (r) priors 
    nSel = length(SL50)
     
    
  
    # Assign Selectivity to abundance indices
    sets.I = select$Selectivity[select$CPUE] 
    
    # Assign Selectivity to catch series
    sets.C = select$Selectivity[select$Catch] # here 1: South, 2: South-East, 3: Trawl
    # Define if index is in numbers: 0 or biomass: 1 
    I.unit = aggregate(CPUE.units~Selectivity,data=select[select$CPUE,],mean)[,2] 
      
    
    
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
  # Optional: Do TAC Projections
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
  Projection = TRUE # Switch on by Projection = TRUE 
  
  # Check final year catch
  apply(catch[,-1],1,sum,na.rm=TRUE)

  # Set range for alternative TAC projections
  TACs = seq(9000,12500,500) #example
  
  # Set year of first TAC implementation
  imp.yr = 2018
  # Set number of projections years
  pyrs = 10
  
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
  # Execute model and produce output
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
  

  
  # MCMC settings
  ni <- 25000 # Number of iterations
  nt <- 2 # Steps saved
  nb <- 3000 # Burn-in
  nc <- 1 # number of chains
  nsaved = (ni-nb)/nt*nc
cat(paste0("\n","- Run Model","\n"))

source(paste0(JABBA.file,"/JABBA_SELECT",version,".r"))

}


#--------------------------------------------
#Compare JABBA-Select base with ss3 base
# JABBA-Select base
#
#Swapping base with hprior0.8 assuming swnbase is the 2017 SWO run with estimated h
#--------------------------------------------



# load
sc =1
output.dir = paste0(File,"/",assessment,"/",Scenarios[sc],"_",Mod.names,"/Output")

load(paste0(output.dir,"/",assessment,"_",Mod.names,"_",Scenarios[sc],"_trajectories.rdata"),verbose=T)

load(paste0(File,"/",assessment,"/swnbase.Rdata"),verbose=TRUE) 
# SS3 base
yr = swnbase$sprseries$Yr[swnbase$sprseries$Yr<2018]
SSB0 = swnbase$sprseries$SSBzero[swnbase$sprseries$Yr<2018]
SSB = swnbase$sprseries$SSB[swnbase$sprseries$Yr<2018]
kobe.ss = swnbase$Kobe[swnbase$Kobe$Year==2017,]
status.ss = swnbase$Kobe[swnbase$Kobe$Year<2018,]
# JABBA output  
ts = jabbaS_out$Stock.trj
Yr = ts$Yr
SB = cbind(ts$SBt.Median,ts$SBt.LCI95.,ts$SBt.UCI95.)
SBtoSB0 = cbind(ts$SBt_SB0.Median,ts$SBt_SB0.LCI95.,ts$SBt_SB0.UCI95.)
SBtoSBmsy= cbind(ts$SBt_SBmsy.Media,ts$SBt_SBmsy.LCI95.,ts$SBt_SBmsy.UCI95.)
kobe.js = kb[kb$year==2017,]

# Compare ss3 wg_base with JABBA-Select
Par = list(mfrow=c(2,2),mai=c(0.5,0.5,0.1,.1),omi = c(0.1,0.1,0.1,0.1) + 0.1,mgp=c(2,0.5,0), tck = -0.02,cex=0.8)
png(file = paste0(File,"/",assessment,"/Status_SS3vsJabbaS",Scenarios[sc],".png"), width = 9, height =8, 
    res = 200, units = "in")
par(Par)
plot(Yr,SB[,1],type="n",ylim=c(0,max(SB,SSB)),ylab="SBB",lwd=2)
polygon(c(Yr,rev(Yr)),c(SB[,2],rev(SB[,3])),border=0,col="grey")
lines(Yr,SB[,1],lwd=2)
lines(yr,SSB,col=4,lwd=2,lty=2)
legend("topright",c("JABBA-Select","wg_ss3"),lwd=2,col=c(1,4),lty=c(1,2),bty="n")
plot(Yr,SBtoSB0[,1],type="n",ylim=c(0,max(SBtoSB0,SSB/SSB0)),ylab="SBB/SSB0",lwd=2)
polygon(c(Yr,rev(Yr)),c(SBtoSB0[,2],rev(SBtoSB0[,3])),border=0,col="grey")
lines(yr,SSB/SSB0,lwd=2,col=4,lty=2)
lines(Yr,SBtoSB0[,1],col=1,lwd=2)
plot(Yr,SBtoSBmsy[,1],type="n",ylim=c(0,max(SBtoSBmsy,status.ss[,2])),ylab="SBB/SSBmsy",lwd=2)
polygon(c(Yr,rev(Yr)),c(SBtoSBmsy[,2],rev(SBtoSBmsy[,3])),border=0,col="grey")
lines(yr,status.ss$B.Bmsy,lwd=2,col=4,lty=2)
lines(Yr,SBtoSBmsy[,1],col=1,lwd=2)
abline(h=1,lty=2)

b = kobe.js$stock; f = kobe.js$harvest
refB ="MSY"
kernelF <- ci2d(b,f,nbins=151,factor=1.5,ci.levels=c(0.50,0.80,0.75,0.90,0.95),show="none",col=1,xlab= ifelse(harvest.label=="Fmsy",expression(paste(F/F[MSY])),expression(paste(H/H[MSY]))),ylab=expression(paste(B/B[MSY])))
#Create plot
plot(1000,1000,type="b", xlim=c(0,2.2), ylim=c(0,2.2),lty=3,ylab=bquote(H/H[.(refB)]),xlab=bquote(SB/SB[.(refB)]),xaxs="i",yaxs="i")
c1 <- c(-1,100)
c2 <- c(1,1)

# extract interval information from ci2d object
# and fill areas using the polygon function
zb2 = c(0,1)
zf2  = c(1,100)
zb1 = c(1,100)
zf1  = c(0,1)
polygon(c(zb1,rev(zb1)),c(0,0,1,1),col="green",border=0)
polygon(c(zb2,rev(zb2)),c(0,0,1,1),col="yellow",border=0)
polygon(c(1,100,100,1),c(1,1,100,100),col=ifelse(KOBE.type=="ICCAT","yellow","orange"),border=0)
polygon(c(0,1,1,0),c(1,1,100,100),col="red",border=0)

polygon(kernelF$contours$"0.95",lty=2,border=NA,col="cornsilk4")
polygon(kernelF$contours$"0.8",border=NA,lty=2,col="grey")
polygon(kernelF$contours$"0.5",border=NA,lty=2,col="cornsilk2")
points(median(b),median(f),cex=2,pch=16)
points(kobe.ss$B.Bmsy,kobe.ss$F.Fmsy,cex=2,col=4,pch=16)
legend("bottomleft",c("JABBA-Select","wg_ss3"),pch=16,col=c(1,4),bty="n",pt.cex = 2)

# Get Propability
Pr.green = sum(ifelse(b>1 & f<1,1,0))/length(b)*100
Pr.red = sum(ifelse(b<1 & f>1,1,0))/length(b)*100

if(KOBE.type=="ICCAT"){               
  Pr.yellow = (sum(ifelse(b<1 & f<1,1,0))+sum(ifelse(b>1 & f>1,1,0)))/length(b)*100} else {
    Pr.yellow = sum(ifelse(b<1 & f<1,1,0))/length(b)*100
    Pr.orange = sum(ifelse(b>1 & f>1,1,0))/length(b)*100
  }
## Add legend
legend('topright', 
         c("50% C.I.","80% C.I.","95% C.I.",paste0(round(c(Pr.red,Pr.yellow,Pr.orange,Pr.green),1),"%")), 
         lty=rep(-1,8),pch=c(rep(22,8)),pt.bg=c("cornsilk2","grey","cornsilk4","red","yellow","orange","green"), 
         col=1,lwd=1.1,cex=0.9,pt.cex=c(rep(1.7,3),rep(2.2,4)),bty="n")  
  

dev.off()




