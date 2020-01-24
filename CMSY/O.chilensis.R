rm(list=ls(all=TRUE))
setwd<-("C:/Users/mauricio.mardones/Documents/IFOP/Prog_Seg_PM_Bentonicos_2017/Stock_assessment_Ancud/Ostra")
dir()
setwd(dir.b)
cdat <- read.csv("catch_ostra_2.csv", sep= ";", header=T)
#cat("catch_raya.csv", File= "catch_ostra_2.csv", "msy_functions")
View(cdat)




                                        #t.span <- 1976:2014
library(lattice)
library(ggplot2)

yr   <- cdat$yr; yr
ct   <- cdat$Ct; ct
nyr  <- length(yr)
cbind(yr,ct)
mean(ct)



png(filename = "C:/Users/mauricio.mardones/Documents/IFOP/Prog_Seg_PM_Bentonicos_2017/Stock_assessment_Ancud/Ostra/Desembarques.png",
    width =14, height =6, units = "in", res = 300)
#x11()
par(mar=c(3,6,1,1))
plot(yr, ct, type="l", ylim = c(0, max(ct)), xlab="", ylab = "Desembarques Ostra (t.)", cex.lab=2, lwd=3, xaxp=c(1976,2017,41))
abline(h ="184", col ="red", lwd=3)#legend('topright', 'A', cex=1.2, bty='n')

dev.off()


#Application of Catch-MSY of S Martell and R Froese
# FISH and FISHERIES 2012 (14)4:504-514.
# source code http://www.fishbase.de/rfroese/


## PARAMETER SECTION

start_r <- c(0.05, 0.5)     #mid resilience
#start_r <- c(0.015, 0.1) # low resilience
## start_r <-c(0.05,0.5) #  mid  resilience will be
## start_r <- c(0.6,1.5) # high resilience
## start_r <- c(0.2,1) # default resilience

start_k     <- c(max(ct),100*max(ct)) ## default for upper k e.g. 100 * max catch
#startbio 	<- c(0.8,1)   ## assumed biomass range at start of time series, as fraction of k
startbio 	<- c(0.5,0.9)   ## assumed biomass range at start of time series, as fraction of k
interyr 	<- yr[2]   ## interim year within time series for which biomass estimate is available; set to yr[2] if no estimates are available
interbio 	<- c(0, 1) ## biomass range for interim year, as fraction of k; set to 0 and 1 if not available
## finalbio 	<- c(0.8, 0.9) ## biomass range after last catches, as fraction of k
#finalbio    <- if(ct[nyr]/max(ct) > 0.5) {c(0.3,0.7)} else {c(0.01,0.4)} ## use for batch processing
finalbio    <- if(ct[nyr]/max(ct) > 0.5) {c(0.3,0.7)} else {c(0.01,0.3)} ## use for batch processing

#n           <- 30000  ## number of iterations, e.g. 100000
n           <- 10000  ## number of iterations, e.g. 100000
sigR        <- 0.0      ## process error; 0 if deterministic model; 0.05 reasonable value? 0.2 is too high

startbt     <- seq(startbio[1], startbio[2], by = 0.05) ## apply range of start biomass in steps of 0.05
parbound <- list(r = start_r, k = start_k, lambda = finalbio, sigR)

cat("Last year =",max(yr),", last catch =",  1000*ct[nyr],"\n")
#cat("Resilience =",res,"\n")
cat("Process error =", sigR,"\n")
cat("Assumed initial biomass (B/k) =", startbio[1],"-", startbio[2], " k","\n")
cat("Assumed intermediate biomass (B/k) in", interyr, " =", interbio[1],"-",interbio[2]," k","\n")
cat("Assumed final biomass (B/k) =", parbound$lambda[1],"-",parbound$lambda[2]," k","\n")
cat("Initial bounds for r =", parbound$r[1], "-", parbound$r[2],"\n")
cat("Initial bounds for k =", format(1000*parbound$k[1], digits=3), "-", format(1000*parbound$k[2],digits=3),"\n")

source("msy_functions.R")

## MAIN
R1 = sraMSY(parbound, n)

#length(r1)
## Get statistics on r, k, MSY and determine new bounds for r and k
r1 	<- R1$r[R1$ell==1]
k1 	<- R1$k[R1$ell==1]
msy1  <- r1*k1/4
mean_msy1 <- exp(mean(log(msy1)))
max_k1a  <- min(k1[r1<1.1*parbound$r[1]]) ## smallest k1 near initial lower bound of r
max_k1b  <- max(k1[r1*k1/4<mean_msy1]) ## largest k1 that gives mean MSY
max_k1 <- if(max_k1a < max_k1b) {max_k1a} else {max_k1b}

if(length(r1)<10) {
cat("Too few (", length(r1), ") possible r-k combinations, check input parameters","\n")
flush.console()
}

if(length(r1)>=10) {

	## set new upper bound of r to 1.2 max r1
	parbound$r[2] <- 1.2*max(r1)
	## set new lower bound for k to 0.9 min k1 and upper bound to max_k1
	parbound$k<- c(0.9 * min(k1), max_k1)


	cat("First MSY =", format(mean_msy1, digits=3),"\n")
	cat("First r =", format(exp(mean(log(r1))), digits=3),"\n")
	cat("New upper bound for r =", format(parbound$r[2],digits=2),"\n")
	cat("New range for k =", format(parbound$k[1], digits=3), "-", format(parbound$k[2],digits=3),"\n")


## Repeat analysis with new r-k bounds
R1 = sraMSY(parbound, n)

## Get statistics on r, k and msy
r = R1$r[R1$ell==1]
k = R1$k[R1$ell==1]
msy = r * k / 4
mean_ln_msy = mean(log(msy))

	## plot MSY over catch data

#png(filename = "C:/Users/mauricio.mardones/Documents/IFOP/FIP 2014_19/Ostrea chilensis/Catch MSY/Output_M&F_ostrea1.png",
#    width =900, height =800, units = "px", pointsize = 20)

x11()

	par(mfcol=c(3,1))
	plot(yr, ct, type="l", ylim = c(0, 1000), xlab = "Año", ylab = "Captura (1000 t)", main = " ")
	abline(h=exp(mean(log(msy))),col="red", lwd=1)
	abline(h=exp(mean_ln_msy - 2 * sd(log(msy))),col="red", lty=3)
	abline(h=exp(mean_ln_msy + 2 * sd(log(msy))),col="red", lty=3)
        legend('topright', 'A', cex=1.2, bty='n')

        
# 	hist(r, freq=F, xlim=c(0, 1.2 * max(r)),  ylab='Densidad', main = "")
# 	abline(v=exp(mean(log(r))),col="red",lwd=1)
# 	abline(v=exp(mean(log(r))-2*sd(log(r))),col="red", lty=3)
# 	abline(v=exp(mean(log(r))+2*sd(log(r))),col="red", lty=3)
#         legend('topright', 'D', cex=1.2, bty='n')

# 	plot(r1, k1, xlim = start_r, ylim = start_k, xlab="r", ylab="k (1000t)")
#         legend('topright', 'B', cex=1.2, bty='n')
# 
# 	hist(k, freq=F, xlim=c(0, 1.2 * max(k)), xlab="k (1000t)", ylab='Densidad', main = "")
# 	abline(v=exp(mean(log(k))),col="red", lwd=1)
# 	abline(v=exp(mean(log(k))-2*sd(log(k))),col="red", lty=3)
# 	abline(v=exp(mean(log(k))+2*sd(log(k))),col="red", lty=3)
#         legend('topright', 'E', cex=1.2, bty='n')

	plot(log(r), log(k),xlab="ln(r)",ylab="ln(k)")
	abline(v=mean(log(r)))
	abline(h=mean(log(k)))
	abline(mean(log(msy))+log(4),-1, col="red", lwd=1)
	abline(mean(log(msy))-2*sd(log(msy))+log(4),-1, col="red", lty=3)
	abline(mean(log(msy))+2*sd(log(msy))+log(4),-1, col="red", lty=3)
        legend('topright', 'B', cex=1.2, bty='n')

	hist(msy, freq=F, xlim=c(0, 1.2 * max(msy)), xlab="RMS (1000t)", ylab='Densidad', main = "")
	abline(v=exp(mean(log(msy))),col="red", lwd=1)
	abline(v=exp(mean_ln_msy - 2 * sd(log(msy))),col="red", lty=3)
	abline(v=exp(mean_ln_msy + 2 * sd(log(msy))),col="red", lty=3)
        legend('topright', 'C', cex=1.2, bty='n')

        #dev.off()
# cat("Possible combinations = ", length(r),"\n")
# cat("geom. mean r =", format(exp(mean(log(r))),digits=3), "\n")
# cat("r +/- 2 SD =", format(exp(mean(log(r))-2*sd(log(r))),digits=3),"-",format(exp(mean(log(r))+2*sd(log(r))),digits=3), "\n")
# cat("geom. mean k =", format(1000*exp(mean(log(k))),digits=3), "\n")
# cat("k +/- 2 SD =", format(1000*exp(mean(log(k))-2*sd(log(k))),digits=3),"-",format(1000*exp(mean(log(k))+2*sd(log(k))),digits=3), "\n")
# cat("geom. mean MSY =", format(1000*exp(mean(log(msy))),digits=3),"\n")
# cat("MSY +/- 2 SD =", format(1000*exp(mean_ln_msy - 2 * sd(log(msy))),digits=3), "-", format(1000*exp(mean_ln_msy + 2 * sd(log(msy))),digits=3), "\n")
#  
# Write results into outfile, in append mode (no header in file, existing files will be continued)
#output = data.frame( sigR, startbio[1], startbio[2], interbio[1], interbio[2], finalbio[1], finalbio[2], min(yr), max(yr), 
#	                    res, max(ct), ct[1], ct[nyr], length(r), exp(mean(log(r))), sd(log(r)), min(r), quantile(r,0.05), quantile(r,0.25), median(r),
#                    quantile(r,0.75), quantile(r,0.95), max(r), exp(mean(log(k))), sd(log(k)), min(k), quantile(k, 0.05), quantile(k, 0.25), 
#	                    median(k), quantile(k, 0.75), quantile(k, 0.95), max(k), exp(mean(log(msy))), sd(log(msy)), min(msy), quantile(msy, 0.05),
#	                    quantile(msy, 0.25), median(msy), quantile(msy, 0.75), quantile(msy, 0.95), max(msy)) 
#	
#write.table((output), file = "outfile_ostra", append = TRUE, sep = ";", dec = ".", row.names = FALSE, col.names = FALSE)
	
}	
	
	    
 
rms <- (R1$r * R1$k * 1/4)
quantile(rms)



#Calcula la Bt por año en base a los valores medianos del análisis
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ct = c(3,2,7,159,85,219,209,321,435,222,343,225,452,119,122,368,75,599,92,110,153,147,119,197,181,170,170,
       152,164,100,111,111,134,105,111,189,181,153,196,186, 284, 284)
length(ct)
k.msy = R1$k[R1$ell==1]
r.msy = R1$r[R1$ell==1]
yrs=seq(1976,2017,1)
bt=matrix(0,length(k.msy),length(yrs))
btem=matrix(0,40,1)
xt=rnorm(1,0, 0)

for (j in 1:length(k.msy)){
  k=k.msy[j]
  r=r.msy[j]
  btem=matrix(0,40,1)
  btem[1]=k
  
  for (i in 1:length(yrs)){
    
    btem[i+1]=(btem[i]+r*btem[i]*(1-btem[i]/k)-ct[i])*exp(xt) #mean
  }
  
  bt[j,]=t(btem[-4])
}

# Obtengo los percentiles para construir los intervalos de confianza

intcon=matrix(0,5,length(yrs))

for (i in 1:length(yrs)){
  intcon[,i]=as.numeric(quantile(bt[,i],probs=c(.05,.25,.50,.75,.95)))
}



png(filename = "C:/Users/mauricio.mardones/Documents/IFOP/Prog_Seg_PM_Bentonicos_2017/Stock_assessment_Ancud/Ostra/Biomasas.png",
    width =14, height =6, units = "in", res = 300)
#x11()
par(mar=c(4,5,1,1),cex=1.0)
plot(1976:2017,intcon[3,],type='l',ylim =c(0,6000), xlim=c(1976,2017),lwd=3,ylab='Biomasa (t)', cex.lab=2, xlab='', xaxp=c(1976,2017,41))
lines(1976:2017,intcon[1,],col='red',lty=2, lwd=3)
lines(1976:2017,intcon[5,],col='red',lty=2, lwd=3)

dev.off()



# Status / Diagrama de Fase de acuerdo a los PBR de un modelo de Schaefer. Lo recomedando
# en el Taller de PBR
# -----------------------------------------------------------------------------------------

k = median(R1$k[R1$ell==1])
r = median(R1$r[R1$ell==1])
yrs=seq(1976,2017,1)

Bmsy=k/2
Fmsy=r/2
Blim=Bmsy/2
Flim=1.5*Fmsy
ct1 = c(3,2,7,159,85,219,209,321,435,222,343,225,452,119,122,368,75,599,92,110,153,147,119,197,181,170,170,
        152,164,100,111,111,134,105,111,189,181,153,196,186,284,284)

Ft=ct1/intcon[3,]

B_Bmrs=intcon[3,]/Bmsy
F_Fmrs<- Ft/Fmsy


Ft=ct1/intcon[3,]

B_Bmrs=intcon[3,]/Bmsy
F_Fmrs<- Ft/Fmsy

#png(filename = "F:/mariella/IFOP_Mlla/IFOP/Proy_2015/2015CBA_Reineta/2ndReport_Reineta/SA_DataPoor/Code_Martell&Froese/DKOBE1.png",
#    width =14, height =10, units = "in", res = 400)

#png("F:/mariella/IFOP_Mlla/IFOP/Proy_2015/2015CBA_Reineta/2ndReport_Reineta/SA_DataPoor/Code_Martell&Froese/DKOBE1.png",
#   width=14, height=10, units="in",res=600)



x11()
par(mfrow=c(1,1),mar=c(4,4.5,2,2),cex=1.5) 
plot(B_Bmrs,F_Fmrs,col='blue',lwd=2, xlim=c(0,2.5),xaxs='i',yaxs='i',
     ylim=c(0,8),main="DIAGRAMA DE FASE",cex.main="1",xlab=expression(paste(B/ B[MRS])),
     ylab=expression(paste(F/ F[MRS])),type='n')
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col ='white')

points(B_Bmrs,F_Fmrs, pch=21, col='black', bg = "black", cex=0.9)
lines(B_Bmrs,F_Fmrs, lty=3, col='darkgrey',lwd=1)

segments(0, 1, 2.5, 1,col='black',lty=2,lwd=1)   #MRS
segments(1, 0,1, 11,col='black',lty=2,lwd=1)      #MRS

segments(0.50, 0, 0.50, 11,col='red',lty=2,lwd=2,cex=1.4) #colapso

mtext('SOBRE-EXPLOTACIÓN', side = 2, line = -6.5, outer = FALSE,adj=0.2,cex=1,
      col="red",font=2)
mtext('SUB-EXPLOTACIÓN', side=1, line=-2,outer=FALSE,adj=1.0,cex=1,col="darkgreen",font=2)
mtext('PLENA', side=1, line=-3.0,outer=FALSE,adj=0.47,cex=1,col="green",font=2)
mtext('EXPLOTACIÓN', side=1, line=-2.5,outer=FALSE,adj=0.55,cex=1,col="green",font=2)
text(2.0,1.1,'SOBREPESCA',col='black',cex=0.5)
text(B_Bmrs[40]+0.16,F_Fmrs[40]+0.1,'2015',col='black',cex=1,font=2)
points(B_Bmrs[40]-0.01,F_Fmrs[40],pch=19,col='red',cex=1.4)

#dev.off()

cx    <-c(0,0.5,0.5,0,0)        #colapso
cy    <-c(-0.3,-0.3,9,9,-0.3)   #colapso
sex   <-c(0.5,1,1,0.5,0.5)      #sobre-explotación
sey   <-c(-0.3,-0.3,9,9,-0.3)   #sobre-explotación
spx   <-c(1,3,3,1,1)            #sobrepesca
spy   <-c(1,1,9,9,1)            #sobrepesca
subex <-c(1,3,3,1,1)            #subexplotado
subey <-c(-0.3,-0.3,1,1,-0.3)   #subexplotado

#x11()

png(filename = "C:/Users/mauricio.mardones/Documents/IFOP/Prog_Seg_PM_Bentonicos_2017/Stock_assessment_Ancud/Ostra/KobePlot.png",
    width =14, height =6, units = "in", res = 300)

par(mar=c(5,5,1,1)+0.5)
plot(B_Bmrs,F_Fmrs,type="o",xlab=expression(B/B[MRS]),ylab=expression(F/F[MRS]),las=1,yaxp=c(0,9,9), xlim=c(0.3,2.5), ylim=c(0.1, 5), cex.lab=2, pch=21)
polygon(cx,cy,col="darkgrey")
polygon(sex,sey,col="dimgray")
polygon(spx,spy,col="gray32")
polygon(subex,subey,col="gray72")
lines(B_Bmrs,F_Fmrs,type="o",pch=16,bg="black",cex=2)
abline(h=1,v=c(0.5,1),lty=2)
points(c(B_Bmrs[1],B_Bmrs[length(yrs)]),c(F_Fmrs[1],F_Fmrs[length(yrs)]),col=c(1,1),bg=c(3,2),pch=21,cex=2.5)
text(B_Bmrs,F_Fmrs+0.2,yrs,cex=1, lwd=2)
text(2.3,0.5,"SUB-EXPLOTADO", cex = 1.5)
text(1.8,3,"SOBREPESCA", cex = 1.5)
text(0.7,3.5,"SOBRE-EXPLOTADO", cex = 1.5)
text(0.35,0.5,"COLAPSADO", cex = 1.5)
box()

dev.off()





