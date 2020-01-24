#############################################################
# PCOM Posterior-focused catch-only method #
# S. Zhou, April 2014 #
#############################################################
# This method requires time series of catch data only.
# However, some life history parameters, such as M, Linf, k, T_max, T_maturation,
# will much improved the performance.
# Also, a rough guess of maximum depletion level D = B_end/K will be helpful.

# This example is for sinle catch series
# model
BDM = function(K, dep, C, r) { ### biomass dynamics model
  B = rep(NA, length(C)) #C
  B[1] = K
  for (i in 2:length(B)) {
    B[i] = max(min(B[i-1] + r*B[i-1]*(1-B[i-1]/K) - C[i-1], K), 0)
  }
  if (all(B > C) & all(B <=K)) abs(B[length(B)]-dep*K) else 10^5
}
# simulation
sim1=function(k25=k25, k75=k75, m25=m25, m75=m75, yr=yr, C=C,
              n.sim=n.sim) {
  Bend.keep=K.keep =r.keep=dep.keep= d.keep = vector()
  nyr=length(yr)
  B= matrix(NA, nyr, n.sim)
  K = r = vector()
  k.low = k25; k.up = k75
  r.low = r25; r.up = r75
  plot(0,0, type='n', xlim=c(min(yr),max(yr)), ylim=c(-1,
                                                      max(k.up)*1.2),xlab='Year', ylab='', las=1)
  nyr = length(C)
  for (j in 1: n.sim) {
    K[1] = runif(1, k.low, k.up )
    r[1] = runif(1, r.low, r.up)
    B[1,j]= K[1]
    for(i in 2 : (nyr)) {
      r[i] = runif(1, r.low, r.up) ;
      K[i] = runif(1, k.low, k.up )
      B[i,j] = B[i-1,j] + r[i]*B[i-1,j]*(1-B[i-1,j]/K[i]) - C[i-
                                                                1]
    }
    lines(yr, B[,j], col=rgb(runif(1,0,j)/n.sim,(n.simrunif(
      1,0,j))/n.sim, (1)/(n.sim+100), alpha=0.6) )
    
    lines(yr, C)
    K.keep[j] = mean(K)
    r.keep[j] = mean(r)
  }
  Bend.keep = B[nyr,]
  d.keep = B[nyr,]/mean(K)
  lines(yr, apply(B, 1, median), lty=2, lwd=2)
  mtext('Year', side=1, line=2, outer=T)
  mtext('Biomass', side=2, line=2, outer=T)
  msy = K.keep*r.keep/4
  return(list(K.keep, r.keep, msy, Bend.keep, d.keep, B))
}
# Jibia Squid
C=c(3476,5589,15191,174768,296575,250639,123726,144001,54545,19714
    4,134171,124376,141597)
yr= seq(2001, 2013)
r.lci=0.5; r.uci=1.5;
# search through K grids
N1 = 200
K1 = exp(seq(log(max(C)), log(max(C)*50), l=N1))
dep = round(seq(0.1, 0.8, 0.05),2)
nd=length(dep)
r1=obj1 = matrix(0, N1,nd)
for (j in 1:nd) {
  for (i in 1:N1) {
    out = optimize(BDM, K=K1[i], C=C, dep=dep[j], interval=c(0.01,1))
    r1[i,j] = out$min
    obj1[i,j] = out$obj }}
#############################
r1.backup=r1 ; r1=r1.backup
r1[obj1> K1*0.01 ]=NA
r1[r1 > r.uci | r1 < r.lci] = NA
kr = as.data.frame(cbind( K1, r1))
colnames(kr) = c('k', dep)
all=cbind(K1, stack(kr[,2:nd]))
colnames(all)=c('k', 'r', 'ind')
all$d=as.numeric(as.character(all$ind))
all=all[, c(1,2,4)]
all=all[!is.na(all[,2]),]
all$msy=all[,1]*all[,2]/4


output=matrix(NA,6,6)
d.lvl=c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
for (i in 1:6) {
  all2= all[!is.na(all$r) & all$d <= d.lvl[i],]
  quan1 = apply(all2, 2, quantile)
  k25=quan1[,1][2] ; k75=quan1[,1][4]
  r25=quan1[,2][2] ; r75=quan1[,2][4]
  msy.median=quan1[,4][3]
  n.sim =100
  out1 = sim1(k25=k25, k75=k75, m25=m25, m75=m75, C=C, yr=yr,
              n.sim=n.sim )
  abline(h=msy.median)
  sp = out1[1:5]
  sp= as.data.frame(sp)
  # colnames(sp)= c('k','r','msy','Bend', 'D')
  BendD= apply(sp, 2, quantile)
  output[i,] = (c(quan1[3,], BendD[3,4:5]) )
}
output[,3] = t(d.lvl)
colnames(output)=c('k','r','d.upper','msy', 'Bend', 'D')
output
# parameters as a function of assumed upper depletion level
par(mfrow=c(2,2))
plot(output[,3], output[,1], ylim=c(0, max(output[,1]*1.1)),
     type='l', xlab='Upper depletion', ylab='K', las=1)
plot(output[,3], output[,2], ylim=c(0, max(output[,2]*1.1)),
     type='l', xlab='Upper depletion', ylab='2', las=1)
plot(output[,3], output[,4], ylim=c(0, max(output[,4]*1.1)),
     type='l', xlab='Upper depletion', ylab='MSY', las=1)
plot(output[,3], output[,6], ylim=c(0, max(output[,6]*1.1)),
     type='l', xlab='Upper depletion', ylab='D', las=1)