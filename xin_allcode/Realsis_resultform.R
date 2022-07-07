setwd("/Users/xinwang/Documents/newto_copy3")
p = 5e3
subsize = 1e3
rounds =100
n = c(1e4,2e4,4e4,8e4,1e5)
k = c(1e3,2e3,3e3,4e3,5e3)
X_distribution = 'fast_t'
datetime = '10-17'
out=list()
ki=1e3
ni=1e5
EPE = data.frame(MSE = numeric(),Method=character(),n_total=numeric())
SENS = data.frame(SENS = numeric(),Method=character(),n_total=numeric(), subdata_size = numeric())
SPEC = data.frame(SPEC = numeric(),Method=character(),n_total=numeric(), subdata_size = numeric())
time = data.frame(lasso_time = numeric(),select_time = numeric(),Method=character(),n_total=numeric())
p_selected = data.frame(p_selected = numeric(),Method=character(),n_total=numeric(), subdata_size = numeric())
methods = c('rand005','sisfull','sislev005','corsis005')
methods_adj = c('SIS-UNIF(250)','SIS-FULL(250)','SIS-LEV(250)','SIS-IBOSS(250)')
for (i in 1:5){
  #ni=n[i]
  ki = k[i]
  out=readRDS(paste0(substr(datetime,1,5),"_n_",ni,"_p_",p,"_k_",ki,"_d_",X_distribution,"_out.rds"))
  records=readRDS(paste0(substr(datetime,1,5),"_n_",ni,"_p_",p,"_k_",ki,"_d_",X_distribution,"_recorsisd.rds"))
  for(j in 1:length(methods)){
    m = methods[j]
    m2 = methods_adj[j]
    EPE = rbind(EPE,data.frame(MSE=out[m][[1]]['EPE'],Method=as.character(methods_adj[j]),n_total=ni,subdata_size =ki))
    EPE[EPE$Method==m2&EPE$n_total==ni&EPE$subdata_size==ki,'MSE'] =mean(sort(records[[m]]$EPE,decreasing = T)[-c(1:4)]) 
    p_selected = rbind(p_selected,data.frame(p_selected=out[m][[1]]['s'],Method=as.character(methods_adj[j]),subdata_size=ki,n_total=ni))
    time = rbind(time,data.frame(lasso_time=out[m][[1]]['lasso_time'],select_time=out[m][[1]]['slct_time'],Method=as.character(methods_adj[j]),subdata_size=ki,n_total=ni))
    SENS = rbind(SENS,data.frame(SENS=out[m][[1]]['TPR'],Method=as.character(methods_adj[j]),subdata_size=ki,n_total=ni))
    SPEC = rbind(SPEC,data.frame(SPEC=out[m][[1]]['SPC'],Method=as.character(methods_adj[j]),subdata_size=ki,n_total=ni))
    # p_selected[p_selected$Method==m2&p_selected$n_total==ni&p_selected$subdata_size==ki,'p_selected'] =mean(records[[m]]$s[records[[m]]$s>0]) 
  }
}


EPE

plot = 
  
  ggplot(EPE,aes(x =n_total,y=log(MSE,10),colours=Method,group=Method,linetype = Method))+
  geom_point(aes(shape=Method))+geom_line(aes())+ggtitle(paste0("MSE Comparison, p = ",p," k = ",subsize," x ~ ",substr(X_distribution,6,7))) +
  theme(text = element_text(size=14),
        axis.text.x = element_text(angle=0, vjust=1)) 
plot
#jpeg(paste0("/Users/xinwang/Google Drive/research reading/2016Spring/BIGDATA/",X_distribution,"p",p,"k",subsize,".jpg"),width = 480, height = 360, quality = 100, pointsize = 10)
#jpeg(paste0("/Users/xinwang/Google Drive/research reading/2016Spring/BIGDATA/",X_distribution,"p",p,"n",ni,".jpg"),width = 480, height = 360, quality = 100, pointsize = 10)
plot
dev.off()
plot = ggplot(EPE,aes(x = subdata_size,y=log(MSE,10),colours=Method,group=Method,linetype = Method))+  geom_point(aes(shape=Method))+geom_line(aes())+ggtitle(paste0("MSE Comparison, p = ",p," n = ",ni," x ~ ",substr(X_distribution,6,7))) +theme(text = element_text(size=14),axis.text.x = element_text(angle=0, vjust=1)) 
plot

time = data.table(time)
time3 = time[,.(total_time = lasso_time+select_time,Method,n_total)]
time1 = time[,.(select_time,Method,n_total)]
time2 = time[,.(lasso_time,Method,n_total)]
#write.csv(reshape(time1, idvar = "n_total", timevar = "Method", direction = "wide"),"/Users/xinwang/Google Drive/research reading/2016Spring/BIGDATA/temp.csv")

reshape(time1, idvar = "n_total", timevar = "Method", direction = "wide")
reshape(time2, idvar = "n_total", timevar = "Method", direction = "wide")
reshape(time3, idvar = "n_total", timevar = "Method", direction = "wide")

time3 = time[,.(total_time = lasso_time+select_time,Method,subdata_size)]
time1 = time[,.(select_time,Method,subdata_size)]
time2 = time[,.(lasso_time,Method,subdata_size)]
reshape(time1, idvar = "subdata_size", timevar = "Method", direction = "wide")
reshape(time2, idvar = "subdata_size", timevar = "Method", direction = "wide")
reshape(time3, idvar = "subdata_size", timevar = "Method", direction = "wide")

(sens_byn = reshape(SENS[,-3], idvar = "n_total", timevar = "Method", direction = "wide"))
(sens_byk = reshape(SENS[,-4], idvar = "subdata_size", timevar = "Method", direction = "wide"))
(spec_byn = reshape(SPEC[,-3], idvar = "n_total", timevar = "Method", direction = "wide"))
(spec_byk = reshape(SPEC[,-4], idvar = "subdata_size", timevar = "Method", direction = "wide"))
(pselect_byn = reshape(p_selected[,-3], idvar = "n_total", timevar = "Method", direction = "wide"))

EPE
A = reshape(EPE[,-3], idvar = "subdata_size", timevar = "Method", direction = "wide")
write.csv(A,file="MSE_sis_byk.csv")
B = reshape(EPE[,-4], idvar = "n_total", timevar = "Method", direction = "wide")
write.csv(B,file="MSE_sis_byn.csv")
(sens_byn = reshape(SENS[,-3], idvar = "n_total", timevar = "Method", direction = "wide"))
write.csv(sens_byn,file="SIS_sens_byn.csv")
(sens_byk = reshape(SENS[,-3], idvar = "subdata_size", timevar = "Method", direction = "wide"))
write.csv(sens_byk,file="SIS_sens_byk.csv")
total_time = reshape(time3, idvar = "subdata_size", timevar = "Method", direction = "wide")
write.csv(total_time,file="SIS_total_time_byk.csv")

total_time = reshape(time3, idvar = "n_total", timevar = "Method", direction = "wide")
write.csv(total_time,file="SIS_total_time_byn.csv")
(spec_byn = reshape(SPEC[,-3], idvar = "n_total", timevar = "Method", direction = "wide"))
write.csv(sens_byn,file="SIS_spec_byn.csv")
(spec_byk = reshape(SPEC[,-3], idvar = "subdata_size", timevar = "Method", direction = "wide"))
write.csv(sens_byn,file="SIS_spec_byk.csv")
# C = reshape(EPE[,-4], idvar = "n_total", timevar = "Method", direction = "wide")
# write.csv(C,file="MSE_p500_byn.csv")
# D = reshape(EPE[,-3], idvar = "subdata_size", timevar = "Method", direction = "wide")
# write.csv(D,file="MSE_p500_byk.csv")


require(tikzDevice)
library(xlsx)
lw <- 4
dat <- read.xlsx(file="effiency_varied_k.xlsx",sheetIndex = 1)
A = reshape(EPE[,-4], idvar = "subdata_size", timevar = "Method", direction = "wide")
dat = A
dat
kss <- c(1000, 2000, 3000, 4000, 5000)
rec.kc <- t(dat[,c(2:6)])
rec.kc <- t(dat[,c(2:5)])
rec.kc
rec.kc <- log10(rec.kc)
pp <- dim(rec.kc)[1]
name <- paste("KLASSO_SIS", ".tex", sep="")
ht <- 3.5
tikz(name, standAlone=TRUE, width=1*ht, height=0.65*ht) 
par(mar=c(2.6, 2.7, 0.6, 0.4))
plot(kss, rec.kc[1,],type="n",ylim=c(min(rec.kc),max(rec.kc)+0.3),lwd=lw)
adl <- c(1,4,5)
cl <- c(1,4,3); pl=c(1,5,3); ll=c(1,3,4);
adl <- c(1:length(methods_adj))
cl <-c(1:length(methods_adj))
pl=c(1:length(methods_adj))
ll=c(1:length(methods_adj))
for(pp in 1:length(adl))
  lines(kss, rec.kc[pp,], type="o", pch=pl[pp], col=cl[pp], lty=ll[pp], lwd=lw)
mtext("$k$", side=1, line=1.6)
mtext("log$_{10}$(Prediction Error)", side=2, line=1.8)
#legend("topright", legend=c("FULL", "COR-DOPT(250)",
#"LEV(250)",'ALEV(250)', "UNIF"), pch=pl, col=cl, lty=ll, lwd=lw, cex=1)

#legend("topright", legend=c("UNIF", "$D$-OPT",
#                           "FULL", "LEV"), pch=pl, col=cl, lty=ll, lwd=lw, cex=1)
legend("topright", legend=methods_adj, pch=pl, col=cl, lty=ll, lwd=lw, cex=1)
dev.off()
tools::texi2dvi(name, pdf=T)
file.remove(paste(gsub(".tex","",name),c(".aux",".log",".tex"),sep=""))


#rm(list=ls())
setwd("/Users/xinwang/Documents/newto_copy3")
require(tikzDevice)
lw <- 4
dat <- read.xlsx(file="effiency_varied_n.xlsx",sheetIndex = 1)
B = reshape(EPE[,-3], idvar = "n_total", timevar = "Method", direction = "wide")
dat = B
dat
nss <- c(10^4, 2*10^4, 4*10^4, 8*10^4)
rec.nc <- t(dat[1:4,c(2:5)])
rec.nc
rec.nc <- log10(rec.nc)
pp <- dim(rec.nc)[1]
name <- paste("NLASSO_SIS", ".tex", sep="")
ht <- 3.5
tikz(name, standAlone=TRUE, width=1*ht, height=0.65*ht) 
par(mar=c(2.6, 2.7, 0.6, 0.4))
plot(nss, rec.nc[1,],type="n",ylim=c(min(rec.nc),max(rec.nc)),lwd=lw)
adl <- c(1,4,5)
cl <- c(1,4,3); 
pl=c(1,5,3); 
ll=c(1,3,4);
adl <- c(1:length(methods_adj))
cl <-c(1:length(methods_adj))
pl=c(1:length(methods_adj))
ll=c(1:length(methods_adj))
for(pp in 1:length(adl))
  lines(nss, rec.nc[pp,], type="o", pch=pl[pp], col=cl[pp], lty=ll[pp], lwd=lw)
mtext("$n$", side=1, line=1.6)
mtext("log$_{10}$(Prediction Error)", side=2, line=1.8)
#legend("topright", legend=c("UNI", "$D$-OPT",
#  "FULL","LEV"), pch=pl, col=cl, lty=ll, lwd=lw, cex=1)
#legend("topright", legend=c("FULL", "COR-DOPT(250)",
#                           "LEV(250)",'ALEV(250)', "UNIF"), pch=pl, col=cl, lty=ll, lwd=lw, cex=1)
legend("topright", legend=methods_adj, pch=pl, col=cl, lty=ll, lwd=lw, cex=1)
dev.off()
tools::texi2dvi(name, pdf=T)
file.remove(paste(gsub(".tex","",name),c(".aux",".log",".tex"),sep=""))

