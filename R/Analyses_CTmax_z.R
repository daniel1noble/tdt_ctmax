# Analyses CTmax and z
data <- read.table("Dataset_CTmax_z.txt",header=TRUE)
filter <- c(which(data$class=="Insecta"),which(data$class=="Bivalvia"),which(data$class=="Actinopterygii"))
data <- data[filter,]
data$bg <- ifelse(data$class=="Insecta","black","white");data$bg <- ifelse(data$class=="Bivalvia","gray",data$bg);

quartz(NA,10,4)
par(oma=c(3,3,0,0),mar=c(2,2,1,0.5))
layout(matrix(c(1,2,3,4),1,4,byrow=TRUE),width=c(1,1,1,1.5))
plot(1,10000,xlim=c(1,1440),ylim=c(0,60),log="x",las=1,xaxt="n")
for (i in 1:nrow(data))
{points(seq(1,1440,100),data$Ctmax[i] - data$z[i]*log10(seq(1,1440,100)),type="l",col=ifelse(data$class[i]=="Insecta","black",0),lty=ifelse(data$LT[i]==50,1,2),lwd=0.7)}
axis(1,at=c(1,10,60,180,720),c("1 min","10 min","1 h","3 h","12 h"))
text(1.5,5,"Insects",adj=c(0,0),cex=1.3)
text(1.5,1,"28 TDT (22 spp)",adj=c(0,0),cex=1.1)
text(900,51,substitute("LT"[100]))
text(180,27,substitute("LT"[50]))

plot(1,10000,xlim=c(1,1440),ylim=c(0,60),log="x",las=1,xaxt="n")
for (i in 1:nrow(data))
{points(seq(1,1440,100),data$Ctmax[i] - data$z[i]*log10(seq(1,1440,100)),type="l",col=ifelse(data$class[i]=="Bivalvia","gray",0),lwd=1)}
axis(1,at=c(1,10,60,180,720),c("1 min","10 min","1 h","3 h","12 h"))
text(1.5,5,"Bivalves",adj=c(0,0),cex=1.3)
text(1.5,1,"14 TDT (13 spp)",adj=c(0,0),cex=1.1)
text(180,41,substitute("LT"[50]))

plot(1,10000,xlim=c(1,1440),ylim=c(0,60),log="x",las=1,xaxt="n")
for (i in 1:nrow(data))
{points(seq(1,1440,100),data$Ctmax[i] - data$z[i]*log10(seq(1,1440,100)),type="l",col=ifelse(data$class[i]=="Actinopterygii","black",0))}
axis(1,at=c(1,10,60,180,720),c("1 min","10 min","1 h","3 h","12 h"))
text(1.5,5,"Fishes",adj=c(0,0),cex=1.3)
text(1.5,1,"14 TDT (14 spp)",adj=c(0,0),cex=1.1)
text(180,41,substitute("LT"[50]))
mtext("Temperature (ºC)",2,line=1,outer=TRUE,cex=0.9)
mtext("Time (log-scale)",1,line=1.2,outer=TRUE,adj=0.32,cex=0.9)


par(mar=c(2,5,4,1))
coef <- coef(lm(data$Ctmax ~ data$z + as.factor(data$LT),subset=data$Ctmax > 20))
data$Ctmax.trans <- ifelse(data$LT==50,data$Ctmax,data$Ctmax - coef[3])
plot(data$z,data$Ctmax.trans,bg=data$bg,cex=2,pch=22,col=ifelse(data$class=="Insecta","white","black"),ylim=c(10,65),lwd=0.9,las=1,ylab=substitute("CT"[max]* " (ºC)"),xlab=expression(paste(italic(z)," (ºC)")), cex.lab=1.4)
abline(coef[1],coef[2],lwd=0.7,lty=2)
axis(1,1:10,lab=NA)
mtext(expression(paste(italic(z)," (ºC)")),1,line=1.2,outer=TRUE,adj=0.87,cex=0.9)




# Analyses CTmin and z
data <- read.table("Dataset_CTmin_z.txt",header=TRUE)

quartz(NA,10,4)
par(oma=c(3,3,0,0),mar=c(2,2,1,0.5))
layout(matrix(c(1,2,3,4),1,4,byrow=TRUE),width=c(1,1.5,1,1))
plot(1,10000,xlim=c(1,1440),ylim=c(-120,20),log="x",las=1,xaxt="n")
for (i in 1:nrow(data))
{points(seq(1,1440,100),data$Ctmin[i] + data$z[i]*log10(seq(1,1440,100)),type="l",col=ifelse(data$class[i]=="Insecta","black",0),lty=ifelse(data$LT[i]==50,1,2),lwd=0.7)}
axis(1,at=c(1,10,60,180,720),c("1 min","10 min","1 h","3 h","12 h"))
text(1.5,15,"Insects",adj=c(0,0),cex=1.3)
text(1.5,6,"22 TDT (13 spp)",adj=c(0,0),cex=1.1)
text(3,-110,substitute("LT"[100]))
text(2,-12,substitute("LT"[50]))
mtext("Temperature (ºC)",2,line=1,outer=TRUE,cex=0.9)
mtext("Time (log-scale)",1,line=1.2,outer=TRUE,adj=0.07,cex=0.9)

# Plot CTmin versus z' 
par(mar=c(2,5,4,1))
coef <- coef(lm(data$Ctmin ~ data$z + as.factor(data$LT)))
data$Ctmin.trans <- ifelse(data$LT==50,data$Ctmin,data$Ctmin - coef[3])
plot(data$z,data$Ctmin.trans,cex=2,pch=22,bg="black",col="white",ylim=c(-100,-20),xlim=c(5,25),lwd=0.9,las=1,ylab=substitute("CT"[min]* " (ºC)"),xlab=expression(paste(italic(z),"' (ºC)")), cex.lab=1.4)
abline(coef[1],coef[2],lwd=0.7,lty=2)
mtext(expression(paste(italic(z),"' (ºC)")),1,line=1.2,outer=TRUE,adj=0.41,cex=0.9)



# -----------------------------------------------------------------------------------

# New figure 5 (suggested by the referee)

data.sub <- data[c(29,31,40),]
data.sub <- data.sub[c(2,1,3),]
par(mar=c(4,4.2,1,8))

layout(matrix(c(1,2,3),3,1,byrow=TRUE))
plot(data.sub$z,data.sub$Ctmax,las=1,ylab=substitute("CT"[max]* " (ºC)"),xlab=expression(paste(italic(z)," (ºC)")),xlim=c(1,7),ylim=c(39,51))
abline(lm(data.sub$Ctmax ~ data.sub$z),lty=2)
points(data.sub$z,data.sub$Ctmax,pch=22,bg="gray",cex=1.5)
segments(4.3,40.5,3.1,43,lwd=0.7)
text(data.sub$z-1,data.sub$Ctmax+0.7,expression(italic("S.corr"),italic("S.nov"),italic("C.gla")),cex=0.9,adj=c(0,0))
text(3.5,39,expression(paste("slope ",italic("b")," = 2.02")),adj=c(0,0))



par(mar=c(4,4,1,1))
plot(1,10000,xlim=c(1,1440),ylim=c(33,52),log="x",las=1,xaxt="n",xlab="Time (log-scale)",ylab="Temperature (ºC)")
for (i in 1:3)
{points(seq(1,2000,100),data.sub$Ctmax[i] - data.sub$z[i]*log10(seq(1,2000,100)),type="l",col="gray",lty=ifelse(data.sub$LT[i]==50,1,2),lwd=1.5)}
axis(1,at=c(1,10,60,180,720),c("1 min","10 min","1 h","3 h","12 h"))
text(100,45,expression(paste(italic("t"),"' = 10"^italic("b")*"x 1 min")),adj=c(0,0))
text(150,39,"105 min",adj=c(0,0),cex=0.9)
segments(180,43.5,105,37.5,lwd=0.7)
#text(250,41.5,"= 105 min",adj=c(0,0))
text(0.9,data.sub$Ctmax+0.5,expression(italic("S.corr"),italic("S.nov"),italic("C.gla")),cex=0.9,adj=c(0,0))


cline <- rbind(data.sub$Ctmax - data.sub$z*log10(10),data.sub$Ctmax - data.sub$z*log10(104.2),data.sub$Ctmax - data.sub$z*log10(720))
plot(1:3,cline[1,],las=1,xlim=c(0.5,3.2),ylim=c(30,46),pch=22,bg="gray",cex=1.5,xaxt="n",type="b",lty=1,ylab="Inferred tolerance (ºC)",xlab="Environmental gradient")
points(1:3,cline[2,],pch=22,bg="gray",cex=1.5,type="b",lty=1)
points(1:3,cline[3,],pch=22,bg="gray",cex=1.5,type="b",lty=1)
axis(1,at=c(1,2,3),lab=expression(italic("S.corr"),italic("S.nov"),italic("C.gla")))
text(rep(0.8,3),cline[,1]+1,c("10 min","105 min","12 h"),adj=c(0,0),cex=0.9)
arrows(2.7,23.7,3.3,23.7,xpd=TRUE,length=0.1)


# -----------------------------------------------------------------------------------
# Data obtained to calculate CTmax and z

# Total number of measurements estimated by Enrico on June 2019 for scaling study with Nacho:
n <- c(4,4,5,5,3,7,7,4,4,4,4,4,4,7,7,7,4,6,5,4,4,4,5,6,5,8,4,3,7,7,7,3,5,6,6,9,3,6,3,3,4,9,4,3,5,5,5,6,5,5,4,3,5,7,7,7,4,3,3,6,6,3)
 

# ArmstrongEA2009_Table_5
t1 <- c(120,30,3,1.5)
t2 <- c(80,25,6,1)
t3 <- c(50,12,3,0.5)
t4 <- c(70,20,4,1)
T <- c(44,46,48,50)

# ArmstrongEA2009_Table_6t1 <- c(60,20,5,1.5)t2 <- c(130,25,4,1.5)
t3 <- c(70,15,2.5,1.0)t4 <- c(55,15,5,1.5)
t5 <- c(150,40,5,2)
T <- c(44,46,48,50)

# AbdelghanyEA2010_Table3
t90 <- c(25,20,3.86,0.18,0.08)*60
t50 <- c(13,10,1.4,0.06,0.05)*60
T <- c(42,45,50,55,60)

# Cerda&Retana2000_Fig3
t1 <- c(57,93,599,1533)/60,#C_rosenhaueri small
T1 <- c(54,52,50,48)
t2 <- c(114,226,597,1704,2691,3604)/60,#C_velox large
T2 <- c(58,56,54,52,50,48)

# Boina&SubramanyamEA2004_Table2
t <- c(213,154,72,57.5,34,14.7)
T <- c(46,48,50,54,58,60)

#HuangEA2006_Table1
t <- c(300,180,60,30,14,8)
T <- c(41,42,43,44,45,46)

#FederEA1997_Fig7
t <- c(31,15,8)
T <- c(39,40,41)
 
#Brett1956 Fig 2
t1 <- c(23.2,72.6,143.2,271.6,853.1),T1 <- c(34.5,34,33.5,33,32.5),
t2 <- c(27.0,54.2,143.5,226.5,547.0),T2 <- c(33,32.5,32,31.5,30.5),t3 <- c(17.9,31.7,170.6,404.6,679.2,1233.1)T3 <- c(33,32.5,31.5,31,31,30.5)t4 <- c(27.1,42.2,106.7,332.7),,T4 <- c(30,29.5,29,28),,t5 <- c(32.7,55.0,86.9,266.7,462.4,1995.3)T5 <- c(29,28.5,28,27,26.5,25.5)t6 <- c(36.6,146.2,226.5,292.4,1425.6,3715.4)T6 <- c(28,27.5,27,26.5,26,25.5)t7 <- c(42.4,92.7,164.8,1002.3),,T7 <- c(27,26.5,26,25),, 
 
#CooperEA2008_Fig2
t <- c(27.5,19.2,8.5,10.2,10.2,5.7,1.4,1.8,1.6)T <- c(38,38.5,39,39.5,40,40.5,41,41.5,42) 

#QuinnEA1994_Table1
t <- c(24,48,96)*60
 T1 <- c(26.8,24.5,22.6)
 T2 <- c(26.9,25.3,23.6)
 T3 <- c(27.8,27,25.9)
 T4 <- c(30.4,26.8,25)
 T5 <- c(32.8,31.8,30.5)
 T6 <- c(27.4,26.5,25.7)
 T7 <- c(27.5,26.3,24.1)
 T8 <- c(30.1,28.9,26.7)
  
#Somero&DeVries1967_Table1
T1 <- c(15,10,8,5)
t1 <- c(6,140,430,12960)
T2 <- c(15,10,5)
t2 <- c(8,60,15840)
T3 <- c(15,10,7)
t3 <- c(7,81,1095) 
 
 
#Cox&Rutherford2000_Fig2
T1 <- c(26.708,25.366,24.537,24.252)
t1 <- c(24,48,72,96)*60

T2 <- c(31.718,31.091,31.004,30.992)
t2 <- c(24,48,72,96)*60 
 
 
#SelongEA2011_Fig1 
T <- c(21,22,23,24,26) 
t <- c(56.501,24.614,10.717,4.735,0.967)*60*24 
 
#Urban1994_Fig2
T1 <- c(36,34,32,30,29)
t1 <- c(1.581,2.291,3.422,11.688,23.748)*60
T2 <- c(36,34,32,29,27)
t2 <- c(2.424,6.268,17.266,35.855,95.832)*60
T3 <- c(36,32,30,29,27,25)
t3 <- c(1.576,11.967,15.164,17.744,24.794,96.307)*60
T4 <- c(36,32,30,29,27,25,23)
t4 <- c(0.376,2.350,8.199,11.568,23.100,48.216,95.895)*60
 
#Doudoroff1945_Fig2
T1 <- c(39.7,39.3,38.5,37.7,37.5,37,36.5)
t1 <- 60*(10^c(-0.311,,-0.008,0.472,0.769,1.074,1.380,1.678)),
T2 <- c(35,34.5,33.5,32.7,32.3,31.5)
t2 <- 60*(10^(c(-0.308,-3.836e-4,0.474,0.775,1.075,1.380))),
T3 <- c(33.4,33.3,32.5,32.4,32.2,32.0,31.8,31.7,31.6)
t3 <- 60*(10^(c(-0.305,-0.002,0.471,0.776,1.075,1.382,1.682,1.858,1.980)))

#AnsellEA1980_Fig1A,1B,1C_Donax
t <- c(3,6,12,24,48,72,96)*60T1 <- c(32.203,30.973,29.578,28.56,28.575,28.468,28.069)T2 <- c(33.701,31.775,30.879,30.306,30.268,30.227,29.386)T3 <- c(36.961,36.212,34.831,33.067,32.545,32.305,32.261)

#AnsellEA1980a_Fig1A,1B,1C_Tellina
t <- c(3,6,12,24,48,72,96)*60T1 <- c(28.865,28.746,27.241,26.724,26.583,26.474,26.356)T2 <- c(34.442,32.895,32.1,30.855,30.845,30.849,30.788)
T3 <- c(35.437,34.732,32.689,32.238,31.374,31.156,30.721)

#AnsellEA1981_Fig1A,B,C_Cardium
t <- c(3,6,12,24,48,72,96)*60
T1 <- c(36.786,36.061,35.085,34.646,34.573,34.423,33.454) 
T2 <- c(34.455,32.921,31.447,31.021,30.851,30.549,29.65)
T3 <- c(36.919,36.609,34.726,33.527,31.163,30.519,29.519) 
 
#Trogoderma granarium larvae Zacher1927 in Strang 1992
T <- c(50.5,52,53,54,55,58)t <- c(300,90,30,20,10,5)  
 
#Tribolium confusum Boina&Subranyam2004 Table 2
T <- c(46,48,50,54,58,60)
t <- c(213.5,153.7,72.2,57.5,34.1,14.7)

 
 
# ----------------------------------------- 
# ----------------------------------------- 
# ----------------------------------------- 

# Data obtained to calculate CTmin and z'

#Fields1992_Fig1
Tcast <- c(-14.105,-4.352,2.315,5.813,6.738,7.01) ,,,#Tribolium castaneum
tcast <- c(12.202,30.476,63.096,122.229,184.579,268.42)*60*24
Tferr <- c(-20.275,-15.734,-3.35,3.729,9.97,13.875) ,,#Cryptolestes ferrugineus tferr <- c(6.372,7.798,40.912,95.141,164.582,242.697)*60*24 
Tgran <- c(-9.168,-5.564,3.59,6.867,8.244,9.8) ,,,,#Sitophilus granariustgran <- c(1.055,3.098,23.426,76.705,157.743,243.335)*60*24 
 
#RenaultEA2004_Fig2
T <- c(-6.989,-4.992,-1.955,-0.002,2.005,5.025,10.033,12.464) #Alphitobius diaperinus_adultt <- c(0.215,0.808,1.408,3.941,3.892,14.649,28.049,28.699)*60*24

#Imai&Harada_Table1_Lasioderma serricorne
T <- c(-15,-10,-5,0,5)
tegg <- c(0.22,1.29,4.92,39.07,137.6)*60
tlarv <- c(0.98,5.2,35.29,105.9,112.9)*60
tpup <- c(1.16,2.44,3,153.5,365.1)*60
tad <- c(0.64,1.21,25.84,126.8,457.6)*60
 
#AbdelghanyEA2010_Table4
#Stegobium paniceum adult 
T <- c(5,0,-5,-10,-15)
t <- c(360,112,25,0.8,0.3)*60
 
#LoganathanEA2011_Tables4_5
#Callosobruchus maculatus
T <- c(0,-5,-10,-15)
tpup <- c(274,122,7,1.8)*60
tegg <- c(44,19,6.4,1)*60

#Oryzaephilus surinamensis adult Mathlein1961_in_Strang1992
T <- c(-2,-5,-7)
t <- c(36000,28800,21600)


#Sitophilus granarius adult Back and Cotton, 1924 in Strang1992
T <- c(-18,-15,-8,-5.5,-2.5,0.5,3)t <- c(300,450,20160,47520,66240,105120,159840)

#Sitophilus granarius Mathlein 1961 in Strang1992
T <- c(4,-1,-4,-6) ,,,,,#eggt <- c(84960,57600,21600,14400)
T1 <- c(-1,-4,-6),,,,,#larvaet1 <- c(72000,43200,21600)

#Sitophilus oryzae adult Back and Cotton1924 in Strang1992
T <- c(-18,-15,-8,-5.5,-2.5,0.5,3,5.5)t <- c(240,270,4320,8640,11520,23040,25920,115200)

#Tribolium castaneum all stages Cotton1950 in Strang1992
T <- c(-1,-4,-7,-10)t <- c(24480,11520,7200,1440)


#Tribolium confusum all stages Cotton1950 in Strang1992
T <- c(-1,-4,-7,-10)t <- c(24480,17280,7200,1440)


#Tineola bisselliella Back&Cotton1927 in Strang1992
T <- c(-16.5,-13.5,-10.5,-5.5,-2.5),#eggt <- c(1440,2880,5760,30240,30240)T1 <- c(-16.5,-13.5,-5.5,-2.5),#larvaet1 <- c(2880,30240,96480,201600),

#Anagasta kuhniella all stages Cotton1950 in Strang1992
T <- c(-4,-7,-10,-12,-15,-18)t <- c(167040,34560,10080,5760,4320,1440)

#Plodia interpunctuella all stages Cotton1950 in Strang1992
T <- c(-4,-7,-9,-12)t <- c(129600,40320,11520,7200)


# ----------------------------------------- 
# ----------------------------------------- 
# ----------------------------------------- 

# Acclimation FryEA1946 Table 3
t3 <- c(39,303,550)
T3 <- c(26,24,23.5)
t11 <- c(24,78,96,136,217,620) 
T11 <- c(28,27.1,26.6,26.3,26,25)
t15 <- c(28,40,51,96,108,151,196,437,859)
T15 <- c(28,28,28,27.5,27,27,26.5,26,25.5)
t20 <- c(34,57,89,136,277,481,2000)
T20 <- c(29,28.5,28,27.6,27.1,26.5,25.5)
t22 <- c(47,77,130,269,519,1010)
T22 <- c(29,28.5,28,27.5,27,26.5)
t24 <- c(15,36,49,86,197,252,745,750,2340,5040)
T24 <- c(30,29.5,29,28.5,28,27.7,27,26.5,26,25.5)
t25 <- c(51,194,2256)
T25 <- c(29,28,26)


reg <- lm(log10(t) ~ T)
ctmax <- reg$coef[1]/reg$coef[2]
ctmax
z <- 1/reg$coef[2]
z
summary(reg)
plot(log10(t) ~ T)
