library(mnonr);library(moments);library(sampling);library(readxl);library(caTools);library(MASS);
set.seed(11)
#############Symmetric Data#######
sigma<-rbind(c(863.7, 59.35, 60.47), c(59.35,  5.91,  3.95), c(60.47,  3.95,  5.82)) # create the mean vector
mu<-c(178, 37,38)
df<-as.data.frame(mvrnorm(n=1000, mu=mu, Sigma=sigma))
 #########Asymmetric Data##########
# sigma<-rbind(c(863.7, 59.35, 60.47), c(59.35,  5.91,  3.95), c(60.47,  3.95,  5.82))
# df=as.data.frame(unonr(1000, c(178, 37,38), sigma, skewness =c(1.19, 1.54, 1.51 ), kurtosis = c(6.5,5.84,7.1),empirical=F))
Y<-df[,1];X<-df[,2];Z<-df[,3];
My<-mean(Y);Mx<-mean(X);Mz<-mean(Z);
N1<-length(Y);
m=300;sr=0.80;r=m*sr;   # Simulation for 20% missing data
fr=(1/r)-(1/N1);fn=(1/m)-(1/N1);frn=(1/r)-(1/m);Sy<-sqrt(var(Y));Sx<-sqrt(var(X));Sz<-sqrt(var(Z));
Cy<-sqrt(var(Y))/My;Cx<-sqrt(var(X))/Mx;Cz<-sqrt(var(Z))/Mz;rho_yx<-cor(Y,X);rho_yz<-cor(Y,Z);rho_xz<-cor(X,Z);
A<-fr*Cy^2;B<-fr*Cx^2;C<-fn*Cx^2;D<-fr*Cz^2;E<-fn*Cz^2;f<-fr*rho_yx*Cy*Cx;G<-fn*rho_yx*Cy*Cx;H<-fr*rho_yz*Cy*Cz;I<-fn*rho_yz*Cy*Cz;J<-fn*Cx^2;K<-fr*rho_xz*Cx*Cz;L<-fn*rho_xz*Cx*Cz;M<-fn*rho_xz*Cx*Cz;N<-fn*rho_xz*Cx*Cz;O<-fn*Cz^2;
Pdelta1_1<-2*(f*E-I*L)/(B*E-L^2);Pdelta1_2<-2*(f*L-I*B)/(B*E-L^2);Pdelta2_1<-2*(f*D-K*H)/(B*D-K^2);Pdelta2_2<-2*(f*K-H*B)/(B*D-K^2);Pdelta3_1<-2*(G*E-N*I)/(C*E-N^2);Pdelta3_2<-2*(G*N-C*I)/(C*E-N^2);
kyz=(rho_yz*(Cy/Cz));syz=(rho_yz*(Sy/Sz));k=(rho_yx*(Cy/Cx));Syx=(rho_yx*(Sy/Sx));alpha_sh<-1-k;
a1<-1/(1+Cy^2*(frn+fn*(1-rho_yx^2))) ;b1<-Syx*a1;a2<-1/(1+Cy^2*(frn*(1-rho_yx^2)));b2<-Syx*a2;a3<-1/(1+Cy^2*(fn+frn*(1-rho_yx^2)));b3<-Syx*a3;b4=b5=b6=k;
a4=(1+((b4^2)/2)*fn*(Cx^2)+(b4/2)*fn*(Cx^2-2*rho_yx*Cy*Cx))/(1+fr*Cy^2+2*b4^2*fn*Cx^2+b4*fn*(Cx^2-4*rho_yx*Cy*Cx));a5=(1+((b5^2)/2)*fr*(Cx^2)+(b5/2)*fr*(Cx^2-2*rho_yx*Cy*Cx))/(1+fr*Cy^2+2*b5^2*fr*Cx^2+b5*fr*(Cx^2-4*rho_yx*Cy*Cx));
a6=(1+((b6^2)/2)*frn*(Cx^2)+(b6/2)*frn*(Cx^2-4*rho_yx*Cy*Cx))/(1+fr*Cy^2+2*b6^2*frn*Cx^2+b6*frn*(Cx^2-4*rho_yx*Cy*Cx));b7=b8=b9=-k;a7=(1+(b7*(b7-2)*fn*Cx^2)/2+(b7*fn*rho_yx*Cx*Cy))/(1+(fr*Cy^2+2*b7*(b7-1)*fn*Cx^2)+(4*b7*fn*rho_yx*Cx*Cy));
a8=(1+(b8*(b8-2)*fr*Cx^2)/2+(b8*fr*rho_yx*Cx*Cy))/(1+(fr*Cy^2+2*b8*(b8-1)*fr*Cx^2)+(4*b8*fr*rho_yx*Cx*Cy));a9=(1+(b9*(b9-2)*frn*Cx^2)/2+(b9*frn*rho_yx*Cx*Cy))/(1+(fr*Cy^2+2*b9*(b9-1)*frn*Cx^2)+(4*b9*frn*rho_yx*Cx*Cy));a10<-Syx; b10=(rho_yx*(Sy/Sz));
a11<-k; b11<-b12<-kyz; a12<-1-k;a14<-k;b14<-kyz;a15<-1-k;b15<-kyz;gamma1<-1/(1+Cy^2*(fr-frn*rho_yx^2-fn*rho_yz^2));gamma2=(1+((a14^2)/2)*frn*Cx^2+((b14^2)/2)*fn*Cz^2+(a14/2)*frn*(Cx^2-2*rho_yx*Cx*Cy)+(b14/2)*fn*(Cz^2-2*rho_yz*Cz*Cy))/(1+fr*Cy^2+2*a14^2*frn*Cx^2+2*b14^2*fn*Cz^2+a14*frn*(Cx^2-4*rho_yx*Cy*Cx)+b14*fn*(Cz^2-4*rho_yz*Cy*Cz));
gamma3=(1+((a15^2)*frn*Cx^2)+((b15^2)*fn*Cz^2)-(a15*frn*(2*Cx^2-2*rho_yx*Cx*Cy))+frn*(Cx^2-2*rho_yx*Cx*Cy)-b15*fn*rho_yz*Cy*Cz);a13=-gamma1*Syx;b13=-gamma1*syz;c1<-1-2*k ; c2<-1-2*((Mx+1)/Mx)*k  ;            
Yn<-numeric(rp);Xn<-numeric(rp);Zn<-numeric(rp);Ynr<-numeric(rp);Xnr<-numeric(rp);Znr<-numeric(rp);Myn<-numeric(rp);Mxn<-numeric(rp);Mzn<-numeric(rp);Mynr<-numeric(rp);Mxnr<-numeric(rp);Mznr<-numeric(rp);YR<-numeric(rp);XR<-numeric(rp);ZR<-numeric(rp);YNR<-numeric(rp)
XNR<-numeric(rp);ZNR<-numeric(rp);yr.icom<-numeric(rp);yr.idp1<-numeric(rp);yr.idp2<-numeric(rp);yr.idp3<-numeric(rp);yr.ibp1<-numeric(rp);yr.ibp2<-numeric(rp);yr.ibp3<-numeric(rp);yr.ibp4<-numeric(rp);yr.ibp5<-numeric(rp)
yr.ibp6<-numeric(rp);yr.ibp10<-numeric(rp);yr.ibp13<-numeric(rp);ynr.im<-numeric(rp);ynr.ir<-numeric(rp);ynr.icom<-numeric(rp);ynr.isd<-numeric(rp);ynr.iss<-numeric(rp);ynr.isd<-numeric(rp);ynr.isf1r<-numeric(rp);
ynr.isf1p<-numeric(rp);ynr.isf1v<-numeric(rp);ynr.isf2r<-numeric(rp);ynr.isf2p<-numeric(rp);ynr.isf2v<-numeric(rp);ynr.idp1<-numeric(rp);ynr.idp2<-numeric(rp);ynr.idp3<-numeric(rp);ynr.ibp1<-numeric(rp);
ynr.ibp2<-numeric(rp);ynr.ibp3<-numeric(rp);ynr.ibp4<-numeric(rp);ynr.ibp5<-numeric(rp);ynr.ibp6<-numeric(rp);ynr.ibp7<-numeric(rp);ynr.ibp8<-numeric(rp);ynr.ibp9<-numeric(rp);ynr.ibp10<-numeric(rp);
ynr.ibp11<-numeric(rp);ynr.ibp12<-numeric(rp);ynr.ibp13<-numeric(rp);ynr.ibp14<-numeric(rp);ynr.ibp15<-numeric(rp);ynr.isps1<-numeric(rp);ynr.isps2<-numeric(rp);ynr.isps3<-numeric(rp);ynr.isps4<-numeric(rp);
ynr.isps5<-numeric(rp);ynr.isps6<-numeric(rp);ynr.isps7<-numeric(rp);ynr.isps8<-numeric(rp);ynr.isps9<-numeric(rp);tnr.iprop1<-numeric(rp);tnr.iprop2<-numeric(rp);tnr.iprop3<-numeric(rp);y.im<-numeric(rp);
y.ir<-numeric(rp);y.icom<-numeric(rp);y.isd<-numeric(rp);y.iss<-numeric(rp);y.isf1r<-numeric(rp);y.isf1p<-numeric(rp);y.isf1v<-numeric(rp);y.isf2r<-numeric(rp);y.isf2p<-numeric(rp);y.isf2v<-numeric(rp);y.idp1<-numeric(rp);
y.idp2<-numeric(rp);y.idp3<-numeric(rp);y.ibp1<-numeric(rp);y.ibp2<-numeric(rp);y.ibp3<-numeric(rp);y.ibp4<-numeric(rp);y.ibp5<-numeric(rp);y.ibp6<-numeric(rp);y.ibp7<-numeric(rp);y.ibp8<-numeric(rp);y.ibp9<-numeric(rp);
y.ibp10<-numeric(rp);y.ibp11<-numeric(rp);y.ibp12<-numeric(rp);y.ibp13<-numeric(rp);y.ibp14<-numeric(rp);y.ibp115<-numeric(rp);y.isps1<-numeric(rp);y.isps2<-numeric(rp);y.isps3<-numeric(rp);y.isps4<-numeric(rp);
y.isps5<-numeric(rp);y.isps6<-numeric(rp);m.im<-numeric(rp) ;m.ir<-numeric(rp) ;m.icom<-numeric(rp) ;m.isd<-numeric(rp);m.iss<-numeric(rp) ;m.isf1r<-numeric(rp);m.isf1p<-numeric(rp);m.isf1v<-numeric(rp);m.isf2r<-numeric(rp);
m.isf2p<-numeric(rp);m.isf2v<-numeric(rp);m.idp1<-numeric(rp);m.idp2<-numeric(rp);m.idp3<-numeric(rp);m.ibp1<-numeric(rp);m.ibp2<-numeric(rp);m.ibp3<-numeric(rp);m.ibp4<-numeric(rp);m.ibp5<-numeric(rp);m.ibp6<-numeric(rp);
m.ibp7<-numeric(rp);m.ibp8<-numeric(rp);m.ibp9<-numeric(rp);m.ibp10<-numeric(rp);m.ibp11<-numeric(rp);m.ibp12<-numeric(rp);m.ibp13<-numeric(rp);m.ibp14<-numeric(rp);m.ibp15<-numeric(rp);m.isps1<-numeric(rp);
m.isps2<-numeric(rp);m.isps3<-numeric(rp);m.isps4<-numeric(rp);m.isps5<-numeric(rp);m.isps6<-numeric(rp);m.multi<-numeric(rp);m.iprop1<-numeric(rp) ;m.iprop2<-numeric(rp) ;m.iprop3<-numeric(rp) ;
rp=50000
for(i in 1:rp){
  smp1<-c(sample(1:N1,m,replace=F));
  mar1<-df[smp1,];
  Yn<-mar1[,1];
  Xn<-mar1[,2];
  Zn=mar1[,3];
  Myn[i]<-mean(Yn);
  Mxn[i]<-mean(Xn);
  Mzn[i]<-mean(Zn);
  d_f <- sample.split(Yn , SplitRatio = sr);
  R <- mar1[d_f,] ;
  NR<- mar1[!d_f,]  ;
  YR<-R[,1];
  XR<-R[,2];
  ZR<-R[,3];
  b<-cov(YR,XR)/var(XR);
  r<-length(YR);
  for(j in 1:length(YR)){
    yr.icom[j]<-((m*alpha_sh*YR[j])/r)+((1-alpha_sh)*(sum(YR)/sum(XR))*XR[j])
    yr.idp1[j]<- (m*YR[j]/r)+b*(Mx-XR[j])  # Diana and Perri
    yr.idp2[j]<- (m*YR[j]/r)-b*(m*XR[j]/r)
    yr.idp3[j]<- (m*YR[j]/r)-b*(m*XR[j]/r)
    yr.ibp1[j]<- a1*YR[j] #bp 2016
    yr.ibp2[j]<- a2*YR[j]
    yr.ibp3[j]<- a3*YR[j]
    yr.ibp4[j]<- a4*YR[j] #bp 2018
    yr.ibp5[j]<- a5*YR[j]
    yr.ibp6[j]<- a6*YR[j]
    yr.ibp10[j]<-(m*YR[j])/r-(a10*m*XR[j])/r+(b10*m*Mz)/r
    yr.ibp13[j]<-(m*XR[j]-a13*m*XR[j]+b13*m*Mz)/r
  }
  YNR<-NR[,1];
  XNR<-NR[,2];
  ZNR<-NR[,3];
    for(j in 1:length(YNR)){
    ynr.im[j]<-mean(YR)
    ynr.ir[j]<- (sum(YR)/sum(XR))*XNR[j]
    ynr.icom[j]<-(1-alpha_sh)*(sum(YR)/sum(XR))*XNR[j]
    ynr.isd[j]<-mean(YR)*(m*(Mxn[i]/mean(XR))^k-r)*(XNR[j]/sum(XNR))
    ynr.iss[j]<-mean(YR)*(((m-r)*Mxn[i]+k*r*(Mxn[i]-mean(XR)))/(k*mean(XR)+(1-k)*Mxn[i]))*(XNR[j]/sum(XNR))
    #Singh et al. (2010)
    ynr.isf1r[j]<-(mean(YR)/(m-r))*(m*(Mx/mean(XR))-r) 
    ynr.isf1p[j]<-(mean(YR)/(m-r))*(m*(mean(XR)/Mx)-r) 
    ynr.isf1v[j]<-(mean(YR)/(m-r))*(m*((N1*Mx-m*mean(XR))/((N1-m)*Mx))-r)
    ynr.isf2r[j]<-(mean(YR)/(m-r))*(m*(Mx/Mxn[i])-r) 
    ynr.isf2p[j]<-(mean(YR)/(m-r))*(m*(Mxn[i]/Mx)-r) 
    ynr.isf2v[j]<-(mean(YR)/(m-r))*(m*((N1*Mx-m*Mxn[i])/((N1-m)*Mx))-r)
    #diana and perri
    ynr.idp1[j]<-b*(Mx-XR[j])
    ynr.idp2[j]<-b*(m*Mx/(m-r))
    ynr.idp3[j]<-b*(m*Mxn[i]/(m-r))
    #bp 2016
    ynr.ibp1[j]<-a1*mean(YR)+((m*b1)/(m-r))*(Mx-Mxn[i])
    ynr.ibp2[j]<-a2*mean(YR)+((m*b2)/(m-r))*(Mx-mean(XR))
    ynr.ibp3[j]<-a3*mean(YR)+b3*(XNR[j]-mean(XR))
    #BP 2018
    ynr.ibp4[j]<-(a4/(m-r))*(m*mean(YR)*(Mx/Mxn[i])^b4-r*mean(YR))
    ynr.ibp5[j]<-(a5/(m-r))*(m*mean(YR)*(Mx/mean(XR))^b5-r*mean(YR))
    ynr.ibp6[j]<-(a6/(m-r))*(m*mean(YR)*(Mxn[i]/mean(XR))^b6-r*mean(YR))
    #Singh et al. (2021)
    ynr.isps1[j]<-m*((c1*mean(YR)+(1-c1)*mean(YR)*exp((Mx-mean(XR))/(Mx+mean(XR))))-(r/m)*mean(YR))*(XNR[j]/sum(XNR))            
    ynr.isps2[j]<-m*((c2*mean(YR)+(1-c2)*mean(YR)*exp((Mx-mean(XR))/(Mx+mean(XR)+2)))-(r/m)*mean(YR))*(XNR[j]/sum(XNR))          
    ynr.isps3[j]<-m*((c1*mean(YR)+(1-c1)*mean(YR)*exp((Mx-Mxn[i])/(Mx+Mxn[i])))-(r/m)*mean(YR))*(XNR[j]/sum(XNR))                
    ynr.isps4[j]<-m*((c2*mean(YR)+(1-c2)*mean(YR)*exp((Mx-Mxn[i])/(Mx+Mxn[i]+2)))-(r/m)*mean(YR))*(XNR[j]/sum(XNR))              
    ynr.isps5[j]<-m*((c1*mean(YR)+(1-c1)*mean(YR)*exp((Mxn[i]-mean(XR))/(Mxn[i]+mean(XR))))-(r/m)*mean(YR))*(XNR[j]/sum(XNR))    
    ynr.isps6[j]<-m*((c2*mean(YR)+(1-c2)*mean(YR)*exp((Mxn[i]-mean(XR))/(Mxn[i]+mean(XR)+2)))-(r/m)*mean(YR))*(XNR[j]/sum(XNR))  
    #bhushan an pandey (2023)
    ynr.ibp7[j]<-(1/(m-r))*(m*a7*mean(YR)*(1+log(Mxn[i]/Mx))^b7-r*mean(YR))
    ynr.ibp8[j]<-(1/(m-r))*(m*a8*mean(YR)*(1+log(mean(XR)/Mx))^b8-r*mean(YR))
    ynr.ibp9[j]<-(1/(m-r))*(m*a9*mean(YR)*(1+log(mean(XR)/Mxn[i]))^b9-r*mean(YR))
    # Bhushan and pandey (2023)
    ynr.ibp10[j]<-(m/(m-r))*(a10*Mxn[i]-b10*Mzn[i])
    ynr.ibp11[j]<-(1/(m-r))*((m*mean(YR)*((Mxn[i]/mean(XR))^a11)*((Mz/Mzn[i])^b11))-r*mean(YR))
    ynr.ibp12[j]<-(1/(m-r))*((m*mean(YR)*(Mxn[i]/(a12*Mxn[i]+(1-a12)*mean(XR)))*(Mz/(b12*Mzn[i]+(1-b12)*Mz)))-r*mean(YR))
    #bhushan and Pandey (2023)
    ynr.ibp13[j]<-(m/(m-r))*((gamma1-1)*mean(YR)+a13*Mxn[i]-b13*Mzn[i])
    ynr.ibp14[j]<-((r*gamma2*mean(YR)*(Mxn[i]/mean(XR))^a14*(Mz/Mzn[i])^b14)-r*mean(YR))/(m-r)
    ynr.ibp15[j]<-((r*gamma3*mean(YR)*(Mxn[i]/(a15*Mxn[i]+(1-a15)*mean(XR)))*(Mz/(b15*Mzn[i]+(1-b15)*Mz)))-r*mean(YR))/(m-r)
    #proposed
    tnr.iprop1[j]<-(1/(m-r))*((m*mean(YR)*exp(Pdelta1_1*((Mx-mean(XR))/(Mx+mean(XR)))+Pdelta1_2*((Mzn[i]-Mz)/(Mzn[i]+Mz))))-r*mean(YR))
    tnr.iprop2[j]<-(1/(m-r))*((m*mean(YR)*exp(Pdelta2_1*((Mx-mean(XR))/(Mx+mean(XR)))+Pdelta2_2*((mean(ZR)-Mz)/(mean(ZR)+Mz))))-r*mean(YR))
    tnr.iprop3[j]<-(1/(m-r))*((m*mean(YR)*exp(Pdelta3_1*((Mx-Mxn[i])/(Mx+Mxn[i]))+Pdelta3_2*((Mzn[i]-Mz)/(Mzn[i]+Mz))))-r*mean(YR))
  }
 
  y.im<-c(YR,ynr.im)
  y.ir<-c(YR,ynr.ir)
  y.icom<-c(yr.icom,ynr.icom)
  y.idp1<-c(yr.idp1,ynr.idp1)
  y.idp2<-c(yr.idp2,ynr.idp2)
  m.im[i]<- (1/m)*(sum(y.im))  
  m.ir[i]<-(1/m)*(sum(y.ir))  
  m.icom[i]<-(1/m)*(sum(y.icom))
  m.idp1[i]<-(1/m)*(sum(y.idp1))
  m.idp2[i]<-(1/m)*(sum(y.idp2))
  #Strategy-I
  y.isd<-c(YR,ynr.isd)
  y.iss<-c(YR,ynr.iss)
  y.ibp6<-c(yr.ibp6,ynr.ibp6)
  y.ibp9<-c(YR,ynr.ibp9)
  y.ibp11<-c(YR,ynr.ibp11)
  y.ibp12<-c(YR,ynr.ibp12)
  y.ibp14<-c(YR,ynr.ibp14)
  y.ibp15<-c(YR,ynr.ibp15)
  y.isps5<-c(YR,ynr.isps5)
  y.isps6<-c(YR,ynr.isps6)
  t.iprop1<-c(YR,tnr.iprop1)
  m.isd[i]<-(1/m)*(sum(y.isd))
  m.iss[i]<-(1/m)*(sum(y.iss))
  m.ibp6[i]<-(1/m)*(sum(y.ibp6))
  m.ibp9[i]<-(1/m)*(sum(y.ibp9))
  m.ibp11[i]<-(1/m)*(sum(y.ibp11))
  m.ibp12[i]<-(1/m)*(sum(y.ibp12))
  m.ibp14[i]<-(1/m)*(sum(y.ibp14))
  m.ibp15[i]<-(1/m)*(sum(y.ibp15))
  m.isps5[i]<-(1/m)*(sum(y.isps5))
  m.isps6[i]<-(1/m)*(sum(y.isps6))
  m.iprop1[i]<-(1/m)*(sum(t.iprop1))
  #Strategy-II
  # y.isf1r<-c(YR,ynr.isf1r)
  # y.isf1p<-c(YR,ynr.isf1p)
  # y.isf1v<-c(YR,ynr.isf1v)
  # y.ibp2<-c(yr.ibp2,ynr.ibp2)
  # y.ibp3<-c(yr.ibp3,ynr.ibp3)
  # y.ibp5<-c(yr.ibp5,ynr.ibp5)
  # y.ibp8<-c(YR,ynr.ibp8)
  # y.isps1<-c(YR,ynr.isps1)
  # y.isps2<-c(YR,ynr.isps2)
  # t.iprop2<-c(YR,tnr.iprop2)
  # m.isf1r[i]<-(1/m)*(sum(y.isf1r))
  # m.isf1p[i]<-(1/m)*(sum(y.isf1p))
  # m.isf1v[i]<-(1/m)*(sum(y.isf1v))
  # m.ibp2[i]<-(1/m)*(sum(y.ibp2))
  # m.ibp3[i]<-(1/m)*(sum(y.ibp3))
  # m.ibp5[i]<-(1/m)*(sum(y.ibp5))
  # m.ibp8[i]<-(1/m)*(sum(y.ibp8))
  # m.isps1[i]<-(1/m)*(sum(y.isps1))
  # m.isps2[i]<-(1/m)*(sum(y.isps2))
  # m.iprop2[i]<-(1/m)*(sum(t.iprop2))
  #Strategy-III
  # y.isf2r<-c(YR,ynr.isf2r)
  # y.isf2p<-c(YR,ynr.isf2p)
  # y.isf2v<-c(YR,ynr.isf2v)
  # y.idp3<-c(yr.idp3,ynr.idp3)
  # y.ibp1<-c(yr.ibp1,ynr.ibp1)
  # y.ibp4<-c(yr.ibp4,ynr.ibp4)
  # y.ibp7<-c(YR,ynr.ibp7)
  # y.ibp10<-c(yr.ibp10,ynr.ibp10)
  # y.ibp13<-c(yr.ibp13,ynr.ibp13)
  # y.isps3<-c(YR,ynr.isps3)
  # y.isps4<-c(YR,ynr.isps4)
  # t.iprop3<-c(YR,tnr.iprop3)
  # m.isf2r[i]<-(1/m)*(sum(y.isf2r))
  # m.isf2p[i]=(1/m)*(sum(y.isf2p))
  # m.isf2v[i]<-(1/m)*(sum(y.isf2v))
  # m.idp3[i]<-(1/m)*(sum(y.idp3))
  # m.ibp1[i]<-(1/m)*(sum(y.ibp1))
  # m.ibp4[i]<-(1/m)*(sum(y.ibp4))
  # m.ibp7[i]<-(1/m)*(sum(y.ibp7))
  # m.ibp10[i]<-(1/m)*(sum(y.ibp10))
  # m.ibp13[i]<-(1/m)*(sum(y.ibp13))
  # m.isps3[i]<-(1/m)*(sum(y.isps3))
  # m.isps4[i]<-(1/m)*(sum(y.isps4))
  # m.iprop3[i]<-(1/m)*(sum(t.iprop3))

}
calculate_mse <- function(estimates, true_value) {
  mse <- mean((estimates - true_value)^2)
  return(round(mse, 4))
}
calculate_PRE <- function(Mse_.im, mse_estimates) {
  if (mse_estimates == 0) {
    stop("Denominator MSE cannot be zero.")
  }
  PRE <- (Mse_.im/ mse_estimates) * 100
  return(round(PRE, 2))
}
calculate_relative_bias <- function(estimates, true_value) {
  rb <- mean(abs(estimates - true_value) / true_value)
  return(round(rb, 4))
}















