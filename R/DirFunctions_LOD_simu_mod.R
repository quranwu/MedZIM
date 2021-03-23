
#------------------------------------------
# data preparation function
#------------------------------------------

dataByZero=function(data,x,y,m,L,commonConfound=NULL,ConfoundX=NULL,ConfoundM=NULL){
  commonConfByMis=ConfoundX[ConfoundX%in%ConfoundM]
  commonConfound=unique(c(commonConfound,commonConfByMis))
  ConfoundX=unique(ConfoundX[!ConfoundX%in%commonConfByMis])
  ConfoundM=unique(ConfoundM[!ConfoundM%in%commonConfByMis])

  allConfoundX=data[,c(commonConfound,ConfoundX),drop=F]
  allConfoundM=data[,c(commonConfound,ConfoundM),drop=F]

  X=data[,x]
  Y=data[,y]
  L=data[,L]
  M.obs=data[,m]
  if((max(M.obs)>1) | (min(M.obs)<0)){
    stop("Mediator variables are out of the range of (0,1)")
    }
  indM=rep(0,length(M.obs))
  indM[which(M.obs>0)]=1

  M.g1=M.obs[indM==1]
  M.g2=M.obs[indM==0]
  Y.g1=Y[indM==1]
  Y.g2=Y[indM==0]
  X.g1=X[indM==1]
  X.g2=X[indM==0]
  L.g1=L[indM==1]
  L.g2=L[indM==0]
  allConfoundX.g1=allConfoundX[indM==1,,drop=F]
  allConfoundX.g2=allConfoundX[indM==0,,drop=F]
  allConfoundM.g1=allConfoundM[indM==1,,drop=F]
  allConfoundM.g2=allConfoundM[indM==0,,drop=F]

  d1.main=data.frame(X.g1=X.g1,Y.g1=Y.g1,M.g1=M.g1,L.g1=L.g1)
  d2.main=data.frame(X.g2=X.g2,Y.g2=Y.g2,M.g2=M.g2,L.g2=L.g2)
  list(X=X,Y=Y,L=L,M.obs=M.obs,indM=indM,d1.main=d1.main,d2.main=d2.main,
       allConfoundX.g1=allConfoundX.g1,allConfoundX.g2=allConfoundX.g2,
       allConfoundM.g1=allConfoundM.g1,allConfoundM.g2=allConfoundM.g2)
}

#--------------------------------------------------
# function for initial values and feasible regions
#-------------------------------------------------
ini.bound=function(X,Y,M.obs,indM,L){
  results=list()
  al=10^(-15)
  perturb=100
  level=1-al
  mui=M.obs[indM==1]
  X.Mobs=X[indM==1]

#--------------------------
# Beta regression work case
#--------------------------
  alpha0.ini=0
  alpha1.ini=0
  #cat("alpha0.ini",alpha0.ini,"\n")
  xi.ini=40
  LBxi.ini=10^(-5)
  UBxi.ini=Inf

  LBalpha0.ini=-Inf
  UBalpha0.ini=Inf
  LBalpha1.ini=-Inf
  UBalpha1.ini=Inf

    # print("LBalpha0.ini:")
    # print(LBalpha0.ini)
    # print("UBalpha0.ini:")
    # print(UBalpha0.ini)

    #-------------------------
    #  logistic regression
    #-------------------------

    logireg=glm((1-indM)~X,family="binomial")
    nulldev=logireg$null.deviance
    # print("summary(logireg):")
    # print(summary(logireg))
    # print("summary(logireg)$coefficients:")
    # print(summary(logireg)$coefficients)

    #-----------------------------------------------
    # logistic regression fails case (all positive)
    #-----------------------------------------------
    # results$logisticFail=0
    # if(nulldev==0){
    #   stop("Logistic regression failed.")
    #   results$logisticFail=1
    #   }
    #-----------------------------------------------
    # logistic regression works (not all M possitive)
    #-----------------------------------------------

    #if(nulldev!=0){
    if(T){

      gamma0.ini=summary(logireg)$coefficients[1,1]
      gamma1.ini=summary(logireg)$coefficients[2,1]

      LBgamma0.ini1=gamma0.ini-abs(qnorm(al/2))*summary(logireg)$coefficients[1,2]-perturb*abs(gamma0.ini)
      UBgamma0.ini1=gamma0.ini+abs(qnorm(al/2))*summary(logireg)$coefficients[1,2]+perturb*abs(gamma0.ini)
      LBgamma1.ini1=gamma1.ini-abs(qnorm(al/2))*summary(logireg)$coefficients[2,2]-perturb*abs(gamma1.ini)
      UBgamma1.ini1=gamma1.ini+abs(qnorm(al/2))*summary(logireg)$coefficients[2,2]+perturb*abs(gamma1.ini)

      # print("LBgamma0.ini1:")
      # print(LBgamma0.ini1)
      # print("UBgamma0.ini1:")
      # print(UBgamma0.ini1)

      # change 10% of the observed zero's to 1 for calculating initial values
      indM2=indM
      indM2[sample(which(indM==0),0.1*length(which(indM==0)))]=1

      #mean(logiM.obs)
      logireg2=glm(indM2~X,family="binomial")
      gamma0.ini2=summary(logireg2)$coefficients[1,1]
      gamma1.ini2=summary(logireg2)$coefficients[2,1]
      # print("summary(logireg2)$coefficients:")
      # print(summary(logireg2)$coefficients)
      # print("gamma0.ini2:")
      # print(gamma0.ini2)
      # print("gamma1.ini2:")
      # print(gamma1.ini2)

      LBgamma0.ini2=gamma0.ini2-abs(qnorm(al/2))*summary(logireg2)$coefficients[1,2]-perturb*abs(gamma0.ini2)
      UBgamma0.ini2=gamma0.ini2+abs(qnorm(al/2))*summary(logireg2)$coefficients[1,2]+perturb*abs(gamma0.ini2)
      LBgamma1.ini2=gamma1.ini2-abs(qnorm(al/2))*summary(logireg2)$coefficients[2,2]-perturb*abs(gamma1.ini2)
      UBgamma1.ini2=gamma1.ini2+abs(qnorm(al/2))*summary(logireg2)$coefficients[2,2]+perturb*abs(gamma1.ini2)

      LBgamma0.ini=min(LBgamma0.ini1,LBgamma0.ini2)
      LBgamma1.ini=min(LBgamma1.ini1,LBgamma1.ini2)
      UBgamma0.ini=max(UBgamma0.ini1,UBgamma0.ini2)
      UBgamma1.ini=max(UBgamma1.ini1,UBgamma1.ini2)

      #---------------------------
      #linear regression
      #--------------------------

      results$glmFail=0
      regY=glm(Y~M.obs+indM+X+X*indM+X*M.obs,family="gaussian")
      # if(sum(is.na(regY$coefficients))>0){
      #   stop("GLM failed.")
      #   results$glmFail=1
      #   }
        print("regY$coefficients:")
        print(regY$coefficients)
        #  print("summary(regY)$coefficients:")
        #  print(summary(regY)$coefficients)

      #if(sum(is.na(regY$coefficients))==0){
      beta0.ini1=regY$coefficients[1]
      beta0.ini=beta0.ini1
      if(is.na(beta0.ini1))beta0.ini=0
      beta1.ini1=regY$coefficients[2]
      beta1.ini=beta1.ini1
      if(is.na(beta1.ini1))beta1.ini=0
      beta2.ini1=regY$coefficients[3]
      beta2.ini=beta2.ini1
      if(is.na(beta2.ini1))beta2.ini=0
      beta3.ini1=regY$coefficients[4]
      beta3.ini=beta3.ini1
      if(is.na(beta3.ini1))beta3.ini=0
      beta4.ini1=regY$coefficients[5]
      beta4.ini=beta4.ini1
      if(is.na(beta4.ini1))beta4.ini=0
      beta5.ini1=regY$coefficients[6]
      beta5.ini=beta5.ini1
      if(is.na(beta5.ini1))beta5.ini=0

        al2=25
        if(!is.na(beta0.ini1)){
         LBbeta0.ini=beta0.ini-al2*summary(regY)$coefficients[1,2]-perturb*abs(beta0.ini)
         }
        if(is.na(beta0.ini1))LBbeta0.ini=-Inf

        if(!is.na(beta0.ini1)){
         UBbeta0.ini=beta0.ini+al2*summary(regY)$coefficients[1,2]+perturb*abs(beta0.ini)
        }
        if(is.na(beta0.ini1))UBbeta0.ini=Inf

        if(!is.na(beta1.ini1)){
         LBbeta1.ini=beta1.ini-al2*summary(regY)$coefficients[2,2]-perturb*abs(beta1.ini)
        }
        if(is.na(beta1.ini1))LBbeta1.ini=-Inf

        if(!is.na(beta1.ini1)){
         UBbeta1.ini=beta1.ini+al2*summary(regY)$coefficients[2,2]+perturb*abs(beta1.ini)
        }
        if(is.na(beta1.ini1))UBbeta1.ini=Inf

        if(!is.na(beta2.ini1)){
         LBbeta2.ini=beta2.ini-al2*summary(regY)$coefficients[3,2]-perturb*abs(beta2.ini)
        }
        if(is.na(beta2.ini1))LBbeta2.ini=-Inf

        if(!is.na(beta2.ini1)){
         UBbeta2.ini=beta2.ini+al2*summary(regY)$coefficients[3,2]+perturb*abs(beta2.ini)
        }
        if(is.na(beta2.ini1))UBbeta2.ini=Inf

        if(!is.na(beta3.ini1)){
         LBbeta3.ini=beta3.ini-al2*summary(regY)$coefficients[4,2]-perturb*abs(beta3.ini)
        }
        if(is.na(beta3.ini1))LBbeta3.ini=-Inf

        if(!is.na(beta3.ini1)){
         UBbeta3.ini=beta3.ini+al2*summary(regY)$coefficients[4,2]+perturb*abs(beta3.ini)
        }
        if(is.na(beta3.ini1))UBbeta3.ini=Inf

        if(!is.na(beta4.ini1)){
         LBbeta4.ini=beta4.ini-al2*summary(regY)$coefficients[5,2]-perturb*abs(beta4.ini)
        }
        if(is.na(beta4.ini1))LBbeta4.ini=-Inf

        if(!is.na(beta4.ini1)){
         UBbeta4.ini=beta4.ini+al2*summary(regY)$coefficients[5,2]+perturb*abs(beta4.ini)
        }
        if(is.na(beta4.ini1))UBbeta4.ini=Inf

        if(!is.na(regY$coefficients[5]) &  !is.na(regY$coefficients[6])){
         LBbeta5.ini=beta5.ini-al2*summary(regY)$coefficients[6,2]-perturb*abs(beta5.ini)
        }else {LBbeta5.ini=-Inf}

        if(!is.na(regY$coefficients[5]) & !is.na(regY$coefficients[6])){
         UBbeta5.ini=beta5.ini+al2*summary(regY)$coefficients[6,2]+perturb*abs(beta5.ini)
        }else {UBbeta5.ini=Inf}

        #}
        # print("summary(regY):")
        # print(summary(regY))

        delta.ini=sqrt(summary(regY)$dispersion)
        if(is.na(delta.ini))delta.ini=1

        LBdelta.ini=10^(-5)

        UBdelta.ini=delta.ini*100+perturb*abs(delta.ini)
        if(is.na(UBdelta.ini))UBdelta.ini=Inf

       }
    out=list(LBbeta0.ini=LBbeta0.ini,beta0.ini=beta0.ini,UBbeta0.ini=UBbeta0.ini,LBbeta1.ini=LBbeta1.ini,
             beta1.ini=beta1.ini, UBbeta1.ini=UBbeta1.ini,LBbeta2.ini=LBbeta2.ini,beta2.ini=beta2.ini,
             UBbeta2.ini=UBbeta2.ini,LBbeta3.ini=LBbeta3.ini,beta3.ini=beta3.ini, UBbeta3.ini=UBbeta3.ini,
             LBbeta4.ini=LBbeta4.ini,beta4.ini=beta4.ini,UBbeta4.ini=UBbeta4.ini,
             LBbeta5.ini=LBbeta5.ini,beta5.ini=beta5.ini,UBbeta5.ini=UBbeta5.ini,
             delta.ini=delta.ini,LBdelta.ini=LBdelta.ini,UBdelta.ini=UBdelta.ini,
             LBalpha0.ini=LBalpha0.ini,alpha0.ini=alpha0.ini,UBalpha0.ini=UBalpha0.ini,
             LBalpha1.ini=LBalpha1.ini,alpha1.ini=alpha1.ini,UBalpha1.ini=UBalpha1.ini,LBxi.ini=LBxi.ini,
             xi.ini=xi.ini, UBxi.ini=UBxi.ini,
             LBgamma0.ini=LBgamma0.ini,gamma0.ini=gamma0.ini,UBgamma0.ini=UBgamma0.ini,
             LBgamma1.ini=LBgamma1.ini,gamma1.ini=gamma1.ini,
             UBgamma1.ini=UBgamma1.ini
             )
    return(out)

}

#--------------------
#likelihood function
#--------------------

hi=function(m,a0,a1,b0,b1,b2,b3,b4,b5,delta,xi,X2,Y2){
  mug2=exp(a0+a1*X2)/(1+exp(a0+a1*X2))
  phig2=xi
  frst=(m^(mug2*phig2-1))*((1-m)^((1-mug2)*phig2-1))
  xbt=b0+b1*m+b2+b3*X2+b4*X2+b5*X2*m
  secd=exp(-(Y2-xbt)^2/(2*delta^2))
  hi_f=frst*secd
  return(hi_f)
  }

hii=function(m,a0,a1,b0,b1,b2,b3,b4,b5,delta,xi,X2,Y2,L2){
  mug2=exp(a0+a1*X2)/(1+exp(a0+a1*X2))
  phig2=xi
  frst=((m*L2)^(mug2*phig2-1))*(((1-m)/(1-1/L2))^((1-mug2)*phig2-1))
  #cat("frst:",frst,"\n")
  xbt=b0+b1*m+b2+b3*X2+b4*X2+b5*X2*m
  xbtL2=b0+b1/L2+b2+b3*X2+b4*X2+b5*X2/L2
  secd=exp(-(Y2-xbt)^2/(2*delta^2)+(Y2-xbtL2)^2/(2*delta^2))
  #cat("secd:",secd,"\n")

  hi_f=frst*secd
  return(hi_f)
  }

hiSimp=function(m,a0,a1,b0,b1,b2,b3,b4,b5,delta,xi,X2,Y2,L2){
  mug2=exp(a0+a1*X2)/(1+exp(a0+a1*X2))
  phig2=xi
  frst=(m^(mug2*phig2-1))*((1-m)^((1-mug2)*phig2-1))

  yxbeta=Y2-b0-b2-(b3+b4)*X2
  b15x=b1+b5*X2
  intgrand1power=(b15x^2*(m^2)-2*yxbeta*b15x*m)/(2*delta^2)
  secd=exp(-intgrand1power)
  hi_f=frst*secd
  return(hi_f)
}

hi3=function(m,a0,a1,b0,b1,b2,b3,b4,b5,delta,xi,X2,Y2,L2){
  mug2=exp(a0+a1*X2)/(1+exp(a0+a1*X2))
  phig2=xi
  frst=(m^(mug2*phig2-1))*((1-m)^((1-mug2)*phig2-1))
  #cat("frst:",frst[1:3],"\n")
  XMbeta=b1*m+b2+b4*X2+b5*X2*m
  yXbeta=Y2-b0-b3*X2
  intgrand1power=XMbeta^2-2*yXbeta*XMbeta

  secd=exp(-intgrand1power)
  #cat("secd:",secd[1:3],"\n")
  #stop()
  hi_f=frst*secd
  return(hi_f)
}

negLogL=function(theta,dg1,dg2){
  beta0=theta[1]
  beta1=theta[2]
  beta2=theta[3]
  beta3=theta[4]
  beta4=theta[5]
  beta5=theta[6]
  delta=theta[7]
  alpha0=theta[8]
  alpha1=theta[9]
  xi=theta[10]
  gamma0=theta[11]
  gamma1=theta[12]

  #cat("theta",theta,"\n\n")
  X.g1=dg1$X.g1
  M.g1=dg1$M.g1
  Y.g1=dg1$Y.g1
  L.g1=dg1$L.g1

  ### no g2
  X.g2=dg2$X.g2
  M.g2=dg2$M.g2
  Y.g2=dg2$Y.g2
  L.g2=dg2$L.g2

  # Likelihood function on group1
  l1_i.1=-log(delta)

  xbt1=beta0+(beta1*M.g1)+beta2+(beta3*X.g1)+beta4*X.g1+beta5*X.g1*M.g1

  l1_i.2=-(Y.g1-xbt1)^2/(2*delta^2)

  l1_i.3=0



  Xgamma=gamma0+(gamma1*X.g1)
  posSlots=which(Xgamma>=0)
  negSlots=which(Xgamma<0)
  XgammaPos=Xgamma[posSlots]
  XgammaNeg=Xgamma[negSlots]
  l1_i.4=rep(NA,length(Xgamma))
  l1_i.4pos=-XgammaPos-log(1+exp(-XgammaPos))
  l1_i.4neg=-log(1+exp(XgammaNeg))
  l1_i.4[posSlots]=l1_i.4pos
  l1_i.4[negSlots]=l1_i.4neg
  rm(Xgamma,posSlots,negSlots,XgammaPos,XgammaNeg,l1_i.4pos,l1_i.4neg)

  Xalpha=alpha0+alpha1*X.g1
  posSlots=which(Xalpha>=0)
  negSlots=which(Xalpha<0)
  XalphaPos=Xalpha[posSlots]
  XalphaNeg=Xalpha[negSlots]
  mui1=rep(NA,length(Xalpha))
  mui1.pos=1/(1+exp(-XalphaPos))
  mui1.neg=exp(XalphaNeg)/(1+exp(XalphaNeg))
  mui1[posSlots]=mui1.pos
  mui1[negSlots]=mui1.neg
  rm(Xalpha,posSlots,negSlots,XalphaPos,XalphaNeg,mui1.pos,mui1.neg)

  phi1=xi
  #cat("mui1",mui1,"\n\n")
  #cat("phi1",phi1,"\n\n")
  a1=mui1*phi1
  b1=(1-mui1)*phi1

  ##B(mu_i phi )
  bet=beta(a1,b1)

  # cat("bet",bet,"\n")
  if(!(is.element(0,bet)) && !(is.element(Inf,bet))){l1_i.5=-log(bet)}
  if(!(is.element(0,bet)) && (is.element(Inf,bet))){
    a1[which(a1==0)]=10^(-300)
    b1[which(b1==0)]=10^(-300)
    l1_i.5=-log(beta(a1,b1))
    }
  if((is.element(0,bet)) && !(is.element(Inf,bet))){
    phi1=10^3
    #cat("phi1",phi1,"\n")
    a1=mui1*phi1
    b1=(1-mui1)*phi1
    #cat("a1",a1,"\n")
    #cat("b1",b1,"\n")
    l1_i.5=-log(beta(a1,b1))
    }
  if((is.element(0,bet)) && (is.element(Inf,bet))){
    a1[which(a1==0)]=10^(-300)
    b1[which(b1==0)]=10^(-300)
    phi1=10^3
    #cat("phi1",phi1,"\n")
    a1=mui1*phi1
    b1=(1-mui1)*phi1
    #cat("a1",a1,"\n")
    #cat("b1",b1,"\n")
    l1_i.5=-log(beta(a1,b1))
    }

  ########
  l1_i.6=(mui1*phi1-1)*log(M.g1)
  l1_i.7=((1-mui1)*phi1-1)*log(1-M.g1)
  # cat("l1_i.1:",l1_i.1,"\n")
  # cat("l1_i.2:",l1_i.2,"\n")
  # cat("l1_i.3:",l1_i.3,"\n")
  # cat("l1_i.4:",l1_i.4,"\n")
  # cat("l1_i.5:",l1_i.5,"\n")
  # cat("l1_i.6:",l1_i.6,"\n")
  # cat("l1_i.7:",l1_i.7,"\n")

  #### delete l1_i.4
  l1_i=l1_i.1+l1_i.2+l1_i.3+l1_i.4+l1_i.5+l1_i.6+l1_i.7

  # Likelihood function on group2
  Xgamma.g2=gamma0+(gamma1*X.g2)
  posSlots=which(Xgamma.g2>=0)
  negSlots=which(Xgamma.g2<0)
  Xgamma.g2Pos=Xgamma.g2[posSlots]
  Xgamma.g2Neg=Xgamma.g2[negSlots]
  expit2=rep(NA,length(Xgamma.g2))
  logExpit2=rep(NA,length(Xgamma.g2))
  expit2.pos=1/(1+exp(-Xgamma.g2Pos))
  logExpit2.pos=-log(1+exp(-Xgamma.g2Pos))

  expit2.neg=exp(Xgamma.g2Neg)/(1+exp(Xgamma.g2Neg))
  logExpit2.neg=Xgamma.g2Neg-log(1+exp(Xgamma.g2Neg))

  expit2[posSlots]=expit2.pos
  logExpit2[posSlots]=logExpit2.pos

  expit2[negSlots]=expit2.neg
  logExpit2[negSlots]=logExpit2.neg

  rm(posSlots,negSlots,Xgamma.g2Pos,Xgamma.g2Neg,expit2.pos,
     expit2.neg,logExpit2.pos,logExpit2.neg)

  Xbeta03=beta0+(beta3*X.g2)
  exp03=exp(-(Y.g2-Xbeta03)^2/(2*delta^2))

  scnd2.1=expit2*exp03

  scnd22.1=logExpit2
  scnd22.2=-(Y.g2-Xbeta03)^2/(2*delta^2)
  expXgamma=exp(-Xgamma.g2)
  rm(Xgamma.g2)

  Xalpha.g2=alpha0+alpha1*X.g2
  posSlots=which(Xalpha.g2>=0)
  negSlots=which(Xalpha.g2<0)
  Xalpha.g2Pos=Xalpha.g2[posSlots]
  Xalpha.g2Neg=Xalpha.g2[negSlots]
  mui2=rep(NA,length(Xalpha.g2))
  mui2.pos=1/(1+exp(-Xalpha.g2Pos))
  mui2.neg=exp(Xalpha.g2Neg)/(1+exp(Xalpha.g2Neg))
  mui2[posSlots]=mui2.pos
  mui2[negSlots]=mui2.neg
  rm(Xalpha.g2,posSlots,negSlots,Xalpha.g2Pos,Xalpha.g2Neg,mui2.pos,mui2.neg)

  phi2=xi

  a2=mui2*phi2
  b2=(1-mui2)*phi2

  #time0=proc.time()[3]

  bet2=beta(a2,b2)
  # #cat("bet2[1:5]:",bet2[1:5],"\n")
  #
  int.approx=c()
  for(i in 1:length(X.g2)){
    fx=function(x)hii(m=x,a0=alpha0,a1=alpha1,b0=beta0,b1=beta1,b2=beta2,b3=beta3,
                     b4=beta4,b5=beta5,
                     delta=delta,xi=xi,X2=X.g2[i],Y2=Y.g2[i],L2=L.g2[i])

    # nPoints=1000
    # gridP=seq(0,(1/L.g2[i]),length=(nPoints+1))
    # gridV=fx(gridP[-1])
    # nonInfGrid=gridV[!is.infinite(gridV)]
    # int.approx[i]=nonInfGrid[1]*(nPoints-length(nonInfGrid)+1)/(nPoints*L.g2[i])+
    #   sum(nonInfGrid[-1])/(nPoints*L.g2[i])
    # cat("gridV[length(gridV)]:",gridV[length(gridV)],"\n")
    # cat("fx(1/L.g2[i]):",fx(1/L.g2[i]),"\n")
    # cat("gridV:",gridV[1],"\n")
    # cat("fx(gridP[2]):",fx(gridP[2]),"\n")
    # cat("gridP[2]:",gridP[2],"\n")
    # cat("gridP[-1]:",gridP[2:5],"\n")
    # cat("a2[i]:",a2[i],"\n")
    # cat("b2[i]:",b2[i],"\n")

    fx_con=function(x)hi(m=x,a0=alpha0,a1=alpha1,b0=beta0,b1=beta1,b2=beta2,b3=beta3,
                      b4=beta4,b5=beta5,
                      delta=delta,xi=xi,X2=X.g2[i],Y2=Y.g2[i])

    # int.approx[i]=fx_con(1/L.g2[i])*quadinf(fx,0,(1/L.g2[i]))$Q

    #cat("fx_con(1/L.g2[i]):",fx_con(1/L.g2[i]),"\n")
    #cat("int.approx[i]:",int.approx[i],"\n")

    #int.approx[i]=quadinfInt(fx,0,1)$Q
    #int.approx[i]=adaptIntegrate(fx,0,1)$integral
    #int.approx[i]=integrate(fx,0,1)$value
    }
  #cat("eta:",eta,"\n")
  #cat("Y.g2[1:5]:",Y.g2[1:5],"\n")
  #cat("L.g2[1:5]:",L.g2[1:5],"\n")
  #cat("int.approx[1:5]:",int.approx[1:5],"\n")

  # betPart3=int.approx/bet2

  #betPart=1
  #betPart1=betPart
  #cat("betPart1[1:5]:",betPart1[1:5],"\n")

  #time1=proc.time()[3]

  #cat("betPart1 time:",time1-time0,"seconds","\n")

  # nMC=10000
  # MCptsMat=matrix(NA,nrow=length(X.g2),ncol=nMC)
  # setSeed(1)
  # u01=torus(nMC)
  # for(i in 1:length(X.g2)){
  #   MCptsMat[i,]=qbeta(u01,a2[i],b2[i])
  #   }
  #
  # betPart<-c()
  #
  # for (i in 1:nrow(MCptsMat)) {
  #   betPart[i]<-MCfunccppvector(m=MCptsMat[i,],b0=beta0,b1=beta1,b2=beta2,b3=beta3,b4=beta4,
  #                                     delta=delta,X2vec=X.g2[i],Y2vec=Y.g2[i],L2vec=L.g2[i])
  # }
  #
  # cat("betPart[1:5]:",betPart[1:5],"\n")
  # time2=proc.time()[3]

  betPart<-c()

  for (i in 1:length(X.g2)) {
    yxbeta=Y.g2[i]-beta0-beta2-(beta3+beta4)*X.g2[i]
    b15x=beta1+beta5*X.g2[i]

    expyx=exp(-yxbeta^2/(2*delta^2))

    # intgrand1up=(b15x^2/((L.g2[i])^2)-2*yxbeta*b15x/(L.g2[i]))/(2*delta^2)
    # # if(intgrand1up>=10^(-2))cat("intgrand1up:",intgrand1up,"\n")
    #
    # if(abs(intgrand1up)<10^(-2)){
    betPart[i]=expyx*pbeta((1/L.g2[i]),a2[i],b2[i])
    # #   }
    # # if(abs(intgrand1up)>=10^(-2)){
    #   nMC=10000
    #   set.seed(1)
    #   MCptsMat=rbeta(nMC,a2[i],b2[i])
    #   MCint<-MCfunccppvector(m=MCptsMat,b0=beta0,b1=beta1,b2=beta2,
    #                          b3=beta3,b4=beta4,
    #                          delta=delta,e=0,X2vec=X.g2[i],
    #                          Y2vec=Y.g2[i],L2vec=L.g2[i])
    #   betPart[i]=MCint
    }
   scnd2.2=(1-expit2)*betPart

  # scnd2.2=(1-expit2)*betPart3

  #cat("scnd2.1",scnd2.1,"\n")
  #cat("scnd2.2",scnd2.2,"\n")

  scnd2=-log(delta)+log(scnd2.1+scnd2.2)
  # print("scnd2.1:")
  # print(scnd2.1[1:10])
  # print("scnd2.2:")
  # print(scnd2.2[1:10])
  # scnd2=-log(delta)+scnd22.1+scnd22.2+log(1+betPart3)

  l2_i=scnd2

  #cat("sum(l2_i):",sum(l2_i),"\n")

  #### delete sum(l2_i)

  l=sum(l1_i)+sum(l2_i)
  #print("theta:")
  #print(theta)
  # cat("logLikehood function value:", l,"\n")
  return(-l)
  }


negLogL_contin=function(theta,dg1){
  beta0=theta[1]
  beta1=theta[2]
  beta2=theta[3]
  beta3=theta[4]
  beta4=theta[5]
  beta5=theta[6]
  delta=theta[7]
  alpha0=theta[8]
  alpha1=theta[9]
  xi=theta[10]
  # gamma0=theta[11]
  # gamma1=theta[12]

  #cat("theta",theta,"\n\n")
  X.g1=dg1$X.g1
  M.g1=dg1$M.g1
  Y.g1=dg1$Y.g1
  L.g1=dg1$L.g1

  ### no g2
  # X.g2=dg2$X.g2
  # M.g2=dg2$M.g2
  # Y.g2=dg2$Y.g2
  # L.g2=dg2$L.g2

  # Likelihood function on group1
  l1_i.1=-log(delta)

  xbt1=beta0+(beta1*M.g1)+beta2+(beta3*X.g1)+beta4*X.g1+beta5*X.g1*M.g1

  l1_i.2=-(Y.g1-xbt1)^2/(2*delta^2)

  l1_i.3=0



  # Xgamma=gamma0+(gamma1*X.g1)
  # posSlots=which(Xgamma>=0)
  # negSlots=which(Xgamma<0)
  # XgammaPos=Xgamma[posSlots]
  # XgammaNeg=Xgamma[negSlots]
  # l1_i.4=rep(NA,length(Xgamma))
  # l1_i.4pos=-XgammaPos-log(1+exp(-XgammaPos))
  # l1_i.4neg=-log(1+exp(XgammaNeg))
  # l1_i.4[posSlots]=l1_i.4pos
  # l1_i.4[negSlots]=l1_i.4neg
  # rm(Xgamma,posSlots,negSlots,XgammaPos,XgammaNeg,l1_i.4pos,l1_i.4neg)

  l1_i.4<-0

  Xalpha=alpha0+alpha1*X.g1
  posSlots=which(Xalpha>=0)
  negSlots=which(Xalpha<0)
  XalphaPos=Xalpha[posSlots]
  XalphaNeg=Xalpha[negSlots]
  mui1=rep(NA,length(Xalpha))
  mui1.pos=1/(1+exp(-XalphaPos))
  mui1.neg=exp(XalphaNeg)/(1+exp(XalphaNeg))
  mui1[posSlots]=mui1.pos
  mui1[negSlots]=mui1.neg
  rm(Xalpha,posSlots,negSlots,XalphaPos,XalphaNeg,mui1.pos,mui1.neg)

  phi1=xi
  #cat("mui1",mui1,"\n\n")
  #cat("phi1",phi1,"\n\n")
  a1=mui1*phi1
  b1=(1-mui1)*phi1

  ##B(mu_i phi )
  bet=beta(a1,b1)

  # cat("bet",bet,"\n")
  if(!(is.element(0,bet)) && !(is.element(Inf,bet))){l1_i.5=-log(bet)}
  if(!(is.element(0,bet)) && (is.element(Inf,bet))){
    a1[which(a1==0)]=10^(-300)
    b1[which(b1==0)]=10^(-300)
    l1_i.5=-log(beta(a1,b1))
  }
  if((is.element(0,bet)) && !(is.element(Inf,bet))){
    phi1=10^3
    #cat("phi1",phi1,"\n")
    a1=mui1*phi1
    b1=(1-mui1)*phi1
    #cat("a1",a1,"\n")
    #cat("b1",b1,"\n")
    l1_i.5=-log(beta(a1,b1))
  }
  if((is.element(0,bet)) && (is.element(Inf,bet))){
    a1[which(a1==0)]=10^(-300)
    b1[which(b1==0)]=10^(-300)
    phi1=10^3
    #cat("phi1",phi1,"\n")
    a1=mui1*phi1
    b1=(1-mui1)*phi1
    #cat("a1",a1,"\n")
    #cat("b1",b1,"\n")
    l1_i.5=-log(beta(a1,b1))
  }

  ########
  l1_i.6=(mui1*phi1-1)*log(M.g1)
  l1_i.7=((1-mui1)*phi1-1)*log(1-M.g1)
  # cat("l1_i.1:",l1_i.1,"\n")
  # cat("l1_i.2:",l1_i.2,"\n")
  # cat("l1_i.3:",l1_i.3,"\n")
  # cat("l1_i.4:",l1_i.4,"\n")
  # cat("l1_i.5:",l1_i.5,"\n")
  # cat("l1_i.6:",l1_i.6,"\n")
  # cat("l1_i.7:",l1_i.7,"\n")

  #### delete l1_i.4
  l1_i=l1_i.1+l1_i.2+l1_i.3+l1_i.4+l1_i.5+l1_i.6+l1_i.7

  # Likelihood function on group2
  # Xgamma.g2=gamma0+(gamma1*X.g2)
  # posSlots=which(Xgamma.g2>=0)
  # negSlots=which(Xgamma.g2<0)
  # Xgamma.g2Pos=Xgamma.g2[posSlots]
  # Xgamma.g2Neg=Xgamma.g2[negSlots]
  # expit2=rep(NA,length(Xgamma.g2))
  # logExpit2=rep(NA,length(Xgamma.g2))
  # expit2.pos=1/(1+exp(-Xgamma.g2Pos))
  # logExpit2.pos=-log(1+exp(-Xgamma.g2Pos))
  #
  # expit2.neg=exp(Xgamma.g2Neg)/(1+exp(Xgamma.g2Neg))
  # logExpit2.neg=Xgamma.g2Neg-log(1+exp(Xgamma.g2Neg))
  #
  # expit2[posSlots]=expit2.pos
  # logExpit2[posSlots]=logExpit2.pos
  #
  # expit2[negSlots]=expit2.neg
  # logExpit2[negSlots]=logExpit2.neg
  #
  # rm(posSlots,negSlots,Xgamma.g2Pos,Xgamma.g2Neg,expit2.pos,
  #    expit2.neg,logExpit2.pos,logExpit2.neg)
  #
  # Xbeta03=beta0+(beta3*X.g2)
  # exp03=exp(-(Y.g2-Xbeta03)^2/(2*delta^2))
  #
  # scnd2.1=expit2*exp03
  #
  # scnd22.1=logExpit2
  # scnd22.2=-(Y.g2-Xbeta03)^2/(2*delta^2)
  # expXgamma=exp(-Xgamma.g2)
  # rm(Xgamma.g2)
  #
  # Xalpha.g2=alpha0+alpha1*X.g2
  # posSlots=which(Xalpha.g2>=0)
  # negSlots=which(Xalpha.g2<0)
  # Xalpha.g2Pos=Xalpha.g2[posSlots]
  # Xalpha.g2Neg=Xalpha.g2[negSlots]
  # mui2=rep(NA,length(Xalpha.g2))
  # mui2.pos=1/(1+exp(-Xalpha.g2Pos))
  # mui2.neg=exp(Xalpha.g2Neg)/(1+exp(Xalpha.g2Neg))
  # mui2[posSlots]=mui2.pos
  # mui2[negSlots]=mui2.neg
  # rm(Xalpha.g2,posSlots,negSlots,Xalpha.g2Pos,Xalpha.g2Neg,mui2.pos,mui2.neg)
  #
  # phi2=xi
  #
  # a2=mui2*phi2
  # b2=(1-mui2)*phi2
  #
  # #time0=proc.time()[3]
  #
  # bet2=beta(a2,b2)
  # # #cat("bet2[1:5]:",bet2[1:5],"\n")
  # #
  # int.approx=c()
  # for(i in 1:length(X.g2)){
  #   fx=function(x)hii(m=x,a0=alpha0,a1=alpha1,b0=beta0,b1=beta1,b2=beta2,b3=beta3,
  #                     b4=beta4,b5=beta5,
  #                     delta=delta,xi=xi,X2=X.g2[i],Y2=Y.g2[i],L2=L.g2[i])
  #
  #   # nPoints=1000
  #   # gridP=seq(0,(1/L.g2[i]),length=(nPoints+1))
  #   # gridV=fx(gridP[-1])
  #   # nonInfGrid=gridV[!is.infinite(gridV)]
  #   # int.approx[i]=nonInfGrid[1]*(nPoints-length(nonInfGrid)+1)/(nPoints*L.g2[i])+
  #   #   sum(nonInfGrid[-1])/(nPoints*L.g2[i])
  #   # cat("gridV[length(gridV)]:",gridV[length(gridV)],"\n")
  #   # cat("fx(1/L.g2[i]):",fx(1/L.g2[i]),"\n")
  #   # cat("gridV:",gridV[1],"\n")
  #   # cat("fx(gridP[2]):",fx(gridP[2]),"\n")
  #   # cat("gridP[2]:",gridP[2],"\n")
  #   # cat("gridP[-1]:",gridP[2:5],"\n")
  #   # cat("a2[i]:",a2[i],"\n")
  #   # cat("b2[i]:",b2[i],"\n")
  #
  #   fx_con=function(x)hi(m=x,a0=alpha0,a1=alpha1,b0=beta0,b1=beta1,b2=beta2,b3=beta3,
  #                        b4=beta4,b5=beta5,
  #                        delta=delta,xi=xi,X2=X.g2[i],Y2=Y.g2[i])
  #
  #   # int.approx[i]=fx_con(1/L.g2[i])*quadinf(fx,0,(1/L.g2[i]))$Q
  #
  #   #cat("fx_con(1/L.g2[i]):",fx_con(1/L.g2[i]),"\n")
  #   #cat("int.approx[i]:",int.approx[i],"\n")
  #
  #   #int.approx[i]=quadinfInt(fx,0,1)$Q
  #   #int.approx[i]=adaptIntegrate(fx,0,1)$integral
  #   #int.approx[i]=integrate(fx,0,1)$value
  # }
  # #cat("eta:",eta,"\n")
  # #cat("Y.g2[1:5]:",Y.g2[1:5],"\n")
  # #cat("L.g2[1:5]:",L.g2[1:5],"\n")
  # #cat("int.approx[1:5]:",int.approx[1:5],"\n")
  #
  # # betPart3=int.approx/bet2
  #
  # #betPart=1
  # #betPart1=betPart
  # #cat("betPart1[1:5]:",betPart1[1:5],"\n")
  #
  # #time1=proc.time()[3]
  #
  # #cat("betPart1 time:",time1-time0,"seconds","\n")
  #
  # # nMC=10000
  # # MCptsMat=matrix(NA,nrow=length(X.g2),ncol=nMC)
  # # setSeed(1)
  # # u01=torus(nMC)
  # # for(i in 1:length(X.g2)){
  # #   MCptsMat[i,]=qbeta(u01,a2[i],b2[i])
  # #   }
  # #
  # # betPart<-c()
  # #
  # # for (i in 1:nrow(MCptsMat)) {
  # #   betPart[i]<-MCfunccppvector(m=MCptsMat[i,],b0=beta0,b1=beta1,b2=beta2,b3=beta3,b4=beta4,
  # #                                     delta=delta,X2vec=X.g2[i],Y2vec=Y.g2[i],L2vec=L.g2[i])
  # # }
  # #
  # # cat("betPart[1:5]:",betPart[1:5],"\n")
  # # time2=proc.time()[3]
  #
  # betPart<-c()
  #
  # for (i in 1:length(X.g2)) {
  #   yxbeta=Y.g2[i]-beta0-beta2-(beta3+beta4)*X.g2[i]
  #   b15x=beta1+beta5*X.g2[i]
  #
  #   expyx=exp(-yxbeta^2/(2*delta^2))
  #
  #   # intgrand1up=(b15x^2/((L.g2[i])^2)-2*yxbeta*b15x/(L.g2[i]))/(2*delta^2)
  #   # # if(intgrand1up>=10^(-2))cat("intgrand1up:",intgrand1up,"\n")
  #   #
  #   # if(abs(intgrand1up)<10^(-2)){
  #   betPart[i]=expyx*pbeta((1/L.g2[i]),a2[i],b2[i])
  #   # #   }
  #   # # if(abs(intgrand1up)>=10^(-2)){
  #   #   nMC=10000
  #   #   set.seed(1)
  #   #   MCptsMat=rbeta(nMC,a2[i],b2[i])
  #   #   MCint<-MCfunccppvector(m=MCptsMat,b0=beta0,b1=beta1,b2=beta2,
  #   #                          b3=beta3,b4=beta4,
  #   #                          delta=delta,e=0,X2vec=X.g2[i],
  #   #                          Y2vec=Y.g2[i],L2vec=L.g2[i])
  #   #   betPart[i]=MCint
  # }
  # scnd2.2=(1-expit2)*betPart
  #
  # # scnd2.2=(1-expit2)*betPart3
  #
  # #cat("scnd2.1",scnd2.1,"\n")
  # #cat("scnd2.2",scnd2.2,"\n")
  #
  # scnd2=-log(delta)+log(scnd2.1+scnd2.2)
  # # print("scnd2.1:")
  # # print(scnd2.1[1:10])
  # # print("scnd2.2:")
  # # print(scnd2.2[1:10])
  # # scnd2=-log(delta)+scnd22.1+scnd22.2+log(1+betPart3)
  #
  # l2_i=scnd2

  #cat("sum(l2_i):",sum(l2_i),"\n")

  #### delete sum(l2_i)

  l=sum(l1_i)
  #print("theta:")
  #print(theta)
  # cat("logLikehood function value:", l,"\n")
  return(-l)
}


# setting interaction parameters to zero
objFunc=function(theta,dg1,dg2,setParam0){
  theta=theta*setParam0
  nL=negLogL(theta=theta,dg1=dg1,dg2=dg2)
  # cat("neg logL value:",nL,"\n")
  return(nL)
}


objFunc_contin=function(theta,dg1,setParam0){
  theta=theta*setParam0
  nL=negLogL_contin(theta=theta,dg1=dg1)
  # cat("neg logL value:",nL,"\n")
  return(nL)
}

#
## fuctions for calculating NIE,NDE,CDE
#


### EFF rewrite for contin


EFF=function(t,x1,x2,ctlM){
  b0=t[1]
  b1=t[2]
  b2=t[3]
  b3=t[4]
  b4=t[5]
  b5=t[6]
  a0=t[7]
  a1=t[8]
  g0=t[9]
  g1=t[10]

  X1a=a0+a1*x1
  X2a=a0+a1*x2
  X2b24=b2+b4*x2
  X2b15=b1+b5*x2

  X1g=g0+g1*x1
  X2g=g0+g1*x2

  if(X1a>=0)expitX1a=1/(1+exp(-X1a))
  if(X1a<0)expitX1a=exp(X1a)/(1+exp(X1a))

  if(X2a>=0)expitX2a=1/(1+exp(-X2a))
  if(X2a<0)expitX2a=exp(X2a)/(1+exp(X2a))

  if(X1g>=0)expitX1g=1/(1+exp(-X1g))
  if(X1g<0)expitX1g=exp(X1g)/(1+exp(X1g))

  if(X2g>=0)expitX2g=1/(1+exp(-X2g))
  if(X2g<0)expitX2g=exp(X2g)/(1+exp(X2g))

  #### delete X2b15 part

  NIE1=as.numeric(X2b15*(expitX2a-expitX1a)-X2b15*(expitX2g*expitX2a-expitX1g*expitX1a))
  NIE2=as.numeric(X2b24*(expitX1g-expitX2g))
  NIE=NIE1+NIE2

  NDE=as.numeric((b3+b4*expitX1g+b5*(1-expitX1g)*expitX1a)*(x2-x1))
  CDE=as.numeric((b3+b4*(ctlM>0)+b5*ctlM)*(x2-x1))

  list(NIE1=NIE1,NIE2=NIE2,NIE=NIE,NDE=NDE,CDE=CDE)
  }

EFF_contin=function(t,x1,x2,ctlM){
  b0=t[1]
  b1=t[2]
  b2=t[3]
  b3=t[4]
  b4=t[5]
  b5=t[6]
  a0=t[7]
  a1=t[8]
  # g0=t[9]
  # g1=t[10]

  X1a=a0+a1*x1
  X2a=a0+a1*x2
  X2b24=b2+b4*x2
  X2b15=b1+b5*x2

  # X1g=g0+g1*x1
  # X2g=g0+g1*x2

  if(X1a>=0)expitX1a=1/(1+exp(-X1a))
  if(X1a<0)expitX1a=exp(X1a)/(1+exp(X1a))

  if(X2a>=0)expitX2a=1/(1+exp(-X2a))
  if(X2a<0)expitX2a=exp(X2a)/(1+exp(X2a))

  # if(X1g>=0)expitX1g=1/(1+exp(-X1g))
  # if(X1g<0)expitX1g=exp(X1g)/(1+exp(X1g))
  #
  # if(X2g>=0)expitX2g=1/(1+exp(-X2g))
  # if(X2g<0)expitX2g=exp(X2g)/(1+exp(X2g))

  #### delete X2b15 part

  NIE1=as.numeric(X2b15*(expitX2a-expitX1a))
  NIE2=0
  NIE=NIE1

  NDE=as.numeric((b3+b5*expitX1a)*(x2-x1))
  CDE=as.numeric((b3+b4*(ctlM>0)+b5*ctlM)*(x2-x1))

  list(NIE1=NIE1,NIE2=NIE2,NIE=NIE,NDE=NDE,CDE=CDE)
}



EFFobject_contin=function(t,x1,x2,ctlM,setParam0,setParam0.a){
  t=t*setParam0
  t=t[setParam0.a==1]
  return(EFF_contin(t,x1=x1,x2=x2,ctlM=ctlM))
}

EFFobject=function(t,x1,x2,ctlM,setParam0,setParam0.a){
  t=t*setParam0
  t=t[setParam0.a==1]
  return(EFF(t,x1=x1,x2=x2,ctlM=ctlM))
}


#------------------------------------------------
#CI and estimators
#------------------------------------------------
est_CI=function(sigmaFull,est,x_from,x_to,M.obs,type1error,setParam0,setParam0.a){
  aveM.star=mean(M.obs[M.obs>0])

  vars=diag(sigmaFull)
  stderror=sqrt(vars)

  ci_lb=est-abs(qnorm(type1error/2))*stderror
  ci_ub=est+abs(qnorm(type1error/2))*stderror

  sigma=sigmaFull
  sigma[setParam0.a==0,setParam0.a==0]=0

  estDIeff=EFFobject(t=est,x1=x_from,x2=x_to,ctlM=aveM.star,
                     setParam0=setParam0,setParam0.a=setParam0.a)

  dNIE1.theta_hat= grad(function(t)EFFobject(t,x1=x_from,x2=x_to,
                                             ctlM=aveM.star,setParam0=setParam0,setParam0.a=setParam0.a)$NIE1,est)
  dNIE2.theta_hat= grad(function(t)EFFobject(t,x1=x_from,x2=x_to,
                                             ctlM=aveM.star,setParam0=setParam0,setParam0.a=setParam0.a)$NIE2,est)
  dNIE.theta_hat= grad(function(t)EFFobject(t,x1=x_from,x2=x_to,
                                            ctlM=aveM.star,setParam0=setParam0,setParam0.a=setParam0.a)$NIE,est)
  dNDE.theta_hat= grad(function(t)EFFobject(t,x1=x_from,x2=x_to,
                                            ctlM=aveM.star,setParam0=setParam0,setParam0.a=setParam0.a)$NDE,est)
  dCDE.theta_hat= grad(function(t)EFFobject(t,x1=x_from,x2=x_to,
                                            ctlM=aveM.star,setParam0=setParam0,setParam0.a=setParam0.a)$CDE,est)

  var.NIE1=dNIE1.theta_hat%*%sigma%*%dNIE1.theta_hat
  var.NIE2=dNIE2.theta_hat%*%sigma%*%dNIE2.theta_hat
  var.NIE=dNIE.theta_hat%*%sigma%*%dNIE.theta_hat
  var.NDE=dNDE.theta_hat%*%sigma%*%dNDE.theta_hat
  var.CDE=dCDE.theta_hat%*%sigma%*%dCDE.theta_hat

  stdiv.NIE1=sqrt(var.NIE1)
  stdiv.NIE2=sqrt(var.NIE2)
  stdiv.NIE=sqrt(var.NIE)
  stdiv.NDE=sqrt(var.NDE)
  stdiv.CDE=sqrt(var.CDE)

  pNIE1=2*(1-pnorm(abs(estDIeff$NIE1/stdiv.NIE1)))
  if(abs(estDIeff$NIE1)<10^(-5))pNIE1=1

  pNIE2=2*(1-pnorm(abs(estDIeff$NIE2/stdiv.NIE2)))
  if(abs(estDIeff$NIE2)<10^(-5))pNIE2=1

  pNIE=2*(1-pnorm(abs(estDIeff$NIE/stdiv.NIE)))
  if(abs(estDIeff$NIE)<10^(-5))pNIE=1

  ciNIE1.lb=estDIeff$NIE1-abs(qnorm(type1error/2))*stdiv.NIE1
  ciNIE1.ub=estDIeff$NIE1+abs(qnorm(type1error/2))*stdiv.NIE1
  ciNIE2.lb=estDIeff$NIE2-abs(qnorm(type1error/2))*stdiv.NIE2
  ciNIE2.ub=estDIeff$NIE2+abs(qnorm(type1error/2))*stdiv.NIE2

  ciNIE.lb=estDIeff$NIE-abs(qnorm(type1error/2))*stdiv.NIE
  ciNIE.ub=estDIeff$NIE+abs(qnorm(type1error/2))*stdiv.NIE
  ciNDE.lb=estDIeff$NDE-abs(qnorm(type1error/2))*stdiv.NDE
  ciNDE.ub=estDIeff$NDE+abs(qnorm(type1error/2))*stdiv.NDE
  ciCDE.lb=estDIeff$CDE-abs(qnorm(type1error/2))*stdiv.CDE
  ciCDE.ub=estDIeff$CDE+abs(qnorm(type1error/2))*stdiv.CDE

  estimators=c(unlist(estDIeff),est)
  names(estimators)=c("NIE1","NIE2","NIE","NDE","CDE","beta0","beta1","beta2",
                      "beta3","beta4","beta5","delta","alpha0","alpha1",
                      "xi","gamma0","gamma1")

  SE=c(stdiv.NIE1,stdiv.NIE2,stdiv.NIE,stdiv.NDE,stdiv.CDE,stderror)
  names(SE)=c("SE.NIE1","SE.NIE2","SE.NIE","SE.NDE","SE.CDE","SE.beta0","SE.beta1",
              "SE.beta2","SE.beta3","SE.beta4","SE.beta5","SE.delta","SE.alpha0",
              "SE.alpha1","SE.xi","SE.gamma0","SE.gamma1")


  CIlbounds=c(ciNIE1.lb, ciNIE2.lb,ciNIE.lb,ciNDE.lb,ciCDE.lb,ci_lb)
  names(CIlbounds)= c("lo.CI.NIE1","lo.CI.NIE2","lo.CI.NIE","lo.CI.NDE",
                      "lo.CI.CDE","lo.CI.beta0","lo.CI.beta1",
                      "lo.CI.beta2","lo.CI.beta3","lo.CI.beta4","lo.CI.beta5",
                      "lo.CI.delta","lo.CI.alpha0",
                      "lo.CI.alpha1","lo.CI.xi","lo.CI.gamma0","lo.CI.gamma1"
  )

  CIubounds=c(ciNIE1.ub, ciNIE2.ub, ciNIE.ub,ciNDE.ub,ciCDE.ub,ci_ub)
  names(CIubounds)= c("up.CI.NIE1","up.CI.NIE2","up.CI.NIE","up.CI.NDE",
                      "up.CI.CDE","up.CI.beta0","up.CI.beta1",
                      "up.CI.beta2","up.CI.beta3","up.CI.beta4","up.CI.beta5",
                      "up.CI.delta","up.CI.alpha0",
                      "up.CI.alpha1","up.CI.xi","up.CI.gamma0","up.CI.gamma1"
  )

  PvalNIE=c(pNIE1,pNIE2,pNIE)
  names(PvalNIE)=c("pNIE1","pNIE2","pNIE")

  ests_resul=c(PvalNIE,estimators,CIlbounds,CIubounds,SE)

  return(ests_resul)
}

est_CI_contin=function(sigmaFull,est,x_from,x_to,M.obs,type1error,setParam0,setParam0.a){
  aveM.star=mean(M.obs[M.obs>0])

  vars=diag(sigmaFull)
  stderror=sqrt(vars)

  ci_lb=est-abs(qnorm(type1error/2))*stderror
  ci_ub=est+abs(qnorm(type1error/2))*stderror

  sigma=sigmaFull
  sigma[setParam0.a==0,setParam0.a==0]=0

  estDIeff=EFFobject_contin(t=est,x1=x_from,x2=x_to,ctlM=aveM.star,
                     setParam0=setParam0,setParam0.a=setParam0.a)

  dNIE1.theta_hat= grad(function(t)EFFobject_contin(t,x1=x_from,x2=x_to,
          ctlM=aveM.star,setParam0=setParam0,setParam0.a=setParam0.a)$NIE1,est)
  dNIE2.theta_hat= grad(function(t)EFFobject_contin(t,x1=x_from,x2=x_to,
          ctlM=aveM.star,setParam0=setParam0,setParam0.a=setParam0.a)$NIE2,est)
  dNIE.theta_hat= grad(function(t)EFFobject_contin(t,x1=x_from,x2=x_to,
          ctlM=aveM.star,setParam0=setParam0,setParam0.a=setParam0.a)$NIE,est)
  dNDE.theta_hat= grad(function(t)EFFobject_contin(t,x1=x_from,x2=x_to,
          ctlM=aveM.star,setParam0=setParam0,setParam0.a=setParam0.a)$NDE,est)
  dCDE.theta_hat= grad(function(t)EFFobject_contin(t,x1=x_from,x2=x_to,
          ctlM=aveM.star,setParam0=setParam0,setParam0.a=setParam0.a)$CDE,est)

  var.NIE1=dNIE1.theta_hat%*%sigma%*%dNIE1.theta_hat
  var.NIE2=dNIE2.theta_hat%*%sigma%*%dNIE2.theta_hat
  var.NIE=dNIE.theta_hat%*%sigma%*%dNIE.theta_hat
  var.NDE=dNDE.theta_hat%*%sigma%*%dNDE.theta_hat
  var.CDE=dCDE.theta_hat%*%sigma%*%dCDE.theta_hat

  stdiv.NIE1=sqrt(var.NIE1)
  stdiv.NIE2=sqrt(var.NIE2)
  stdiv.NIE=sqrt(var.NIE)
  stdiv.NDE=sqrt(var.NDE)
  stdiv.CDE=sqrt(var.CDE)

  pNIE1=2*(1-pnorm(abs(estDIeff$NIE1/stdiv.NIE1)))
  if(abs(estDIeff$NIE1)<10^(-5))pNIE1=1

  pNIE2=2*(1-pnorm(abs(estDIeff$NIE2/stdiv.NIE2)))
  if(abs(estDIeff$NIE2)<10^(-5))pNIE2=1

  pNIE=2*(1-pnorm(abs(estDIeff$NIE/stdiv.NIE)))
  if(abs(estDIeff$NIE)<10^(-5))pNIE=1

  ciNIE1.lb=estDIeff$NIE1-abs(qnorm(type1error/2))*stdiv.NIE1
  ciNIE1.ub=estDIeff$NIE1+abs(qnorm(type1error/2))*stdiv.NIE1
  ciNIE2.lb=estDIeff$NIE2-abs(qnorm(type1error/2))*stdiv.NIE2
  ciNIE2.ub=estDIeff$NIE2+abs(qnorm(type1error/2))*stdiv.NIE2

  ciNIE.lb=estDIeff$NIE-abs(qnorm(type1error/2))*stdiv.NIE
  ciNIE.ub=estDIeff$NIE+abs(qnorm(type1error/2))*stdiv.NIE
  ciNDE.lb=estDIeff$NDE-abs(qnorm(type1error/2))*stdiv.NDE
  ciNDE.ub=estDIeff$NDE+abs(qnorm(type1error/2))*stdiv.NDE
  ciCDE.lb=estDIeff$CDE-abs(qnorm(type1error/2))*stdiv.CDE
  ciCDE.ub=estDIeff$CDE+abs(qnorm(type1error/2))*stdiv.CDE

  estimators=c(unlist(estDIeff),est)
  names(estimators)=c("NIE1","NIE2","NIE","NDE","CDE","beta0","beta1","beta2",
                      "beta3","beta4","beta5","delta","alpha0","alpha1",
                      "xi")

  SE=c(stdiv.NIE1,stdiv.NIE2,stdiv.NIE,stdiv.NDE,stdiv.CDE,stderror)
  names(SE)=c("SE.NIE1","SE.NIE2","SE.NIE","SE.NDE","SE.CDE","SE.beta0","SE.beta1",
              "SE.beta2","SE.beta3","SE.beta4","SE.beta5","SE.delta","SE.alpha0",
              "SE.alpha1","SE.xi")


  CIlbounds=c(ciNIE1.lb, ciNIE2.lb,ciNIE.lb,ciNDE.lb,ciCDE.lb,ci_lb)
  names(CIlbounds)= c("lo.CI.NIE1","lo.CI.NIE2","lo.CI.NIE","lo.CI.NDE",
                    "lo.CI.CDE","lo.CI.beta0","lo.CI.beta1",
                    "lo.CI.beta2","lo.CI.beta3","lo.CI.beta4","lo.CI.beta5",
                    "lo.CI.delta","lo.CI.alpha0",
                    "lo.CI.alpha1","lo.CI.xi"
                    )

  CIubounds=c(ciNIE1.ub, ciNIE2.ub, ciNIE.ub,ciNDE.ub,ciCDE.ub,ci_ub)
  names(CIubounds)= c("up.CI.NIE1","up.CI.NIE2","up.CI.NIE","up.CI.NDE",
                    "up.CI.CDE","up.CI.beta0","up.CI.beta1",
                    "up.CI.beta2","up.CI.beta3","up.CI.beta4","up.CI.beta5",
                    "up.CI.delta","up.CI.alpha0",
                    "up.CI.alpha1","up.CI.xi"
                    )

  PvalNIE=c(pNIE1,pNIE2,pNIE)
  names(PvalNIE)=c("pNIE1","pNIE2","pNIE")

  ests_resul=c(PvalNIE,estimators,CIlbounds,CIubounds,SE)

  return(ests_resul)
}

