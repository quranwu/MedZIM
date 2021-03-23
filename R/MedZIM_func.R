##' @description {
##' The Outcome model is:
##' \deqn{Y=\beta0+\beta1M+\beta2 1_{(M>0)}+\beta3X+\beta4X1_{(M>0)}+\beta5XM+\epsilon}
##' }
##' @title Perform MedZIM model
##' @param dat The dataset used for the model.
##' @param xVar The name of X (independent variable) in the data.
##' @param yVar The name of Y (outcome variable) in the data.
##' @param taxon_name The string that can identify the taxon variables.
##' @param libSize_name The name of library size variable in the data.
##' @param obs_gt_0 The at least number of non-zero observations that a taxon should have to be included into the analysis.
##' @param obs_eq_0 The at least number of zero observations that a taxon should have to be included in group 1.
##' @param inter_x_mg0 Whether to include the \eqn{\beta4} term into the model.
##' @param inter_x_m Whether to include the \eqn{\beta5} term into the model.
##' @param eval.max Maximum number of evaluations of objective functions allowed in nlminb function.
##' @param iter.max Maximum number of iterations allowed in nlminb function.
##' @param x_from The value X starts from in evaluating the effect of X.
##' @param x_to The value X ends to in evaluating the effect of X.
##' @param type1error The type I error.
##' @param paraJobs The number of cores for parallel computing.
##' @return A list contains all results.
##' \item{fullList}{A list contains resutls for each taxon. See details.}
##' \item{validTaxaNames}{The taxons that were in group 1, which had enough zero and non-zero observations.}
##' \item{continuTaxa}{The taxons that were in group 2, which only had non-zero observations.}
##'
##' @details The output list contains information of all taxons to be analyzed for mediation effects. Each taxon result is also
##' a list, with first object contains p value, estimates, and CI of NIE1, NIE2, NIE, and each parameters.
##'
##' @examples
##' MedZIM_func(dat=data_ZIM, xVar = "x",yVar = "y1",
##' taxon_name = "taxon",libSize_name = "libSize",paraJobs=2)
##'
##'
##' @import parallel doParallel pracma betareg foreach stats
##' @export
MedZIM_func<-function(dat,xVar,yVar,taxon_name,libSize_name,obs_gt_0=2,obs_eq_0=2,
                   inter_x_mg0=T,inter_x_m=F,eval.max=200,iter.max=200,x_from=0,
                   x_to=1,type1error=0.05,paraJobs=4) {
  # xVar<-xVar
  # yVar<-yVar
  data<-dat
  libSize<-libSize_name
  Mprefix<-taxon_name

  setParam0=c(1,1,rep(1,2),inter_x_mg0,inter_x_m,rep(1,6)) # for model with continuous M

  # setParam_hess_contin<-c(1,1,rep(1,3),0,rep(1,4),0,0)

  setParam0.contin=c(rep(1,2),0,1,0,inter_x_m,rep(1,4)) # identify non-zero beta in Y model for M without 0

  setParam0.a=c(rep(1,6),0,rep(1,2),0,rep(1,2)) # identify elements for NIE calculations

  setParam0.a_contin<-c(rep(1,6),0,rep(1,2),0) # for NIE calculation with no gamma
  # environment(ZIM_execute)<-environment()
  #
  # res<-ZIM_execute()
  # cat(xVar)
  time0=proc.time()[3]
  # Mprefix<-taxon_name
  # get taxa variable names
  MVarNamLength=nchar(Mprefix)
  micros = sapply(substr(colnames(data),1,MVarNamLength), function(x) {grep(Mprefix, x)})
  microPositions=which(micros == 1)
  rm(micros)
  sampleSiz=nrow(data)
  nTaxa=length(microPositions)

  allTaxaNames=colnames(data)[microPositions]

  # convert relative abundance
  data[,allTaxaNames]=data[,allTaxaNames]/rowSums(data[,allTaxaNames])
  # print(rowSums(data[,allTaxaNames]))

  allTaxaMat=data[,allTaxaNames]

  data[,allTaxaNames]=allTaxaMat/rowSums(allTaxaMat)
  #sum(allTaxaMat[1,])
  #save(data,file="real_data_NH.RData")

  validTaxaNames=allTaxaNames[(colSums(allTaxaMat>0)<(nrow(data)-obs_eq_0))&(colSums(allTaxaMat>0)>=obs_gt_0)]

  #validTaxaNames=validTaxaNames[1:4]
  # print(validTaxaNames)

  continuTaxa=allTaxaNames[colSums(allTaxaMat>0)>=(nrow(data)-obs_eq_0)]
  #continuTaxa=c()

  allAnalyTaxaNames=c(validTaxaNames,continuTaxa)

  notAnalyTaxa=allTaxaNames[colSums(allTaxaMat>0)<obs_gt_0]

  allFunc=c("dataByZero","ini.bound","hi","hii","hiSimp","hi3",
            "negLogL","negLogL_contin","objFunc","objFunc_contin","EFF","EFF_contin","EFFobject",
            "EFFobject_contin","est_CI","est_CI_contin")

  if(length(validTaxaNames)>=1){
    c1<-makeCluster(paraJobs)

    clusterExport(c1, allFunc, envir = environment())
    registerDoParallel(c1)

    # dataForOneM=dataByZero(data=data,x=xVar,y=yVar,L=libSize,m="taxon1")
    # #dataForOneM
    #
    # X=dataForOneM$X;Y=dataForOneM$Y;M.obs=dataForOneM$M.obs;indM=dataForOneM$indM;
    # L=dataForOneM$L;d1=dataForOneM$d1;d2=dataForOneM$d2
    #
    # iniV=ini.bound(X,Y,M.obs,indM,L)

    # cat("validTaxaNames:",validTaxaNames,"\n")

    binaryInter=foreach(jj=1:length(validTaxaNames),.multicombine=T,
                        .packages=c("pracma","betareg"),
                        .errorhandling="pass") %dopar% {

                          # print(xVar)
                          taxon=validTaxaNames[jj]
                          dataForOneM=dataByZero(data=data,x=xVar,y=yVar,L=libSize,m=taxon)
                          #dataForOneM

                          X=dataForOneM$X;Y=dataForOneM$Y;M.obs=dataForOneM$M.obs;indM=dataForOneM$indM;
                          L=dataForOneM$L;d1=dataForOneM$d1;d2=dataForOneM$d2

                          iniV=ini.bound(X,Y,M.obs,indM,L)

                          initial=c(iniV$beta0.ini,iniV$beta1.ini,iniV$beta2.ini,iniV$beta3.ini,
                                    iniV$beta4.ini,iniV$beta5.ini,
                                    iniV$delta.ini,iniV$alpha0.ini,iniV$alpha1.ini,iniV$xi.ini,
                                    iniV$gamma0.ini,iniV$gamma1.ini)
                          LBsnlminb=c(iniV$LBbeta0.ini,iniV$LBbeta1.ini,iniV$LBbeta2.ini,iniV$LBbeta3.ini,
                                      iniV$LBbeta4.ini,iniV$LBbeta5.ini,iniV$LBdelta.ini,
                                      iniV$LBalpha0.ini,
                                      iniV$LBalpha1.ini,iniV$LBxi.ini,iniV$LBgamma0.ini,
                                      iniV$LBgamma1.ini)
                          UBsnlminb=c(iniV$UBbeta0.ini,iniV$UBbeta1.ini,iniV$UBbeta2.ini,iniV$UBbeta3.ini,
                                      iniV$UBbeta4.ini,iniV$UBbeta5.ini,iniV$UBdelta.ini,
                                      iniV$UBalpha0.ini,iniV$UBalpha1.ini,iniV$UBxi.ini,
                                      iniV$UBgamma0.ini,iniV$UBgamma1.ini)

                          time1=proc.time()[3]

                          initial=initial*setParam0
                          LBsnlminb=LBsnlminb*setParam0
                          UBsnlminb=UBsnlminb*setParam0

                          # source("VSLfunctions_LOD_simu.R")
                          # print(initial)

                          est.val=nlminb(initial,objFunc,
                                         lower=LBsnlminb,
                                         upper=UBsnlminb,
                                         control=list(eval.max=eval.max,iter.max=iter.max),
                                         dg1=d1,dg2=d2,setParam0=setParam0)

                          print(est.val)

                          time2=proc.time()[3]
                          # cat("Optimization time:",(time2-time1)/60,"mins","\n")
                          iter=est.val$iterations
                          # iter
                          evalu=est.val$evaluations
                          # evalu
                          if(est.val$convergence==0)success="MLE successfully converged."
                          if(est.val$convergence==1)success="MLE Not successfully converged."

                          est=est.val$par
                          names(est)=c("beta0","beta1","beta2","beta3","beta4","beta5","delta",
                                       "alpha0","alpha1","xi","gamma0","gamma1")

                          # cat("est=",est,"\n")

                          names(initial)= names(est)
                          # cat("initial=",initial,"\n")

                          est_ini=cbind(est,initial,LBsnlminb,UBsnlminb)
                          colnames(est_ini)=c("est","initial","lowBound","upBound")

                          # save a csv file for the estimates and initial values
                          # write.csv(est_ini, paste0("./results/est_LOD_",start.sample,"_",taxon,".csv"))

                          # start calculae hessian marix
                          hessFull <- hessian(objFunc,est,dg1=d1,dg2=d2,setParam0=setParam0)

                          # hessFull=diag(rep(1,12))+rep(1,12)%*%t(rep(1,12))

                          colnames(hessFull)=names(est)
                          rownames(hessFull)=colnames(hessFull)

                          time3=proc.time()[3]
                          # cat("Hessian matrix calculation time:",(time3-time2)/60,"mins","\n")

                          #  write.csv(hessFull, paste0("./results/hess_LOD_",start.sample,"_",taxon,".csv"))

                          hess=hessFull[setParam0==1,setParam0==1]
                          ei=eigen(hess)
                          hessianPosDef=1
                          if(length(which(ei$values<=0))>0){
                            hessianPosDef=0
                            ei$values[ei$values<=0]=10^(-3)
                          }
                          sigma=ei$vectors%*%diag(1/ei$values)%*%solve(ei$vectors)

                          #
                          ## transform the sigma to full sigma matrix when there parameters set to be 0
                          #
                          sigmaFull=sigma
                          if(is.element(0,setParam0)){
                            zeroSet=which(setParam0==0)
                            for (i in 1:length(zeroSet)){
                              sigmaTemp=matrix(NA,nrow=(nrow(sigmaFull)+1),ncol=(ncol(sigmaFull)+1))
                              nrowTemp=nrow(sigmaFull)+1
                              nrowFull=nrow(sigmaFull)

                              if(zeroSet[i]==1){
                                sigmaTemp[zeroSet[i],]=0
                                sigmaTemp[,zeroSet[i]]=0
                                sigmaTemp[(zeroSet[i]+1):nrowTemp,(zeroSet[i]+1):nrowTemp]=sigmaFull[zeroSet[i]:nrowFull,zeroSet[i]:nrowFull]
                                sigmaFull=sigmaTemp
                              }
                              if(zeroSet[i]==nrowTemp){
                                sigmaTemp[zeroSet[i],]=0
                                sigmaTemp[,zeroSet[i]]=0
                                sigmaTemp[1:nrowFull,1:nrowFull]=sigmaFull[1:nrowFull,1:nrowFull]
                                sigmaFull=sigmaTemp
                              }
                              if(zeroSet[i]>1 & zeroSet[i]<nrowTemp){
                                sigmaTemp[zeroSet[i],]=0
                                sigmaTemp[,zeroSet[i]]=0
                                sigmaTemp[(zeroSet[i]+1):nrowTemp,(zeroSet[i]+1):nrowTemp]=sigmaFull[zeroSet[i]:nrowFull,zeroSet[i]:nrowFull]
                                sigmaTemp[1:(zeroSet[i]-1),(zeroSet[i]+1):nrowTemp]=sigmaFull[1:(zeroSet[i]-1),zeroSet[i]:nrowFull]
                                sigmaTemp[(zeroSet[i]+1):nrowTemp,1:(zeroSet[i]-1)]=sigmaFull[zeroSet[i]:nrowFull,1:(zeroSet[i]-1)]
                                sigmaTemp[1:(zeroSet[i]-1),1:(zeroSet[i]-1)]=sigmaFull[1:(zeroSet[i]-1),1:(zeroSet[i]-1)]
                                sigmaFull=sigmaTemp
                              }
                            }
                          }

                          colnames(sigmaFull)=colnames(hessFull)
                          rownames(sigmaFull)=colnames(sigmaFull)

                          # write.csv(sigmaFull, paste0("./results/sigma_LOD_",start.sample,"_",taxon,".csv"))

                          ests_ci=est_CI(sigma=sigmaFull,est=est,x_from=x_from,x_to=x_to,
                                         M.obs=M.obs,type1error=type1error,
                                         setParam0=setParam0,setParam0.a=setParam0.a)

                          # save a csv file for the estimates and CIs of the NIE,NDE and CDE
                          # write.csv(as.matrix(ests_ci),paste0("./results/NIE_LOD_",start.sample,"_",taxon,".csv"))

                          # fullList[[taxon]]=list(ests_ci,est_ini,success,hessianPosDef,taxon,iter,
                          #                    evalu,hessFull,sigmaFull)

                          return(list(ests_ci,est_ini,success,hessianPosDef,taxon,iter,evalu,
                                      hessFull,sigmaFull))
                        }

    stopCluster(c1)
    names(binaryInter)=validTaxaNames
    fullList=binaryInter
  }

  cat("continuTaxa:",continuTaxa,"\n")

  if(length(continuTaxa)>=1){

    c2<-makeCluster(paraJobs)

    clusterExport(c2,allFunc,envir = environment())
    registerDoParallel(c2)

    continuInter=foreach(ii=1:length(continuTaxa),.multicombine=T,
                         .packages=c("pracma","betareg"),
                         .errorhandling="pass") %dopar% {
                           # for(ii in 2){
                           #   ii=2
                           taxon=continuTaxa[ii]
                           dataForOneM=dataByZero(data=data,x=xVar,y=yVar,L=libSize,m=taxon)
                           #dataForOneM

                           X=dataForOneM$X;Y=dataForOneM$Y;M.obs=dataForOneM$M.obs;indM=dataForOneM$indM;
                           L=dataForOneM$L;d1=dataForOneM$d1;d2=dataForOneM$d2

                           iniV=ini.bound(X,Y,M.obs,indM,L)

                           initial=c(iniV$beta0.ini,iniV$beta1.ini,iniV$beta2.ini,iniV$beta3.ini,
                                     iniV$beta4.ini,iniV$beta5.ini,
                                     iniV$delta.ini,iniV$alpha0.ini,iniV$alpha1.ini,iniV$xi.ini,
                                     iniV$gamma0.ini,iniV$gamma1.ini)
                           LBsnlminb=c(iniV$LBbeta0.ini,iniV$LBbeta1.ini,iniV$LBbeta2.ini,iniV$LBbeta3.ini,
                                       iniV$LBbeta4.ini,iniV$LBbeta5.ini,iniV$LBdelta.ini,
                                       iniV$LBalpha0.ini,
                                       iniV$LBalpha1.ini,iniV$LBxi.ini,iniV$LBgamma0.ini,
                                       iniV$LBgamma1.ini)
                           UBsnlminb=c(iniV$UBbeta0.ini,iniV$UBbeta1.ini,iniV$UBbeta2.ini,iniV$UBbeta3.ini,
                                       iniV$UBbeta4.ini,iniV$UBbeta5.ini,iniV$UBdelta.ini,
                                       iniV$UBalpha0.ini,iniV$UBalpha1.ini,iniV$UBxi.ini,
                                       iniV$UBgamma0.ini,iniV$UBgamma1.ini)

                           time1=proc.time()[3]

                           initial=(initial*setParam0.contin)[1:10]
                           LBsnlminb=(LBsnlminb*setParam0.contin)[1:10]
                           UBsnlminb=(UBsnlminb*setParam0.contin)[1:10]

                           LBsnlminb[is.nan(LBsnlminb)]<-0
                           UBsnlminb[is.nan(UBsnlminb)]<-0
                           # source("VSLfunctions_LOD_simu.R")
                           # print(initial)

                           est.val=nlminb(initial,objFunc_contin,
                                          lower=LBsnlminb,
                                          upper=UBsnlminb,
                                          control=list(eval.max=eval.max,iter.max=iter.max),
                                          dg1=d1,setParam0=setParam0.contin[1:10])

                           print(est.val)

                           time2=proc.time()[3]
                           # cat("Optimization time:",(time2-time1)/60,"mins","\n")
                           iter=est.val$iterations
                           # iter
                           evalu=est.val$evaluations
                           # evalu
                           if(est.val$convergence==0)success="MLE successfully converged."
                           if(est.val$convergence==1)success="MLE Not successfully converged."

                           est=est.val$par
                           names(est)=c("beta0","beta1","beta2","beta3","beta4","beta5","delta",
                                        "alpha0","alpha1","xi")

                           # cat("est=",est,"\n")

                           names(initial)= names(est)
                           # cat("initial=",initial,"\n")

                           est_ini=cbind(est,initial,LBsnlminb,UBsnlminb)
                           colnames(est_ini)=c("est","initial","lowBound","upBound")

                           # save a csv file for the estimates and initial values
                           # write.csv(est_ini, paste0("./results/est_LOD_",start.sample,"_",taxon,".csv"))

                           # start calculae hessian marix
                           hessFull <- hessian(objFunc_contin,est,dg1=d1,setParam0=setParam0.contin[1:10])

                           # hessFull=diag(rep(1,12))+rep(1,12)%*%t(rep(1,12))

                           colnames(hessFull)=names(est)
                           rownames(hessFull)=colnames(hessFull)

                           time3=proc.time()[3]
                           # cat("Hessian matrix calculation time:",(time3-time2)/60,"mins","\n")

                           #  write.csv(hessFull, paste0("./results/hess_LOD_",start.sample,"_",taxon,".csv"))

                           hess=hessFull[setParam0.contin[1:10]==1,setParam0.contin[1:10]==1]
                           ei=eigen(hess)
                           hessianPosDef=1
                           if(length(which(ei$values<=0))>0){
                             hessianPosDef=0
                             ei$values[ei$values<=0]=10^(-3)
                           }
                           sigma=ei$vectors%*%diag(1/ei$values)%*%solve(ei$vectors)

                           #
                           ## transform the sigma to full sigma matrix when there parameters set to be 0
                           #
                           sigmaFull=sigma
                           if(is.element(0,setParam0.contin)){
                             zeroSet=which(setParam0.contin==0)
                             for (i in 1:length(zeroSet)){
                               sigmaTemp=matrix(NA,nrow=(nrow(sigmaFull)+1),ncol=(ncol(sigmaFull)+1))
                               nrowTemp=nrow(sigmaFull)+1
                               nrowFull=nrow(sigmaFull)

                               if(zeroSet[i]==1){
                                 sigmaTemp[zeroSet[i],]=0
                                 sigmaTemp[,zeroSet[i]]=0
                                 sigmaTemp[(zeroSet[i]+1):nrowTemp,(zeroSet[i]+1):nrowTemp]=sigmaFull[zeroSet[i]:nrowFull,zeroSet[i]:nrowFull]
                                 sigmaFull=sigmaTemp
                               }
                               if(zeroSet[i]==nrowTemp){
                                 sigmaTemp[zeroSet[i],]=0
                                 sigmaTemp[,zeroSet[i]]=0
                                 sigmaTemp[1:nrowFull,1:nrowFull]=sigmaFull[1:nrowFull,1:nrowFull]
                                 sigmaFull=sigmaTemp
                               }
                               if(zeroSet[i]>1 & zeroSet[i]<nrowTemp){
                                 sigmaTemp[zeroSet[i],]=0
                                 sigmaTemp[,zeroSet[i]]=0
                                 sigmaTemp[(zeroSet[i]+1):nrowTemp,(zeroSet[i]+1):nrowTemp]=sigmaFull[zeroSet[i]:nrowFull,zeroSet[i]:nrowFull]
                                 sigmaTemp[1:(zeroSet[i]-1),(zeroSet[i]+1):nrowTemp]=sigmaFull[1:(zeroSet[i]-1),zeroSet[i]:nrowFull]
                                 sigmaTemp[(zeroSet[i]+1):nrowTemp,1:(zeroSet[i]-1)]=sigmaFull[zeroSet[i]:nrowFull,1:(zeroSet[i]-1)]
                                 sigmaTemp[1:(zeroSet[i]-1),1:(zeroSet[i]-1)]=sigmaFull[1:(zeroSet[i]-1),1:(zeroSet[i]-1)]
                                 sigmaFull=sigmaTemp
                               }
                             }
                           }

                           colnames(sigmaFull)=colnames(hessFull)
                           rownames(sigmaFull)=colnames(sigmaFull)

                           # write.csv(sigmaFull, paste0("./results/sigma_LOD_",start.sample,"_",taxon,".csv"))

                           ests_ci=est_CI_contin(sigma=sigmaFull,est=est,x_from=x_from,x_to=x_to,
                                                 M.obs=M.obs,type1error=type1error,
                                                 setParam0=setParam0.contin,setParam0.a=setParam0.a_contin)

                           # save a csv file for the estimates and CIs of the NIE,NDE and CDE
                           # write.csv(as.matrix(ests_ci),paste0("./results/NIE_LOD_",start.sample,"_",taxon,".csv"))

                           # fullList[[taxon]]=list(ests_ci,est_ini,success,hessianPosDef,taxon,iter,
                           #                    evalu,hessFull,sigmaFull)

                           return(list(ests_ci,est_ini,success,hessianPosDef,taxon,iter,evalu,
                                       hessFull,sigmaFull))
                         }

    stopCluster(c2)
    names(continuInter)=continuTaxa
    fullList=do.call(c, list(fullList, continuInter))
  }

  time4=proc.time()[3]
  totTime=(time4-time0)/60
  cat("Entire analysis took:",(time4-time0)/60,"mins","\n")

  # save a .Rdata file including all results
  # output.data=paste("./realResultsNH/Realresults_NH_Dir",".RData",sep="")

  finalList=list(fullList=fullList,
                 totTime=totTime,
                 setParam0=setParam0,
                 setParam0.a=setParam0.a,
                 setParam0.contin=setParam0.contin,
                 type1error=type1error,
                 allTaxaNames=allTaxaNames,
                 allAnalyTaxaNames=allAnalyTaxaNames,
                 validTaxaNames=validTaxaNames,
                 continuTaxa=continuTaxa,
                 notAnalyTaxa=notAnalyTaxa,
                 nTaxa=nTaxa,
                 sampleSiz=sampleSiz)
  return(finalList)
}
