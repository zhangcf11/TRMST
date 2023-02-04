#' Title: the function of time-dependent RMST regression
#' @importFrom survival coxph survfit Surv
#' @importFrom stats coef as.formula
#' @importFrom geepack geeglm
#' @param data Raw data set
#' @param time the name of the variable in data that indicates survival time
#' @param tstop the name of the variable in data that indicates the right  endpoint of the covariate observation time interval
#' @param id the name of the variable in data that identifies the different subjects
#' @param status the name of the variable in data that indicates  survival  status
#' @param event the name of the variable in date that indicates survival status in the covariate observation time interval
#' @param fixed the fixed covariates
#' @param varing the time-dependent covariates
#' @param tau the value to specify the truncation time point for the RMST,if it is null, indicates it is the largest follow-up time
#' @usage T_RMST(data,time,tstop,id,status,event,fixed,varing,tau)
#' @return An object of type "geeglm"
#' @export
#' @references Chengfeng Zhang,Baoyi Huang,Hongji Wu,Hao Yuan,Yawen Hou*,
#' Zheng Chen*.Restricted mean survival time regression model with time-dependent covariates.
#' Statistics in Medicine.2022;41(21):4081–4090.
#' @examples
#' library(TRMST)
#' data(heartdata)
#' fixed<-c("age1","age2","year","surgery")
#' varing<-c("transplant")
#' TRMST<-T_RMST(data=heartdata,time="time",tstop="tstop",id="id",status="status",event="event",fixed,varing,tau=NULL)
#' summary(TRMST)
T_RMST<-function(data,time,tstop,id,status,event,fixed,varing,tau=NULL){
  #check to see if the package is already installed
  if (is.element("survival", installed.packages()[,1])==FALSE){
    install.packages("survival")
  }
  require(survival) #same as library statement
  if (is.element("geepack", installed.packages()[,1])==FALSE){
    install.packages("geepack")
  }
  require(geepack)
  if (is.null(tau)){
    data$time0<-data[[time]]
    data$status0<-data[[status]]
    data$tstop0<-data[[tstop]]
    data$event0<-data[[event]]
  }
  else{
    data$time0<-pmin(data[[time]],tau)
    data$status0<-ifelse(data$time0==tau,0,data[[status]])
    data$tstop0<-pmin(data[[tstop]],tau)
    data$event0<-ifelse(data$tstop0==tau,0,data[[event]])
  }

  ##Adds a column of censoring status
  Censor<-function(data,id,status)
  {
    data$z<-1:length(data[[id]])
    data1<-data[!duplicated(data[[id]],fromLast = TRUE),]
    data1$cens<-ifelse(data1[[status]]==0,1,0)
    data$cens<-0
    data2<-rbind(data,data1)
    data2<-data2[order(data2[[id]]),]
    data2<-data2[!duplicated(data2$z,fromLast = TRUE),]
    return(data2)
  }
  Dat<-Censor(data,id,"status0")
  basecov1<-paste(fixed,collapse="+")
  formular1<-paste("Surv(","time0",",","cens",")","~",basecov1)
  formular1<-as.formula(formular1)
  CZ <- coxph(formular1,data=Dat)
  betaCHat <- matrix(coef(CZ), ncol=1)
  fc<-as.data.frame(matrix(0,1,length(fixed)))
  colnames(fc)<-fixed
  baseline <- summary(survfit(CZ,newdata=fc))
  BaseCumHaz <- data.frame(cbind(baseline$time, baseline$cumhaz))
  colnames(BaseCumHaz) <- c("time", "cumhaz")
  BaseCumHaz <- BaseCumHaz[order(BaseCumHaz$time, decreasing = T),] # needs descending order
  if(min(BaseCumHaz[,"time"])>0) {BaseCumHaz <- rbind(BaseCumHaz, c(0,0))}
  WYHat <- sapply(Dat$time0,  function(t){BaseCumHaz[which(BaseCumHaz[,"time"] <= t)[1], "cumhaz"]})
  WYHat <- exp(WYHat*exp(as.matrix(Dat[,fixed]) %*% betaCHat))

  basecov2<-paste(varing,collapse="+")
  formular2<-paste("Surv(","time0",",","cens",")","~",basecov2)
  formular2<-as.formula(formular2)
  CZ1<- coxph(formular2, data=Dat)
  betaCHat1 <- matrix(coef(CZ1), ncol=1)
  vc<-as.data.frame(matrix(0,1,length(varing)))
  colnames(vc)<-varing
  baseline1 <- summary(survfit(CZ1, newdata=vc))
  BaseCumHaz1 <- data.frame(cbind(baseline1$time, baseline1$cumhaz))
  colnames(BaseCumHaz1) <- c("time", "cumhaz")
  BaseCumHaz1 <- BaseCumHaz1[order(BaseCumHaz1$time, decreasing = T),] # needs descending order
  if(min(BaseCumHaz1[,"time"])>0) {BaseCumHaz1 <- rbind(BaseCumHaz1, c(0,0))}
  WYHat1 <- sapply(Dat$tstop0, function(t){BaseCumHaz1[which(BaseCumHaz1[,"time"] <= t)[1], "cumhaz"]})
  WYHat1 <- exp(WYHat1*exp(as.matrix(Dat[,varing]) %*% betaCHat1))

  Weights<-WYHat*WYHat1
  Dat$weight<-Weights
  Dat$id<-Dat[[id]]
  data1<-Dat[which(Dat$event0==1),]
  formular3<-paste("time0","~",paste(c(fixed,varing),collapse="+"))
  formular3<-as.formula(formular3)
  TRMST<-geeglm(formular3,data=data1,weights=weight,id=id)
  return(TRMST)
}
