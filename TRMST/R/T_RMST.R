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
#' @param fixed the fixed covariates parameters (including the name of fixed covariates and reference group)
#' @param varing the time-dependent covariates parameters
#'
#' @return the result of time-dependent RMST
#' @export
#'
#' @examples
#' library(TRMST)
#' library(survival)
#' library(geepack)
#' data(heartdata)
#' fixed<-list(ZC=c("age1","age2","year","surgery"),ref=c(age1=0,age2=0,year=0,surgery=0))
#' varing<-list(ZC=c("transplant"),ref=c(transplant=0))
#' T_RMST(heartdata,time="time",tstop="tstop",id="id",status="status",event="event",fixed,varing)
T_RMST<-function(data,time,tstop,id,status,event,fixed,varing){
  #check to see if the package is already installed
  if (is.element("survival", installed.packages()[,1])==FALSE){
    install.packages("survival")
  }
  require(survival) #same as library statement
  if (is.element("geepack", installed.packages()[,1])==FALSE){
    install.packages("geepack")
  }
  require(geepack)
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
  Dat<-Censor(data,id,status)
  basecov1<-paste(fixed$ZC,collapse="+")
  formular1<-paste("Surv(",time,",","cens",")","~",basecov1)
  formular1<-as.formula(formular1)
  CZ <- coxph(formular1,data=Dat)
  betaCHat <- matrix(coef(CZ), ncol=1)
  baseline <- summary(survfit(CZ,newdata=fixed$ref))
  BaseCumHaz <- data.frame(cbind(baseline$time, baseline$cumhaz))
  colnames(BaseCumHaz) <- c("time", "cumhaz")
  BaseCumHaz <- BaseCumHaz[order(BaseCumHaz$time, decreasing = T),] # needs descending order
  if(min(BaseCumHaz[,"time"])>0) {BaseCumHaz <- rbind(BaseCumHaz, c(0,0))}
  WYHat <- sapply(Dat[[time]],  function(t){BaseCumHaz[which(BaseCumHaz[,"time"] <= t)[1], "cumhaz"]})
  WYHat <- exp(WYHat*exp(as.matrix(Dat[,fixed$ZC]) %*% betaCHat))

  basecov2<-paste(varing$ZC,collapse="+")
  formular2<-paste("Surv(",time,",","cens",")","~",basecov2)
  formular2<-as.formula(formular2)
  CZ1<- coxph(formular2, data=Dat)
  betaCHat1 <- matrix(coef(CZ1), ncol=1)
  baseline1 <- summary(survfit(CZ1, newdata=varing$ref))
  BaseCumHaz1 <- data.frame(cbind(baseline1$time, baseline1$cumhaz))
  colnames(BaseCumHaz1) <- c("time", "cumhaz")
  BaseCumHaz1 <- BaseCumHaz1[order(BaseCumHaz1$time, decreasing = T),] # needs descending order
  if(min(BaseCumHaz1[,"time"])>0) {BaseCumHaz1 <- rbind(BaseCumHaz1, c(0,0))}
  WYHat1 <- sapply(Dat[[tstop]], function(t){BaseCumHaz1[which(BaseCumHaz1[,"time"] <= t)[1], "cumhaz"]})
  WYHat1 <- exp(WYHat1*exp(as.matrix(Dat[,varing$ZC]) %*% betaCHat1))

  Weights<-WYHat*WYHat1
  Dat$weight<-Weights
  data1<-Dat[which(Dat[[event]]==1),]
  formular3<-paste(time,"~",paste(c(fixed$ZC,varing$ZC),collapse="+"))
  formular3<-as.formula(formular3)
  fit<-geeglm(formular3,data=data1,weights=weight,id=id)
  return(fit)
}
