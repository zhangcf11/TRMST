#' @title Data from the Stanford Heart Transplant Center
#' @docType data
#' @keywords datasets
#' @name heartdata
#' @usage heartdata
#' @format A numeric dataframe
#' \describe{
#' \item{status}{the name of the variable in data that indicates  survival  status}
#' \item{transplant}{the time-dependent covariate,coded as 1 for undergo heart transplant in follow-up time and 0 for none}
#' \item{id}{the name of the variable in data that identifies the different subjects}
#' \item{tstart}{the name of the variable in data that indicates the left  endpoint of the covariate observation time interval}
#' \item{tstop}{the name of the variable in data that indicates the right  endpoint of the covariate observation time interval}
#' \item{event}{the name of the variable in date that indicates survival status in the covariate observation time interval}
#' \item{age1}{the Dummy variable of age, coded as 1 for 45<=age <60 and 0 for else}
#' \item{age2}{the Dummy variable of age, coded as 1 for age>=60 and 0 for else}
#' \item{year}{enrollment time of subjects}
#' \item{surgery}{whether the patient had had undergone heart bypass surgery, coded as 1 for undergone and 0 for none}
#' }
#' @source the data is originally from the package named survival and transformed the data into a time-dependent covariate
#' @examples
#' library(TRMST)
#' data(heartdata)
#' str(heartdata)
NULL
