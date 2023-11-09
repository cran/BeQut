#' dataLong
#'
#' 'dataLong' is a dataset simulated from a joint model for longitudinal and time-to-event data.
#' This dataset is used to illustrate both 'lqmm' and 'qrjm' functions.
#'
#' @usage data(dataLong)
#'
#' @name dataLong
#' @docType data
#'
#' @format A \code{data.frame} with 1562 observations from 300 subjects. The columns are:
#' \describe{
#'   \item{ID}{integer: number for patient identification.}
#'   \item{visit}{numeric: measurement times for the repeated blood pressure measurements.}
#'   \item{y}{numeric: longitudinal measurements.}
#'   \item{time}{numeric: time to event (or censoring).}
#'   \item{event}{integer: event indicator. Coded as 0 = right-censoring, and 1 = event.}
#'   \item{X1}{integer: time-independent binary explanatory variable.}
#'   \item{X2}{numeric: time-independent continuous explanatory variable.}
#' }
#'
#'
#' @examples
#' data(dataLong)
NULL
