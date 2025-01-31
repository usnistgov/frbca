#' Input parameters.
#'
#' Nested list of analysis parameters used for FR-BCA.
#'
#' @format A list with with 2 variables:
#' \describe{
#' \itemize{
#'   \item{model}{model name (NOT USED)}
#'   \item{parameters}{list of analysis parameters}
#'   \itemize{
#'     \item{base}{baseline analysis parameters}
#'     \itemize{
#'       \item{floor_area}{total building floor area}
#'       \item{delta}{discount rate}
#'       \item{T}{time horizon}
#'       \item{loss}{list of economic losses}
#'       \item{tenant}{tenants per square foot}
#'       \item{recapture}{loss recapture rate}
#' }
#'     \item{sensitivity}{sensitivity analysis parameters}
#'     \itemize{
#'       \item{loss}{list of economic losses}
#'       \item{delta}{list of discount rates for sensitivity analysis}
#'       \item{T}{list of time horizons for sensitivity analysis}
#' }
#' }
#' }
#' }
#' @source Fung et al. (2025)
"input_param"
