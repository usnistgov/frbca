#' Construction costs for FR-BCA
#'
#' Structural, nonstructural, total construction, and total project costs.
#'
#' @format A data frame with 48 rows and 10 variables:
#' \describe{
#'   \item{system}{Structural system (Lateral Force Resisting System)}
#'   \item{num_stories}{Number of stories}
#'   \item{design_s}{Structural design intervention (baseline/RC IV/backup-frame)}
#'   \item{design_ns}{Nonstructural design intervention (baseline/nsfr)}
#'   \item{model}{Model name (system-num_stories-design_s-design_ns)}
#'   \item{intervention}{Recovery-based design intervention (0/1)}
#'   \item{c_s}{Structural construction cost}
#'   \item{c_ns}{Nonstructural construction cost}
#'   \item{construction}{Total construction cost}
#'   \item{project}{Total project cost}
#' }
#' @source Fung et al. (2025)
"input_cost"
