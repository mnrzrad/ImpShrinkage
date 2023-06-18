#' Hald's Cement Data
#'
#' Heat evolved in setting of cement, as a function of its chemical composition.
#'
#' \code{cement} is taken from an example in Little & Rubin's book on
#' Statistical Analysis with Missing Data (2002), pp.~154,
#' for demonstrating estimation of multivariate means and variances
#' when the missing data pattern is monotone.
#'
#' @name cement
#' @format
#' A data.frame with 13 observations on the following 5 variables.
#' \describe{
#'   \item{x1}{percentage weight in clinkers of 3CaO.Al2O3}
#'   \item{x2}{percentage weight in clinkers of 3CaO.SiO2}
#'   \item{x3}{percentage weight in clinkers of 4CaO.Al2O3.Fe2O3}
#'   \item{x4}{percentage weight in clinkers of 2CaO.SiO2}
#'   \item{y}{heat evolved (calories/gram)}
#' }
#'
#' @source Woods, H., Steinour, H. H. and Starke, H. R. (1932) Effect of composition of Portland cement on heat evolved during hardening. Industrial Engineering and Chemistry, 24, 1207–1214.
#'
#' @references
#' Davison, A. C. (2003) Statistical Models. Cambridge University Press. Page 355.
#'
#' Draper, N.R. and Smith, H. (1998) Applied Regression Analysis. Wiley. Page 630.
#'
#' Roderick J.A. Little and Donald B. Rubin (2002). Statistical Analysis with Missing Data, Second Edition. Wilely. Page 154.
#' @examples
#' data("cement")
NULL