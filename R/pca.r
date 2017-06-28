#' Biplot of a prcomp object
#'
#' @param r prcomp object
#' @param d1,d2 An integer for dimension
#' @param show.var.name Show variable names. Defaulr is \code{TRUE}.
#'
#' @export
biplot <- function(r, d1 = 1, d2 = 2, show.var.name = TRUE) {
  r <- r$rotation
  r <- r[, c(d1, d2)]
  mx <- max(abs(range(r[,1])))
  plot(r, type = "n", ylim = c(-1, 1), xlim = c(-mx*1.2, mx*1.2),
       xlab = paste0("PC", d1), ylab = paste0("PC", d2))
  abline(v = 0, h = 0, lty = 3)
  arrows(0, 0, r[,1], r[,2], len = 0.1, col = "darkred")
  if (show.var.name) text(1.1 * r, rownames(r), col = "darkred",
                          xpd = TRUE, cex = 1)
}

#' Variable Coordinates for prcomp object
#'
#' @param pca prcomp object
#'
#' @export
var.coord <- function(pca) {
  var_cor_func <- function(var.loadings, comp.sdev) var.loadings * comp.sdev[1:length(var.loadings)]
  # Variable correlation/coordinates
  t(apply(pca$rotation, 1, var_cor_func, pca$sdev))
}

#' Contribution to PC
#'
#' @param pca prcomp object
#'
#' @export
contrib <- function(pca) {
  var.cos2 <- var.coord(pca)^2
  comp.cos2 <- apply(var.cos2, 2, sum)
  cntrb <- function(var.cos2, comp.cos2){ var.cos2 * 100/comp.cos2 }
  t(apply(var.cos2,1, cntrb, comp.cos2))
}
