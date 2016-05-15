#' Bicubic interpolation on an array
#'
#' First two dimensions of the array is interpolated by bicubic grid.
#' Other dimensions are repeated. If \code{xlim} and \code{ylim} are
#' not set, range of first two dimensions are used. If \code{dx} and
#' \code{dy} are not set, half of the grid spacing is used. If \code{x}
#' matrix/array has \code{\link{dimnames}}, they are used as coordinates
#' of the rectangular grid.
#'
#' @param x A matrix or higher dimensonal array  containing the
#'          z[i,j] data values for the grid points (x[i],y[j]).
#' @param xlim Lower and upper limits (x1, x2) for x-coordinate used for output.
#' @param ylim Lower and upper limits (y1, y2) for y-coordinate used for output.
#' @param dx Output grid spacing in x direction.
#' @param dy Output grid spacing in y direction.
#' @param par Use \link{parallel} package to interpolate.
#' @return Interpolated array
#'
#' @examples
#' x <- array(rnorm(1000), dim = c(10, 10, 10))
#' xi <- ibicubic(x)
#' str(xi) # should have dim = c(19, 19, 10)
#'
#' x <- array(rnorm(1000), dim = c(10, 10, 10),
#'            dimnames = list(seq(2, 20, 2), seq(2, 20, 2), NULL))
#' xi <- ibicubic(x)
#' str(xi) # note dimnames
#'
#' xi <- ibicubic(x, xlim = c(4, 6), ylim = c(6, 8))
#' str(xi) # note dimnames
#'
#' @export
ibicubic <- function(x, xlim = NULL, ylim = xlim, dx = NULL, dy = dx, par = T) {
  stopifnot( length(dim(x)) > 1 )

  Rev <- function(x, margin, ...) {
    if (!is.array(x)) stop("'x' is not an array")
    newdim <- rep("", length(dim(x)))
    newdim[margin] <- paste(dim(x), ":1", sep = "")[margin]
    z <- eval(parse(text = gettextf("x[%s,drop = F]", paste(newdim, sep = "",
                                                            collapse = ","))))
    class(z) <- oldClass(x)
    return(z)
  }

  set_dim_monoton_inc <- function(x) {
    dims <- vector()
    for (i in 1:2) { # check only first two dims
      n <- as.numeric(dimnames(x)[[i]])
      if (length(n)) {
        if (!all(n == cummax(n))) {
          x <- Rev(x, i)
          dims <- c(dims, i)
        }
      } else {
        dimnames(x)[i] <- list(1:dim(x)[i])
      }
    }
    return(list(x = x, dims = dims))
  }

  new_xy <- function(x, xlim, dx) {
    rng <- range(x)
    new_r <- seq(rng[1], rng[2], dx)
    return(new_r[new_r >= xlim[1] & new_r <= xlim[2]])
  }

  x <- set_dim_monoton_inc(x)
  dims <- x$dims; x <- x$x

  r <- as.numeric(dimnames(x)[[1]])
  c <- as.numeric(dimnames(x)[[2]])
  if (is.null(dx)) dx <- abs((r[2] - r[1]) / 2)
  if (is.null(dy)) dy <- abs((c[2] - c[1]) / 2)
  if (is.null(xlim)) xlim <- range(r)
  if (is.null(ylim)) ylim <- range(c)
  new_r <- new_xy(r, xlim, dx)
  new_c <- new_xy(c, ylim, dy)

  l <- length(dim(x))
  if (require(akima)) {
    bicg <- function(z) akima::bicubic.grid(r, c, z, xlim, ylim, dx, dy)$z
    if (l > 2) {
      if (require(parallel) && par) {
        cl <- parallel::makeCluster(getOption("cl.cores", detectCores()/2))
        parallel::clusterExport(cl, c("r", "c", "xlim", "ylim", "dx", "dy"),
                                envir = environment())
        parallel::clusterEvalQ(cl, library(akima))
        xi <- parApply(cl, x, (1:l)[-c(1, 2)], bicg)
        parallel::stopCluster(cl)
      } else {
        xi <- apply(x, (1:l)[-c(1,2)], bicg)
      }
    } else {
      xi <- bicg(x)
    }
  } else {
    stop("Install akima package")
  }

  xi <- array(xi, dim = c(length(new_r), length(new_c), dim(x)[-c(1,2)]))
  dimnames(xi)[[1]] <- new_r
  dimnames(xi)[[2]] <- new_c
  dimnames(xi)[-c(1,2)] <- dimnames(x)[-c(1,2)]
  names(dimnames(xi)) <- names(dimnames(x))
  if (length(dims)) xi <- Rev(xi, dims)

  # copy attributes
  attrs <- names(attributes(x))
  attrs <- attrs[!(attrs %in% c("dim", "dimnames"))]
  for (a in attrs) attr(xi, a) <- attr(x, a)
  return(xi)
}

