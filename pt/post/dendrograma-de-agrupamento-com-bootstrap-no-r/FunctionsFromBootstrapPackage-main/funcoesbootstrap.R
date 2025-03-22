#Bootstrap functions adapted from 'booststrap' package (Gibb S. 2013) - See http://www.github.com/sgibb/bootstrap/
#--------------------------------------------------------------------------------------------------------------------------#
bootstrap <- function(x, fun, n=1000L, mc.cores=getOption("mc.cores", 1L)) {
  fun <- match.fun(fun)

  origin <- .clust(x, fun=fun)

  v <- mclapply(seq_len(n), function(y, origin, size, fun, nco) {
    current <- .clust(.resample(x, size=size), fun=fun)
    return(.calculateMatches(origin, current, nco))
  }, mc.cores=mc.cores, origin=origin, size=ncol(x), fun=fun, nco=ncol(origin))

  return(colSums(do.call(rbind, v))/n)
}
#--------------------------------------------------------------------------------------------------------------------------#
as.binary.matrix.hclust <- function(x) {
  nr <- as.integer(nrow(x$merge))

  m <- matrix(0L, nrow=nr, ncol=nr+1L)

  for (i in seq_len(nr)) {
    left <- x$merge[i, 1L]

    if (left < 0L) {
      ## negative values correspond to observations
      m[i, -left] <- 1L
    } else {
      ## positive values correspond to childcluster
      m[i, ] <- m[left, ]
    }

    right <- x$merge[i, 2L]

    if (right < 0L) {
      ## negative values correspond to observations
      m[i, -right] <- 1L
    } else {
      ## positive values correspond to childcluster
      m[i, ] <- m[i,] | m[right, ]
    }
  }

  return(m)
}
#--------------------------------------------------------------------------------------------------------------------------#
.text.coord.hclust <- function(x) {
  nr <- as.integer(nrow(x$merge))

  p <- matrix(c(rep(0L, nr), x$height), nrow=nr, ncol=2, byrow=FALSE,
              dimnames=list(c(), c("x", "y")))
  o <- order(x$order)
  tmp <- double(2)

  for (i in seq_len(nr)) {
    left <- x$merge[i, 1L]

    if (left < 0L) {
      ## negative values correspond to observations
      tmp[1L] <- o[-left]
    } else {
      ## positive values correspond to childcluster
      tmp[1L] <- p[left, 1L]
    }

    right <- x$merge[i, 2L]

    if (right < 0L) {
      ## negative values correspond to observations
      tmp[2L] <- o[-right]
    } else {
      ## positive values correspond to childcluster
      tmp[2L] <- p[right, 1L]
    }

    p[i, 1L] <- mean(tmp)
  }

  return(p)
}
#--------------------------------------------------------------------------------------------------------------------------#
bootlabels.hclust <- function(x, bootstrapValues, horiz=FALSE, ...) {
  p <- .text.coord.hclust(x)
  if (horiz) {
    p[, c(2,1)] <- p
  }
  labels <- sprintf("%.2f", bootstrapValues)
  text(p, labels=labels, ...)
  invisible(NULL)
}
#--------------------------------------------------------------------------------------------------------------------------#
.clust <- function(x, fun) {
  hc <- fun(x)
  return(as.binary.matrix.hclust(hc))
}
#--------------------------------------------------------------------------------------------------------------------------#
.calculateMatches <- function(origin, current, nc=ncol(origin)) {
  ## both 1
  one <- tcrossprod(origin, current)
  ## both 0
  zero <- tcrossprod(1-origin, 1-current)

  ## calc matches
  return(rowSums(one + zero == nc))
}
#--------------------------------------------------------------------------------------------------------------------------#
.resample <- function(x, size=ncol(x)) {
  sel <- sample.int(ncol(x), size=size, replace=TRUE)
  return(x[, sel])
}
#--------------------------------------------------------------------------------------------------------------------------#