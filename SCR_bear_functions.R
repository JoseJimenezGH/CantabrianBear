# Create a new mask object including a specified subset of good-habitat pixels

addCore <- function(JAGSmask, type = c("perimeter", "traps", "polygon"),
    buffer = 0, poly=NULL, cell.overlap = c("centre","any","all"), plot=TRUE)  {

  if(!inherits(JAGSmask, "JAGSmask"))
    stop("JAGSmask is not a valid 'JAGSmask' object.")

  type <- match.arg(type)
  olap <- match.arg(cell.overlap)

  pixCoord <- cbind(x=as.vector(row(JAGSmask$habMat))+0.5,
                    y=as.vector(col(JAGSmask$habMat))+0.5)
  pixOriginal <- convertOutput(pixCoord, JAGSmask)
  if(olap == "any" || olap == "all") {
    corns <- cbind(c(-0.5,-0.5,+0.5,+0.5), c(-0.5,+0.5,-0.5,+0.5))
    corners <- array(0, c(nrow(pixCoord), 4, 2))
    for(i in 1:4)
      corners[,i,] <- sweep(pixCoord, 1:2, corns[i, ], "+")
    cornOriginal <- convertOutput(corners, JAGSmask)
  }

  if(type == "traps") {
    if(buffer == 0)
      warning("type='traps' with buffer = 0 will select zero pixels!")
    bufferPix <- buffer / pixelWidth(JAGSmask)
    if(olap == "centre") {
      isCore <- secr::distancetotrap(pixCoord, JAGSmask$trapMat) <= bufferPix
    } else {
      isCornerCore <- matrix(FALSE, nrow(pixCoord), 4)
      for(i in 1:4) {
        isCornerCore[,i] <- secr::distancetotrap(corners[,i,], JAGSmask$trapMat) <= bufferPix
      }
      sumcorn <- rowSums(isCornerCore)
      if(olap=='any') {
        isCore <- sumcorn > 0
      } else {
        isCore <- sumcorn == 4
      }
    }
  }

  if(type == "perimeter") {
    tmp <- JAGSmask$trapMat[chull(JAGSmask$trapMat), ]
    poly <- convertOutput(tmp, JAGSmask)
    type <- "polygon"
  }

  if(type == "polygon") {
    if(is.null(poly))
      stop("'poly' must be supplied for type='polygon'")
    if(!inherits(poly, "SpatialPolygons"))
      poly <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(poly)), 1)))
    if(buffer > 0)
      poly <- rgeos::gBuffer(poly, width = buffer)
    if(olap == "centre") {
      isCore <- secr::pointsInPolygon(pixOriginal, poly)
    } else {
      isCornerCore <- matrix(FALSE, nrow(pixCoord), 4)
      for(i in 1:4) {
        isCornerCore[,i] <- secr::pointsInPolygon(cornOriginal[,i,], poly)
      }
      sumcorn <- rowSums(isCornerCore)
      if(olap=='any') {
        isCore <- sumcorn > 0
      } else {
        isCore <- sumcorn == 4
      }
    }
  }

  JAGSmask$coreMat <- JAGSmask$habMat * isCore
  if(plot)
    plot(JAGSmask)
  attr(JAGSmask, "poly") <- poly
  return(JAGSmask)
}

# Convert secr mask and traps objects for use with JAGS
# Output: list with habitat matrix, new traps coords as a matrix,
#   extent, pixel width, and area in original units;
#   attributes: original bounding box as a matrix, false origin, and pixel width.

convertMask <- function(secrmask, secrtraps, plot=TRUE) {

  # Sanity checks
  if(!inherits(secrmask, "mask"))
    stop("'", deparse(substitute(secrmask)), "' is not a valid 'mask' object.")
  if(!inherits(secrtraps, "traps"))
    stop("'", deparse(substitute(secrtraps)), "' is not a valid 'traps' object.")
  # Get point spacing, which will be our pixelWidth
  pixWidth <- attr(secrmask, "spacing")
  # Do we need to deal with masks without spacing?
  # min(abs(diff(secrmask$x))) # should be same

  bbox <- as.matrix(attr(secrmask, "boundingbox"))
  # Create 'false origin' so that SW corner of matrix is at [1, 1]
  origin <- bbox[1, ] - pixWidth

  # Get dimensions of the matrix
  nrows <- round((bbox[2, 1] - bbox[1, 1]) / pixWidth)
  ncols <- round((bbox[3, 2] - bbox[2, 2]) / pixWidth)
  habMat <- matrix(0L, nrow=nrows, ncol=ncols)
  # Convert mask x and y to col/row numbers
  dex <- as.matrix(floor(sweep(secrmask, 2, origin) / pixWidth))
  for(i in 1:nrow(dex))
    habMat[dex[i,1], dex[i,2]] <- 1L

  # Convert trap coordinates to the new units:
  newtraps <- sweep(secrtraps, 2, origin) / pixWidth

  out <- list(habMat = habMat, 
              trapMat = as.matrix(newtraps), 
              upperLimit = c(x=nrows+1, y=ncols+1),
              pixelWidth = pixWidth,
              area = sum(habMat) * pixWidth^2)
  attr(out, "boundingbox") <- bbox
  attr(out, "origin") <- origin
  attr(out, "pixelWidth") <- pixWidth
  class(out) <- "JAGSmask"
  
  if(plot)
    plot.JAGSmask(out)
  
  return(out)
}

# Convert JAGS AC location output back to the original
#   coordinate reference system
# Output...

# Helper function for x/y names
pmatch_xy <- function(names) {
  xy <- c(pmatch("x", tolower(names)),
          pmatch("y", tolower(names)))
  if(is.na(diff(xy)) || diff(xy) == 0)
    return(1:2)
  return(xy)
}

convertOutput <- function(ACs, JAGSmask) {
  if(!inherits(JAGSmask, "JAGSmask"))
    stop("'", deparse(substitute(JAGSmask)), "' is not a valid 'JAGSmask' object.")
  classOut <- class(ACs)[1]
  if(is.list(ACs) && length(ACs) == 2) {
    xy <- pmatch_xy(names(ACs))
    x <- ACs[[xy[1]]]
    y <- ACs[[xy[2]]]
  } else if(is.matrix(ACs) && ncol(ACs) == 2) {
    xy <- pmatch_xy(colnames(ACs))
    x <- ACs[, xy[1]]
    y <- ACs[, xy[2]]
  } else if(is.array(ACs) && length(dim(ACs)) == 3 && dim(ACs)[3] == 2) {
    x <- ACs[, , 1]
    y <- ACs[, , 2]
  } else {
    stop("invalid input")
  }
  # Get pixel width and original false origin
  pixWidth <- pixelWidth(JAGSmask)
  origin <- attr(JAGSmask, "origin")

  x1 <- x * pixWidth + origin[1]
  y1 <- y * pixWidth + origin[2]

  out <- switch(classOut,
    matrix = cbind(x = x1, y = y1),
    array = abind(x = x1, y = y1, along=3),
    data.frame = data.frame(x = x1, y = y1),
    list(x = x1, y = y1)
  )
  return(out)
}
  
# Convert a habitat raster and trap coordinates for use with JAGS
# Output: list with habitat matrix, new traps coords as a matrix,
#   extent, pixel width, and area in original units;
#   attributes: original bounding box as a matrix, false origin, and pixel width.

convertRaster <- function(raster, traps, plot=TRUE) {

  # Sanity checks
  if(!inherits(raster, "Raster"))
    stop("'", deparse(substitute(raster)), "' is not a valid 'raster' object.")
  # raster <- raster::trim(raster)
  if(!inherits(traps, "data.frame"))
    stop("'", deparse(substitute(traps)), "' is not a valid 'data.frame' object.")
  # Get point spacing, which will be our pixelWidth
  pixWidth <- xres(raster)
  stopifnot(all.equal(pixWidth, yres(raster)))

  bbox <- matrix(extent(raster)[c(1,2,2,1,3,3,4,4)], 4, 2)
  origin <- bbox[1, ] - pixWidth

  # Generate habMat
  r <- raster[[1]]
  values(r) <- !is.na(values(r))
  tmp <- raster::as.matrix(r) * 1
  habMat <- t(tmp[nrow(tmp):1, ])

  # Convert raster or stack to a list of matrices # new 2019-07-2
  nl <- nlayers(raster)
  covs <- vector('list', nl)
  if(nl == 1) {
    names(covs) <- "covMat"
  } else {
    names(covs) <- names(raster)
  }
  for(i in 1:nl) {
    tmp <- raster::as.matrix(raster[[i]])
    tmp[is.na(tmp)] <- 0
    covs[[i]] <- t(tmp[nrow(tmp):1, ])
  }
  # str(covs)
  
  # Convert trap coordinates to the new units:
  newtraps <- sweep(traps, 2, origin) / pixWidth

  out <- c(covs, list(habMat = habMat,
              trapMat = as.matrix(newtraps),
              upperLimit = c(x=nrow(habMat)+1, y=ncol(habMat)+1),
              pixelWidth = pixWidth,
              area = sum(habMat) * pixWidth^2))
  attr(out, "boundingbox") <- bbox
  attr(out, "origin") <- origin
  attr(out, "pixelWidth") <- pixWidth
  class(out) <- "JAGSmask"

  if(plot)
    plot.JAGSmask(out)

  return(out)
}


if(FALSE) {

str(tmp <- convertRaster(stack, traps, plot=TRUE))
str(tmp <- convertRaster(patchR, traps, plot=TRUE))



}


# Get the number of activity centres in the core

# Returns an MCMC chain with posterior probability of number of ACs in the core.

getACinCore <- function(S,          # iters x animals x 2 array with AC locations on the pixel scale
                        w,          # iters x animals 1/0 matrix, 1 if the animal is real
                        JAGSmask){  # the JAGSmask used for the analysis with a component coreMat
  if(!inherits(JAGSmask, "JAGSmask"))
    stop(deparse(substitute(JAGSmask)), " is not a valid JAGSmask object.", call.=FALSE)
  if(is.null(JAGSmask$coreMat))
    stop(deparse(substitute(JAGSmask)), " does not contain core information.", call.=FALSE)
  if(length(dim(S)) != 3 || dim(S)[3] != 2)
    stop(deparse(substitute(JAGSmask)), " is not a proper activity centre array.", call.=FALSE)
  dimMat <- dim(JAGSmask$coreMat)
  if(max(S[,,1], na.rm=TRUE) > dimMat[1] + 1 ||
      min(S[,,1], na.rm=TRUE) < 1 ||
      max(S[,,2], na.rm=TRUE) > dimMat[2] + 1 ||
      min(S[,,2], na.rm=TRUE) < 1)
  stop(deparse(substitute(S)), " has values outside the range of the matrix. Are these on the pixel scale?", call.=FALSE)
      
  if(!missing(w))
    S[w==0] <- NA
  inCore <- array(NA, dim(S)[1:2])
  for(j in 1:ncol(inCore))
    inCore[, j] <- JAGSmask$coreMat[S[, j, ]]
  return(rowSums(inCore, na.rm=TRUE))
}


.onAttach <- function(libname, pkgname) {
  version <- try(utils::packageVersion('makeJAGSmask'), silent=TRUE)
  if(!inherits(version, "try-error"))
    packageStartupMessage("This is makeJAGSmask ", version,
      ". For overview type ?makeJAGSmask; for changes do news(p='makeJAGSmask').")
}

# Extract pixel width from a JAGSmask object.

pixelWidth <- function(JAGSmask) {
  pixWidth <- JAGSmask$pixelWidth
  if(is.null(pixWidth))
    pixWidth <- attr(JAGSmask, "pixelWidth")
  return(pixWidth)
}


# Function to plot a JAGSmask object

plot.JAGSmask <- function(x, col, verify=TRUE, ...) {

  if(names(x)[1] == "habMat") {  # no covariate matrices
    if(missing(col))
      col <- c("grey", "white", "yellow")
    if(is.null(x$coreMat)) {
      toPlot <- x$habMat
    } else {
      toPlot <- x$habMat + x$coreMat
    }
    image(x=1:x$upperLimit[1],
          y=1:x$upperLimit[2],
          z=toPlot,
          ann=FALSE, axes=FALSE, col=col, asp=1, ...)
  } else {
    if(missing(col))
      col <- terrain.colors(100)
    image(x=1:x$upperLimit[1],
          y=1:x$upperLimit[2],
          z=x[[1]],
          ann=FALSE, axes=FALSE, col='black', asp=1)#, ...)
    tmp1 <- x[[1]]
    tmp1[x$habMat==0] <- NA
    image(x=1:x$upperLimit[1],
          y=1:x$upperLimit[2],
          z=tmp1, col=col, add=TRUE)#, ...)
    if(!is.null(x$coreMat)) {
      tmp2 <- x$coreMat
      tmp2 <- 1-tmp2
      tmp2[tmp2==0] <- NA
      tmp2[x$habMat==0] <- NA
      image(x=1:x$upperLimit[1],
            y=1:x$upperLimit[2],
            z=tmp2, col=adjustcolor('grey', 0.7), add=TRUE)
    }
  }

  bbox <- attr(x, "boundingbox")
  xlabels <- pretty(bbox[1:2, 1], n=5)
  xpos <- (xlabels - bbox[1, 1]) / diff(bbox[1:2, 1]) * nrow(x$habMat) + 1
  ylabels <- pretty(bbox[2:3, 2], n=5)
  ypos <- (ylabels - bbox[2, 2]) / diff(bbox[2:3, 2]) * ncol(x$habMat) + 1

  title(xlab="Easting", ylab="Northing")
  axis(1, at=xpos, labels = xlabels)
  axis(2, at=ypos, labels = ylabels)
  box()
  points(x$trapMat, pch=3, col='red', xpd=TRUE)

  if(verify) {
    # Check locations of traps
    trapcells <- floor(x$trapMat)
    ok <- numeric(nrow(trapcells))
    for(i in 1:nrow(trapcells))
      ok[i] <- x$habMat[trapcells[i,1], trapcells[i,2]]
    if(!all(ok == 1)) {
      cat("The following traps appear to be in bad habitat:\n", which(!ok), "\n")
      cat("They are circled in the plot.\n")
      cat("If on the edge, this is probably due to rasterization of the habitat polygon.\n")
    }
    if(!all(ok == 1))
      points(x=x$trapMat[!ok, 1], y=x$trapMat[!ok, 2], col='red', cex=2, xpd=TRUE)
  }
}


# Function to plot activity centres

plotACs <- function(
    which=NA,     # which ACs to plot (don't usually want to do all in one plot)
    ACs,          # iters x animals x 2 array, as produced by convertOutput; NA out phantoms
    traps,        # 2-column matrix or data frame with trap coordinates
    Y,            # animals x traps matrix with capture histories; rownames assumed to be animal IDs
    hab,          # spatialPolygons object with the extent of the habitat
    howMany=3000, # number of points to plot for each animal
    show.labels=TRUE, # whether to label plot with animal IDs
    rad=50,       # amount of jitter to add to capture locations
    link=TRUE,    # if TRUE, link capture locations with dotted line
    colors        # vector of colors to use for plotting
  )  {

  # Reduce number of iterations
  if(dim(ACs)[1] > howMany) {
    keep <- seq(1, dim(ACs)[1], length = howMany)
    ACs <- ACs[keep,,]
  }
  # Get posterior means for locations
  x <- colMeans(ACs[, ,1], na.rm=TRUE)
  y <- colMeans(ACs[, ,2], na.rm=TRUE)

  # Recover animal IDs if possible
  M <- dim(ACs)[2] # total, incl. uncaptures
  animalIDs <- sprintf("id%03d", 1:M)
  captLocList <- NULL
  if(!missing(Y)) {
    ncap <- nrow(Y) # number of animals captured
    if(!is.null(rownames(Y)))
      animalIDs[1:nrow(Y)] <- rownames(Y)
    # which to plot
    if(any(is.na(which)))
      which <- 1:ncap

    # Get capture locations
    if(!missing(traps)) {
      captLocList <- vector('list', ncap)
      for(i in 1:ncap) {
        captTraps <- which(Y[i, ] > 0) # Which traps caught the animal
        tmp <- traps[captTraps, , drop=FALSE] # Locations of the traps
        jitangle <- runif(nrow(tmp), -pi, pi)
        jit <- cbind(rad*sin(jitangle), rad*cos(jitangle))
        captLocList[[i]] <- tmp + jit
      }
    }
  }
  # if(any(is.na(which)))
    # which <- 1:M

  # do the plot
  # -----------
  MASS::eqscplot(x, y, type='n', xlim=range(ACs[, , 1], na.rm=TRUE),
    ylim=range(ACs[, , 2], na.rm=TRUE),
    ann=FALSE, axes=FALSE)
  if(!missing(hab)) {
    if(inherits(hab, "SpatialPolygons")) {
    sp::plot(hab, add=TRUE)
    } else if((is.matrix(hab) || is.data.frame(hab)) && ncol(hab) == 2) {
      lines(hab)
    }
  }
  if(!missing(traps))
    points(traps, pch=3, col='red')
  if(missing(colors))
    colors <- palette()[-1]
  colno <- 1
  for(i in which) {
    col <- colors[colno]
    points(ACs[, i, ], cex=0.1, col=adjustcolor(col, 0.3))
    if(!is.null(captLocList) && i <= ncap){
      points(captLocList[[i]], pch=21, col='black', bg=col, cex=1.2)
      if(link && nrow(captLocList[[i]]) > 1)
        lines(captLocList[[i]], col=col)
    }
    colno <- colno+1
    if(colno > length(colors))
    colno <- 1
  }
  if(show.labels)
    plotrix::boxed.labels(x[which], y[which], labels=animalIDs[which])

  # data frame to return
  # --------------------
  niter <- dim(ACs)[1]
  tmp <- matrix(ACs[, which, ], ncol=2)
  IDs <- make.unique(animalIDs[which])
  out <- data.frame(ID=factor(rep(IDs, each=niter)), x=tmp[,1], y=tmp[,2])
  out <- out[complete.cases(out), ]
  return(invisible(out))
}

# Function to print a JAGSmask object

print.JAGSmask <- function(x, ...) {
  cat("An object of class 'JAGSmask'\n")
  size <- dim(x$habMat)
  ntraps <- nrow(x$trapMat)
  cat("The habitat mask has", size[1], "rows and", size[2], "columns.\n")
  cat("The trap matrix has coordinates for", ntraps, "traps.\n")
  cat("Original bounding box is:\n")
  print(attr(x, "boundingbox"))
}


# Generate the coordinates of randomly selected points in good habitat.

# This is intended to provide starting values. The starting locations for animals
#   caught can be fixed.

randomPoints <- function(n, JAGSmask, fixed) {

  # generate 2 x n random locations:
  XY <- cbind(x=runif(2*n, 1, JAGSmask$upperLimit[1]),
              y=runif(2*n, 1, JAGSmask$upperLimit[2]))
  # check habitat
  ok <- JAGSmask$habMat[XY] == 1
  # reject locations in bad habitat
  XY <- XY[ok, , drop=FALSE]
  # if not enough values, do more
  while(nrow(XY) < n) {
    more <- cbind(runif(n, 1, JAGSmask$upperLimit[1]),
                  runif(n, 1, JAGSmask$upperLimit[2]))
    ok <- JAGSmask$habMat[more] == 1
    XY <- rbind(XY, more[ok, , drop=FALSE])
  }
  # Truncate to length n
  XY <- XY[1:n, ]
  # Insert fixed values, ignoring NAs
  if(!missing(fixed)) {
    fix <- which(!is.na(rowSums(fixed)))
    fixed0 <- fixed[fix, , drop=FALSE]
    if(any(fixed0 < 1) ||
        any(fixed0[, 1] > JAGSmask$upperLimit[1]) ||
        any(fixed0[, 2] > JAGSmask$upperLimit[2]))
      stop("At least 1 fixed location is outside the habitat mask.")
    ok <- JAGSmask$habMat[fixed0] == 1
    if(any(!ok))
      warning("Fixed locations ", which(!ok), " are in bad habitat.")
    XY[fix, ] <- fixed[fix, ]
  }
  return(XY)
}


GetDetectionRate <- nimbleFunction(
  run = function(s = double(1), p0=double(1), alpha1=double(0), 
                 X=double(2), J=double(0), z=double(0)){ 
    returnType(double(1))
    if(z==0) return(rep(0,J))
    if(z==1){
     d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
     ans <- p0*exp(-alpha1*d2)
     return(ans)
    }
  }
)


# Vectorized Poisson
dPoissonVector <- nimbleFunction(
  run = function(x = double(1), lambda = double(1),
  log = integer(0, default = 0)) {
    J <- length(x)
    ans <- 0.0
    for(j in 1:J)
      ans <- ans + dpois(x[j], lambda[j], 1)
    returnType(double())
    if(log) return(ans)
    else return(exp(ans))
  })

rPoissonVector  <- nimbleFunction(
  run = function(n = integer(), lambda = double(1)) {
    J <- length(lambda)
    ans<- numeric(J)
    for(j in 1:J)
      ans[j] <- rpois(1, lambda[j])
    returnType(double(1))
    return(ans)
  })

registerDistributions(list(
  dPoissonVector = list(
    BUGSdist = "dPoissonVector(lambda)",
    Rdist = "dPoissonVector(lambda)",
    discrete = TRUE,
    range = c(0, Inf),
    types = c('value = double(1)', 'lambda = double(1)'))
))


# Vectorized binomial
dBinomVector <- nimbleFunction(
  run = function(x = double(1), P = double(1), K = double(0), log = integer(0, default = 0)) {
    J <- length(x)
    ans <- 0.0
    for(j in 1:J)
      ans <- ans + dbinom(x[j], K, P[j], 1)
    returnType(double())
    if(log) return(ans)
    else return(exp(ans))
  })

rBinomVector  <- nimbleFunction(
  run = function(n = integer(), P = double(1), K = double(0)) {
    J <- length(P)
    ans<- numeric(J)
    for(j in 1:J)
      ans[j] <- rbinom(1, K, P[j])
    returnType(double(1))
    return(ans)
  })

registerDistributions(list(
  dBinomVector = list(
    BUGSdist = "dBinomVector(P, K)",
    Rdist = "dBinomVector(P, K)",
    discrete = TRUE,
    range = c(0, Inf),
    types = c('value = double(1)', 'P = double(1)', 'K = double(0)'))
))


# Function to initialize complex nodes
InitNod<-function(simNodes){
  simNodeScalar <- Rmodel$expandNodeNames(simNodes)
  allNodes <- Rmodel$getNodeNames()
  nodesSorted <- allNodes[allNodes %in% simNodeScalar]
  set.seed(1) # to fix simulations
  for(n in nodesSorted) {
    Rmodel$simulate(n)
    depNodes <- Rmodel$getDependencies(n)
    Rmodel$calculate(depNodes)
  }
}