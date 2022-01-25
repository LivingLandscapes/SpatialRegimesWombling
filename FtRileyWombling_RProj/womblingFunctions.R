
###### read_project()
# Function to read in rasters from a folder and reproject them
read_project <- function(folder,
                         all.files = FALSE,
                         specifyFiles = c(1991:2015, 2017),
                         proj_str = "+proj=utm +zone=14 +datum=WGS84 +units=m",
                         normalize = TRUE,
                         mean.center = FALSE) {
  
  if (all.files) {
    # select all files
    selected_files <- list.files(folder)
  } else {
    # Select only files specified
    selected_files <- list.files(folder)[unlist(sapply(specifyFiles,
                                                       function(x) which(str_detect(list.files(folder),
                                                                                    as.character(x)))))]
  }
  
  # Paste together folder and filenames
  Files <- sapply(selected_files,
                  function(x) paste0(folder, "/", x))
  
  # Read in raster for each year
  rasters <- lapply(Files, raster)
  
  #
  if (normalize) {
    # re-project rasters
    proj_rasts <- lapply(rasters,
                         function(x) scale(projectRaster(x,
                                                         crs = proj_str),
                                           center = mean.center))
  } else {
    # re-project rasters
    proj_rasts <- lapply(rasters,
                         function(x) projectRaster(x,
                                                   crs = proj_str))
  }
  
  
  
  return(proj_rasts)
}
# debug(read_project)
# undebug(read_project)

#=============================================================================
## Wombling functions
#=============================================================================

GeoGWR <-function(z,
                  coords,
                  nsim,
                  minalpha,
                  maxalpha,
                  confInt,
                  PLOT_alpha = FALSE){
  
  x <-coords[,1] 
  y <-coords[,2]
  stGWR <-matrix(0,nsim,3)
  
  
  for(k in 1:nsim){
    
    alpha <-runif(1,minalpha,maxalpha)
    D <-as.matrix(dist(coords,up=T,diag=T))
    D <-D/max(D)
    W <-exp(-alpha*D)
    W <-as.matrix(1/(D^alpha))	
    diag(W)<-1
    
    ze <-numeric()
    for(i in 1:nrow(coords)){
      ze[i] <-lm(z~x+y,weight=W[,i])$fitted.values[i]
    }
    # ze <- sapply( ,
    #              function(X)  lm(z ~ x + y,
    #                              weight = W[ , X]))
    
    stGWR[k,1] <-alpha
    stGWR[k,2] <-cor(z,ze)^2
    stGWR[k,3] <-AIC(lm(ze~z))
    
  }
  
  if (PLOT_alpha){
    plot(stGWR[,1],stGWR[,2],xlab="Alpha",ylab="GWR R2",cex=1,pch=16)
    
  }
  
  best <-which(stGWR[,2]==max(stGWR[,2]))
  alpha_best <-stGWR[best,1]
  W <-exp(-alpha_best*D)
  W <-as.matrix(1/(D^alpha_best))	
  diag(W)<-1
  
  
  
  
  r2_coords <-numeric()
  bx <-numeric()
  by <-numeric()
  # pbx <-numeric()
  # pby <-numeric()
  ySE <- numeric()
  xSE <- numeric()
  ze <-numeric()
  
  for(i in 1:nrow(coords)){
    
    l_reg <-summary(lm(z~x+y,weight=W[,i]))
    
    if ( !(identical(as.numeric(dim(l_reg$coefficients)),
                     c(3, 4))) ) {
      next
    }
    
    r2_coords[i] <-l_reg$r.squared
    bx[i] <-l_reg$coefficients[2,1]
    by[i] <-l_reg$coefficients[3,1]
    # pbx[i] <-l_reg$coefficients[2,4]
    # pby[i] <-l_reg$coefficients[3,4]
    xSE[i] <-l_reg$coefficients[2,2]
    ySE[i] <-l_reg$coefficients[3,2]
    
    ze[i] <-lm(z~x+y,weight=W[,i])$fitted.values[i]
    
    # Print progress
    Sys.sleep(.1)
    cat(paste0(i, " of ", nrow(coords), "\r"))
    flush.console()
  }
  
  r2gwr <-cor(z,ze)^2
  AICgwr <-AIC(lm(ze~z),k=4)
  r2ols <-summary(lm(z~x+y))$r.squared
  AICols <-AIC(lm(z~x+y))
  
  statsGWR <-cbind(alpha_best,r2ols,r2gwr,AICols,AICgwr)
  
  output <-list("Stats"=statsGWR,
                "GWRr2"=r2_coords, 
                "GWRbx"=bx,
                "GWRby"=by,
                # "GWR_Px"=pbx,
                # "GWR_Py"=pby,
                "GWR_xSE" = xSE,
                "GWR_ySE" = ySE
                # ,
                # "GWR_yCI" = ySE * confInt,
                # "GWR_xCI" = xSE * confInt
  )
  
  return(output)
  
}

Run_GeoGWR <- function(DATA,
                       SpatialObject = TRUE,
                       YRcol = "Year",
                       confINT = 2.576,
                       spCols,
                       IDcols,
                       coordCols,
                       projFOURstr,
                       Nsim = 100,
                       minA = 5,
                       maxA = 15) {
  
  # reset to data.frame
  DATA <- as.data.frame(DATA)
  
  # Make sure Year is numeric
  DATA[[YRcol]] <- as.numeric(as.character(DATA[[YRcol]]))
  
  # Empty list  
  outGWR <- vector("list", length(unique(DATA[[YRcol]])))
  
  # Loop through years
  for (i in 1:length(unique(DATA[[YRcol]]))) {
    
    thisYear <- sort(unique(DATA[[YRcol]]))[i]
    
    DATA_sub <- DATA[which(DATA[[YRcol]] == thisYear), ]
    # pca1_sub <- DATA[which(DATA[[YRcol]] == thisYear), "PC1"]
    pca1_sub <- prcomp(decostand(DATA_sub[ , spCols],
                                 method = "hellinger"))$x[,1]
    
    # Test GWR function
    outGWR[[i]] <- GeoGWR(pca1_sub,
                          DATA_sub[ , coordCols],
                          Nsim,
                          minA,
                          maxA,
                          confInt = confINT)
    
    # # Clear the Console
    # cat("\014")
    
  }
  
  # Collapse list into a data.frame (while omitting any NAs dropped
  # from wombling loop)
  outGWR <- na.omit(do.call(rbind,
                    lapply(outGWR,
                           function(x) data.frame(GWRr2 = x$GWRr2,
                                                  GWRbx = x$GWRbx,
                                                  GWRby = x$GWRby,
                                                  GWR_xSE = x$GWR_xSE,
                                                  GWR_ySE = x$GWR_ySE))))
  
  
  # Create a data.frame with output 
  output <- data.frame(outGWR,
                          DATA[row.names(DATA) %in% row.names(outGWR),
                               c(coordCols, YRcol, IDcols)])
  
  if (SpatialObject) {
    
    # Convert to spatialpointsdataframe
    output <- SpatialPointsDataFrame(coords = output[ , coordCols],
                                     data = output,
                                     proj4string = CRS(projFOURstr))

  } 
  
  return(output)
  
}
# debug(Run_GeoGWR)

# The modified draw.gam function
mydraw.gam <- function (object, parametric = TRUE, select = NULL, scales = c("free", 
                                                                             "fixed"), align = "hv", axis = "lrtb", n = 100, unconditional = FALSE, 
                        overall_uncertainty = TRUE, dist = 0.1, ...) 
{
  scales <- match.arg(scales)
  S <- smooths(object)
  select <- gratia:::check_user_select_smooths(smooths = S, select = select)
  d <- gratia:::smooth_dim(object)
  take <- d <= 2L
  select <- select[take]
  S <- S[take]
  d <- d[take]
  is_re <- vapply(object[["smooth"]], gratia:::is_re_smooth, logical(1L))
  is_by <- vapply(object[["smooth"]], gratia:::is_by_smooth, logical(1L))
  if (any(is_by)) {
    S <- vapply(strsplit(S, ":"), `[[`, character(1L), 1L)
  }
  npara <- 0
  nsmooth <- length(S)
  if (isTRUE(parametric)) {
    terms <- parametric_terms(object)
    npara <- length(terms)
    p <- vector("list", length = npara)
  }
  g <- l <- vector("list", length = nsmooth)
  for (i in unique(S)) {
    eS <- evaluate_smooth(object, smooth = i, n = n, unconditional = unconditional, 
                          overall_uncertainty = overall_uncertainty, dist = dist)
    l[S == i] <- split(eS, eS[["smooth"]])
  }
  l <- l[select]
  d <- d[select]
  g <- g[select]
  if (length(g) == 0L) {
    message("Unable to draw any of the model terms.")
    return(invisible(g))
  }
  for (i in seq_along(l)) {
    g[[i]] <- draw(l[[i]])
  }
  if (isTRUE(parametric)) {
    for (i in seq_along(terms)) {
      p[[i]] <- evaluate_parametric_term(object, term = terms[i])
      g[[i + length(g)]] <- draw(p[[i]])
    }
  }
  if (isTRUE(identical(scales, "fixed"))) {
    wrapper <- function(x) {
      range(x[["est"]] + (2 * x[["se"]]), x[["est"]] - 
              (2 * x[["se"]]))
    }
    ylims <- range(unlist(lapply(l, wrapper)))
    if (isTRUE(parametric)) {
      ylims <- range(ylims, unlist(lapply(p, function(x) range(x[["upper"]], 
                                                               x[["lower"]]))))
    }
    gg <- seq_along(g)[c(d == 1L, rep(TRUE, npara))]
    for (i in gg) {
      g[[i]] <- g[[i]] + lims(y = ylims)
    }
  }
  g
}
