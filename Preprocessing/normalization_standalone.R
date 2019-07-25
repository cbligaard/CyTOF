### Normalization tool strongly inspired from the premessa implementation (some of the functions are even taken from there) and the CATALYST version

#### Helper functions ####
# A plotting function to help users identify beads in their files
plot_beads <- function(fcs, save = F, xleft = 4, xright = 9, ybottom = -0.2, ytop = 4, wd = '.') {
  exprs <- asinh(fcs@exprs / 5)
  
  # Downampling if necessary
  if (nrow(exprs) > 500000) {
    sample <- sample(nrow(exprs), 500000, replace = F)
  } else {
    sample <- 1:nrow(exprs)
  }
  
  if (save) {
    png(paste0(wd, '/bead_plot_', gsub('\\.fcs$', '\\.png', fcs@description$GUID)), width = 20, height = 4, units = 'in', res = 300)
  }
  
  par(mfrow = c(1,5))
  
  # Plot for all bead channels (assume EQ beads)
  for (bead_channel in c('Ce140', 'Eu151', 'Eu153', 'Ho165', 'Lu175')) {
    plot(x = exprs[sample,grep(bead_channel, colnames(exprs))], y = exprs[sample,grep("Ir193", colnames(exprs))], col = scales::alpha('grey30', 0.5), pch = '.', xlab = bead_channel, ylab = 'Ir193')
    rect(xleft = xleft, xright = xright, ybottom = ybottom, ytop = ytop, border="red") # Adding a square to make a reference for the bead gating
  }
  
  
  if (save) {
    dev.off()
  }
  
}



find_bead_channels <- function(fcs, bead.type = c("Fluidigm", "Beta")) {
  bead.type <- match.arg(bead.type)
  grep.string <- switch(bead.type,
                        Fluidigm = "Ce140|Eu151|Eu153|Ho165|Lu175",
                        Beta = "La139|Pr141|Tb159|Tm169|Lu175"
  )
  
  return(grep(grep.string, flowCore::parameters(fcs)$name))
}

find_dna_channel <- function(fcs) {
  return(grep("Ir193", flowCore::parameters(fcs)$name))
}

find_beads_channels_names <-  function(fcs, bead.type = c("Fluidigm", "Beta")) {
  beads.cols <- find_bead_channels(fcs, bead.type)
  return(get_parameter_name(fcs, beads.cols))
}

get_parameter_name <- function(fcs, i) {
  return(as.vector(unname(flowCore::parameters(fcs)$name[i])))
}


# Function that applies the beads gates to the data matrix and idenfies the beads
identify_beads_cbp <- function(m, gates, beads.cols.names, dna.col) {
  sel <- lapply(beads.cols.names, function(n) {
    g <- gates[[n]]
    ret <- (m[,n] > g$x[1]) & (m[,n] < g$x[2]) & (m[,dna.col] > g$y[1]) & (m[,dna.col] < g$y[2])
    
  })
  
  return(Reduce("&", sel))
}


# Calculate baseline for normalization
calculate_baseline_cbp <- function(fcs_list, bead_type, beads_gates) {

  message('Calculating baseline using the ', length(fcs_list), ' fcs files supplied.\n')
  
  ret <- lapply(1:length(fcs_list), function(batch) {
    
    # Get basic information
    exprs <- flowCore::exprs(fcs_list[[batch]])
    beads_cols_names <- find_beads_channels_names(fcs_list[[batch]], bead_type)
    dna_col <- find_dna_channel(fcs_list[[batch]])
    
    
    m.transformed <- asinh(exprs / 5)
    sel <- identify_beads_cbp(m.transformed, beads_gates[[batch]], beads_cols_names, dna_col)
    
    message('Identified beads in file.\n')
    
    # The gating is done on transformed data, but we return the untransformed
    return(exprs[sel, beads_cols_names])
  })
  
  # Getting column medians
  ret <- Reduce("rbind", ret)
  ret <- apply(ret, 2, median)
  return(ret)
}

smooth_beads <- function(beads.data, window.size = 201) {
  ret <- apply(beads.data, 2, runmed, k = window.size)
  return(ret)
  
}

compute_bead_slopes <- function(beads.data, baseline, beads.cols.names) {
  beads.data <- beads.data[, beads.cols.names]
  baseline <- baseline[beads.cols.names]
  
  stopifnot(length(baseline) == dim(beads.data)[2])
  m <- t(t(beads.data) * baseline)
  
  #Minimizing residual sum of squares, assuming a no-intercept model
  ret <- rowSums(m) / rowSums(beads.data ^ 2)
  return(ret)
}


# Performing correction in data
correct_data_channels_cbp <- function(m, beads.data, baseline, beads.cols.names, time.col.name = "Time") {
  sm.beads <- smooth_beads(beads.data)
  slopes <- compute_bead_slopes(sm.beads, baseline, beads.cols.names)
  
  tt <- beads.data[, time.col.name]
  y <- slopes
  int.slopes <- approx(tt, y, m[, time.col.name], method = "linear", rule = 2) # If rule is 2, the value at the closest data extreme is used (extrapolation)
  
  include.channels <- grepl("Di$|Dd$", colnames(m))
  
  
  message(sprintf("Modfying channels: %s", paste(colnames(m)[include.channels], collapse = ", ")))
  message(sprintf("Leaving alone channels: %s", paste(colnames(m)[!include.channels], collapse = ", ")))
  
  m[, include.channels] <- m[, include.channels] * int.slopes$y
  beads.slopes <- cbind(time = tt, slope = slopes)
  beads.normed <- beads.data
  beads.normed[, beads.cols.names] <- beads.normed[, beads.cols.names] * beads.slopes[, "slope"]
  
  return(list(m.normed = m, beads.smoothed = sm.beads,
              beads.normed = beads.normed, beads.slopes = beads.slopes))
}


# Calculates the Mahalanobis distance of each event from the centroid of the beads population
get_mahalanobis_distance_from_beads <- function(m, beads.events, beads.cols.names) {
  m <- m[, beads.cols.names]
  m <- asinh(m / 5)
  beads.data <- m[beads.events,]
  
  cov.m <- cov(beads.data)
  ret <- sqrt(mahalanobis(m, colMeans(beads.data), cov.m))
  return(ret)
}


# Plotting bead medians across files
plot_beads_medians <- function(tab, x.var = "sample") {
  
  m <- reshape::melt.data.frame(tab, id.vars = c(x.var, "type"))
  
  m$value <- asinh(m$value / 5)
  
  (p <- ggplot2::ggplot(ggplot2::aes_string(x = x.var, y = "value", color = "variable", group = "variable"), data = m)
    + ggplot2::facet_wrap(~type, ncol = 1)
    + ggplot2::geom_line()
    + ggplot2::geom_point(alpha = as.numeric(x.var == 'sample'))
    + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    + ggplot2::scale_y_continuous("Median intensity (asinh transformed)")
  )
}



# Plotting beads over time
plot_beads_over_time <- function(beads.data, beads.normed, beads.cols) {
  beads.data <- data.frame(beads.data, type = "Before", check.names = F, stringsAsFactors =  F)
  beads.normed <- data.frame(beads.normed, type = "After", check.names =  F, stringsAsFactors =  F)
  m <- rbind(beads.data, beads.normed)
  m <- data.frame(m[, beads.cols], Time = m$Time, type = factor(m$type, levels = c('Before', 'After')), check.names = F, stringsAsFactors = F)
  
  plot_beads_medians(m, "Time")
}

# Normalization main function - only allows normalization to the files in the same run - for normalization to another run, see premessa
normalizeFCS_cbp <- function(fcs_list, bead_type, beads_gates, wd) {
  
  # Calculate baseline
  baseline <- calculate_baseline_cbp(fcs_list, bead_type, beads_gates)
  message('These are your specified gates.')
  message(beads_gates)
  
  
  # Going through the different files to do correction
  ll <- lapply(1:length(fcs_list), function(batch) {
    
    fcs <- fcs_list[[batch]]
    message('Starting to process file ', fcs@description$GUID, '\n')
    
    beads.cols <- find_bead_channels(fcs, bead_type)
    beads.cols.names <- get_parameter_name(fcs, beads.cols)
    dna.col <- find_dna_channel(fcs)
    
    m <- flowCore::exprs(fcs)
    beads.events <- identify_beads_cbp(asinh(m / 5), beads_gates[[batch]], beads.cols.names, dna.col)
    beads.data <- m[beads.events,]
    
    norm.res <- correct_data_channels_cbp(m, beads.data, baseline, beads.cols.names)
    
    # Do some cleanup to save memory
    rm(m)
    
    # get final data for normalized file
    m.normed <- norm.res$m.normed
    m.normed <- cbind(m.normed,
                      beadDist = get_mahalanobis_distance_from_beads(m.normed, beads.events, beads.cols.names))
    
    m.normed <- premessa::as_flowFrame(m.normed, fcs)

    # Calculating normed and smoothed beads    
    message('Calculating normed and smoothed beads.\n')
    beads.normed <- apply(norm.res$beads.normed[, beads.cols], 2, median)
    beads.smoothed <- apply(norm.res$beads.smoothed[, beads.cols], 2, median)
    
    # Getting a plot for beads over time
    p <- plot_beads_over_time(norm.res$beads.smoothed, smooth_beads(norm.res$beads.normed), beads.cols)
    ggplot2::ggsave(paste0(wd, '/', gsub('\\.fcs', '', fcs@description$GUID), '_beads.pdf'), plot = p, width = 11, height = 8.5, units = "in")
    
    return(list(beads.normed = beads.normed, beads.smoothed = beads.smoothed, m.normed = m.normed))
  })
  
  message('Processed all the files. Now plotting overview.\n')
  beads.medians <- t(sapply(ll, function(x) {return(x[["beads.normed"]])}))
  beads.medians <- rbind(beads.medians,
                         t(sapply(ll, function(x) {return(x[["beads.smoothed"]])})))
  
  beads.medians <- data.frame(beads.medians, row.names = NULL)
  
  beads.medians$sample <- rep(1:length(fcs_list), 2)
  beads.medians$type <- c(rep("After", length(ll)), rep("Before", length(ll)))
  beads.medians$type <- factor(beads.medians$type, levels = c("Before", "After"))
  p <- plot_beads_medians(beads.medians)
  ggplot2::ggsave(paste0(wd, "/beads_before_and_after.pdf"), plot = p, width = 11, height = 8.5, units = "in")
  
  
  return(lapply(ll, function(x) {return(x[["m.normed"]])}))
}
