### Removal of debris, beads, doublets and deads 

# Too for visualizing a multi-step gating
visualize_strategy <- function(fcs, name, wd, channels, gates) {
  
  exprs <- fcs@exprs
  start_rows <- nrow(exprs)
  
  png(paste0(wd, '/pre_gating_strategy_', name, '.png'), width = 4*length(channels), height = 5, units = 'in', res = 300)
  par(mfrow = c(1,length(channels)), oma = c(0,0,5,0))
  
  for (g in 1:length(channels)) {
    message('Generating plot ', g)
    plot_biaxial(exprs, channels[[g]][1], channels[[g]][2], xleft = gates[[g]]$x[1], xright = gates[[g]]$x[2], ybottom = gates[[g]]$y[1], ytop = gates[[g]]$y[2])
    message('Gating for gate ', g)
    exprs <- remove_from_fcs(exprs, x.var = channels[[g]][1], y.var = channels[[g]][2], gates[[g]])
  }
  
  title(paste0('Pre-gating plot for ', name, ' - recovery ', round((nrow(exprs)/start_rows)*100, 2), ' %'), outer = T)
  dev.off()
  
  message('With the current settings, the recovery for this run is ', round((nrow(exprs)/start_rows)*100, 2), ' %.\n')
  
}


# Too for visualizing and APPLYING a multi-step gating
apply_strategy <- function(fcs, name, wd, channels, gates) {
  
  exprs <- fcs@exprs
  start_rows <- nrow(exprs)
  
  png(paste0(wd, '/pre_gating_strategy_', name, '.png'), width = 4*length(channels), height = 5, units = 'in', res = 300)
  par(mfrow = c(1,length(channels)), oma = c(0,0,5,0))
  
  for (g in 1:length(channels)) {
    
    plot_biaxial(exprs, channels[[g]][1], channels[[g]][2], xleft = gates[[g]]$x[1], xright = gates[[g]]$x[2], ybottom = gates[[g]]$y[1], ytop = gates[[g]]$y[2])
    exprs <- remove_from_fcs(exprs, x.var = channels[[g]][1], y.var = channels[[g]][2], gates[[g]])
  }
  
  title(paste0('Pre-gating plot for ', name, ' - recovery ', round((nrow(exprs)/start_rows)*100, 2), ' %'), outer = T)
  dev.off()
  
  
  fcs <- premessa::as_flowFrame(exprs, fcs)
  
  message('With the applied settings, the recovery for this run is ', round((nrow(exprs)/start_rows)*100, 2), ' %.\n')
  
  return(fcs)
  
}



# Biaxial plot
plot_biaxial <- function(m, x.var, y.var, xleft, xright, ybottom, ytop) {
  
  if (grepl("Di$|Dd$", x.var)) {
    m[, x.var] <- asinh(m[, x.var] / 5)
  } 
  
  if (grepl("Di$|Dd$", y.var)) {
    m[, y.var] <- asinh(m[, y.var] / 5)
  }  
  
  
  if (nrow(m) > 100000) {
   m <- m[sample(nrow(m), 100000, replace = F),]
  }
  
  # A little hack to avoid plots becoming to ugly for the event length variable
  if (x.var == 'Event_length') {
    m <- m[m[,'Event_length'] <= 150,]
  }
  
  k <- 11; my.cols <- rev(RColorBrewer::brewer.pal(k, "RdYlBu"))
  z <- MASS::kde2d(x=m[, x.var], y=m[, y.var], n=50, h=max(MASS::bandwidth.nrd(m[, x.var]), MASS::bandwidth.nrd(m[, y.var])))
  
  
  plot(x = m[, x.var], y = m[, y.var], col = scales::alpha('grey30', 0.5), pch = '.', xlab = x.var, ylab = y.var, cex = 1.2)
  contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)
  rect(xleft = xleft, xright = xright, ybottom = ybottom, ytop = ytop, border="red") # Adding a square to make a reference for the bead gating
  
}



# Function for removal from fcs
remove_from_fcs <- function(exprs, x.var, y.var, g) {
  m <- cbind(exprs[,c('Time', 'Event_length')], asinh(exprs[,grepl("Di$|Dd$", colnames(exprs))] / 5))
  
  keep <- (m[,x.var] > g$x[1]) & (m[,x.var] < g$x[2]) & (m[,y.var] > g$y[1]) & (m[,y.var] < g$y[2])
  
  keep_events <- exprs[keep,]
  
  return(keep_events)
}
