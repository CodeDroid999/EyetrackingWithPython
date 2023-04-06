
#' Shows the raw samples and the detected fixations in an interactive
#' plot.  
#' 
#' @title Interactive Diagnostic Plot of Samples and Fixations
#' @param samples a data frame containing the raw samples as recorded
#' by the eye-tracker.  This data frame has to have the following
#' columns:
#' \describe{
#'  \item{time:}{the time at which the sample was recorded}
#'  \item{x:}{the x-coordinate of the sample}
#'  \item{y:}{the y-coordinate of the sample}
#' }
#' 
#' 
diagnostic.plot <- function(samples, fixations, start.event=1, start.time=NULL, duration=2000, interactive=TRUE, ...) {

  stopifnot(start.event >= 1)
  stopifnot(start.event <= nrow(fixations))
  stopifnot(start.time  >= min(samples$time))
  stopifnot(start.time  <= max(samples$time))

  if (is.null(start.time))
    start.time <- fixations$start[start.event]

  if (interactive) grDevices::dev.new()

  graphics::par(mar=c(2,2,0,0))
  if ("ylim" %in% list(...))
    with(samples,   plot(time, x, pch=20, cex=0.3, col="red",
                         ylim=c(min(fixations$x, fixations$y),
                                max(fixations$x, fixations$y)),
                         xlim=c(start.time, start.time+duration)))
  else
    with(samples,   plot(time, x, pch=20, cex=0.3, col="red",
                         xlim=c(start.time, start.time+duration), ...))
  with(samples,   points(time, y, pch=20, cex=0.3, col="orange"))
  with(fixations, lines(zip(start, end, NA), rep(x, each=3)))
  with(fixations, lines(zip(start, end, NA), rep(y, each=3)))
  with(fixations, abline(v=c(start, end), col="lightgrey"))
  
  if (interactive) {
    p <- grDevices::recordPlot()
    grDevices::dev.off()
    zoom::zm(rp=p)
  }

}

diagnostic.plot.event.types <- function(fixations) {
  graphics::pairs( ~ log10(mad.x+0.001) + log10(mad.y+0.001) + log10(dur),
        data=fixations,
        col=fixations$event)
}

#' Calculates summary statistics about the trials and fixations in the
#' given data frame.

calculate.summary <- function(fixations) {
  
  stats <- data.frame(mean=double(8), sd=double(8), row.names=c("Number of trials", "Duration of trials", "No. of fixations per trial", "Duration of fixations", "Dispersion horizontal", "Dispersion vertical", "Peak velocity horizontal", "Peak velocity vertical"))

  stats["Number of trials",] <- c(length(unique(fixations$trial)), NA)

  s <- tapply(fixations$start, fixations$trial, min)
  e <- tapply(fixations$end, fixations$trial, max)
  tdur <- e - s
  stats["Duration of trials",] <- c(mean(tdur), stats::sd(tdur))
    
  n <- tapply(fixations$start, fixations$trial, length)
  stats["No. of fixations per trial",] <- c(mean(n), stats::sd(n))
  stats["Duration of fixations",] <- c(mean(fixations$dur), stats::sd(fixations$dur))

  stats["Dispersion horizontal",] <- c(mean(fixations$mad.x, na.rm=TRUE), stats::sd(fixations$mad.x, na.rm=TRUE))
  stats["Dispersion vertical",]   <- c(mean(fixations$mad.y, na.rm=TRUE), stats::sd(fixations$mad.y, na.rm=TRUE))
          
  stats["Peak velocity horizontal",] <- c(mean(fixations$peak.vx, na.rm=T), stats::sd(fixations$peak.vx, na.rm=T))
  stats["Peak velocity vertical",]   <- c(mean(fixations$peak.vy, na.rm=T), stats::sd(fixations$peak.vy, na.rm=T))

  stats

}

zip <- function(...) {
  x <- cbind(...)
  t(x)[1:length(x)]
}