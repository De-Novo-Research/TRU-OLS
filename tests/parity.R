script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
root <- dirname(dirname(normalizePath(script)))
source(file.path(root, "TRU-OLS.R"))

args <- commandArgs(TRUE)
if (length(args) != 1L) stop("usage: Rscript tests/parity.R REFERENCE")

namevec <- c("A dye", "B-dye", "AF")
mixmat <- matrix(c(
  1, .2, .1,
  .1, 1, .2,
  .2, .1, 1,
  .3, .2, .4,
  .5, .3, .2
), 5, 3, byrow = TRUE)
unstained_abundances <- rbind(
  c(.1, .3, 1), c(.2, .4, 1.2), c(.3, .5, 1.4),
  c(.4, .6, 1.6), c(.5, .7, 1.8), c(.6, .8, 2)
)
event_abundances <- rbind(
  c(8, .2, 2), c(.1, 7, 1.5), c(6, 5, 2.5),
  c(.2, .1, 1), c(3, .3, 2), c(.4, 4, 2)
)
unstained <- unstained_abundances %*% t(mixmat)
dataset <- event_abundances %*% t(mixmat)

fit <- TRU_OLS(mixmat, dataset, c(1, 1, 1e6), namevec)
core <- matrix(0, nrow(dataset), length(namevec), dimnames = list(NULL, namevec))
survivors <- matrix(
  FALSE, nrow(dataset), length(namevec), dimnames = list(NULL, namevec)
)
for (i in seq_along(fit$names)) {
  columns <- match(fit$names[[i]], namevec)
  core[i, columns] <- fit$coefficients[[i]]
  survivors[i, columns] <- TRUE
}
calibration <- mean_unmix(mixmat, unstained)
complete <- create_complete_dataframe(mixmat, namevec, dataset, unstained)

read_reference <- function(file) read.csv(
  file.path(args[1], file), check.names = FALSE
)
compare <- function(x, y) {
  x <- as.double(unlist(x, use.names = FALSE))
  y <- as.double(unlist(y, use.names = FALSE))
  if (length(x) != length(y) || any(!is.finite(c(x, y))) ||
      any(abs(x - y) > 1e-9 + 1e-12 * pmax(abs(x), abs(y)))) {
    stop("Julia and R results differ", call. = FALSE)
  }
}

julia_calibration <- read_reference("calibration.csv")
compare(calibration$means, julia_calibration$mean)
compare(calibration$cutoffs, julia_calibration$cutoff)
compare(core, as.matrix(read_reference("core.csv")))
compare(complete, read_reference("complete.csv"))
reference_survivors <- as.matrix(read_reference("survivors.csv"))
if (!identical(dim(survivors), dim(reference_survivors)) ||
    !all(survivors == reference_survivors)) {
  stop("survivor decisions differ", call. = FALSE)
}
if (!identical(names(complete), names(read_reference("complete.csv")))) {
  stop("column order differs", call. = FALSE)
}

cat("PASS: Julia/R parity\n")
