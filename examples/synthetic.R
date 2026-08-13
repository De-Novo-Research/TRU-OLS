script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
root <- dirname(dirname(normalizePath(script)))
source(file.path(root, "TRU-OLS.R"))

mixmat <- matrix(c(
  1, .2, .1,
  .1, 1, .2,
  .2, .1, 1,
  .3, .2, .4,
  .5, .3, .2
), 5, 3, byrow = TRUE)
namevec <- c("dye A", "dye B", "autofluorescence")
unstained <- rbind(
  c(.1, .3, 1), c(.2, .4, 1.2), c(.3, .5, 1.4),
  c(.4, .6, 1.6), c(.5, .7, 1.8), c(.6, .8, 2)
) %*% t(mixmat)
dataset <- rbind(
  c(8, .2, 2), c(.1, 7, 1.5), c(6, 5, 2.5), c(.2, .1, 1)
) %*% t(mixmat)

print(create_complete_dataframe(mixmat, namevec, dataset, unstained))
