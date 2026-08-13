script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
root <- dirname(dirname(normalizePath(script)))
source(file.path(root, "TRU-OLS.R"))

tests <- 0L
expect <- function(value) {
  tests <<- tests + 1L
  if (!isTRUE(value)) stop("test failed", call. = FALSE)
}
close <- function(x, y, tolerance = 1e-10) {
  isTRUE(all.equal(x, y, tolerance = tolerance, check.attributes = FALSE))
}
errors <- function(expr, pattern) {
  message <- tryCatch({ force(expr); "" }, error = conditionMessage)
  grepl(pattern, message)
}

identity_mix <- diag(3)
control <- rbind(c(1, 10, 2), c(2, 20, 3), c(3, 30, 4), c(4, 40, 5))
expect(close(
  mean_unmix(identity_mix, control),
  mean_unmix(identity_mix, control, 99)
))
expect(close(
  mean_unmix(identity_mix, control)$cutoffs,
  apply(control, 2, quantile, probs = 0.99, type = 7, names = FALSE)
))
expect(errors(
  TRU_OLS(identity_mix, matrix(c(1, 1, 1), 1), c(1, 1), c("A", "B", "AF")),
  "cutoffs"
))

fit <- TRU_OLS(
  identity_mix,
  rbind(c(0.5, 0.2, 4), c(1, 1, 4)),
  c(1, 1, 1e6),
  c("A", "B", "AF")
)
expect(identical(fit$names[[1]], "AF"))
expect(close(fit$coefficients[[1]], 4))
expect(identical(fit$names[[2]], c("A", "B", "AF")))

multi_mix <- cbind(c(2, 0, 1), c(1, 0, 0), c(0, 1, 0))
multi_event <- matrix(as.double(multi_mix %*% c(-1, 2, 5)), 1)
multi <- TRU_OLS(
  multi_mix, multi_event, c(0, 1, 1e6), c("A", "B", "AF")
)
expect(identical(multi$names[[1]], "AF"))
expect(identical(names(multi$removed[[1]]), c("A", "B")))
expect(close(multi$coefficients[[1]], 5))

expect(close(.map_distribution(c(30, 10, 20), c(-2, -1, 0, 1, 2)),
             c(2, -2, 0)))
expect(close(.map_distribution(100, 1:4), 3))
expect(close(.map_distribution(c(2, 1, 1), c(10, 20, 30)), c(30, 10, 20)))
expect(close(.map_distribution(c(3, 1, 2), 7), c(7, 7, 7)))

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
  c(8, .2, 2), c(.1, 7, 1.5), c(6, 5, 2.5), c(.2, .1, 1)
)
unstained <- unstained_abundances %*% t(mixmat)
dataset <- event_abundances %*% t(mixmat)
complete <- create_complete_dataframe(
  mixmat, c("A dye", "B-dye", "AF"), dataset, unstained
)
expect(identical(names(complete), c("A dye", "B-dye", "AF")))
expect(all(is.finite(as.matrix(complete))))
expect(nrow(create_complete_dataframe(
  mixmat, c("A dye", "B-dye", "AF"), dataset[FALSE, ], unstained
)) == 0L)
expect(errors(
  create_complete_dataframe(mixmat, c("A", "B", "AF"), dataset,
                            unstained[FALSE, ]),
  "at least one"
))
expect(errors(
  TRU_OLS(cbind(1:3, 2 * (1:3)), matrix(1:3, 1), c(0, 0), c("A", "AF")),
  "full column rank"
))

cat("PASS:", tests, "tests\n")
