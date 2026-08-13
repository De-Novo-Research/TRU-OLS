# SPDX-License-Identifier: GPL-3.0-only

.numeric_matrix <- function(x, name) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x) || !is.numeric(x) || any(!is.finite(x))) {
    stop(name, " must be a finite numeric matrix", call. = FALSE)
  }
  storage.mode(x) <- "double"
  x
}

.check_inputs <- function(mixmat, dataset, namevec, cutoffs = NULL) {
  if (!ncol(mixmat)) stop("mixmat must not be empty", call. = FALSE)
  if (nrow(mixmat) != ncol(dataset)) {
    stop("mixmat rows must match dataset columns", call. = FALSE)
  }
  if (nrow(mixmat) < ncol(mixmat)) {
    stop("mixmat must have at least as many rows as columns", call. = FALSE)
  }
  singular_values <- svd(mixmat, nu = 0, nv = 0)$d
  tolerance <- min(dim(mixmat)) * .Machine$double.eps * max(singular_values)
  if (sum(singular_values > tolerance) < ncol(mixmat)) {
    stop("mixmat must have full column rank", call. = FALSE)
  }
  if (length(namevec) != ncol(mixmat) || anyNA(namevec) ||
      any(!nzchar(namevec)) || anyDuplicated(namevec)) {
    stop("namevec must be unique and match mixmat columns", call. = FALSE)
  }
  if (!is.null(cutoffs) && length(cutoffs) != ncol(mixmat)) {
    stop("cutoffs must match mixmat columns", call. = FALSE)
  }
}

.left_divide <- function(a, b) {
  answer <- if (nrow(a) == ncol(a)) {
    solve(a, b, tol = 0)
  } else {
    qr.coef(qr(a, LAPACK = TRUE), b)
  }
  if (any(!is.finite(answer))) stop("mixmat is rank deficient", call. = FALSE)
  answer
}

mean_unmix <- function(mixmat, observations, percentile_cutoff = 99) {
  # This function takes in a mixing matrix, a matrix of observations
  # (observations as rows), and a percentile.  It returns the mean unmixed
  # abundances and the percentile cutoff for each endmember.  Run this with a
  # full mixing matrix and unstained control observations in order to estimate
  # endmember means used for baseline subtraction and the cutoffs for TRU-OLS.
  #
  # Parameters:
  # - mixmat: The full mixing matrix. (detectors x endmembers)
  # - observations: The unstained data matrix. (events x detectors)
  # - percentile_cutoff: Cutoff percentile on a 0-100 scale.  Defaults to 99.
  #
  # Return:
  # - means: Mean unmixed abundance for each endmember
  # - cutoffs: Cutoff values for each endmember
  mixmat <- .numeric_matrix(mixmat, "mixmat")
  observations <- .numeric_matrix(observations, "observations")
  .check_inputs(mixmat, observations, as.character(seq_len(ncol(mixmat))))
  if (!nrow(observations)) {
    stop("observations must contain at least one row", call. = FALSE)
  }
  if (length(percentile_cutoff) != 1L || !is.finite(percentile_cutoff) ||
      percentile_cutoff < 0 || percentile_cutoff > 100) {
    stop("percentile_cutoff must be between 0 and 100", call. = FALSE)
  }

  abundances <- t(.left_divide(mixmat, t(observations)))
  means <- colSums(abundances) / nrow(abundances)
  cutoffs <- vapply(
    seq_len(ncol(abundances)),
    function(j) quantile(
      abundances[, j], percentile_cutoff / 100,
      type = 7, names = FALSE
    ),
    numeric(1)
  )
  list(means = means, cutoffs = cutoffs)
}

TRU_OLS <- function(mixmat, dataset, cutoffs, namevec,
                    autofluorescence = length(namevec)) {
  # This function takes in a mixing matrix, a dataset matrix, a vector of
  # cutoffs, and a vector of names.  The vector of names are the names of the
  # endmembers in the mixing matrix.  Must be in same order as endmember
  # columns.  This function runs TRU-OLS regression.  It returns the abundances
  # that are relevant for each event in the dataset as a list of lists,
  # corresponding names as a separate list of lists, and removed endmembers with
  # their unmixed values.
  #
  # Parameters:
  # - mixmat: The mixing matrix (detectors x endmembers)
  # - dataset: The dataset matrix (events x detectors)
  # - cutoffs: Threshold vector for each endmember
  # - namevec: Names of endmembers
  # - autofluorescence: Endmember column retained during threshold refinement
  #
  # Return:
  # - coefficients: A list of relevant unmixed values for each event
  # - names: A list of names of endmembers with relevant unmixed values
  # - removed: A list associating events with removed endmembers and values
  mixmat <- .numeric_matrix(mixmat, "mixmat")
  dataset <- .numeric_matrix(dataset, "dataset")
  namevec <- as.character(namevec)
  cutoffs <- as.double(cutoffs)
  .check_inputs(mixmat, dataset, namevec, cutoffs)
  if (any(!is.finite(cutoffs))) stop("cutoffs must be finite", call. = FALSE)
  if (length(autofluorescence) != 1L || !autofluorescence %in% seq_along(namevec)) {
    stop("autofluorescence must index a mixmat column", call. = FALSE)
  }

  coefficients <- vector("list", nrow(dataset))
  retained_names <- vector("list", nrow(dataset))
  removed <- vector("list", nrow(dataset))

  for (i in seq_len(nrow(dataset))) {
    active <- seq_len(ncol(mixmat))
    event_removed <- setNames(numeric(0), character(0))
    repeat {
      estimate <- as.double(
        .left_divide(mixmat[, active, drop = FALSE], dataset[i, ])
      )
      excluded <- which(
        active != autofluorescence & estimate < cutoffs[active]
      )
      if (!length(excluded)) {
        coefficients[[i]] <- estimate
        retained_names[[i]] <- namevec[active]
        break
      }
      event_removed[namevec[active[excluded]]] <- estimate[excluded]
      active <- active[-excluded]
    }
    removed[[i]] <- event_removed
  }

  list(coefficients = coefficients, names = retained_names, removed = removed)
}

.map_distribution <- function(values, control) {
  # This function takes in two vectors of values.  The first is a vector of all
  # irrelevant abundances for a single endmember over a dataset.  The second is
  # a vector of unmixed control abundances.  The output is a vector with the
  # irrelevant abundances replaced with their percentile match from the control.
  #
  # Parameters:
  # - values: Irrelevant unmixed values for a single endmember over a dataset
  # - control: Unmixed control abundances for the same endmember
  #
  # Return:
  # - values: Percentile-matched control values to replace irrelevant data
  values <- as.double(values)
  control <- as.double(control)
  if (!length(values)) return(values)
  if (!length(control)) stop("control must not be empty", call. = FALSE)

  sorted_control <- sort(control, method = "radix")
  order <- order(values, method = "radix")
  if (length(values) == 1L) {
    values[1L] <- sorted_control[round((length(control) - 1L) / 2) + 1L]
    return(values)
  }

  for (rank in seq_along(order)) {
    p <- (rank - 1L) / (length(values) - 1L)
    control_index <- round(p * (length(control) - 1L)) + 1L
    values[order[rank]] <- sorted_control[control_index]
  }
  values
}

create_complete_dataframe <- function(mixmat, namevec, dataset, unstained,
                                      match = TRUE,
                                      percentile_cutoff = 99,
                                      autofluorescence = length(namevec)) {
  # This function creates a complete dataframe of unmixed values where columns
  # are the original endmembers.  For endmembers that survived TRU_OLS, the
  # refitted coefficients are used.  For endmembers removed during TRU_OLS,
  # the values are matched to the unstained control or set to zero.
  #
  # Parameters:
  # - mixmat: The mixing matrix (detectors x endmembers)
  # - namevec: List of endmember names
  # - dataset: The dataset matrix (events x detectors)
  # - unstained: Unstained data (events x detectors)
  # - match: Boolean to either match unstained control or not
  # - percentile_cutoff: Cutoff percentile on a 0-100 scale.  Defaults to 99.
  # - autofluorescence: Endmember column excluded from baseline subtraction and
  #   retained during threshold refinement
  #
  # Returns:
  # - A dataframe containing retained unmixed values and matched or zero values
  mixmat <- .numeric_matrix(mixmat, "mixmat")
  dataset <- .numeric_matrix(dataset, "dataset")
  unstained <- .numeric_matrix(unstained, "unstained")
  namevec <- as.character(namevec)
  .check_inputs(mixmat, dataset, namevec)
  .check_inputs(mixmat, unstained, namevec)
  if (!is.logical(match) || length(match) != 1L || is.na(match)) {
    stop("match must be TRUE or FALSE", call. = FALSE)
  }
  if (length(autofluorescence) != 1L || !autofluorescence %in% seq_along(namevec)) {
    stop("autofluorescence must index a mixmat column", call. = FALSE)
  }

  calibration <- mean_unmix(mixmat, unstained, percentile_cutoff)
  baseline_matrix <- mixmat
  baseline_matrix[, autofluorescence] <- 0
  adjusted <- sweep(
    dataset, 2L, as.double(baseline_matrix %*% calibration$means), "-"
  )
  fit <- TRU_OLS(
    mixmat, adjusted, calibration$cutoffs, namevec,
    autofluorescence = autofluorescence
  )

  result <- matrix(
    0, nrow(dataset), length(namevec), dimnames = list(NULL, namevec)
  )
  name_index <- setNames(seq_along(namevec), namevec)
  for (i in seq_len(nrow(dataset))) {
    if (!length(fit$names[[i]])) next
    result[i, name_index[fit$names[[i]]]] <- fit$coefficients[[i]]
  }

  if (match) {
    control_abundances <- t(.left_divide(mixmat, t(unstained)))
    for (j in seq_along(namevec)) {
      events <- which(vapply(
        fit$removed, function(x) namevec[j] %in% names(x), logical(1)
      ))
      if (!length(events)) next
      values <- vapply(
        fit$removed[events], function(x) unname(x[[namevec[j]]]), numeric(1)
      )
      result[events, j] <- .map_distribution(values, control_abundances[, j])
    }
  }

  as.data.frame(result, check.names = FALSE)
}
