script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
root <- dirname(dirname(normalizePath(script)))
source(file.path(root, "TRU-OLS.R"))

args <- commandArgs(TRUE)
if (length(args) != 4L) {
  stop("usage: Rscript examples/complete_workflow.R MIX MULTICOLOR UNSTAINED OUTPUT")
}

mix_file <- read.csv(args[1], check.names = FALSE)
detectors <- as.character(mix_file[[1]])
mixmat <- as.matrix(mix_file[, 3:ncol(mix_file), drop = FALSE])
storage.mode(mixmat) <- "double"

read_events <- function(path) {
  frame <- read.csv(path, check.names = FALSE)
  value <- as.matrix(frame[, detectors, drop = FALSE])
  storage.mode(value) <- "double"
  value
}

result <- create_complete_dataframe(
  mixmat,
  colnames(mixmat),
  read_events(args[2]),
  read_events(args[3])
)
write.csv(result, args[4], row.names = FALSE)
