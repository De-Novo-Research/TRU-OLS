# TRU-OLS

Julia and R implementations of threshold-refined ordinary least squares for
spectral flow-cytometry unmixing.  The Julia version should be considered more reliable, and has recently
been updated to fix small bugs.
Please note that the R implementation was AI
generated using ChatGPT-5.6 Sol.  We are providing it for ease of use due to requests
but with very minimal testing. 
If you have any bugs, please send them to ryan.kmet@denovoresearch.org, though
we cannot promise that we will resolve them.

`mixmat` is detectors by endmembers. Event and unstained-control matrices are
events by detectors. The last mixing-matrix column is treated as
autofluorescence by default: it is excluded from baseline subtraction and is
retained during threshold refinement. A different column can be selected with
`autofluorescence`.

The cutoff percentile uses a 0-100 scale and defaults to the 99th percentile.

## R

```r
source("TRU-OLS.R")

# Read in a mixing matrix, an unstained dataset, and a multicolor dataset.
mixmat_path <- "mix.csv"
unstained_path <- "unstained.csv"
multicolor_path <- "multicolor.csv"

mix_df <- read.csv(mixmat_path, check.names = FALSE)

# Get detector and endmember names, then cast the mixing matrix.
detectors <- as.character(mix_df[[1]])
endmember_names <- names(mix_df)[3:ncol(mix_df)]
mixmat <- as.matrix(mix_df[, 3:ncol(mix_df), drop = FALSE])
storage.mode(mixmat) <- "double"

# Read event data in mixing-matrix detector order.
read_events <- function(path) {
  frame <- read.csv(path, check.names = FALSE)
  values <- as.matrix(frame[, detectors, drop = FALSE])
  storage.mode(values) <- "double"
  values
}

unstained <- read_events(unstained_path)
multicolor <- read_events(multicolor_path)

# Invoke TRU-OLS using control matching and a 99th-percentile cutoff.
result <- create_complete_dataframe(
  mixmat,
  endmember_names,
  multicolor,
  unstained,
  match = TRUE,
  percentile_cutoff = 99
)

# Write the complete unmixed result.
write.csv(result, "unmixed.csv", row.names = FALSE)
```

No contributed R packages are required. The CSV command-line workflow expects
the mixing file to contain detector names in column 1, metadata in column 2,
and endmembers in the remaining columns. The event files must contain columns
with those detector names.

```sh
Rscript examples/complete_workflow.R mix.csv multicolor.csv unstained.csv output.csv
```

## Julia

Start Julia from the repository root:

```sh
julia --project=.
```

```julia
using Pkg
Pkg.instantiate()
using CSV, DataFrames
include("TRU-OLS.jl")

# Read in a mixing matrix, an unstained dataset, and a multicolor dataset.
mixmat_path = "mix.csv"
unstained_path = "unstained.csv"
multicolor_path = "multicolor.csv"

mix_df = CSV.read(mixmat_path, DataFrame)

# Get detector and endmember names, then cast the mixing matrix.
detectors = String.(mix_df[:, 1])
endmember_names = names(mix_df)[3:end]
mixmat = Matrix{Float64}(mix_df[:, 3:end])

# Read event data in mixing-matrix detector order.
unstained_df = CSV.read(unstained_path, DataFrame)
multicolor_df = CSV.read(multicolor_path, DataFrame)
unstained = Matrix{Float64}(unstained_df[:, detectors])
multicolor = Matrix{Float64}(multicolor_df[:, detectors])

# Invoke TRU-OLS using control matching and a 99th-percentile cutoff.
result = create_complete_dataframe(
    mixmat,
    endmember_names,
    multicolor,
    unstained;
    match=true,
    percentile_cutoff=99.0
)

# Write the complete unmixed result.
CSV.write("unmixed.csv", result)
```

The included synthetic example can be run with
`julia --project=. examples/synthetic.jl`.

## Tests

```sh
Rscript tests/test.R
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. tests/reference.jl tests/reference
Rscript tests/parity.R tests/reference
```

The parity test compares the full workflow, including 99th-percentile
cutoffs, iterative refitting, autofluorescence retention, and distribution
matching.

## Reference

Kmet R, Novo D. Reducing Spreading: Removing the Impact of Irrelevant Dyes
Improves Unmixed Flow Cytometry Data. *Cytometry Part A*. 2025;107(9):573-586.
[doi:10.1002/cyto.a.24957](https://doi.org/10.1002/cyto.a.24957)

## License

GNU General Public License v3.0. See `LICENSE`.
