using CSV
using DataFrames

include(joinpath(@__DIR__, "..", "TRU-OLS.jl"))

length(ARGS) == 1 || error("usage: julia --project=. tests/reference.jl OUTPUT")
output = ARGS[1]
mkpath(output)

namevec = ["A dye", "B-dye", "AF"]
mixmat = [
    1.0 .2 .1;
    .1 1.0 .2;
    .2 .1 1.0;
    .3 .2 .4;
    .5 .3 .2
]
unstained_abundances = [
    .1 .3 1.0;
    .2 .4 1.2;
    .3 .5 1.4;
    .4 .6 1.6;
    .5 .7 1.8;
    .6 .8 2.0
]
event_abundances = [
    8.0 .2 2.0;
    .1 7.0 1.5;
    6.0 5.0 2.5;
    .2 .1 1.0;
    3.0 .3 2.0;
    .4 4.0 2.0
]
unstained = unstained_abundances * transpose(mixmat)
dataset = event_abundances * transpose(mixmat)

default_calibration = mean_unmix(mixmat, unstained)
explicit_calibration = mean_unmix(
    mixmat, unstained; percentile_cutoff=99.0
)
@assert default_calibration == explicit_calibration
@assert _map_distribution!([100.0], [1.0, 2.0, 3.0, 4.0]) == [3.0]

fit = TRU_OLS(mixmat, dataset, [1.0, 1.0, 1e6], namevec)
@assert all(names -> "AF" in names, fit.names)

core = zeros(size(dataset, 1), length(namevec))
survivors = falses(size(dataset, 1), length(namevec))
for i in axes(core, 1), j in eachindex(fit.names[i])
    column = findfirst(==(fit.names[i][j]), namevec)
    core[i, column] = fit.coefficients[i][j]
    survivors[i, column] = true
end

complete = create_complete_dataframe(mixmat, namevec, dataset, unstained)
CSV.write(joinpath(output, "calibration.csv"), DataFrame(
    endmember=namevec,
    mean=default_calibration.means,
    cutoff=default_calibration.cutoffs
))
CSV.write(joinpath(output, "core.csv"), DataFrame(core, namevec))
CSV.write(joinpath(output, "survivors.csv"), DataFrame(Int.(survivors), namevec))
CSV.write(joinpath(output, "complete.csv"), complete)
