include(joinpath(@__DIR__, "..", "TRU-OLS.jl"))

mixmat = [
    1.0 .2 .1;
    .1 1.0 .2;
    .2 .1 1.0;
    .3 .2 .4;
    .5 .3 .2
]
namevec = ["dye A", "dye B", "autofluorescence"]
unstained = [
    .1 .3 1.0;
    .2 .4 1.2;
    .3 .5 1.4;
    .4 .6 1.6;
    .5 .7 1.8;
    .6 .8 2.0
] * transpose(mixmat)
dataset = [
    8.0 .2 2.0;
    .1 7.0 1.5;
    6.0 5.0 2.5;
    .2 .1 1.0
] * transpose(mixmat)

println(create_complete_dataframe(mixmat, namevec, dataset, unstained))
