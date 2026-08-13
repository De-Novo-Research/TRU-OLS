# SPDX-License-Identifier: GPL-3.0-only

using DataFrames
using LinearAlgebra
using StatsBase

function _check_inputs(mixmat, dataset, namevec; cutoffs=nothing)
    size(mixmat, 2) > 0 || throw(ArgumentError("mixmat must not be empty"))
    size(mixmat, 1) == size(dataset, 2) || throw(DimensionMismatch(
        "mixmat rows must match dataset columns"
    ))
    size(mixmat, 1) >= size(mixmat, 2) || throw(DimensionMismatch(
        "mixmat must have at least as many rows as columns"
    ))
    rank(mixmat) == size(mixmat, 2) || throw(ArgumentError(
        "mixmat must have full column rank"
    ))
    length(namevec) == size(mixmat, 2) || throw(DimensionMismatch(
        "namevec must match mixmat columns"
    ))
    all(name -> !isempty(name), namevec) || throw(ArgumentError(
        "namevec must not be empty"
    ))
    length(unique(namevec)) == length(namevec) || throw(ArgumentError(
        "namevec must be unique"
    ))
    all(isfinite, mixmat) && all(isfinite, dataset) || throw(ArgumentError(
        "inputs must be finite"
    ))
    if cutoffs !== nothing
        length(cutoffs) == size(mixmat, 2) || throw(DimensionMismatch(
            "cutoffs must match mixmat columns"
        ))
        all(isfinite, cutoffs) || throw(ArgumentError("cutoffs must be finite"))
    end
end

function mean_unmix(mixmat::AbstractMatrix{<:Real},
                    observations::AbstractMatrix{<:Real};
                    percentile_cutoff::Real=99.0)
    """
    This function takes in a mixing matrix, a matrix of observations (observations as rows),
    and a percentile.  It returns the mean unmixed abundances and the percentile cutoff for each endmember.
    Run this with a full mixing matrix and unstained control observations in order to estimate endmember means used for baseline subtraction
    and the cutoffs for TRU-OLS.

    Parameters:
    - mixmat: The full mixing matrix. (detectors x endmembers)
    - observations: The unstained data matrix. (events x detectors)
    - percentile_cutoff: Cutoff percentile on a 0-100 scale.  Defaults to 99.

    Return:
    - means: Mean unmixed abundance for each endmember
    - cutoffs: Cutoff values for each endmember
    """
    mixmat = Matrix{Float64}(mixmat)
    observations = Matrix{Float64}(observations)
    _check_inputs(mixmat, observations, string.(1:size(mixmat, 2)))
    size(observations, 1) > 0 || throw(ArgumentError(
        "observations must contain at least one row"
    ))
    0 <= percentile_cutoff <= 100 || throw(ArgumentError(
        "percentile_cutoff must be between 0 and 100"
    ))

    abundances = transpose(mixmat \ transpose(observations))
    means = vec(sum(abundances, dims=1)) ./ size(abundances, 1)
    cutoffs = [percentile(view(abundances, :, j), percentile_cutoff)
               for j in axes(abundances, 2)]
    (means=means, cutoffs=cutoffs)
end

function TRU_OLS(mixmat::AbstractMatrix{<:Real},
                 dataset::AbstractMatrix{<:Real},
                 cutoffs::AbstractVector{<:Real},
                 namevec::AbstractVector{<:AbstractString};
                 autofluorescence::Integer=length(namevec))
    """
    This function takes in a mixing matrix, a dataset matrix, a vector of cutoffs, and a vector of names.  The vector
    of names are the names of the endmembers in the mixing matrix.  Must be in same order as endmember columns.
    This function runs TRU-OLS regression.  It returns the abundances that are relevant for each event in the dataset as a list of lists,
    corresponding names as a separate list of lists, and removed endmembers with their unmixed values.

    Parameters:
    - mixmat: The mixing matrix (detectors x endmembers)
    - dataset: The dataset matrix (events x detectors)
    - cutoffs: Threshold vector for each endmember
    - namevec: Names of endmembers
    - autofluorescence: Endmember column retained during threshold refinement

    Return:
    - coefficients: A list of relevant unmixed values for each event
    - names: A list of names of endmembers with relevant unmixed values for each event
    - removed: A list associating each event with irrelevant endmembers and their unmixed values
    """
    mixmat = Matrix{Float64}(mixmat)
    dataset = Matrix{Float64}(dataset)
    cutoffs = Vector{Float64}(cutoffs)
    _check_inputs(mixmat, dataset, namevec; cutoffs=cutoffs)
    autofluorescence in eachindex(namevec) || throw(ArgumentError(
        "autofluorescence must index a mixmat column"
    ))

    coefficients = Vector{Vector{Float64}}(undef, size(dataset, 1))
    retained_names = Vector{Vector{String}}(undef, size(dataset, 1))
    removed = Vector{Dict{String, Float64}}(undef, size(dataset, 1))

    for i in axes(dataset, 1)
        active = collect(axes(mixmat, 2))
        event_removed = Dict{String, Float64}()
        observation = Vector(dataset[i, :])

        while true
            estimate = Vector{Float64}(mixmat[:, active] \ observation)
            excluded = findall(j ->
                active[j] != autofluorescence &&
                estimate[j] < cutoffs[active[j]],
                eachindex(active)
            )
            if isempty(excluded)
                coefficients[i] = estimate
                retained_names[i] = String.(namevec[active])
                break
            end
            for j in excluded
                event_removed[String(namevec[active[j]])] = estimate[j]
            end
            deleteat!(active, excluded)
        end
        removed[i] = event_removed
    end

    (coefficients=coefficients, names=retained_names, removed=removed)
end

function _map_distribution!(values::Vector{Float64}, control::Vector{Float64})
    """
    This function takes in two vectors of values.  The first is a vector of all irrelevant abundances for a single endmember over a dataset.
    The second is a vector of unmixed control abundances.  The output is a vector with the irrelevant abundances replaced with
    their percentile match from the control.

    Parameters:
    - values: Irrelevant unmixed values for a single endmember over a whole dataset
    - control: Unmixed control abundances for the same endmember

    Return:
    - values: Percentile-matched control values to replace irrelevant data
    """
    isempty(values) && return values
    isempty(control) && throw(ArgumentError("control must not be empty"))

    sorted_control = sort(control)
    order = sortperm(values; alg=MergeSort)
    if length(values) == 1
        values[1] = sorted_control[round(Int, (length(control) - 1) / 2) + 1]
        return values
    end

    for rank in eachindex(order)
        p = (rank - 1) / (length(values) - 1)
        control_index = round(Int, p * (length(control) - 1)) + 1
        values[order[rank]] = sorted_control[control_index]
    end
    values
end

function create_complete_dataframe(mixmat::AbstractMatrix{<:Real},
                                   namevec::AbstractVector{<:AbstractString},
                                   dataset::AbstractMatrix{<:Real},
                                   unstained::AbstractMatrix{<:Real};
                                   match::Bool=true,
                                   percentile_cutoff::Real=99.0,
                                   autofluorescence::Integer=length(namevec))
    """
    This function creates a complete dataframe of unmixed values where columns are the original endmembers.
    For endmembers that survived TRU_OLS, the refitted coefficients are used.
    For endmembers removed during TRU_OLS, the values are matched to the unstained control or set to zero.

    Parameters:
    - mixmat: The mixing matrix (detectors x endmembers)
    - namevec: List of endmember names
    - dataset: The dataset matrix (events x detectors)
    - unstained: Unstained data (events x detectors)
    - match: Boolean to either match unstained control or not
    - percentile_cutoff: Cutoff percentile on a 0-100 scale.  Defaults to 99.
    - autofluorescence: Endmember column excluded from baseline subtraction and retained during threshold refinement

    Returns:
    - A dataframe with all original endmembers, containing retained unmixed values and matched or zero replacement values
    """
    mixmat = Matrix{Float64}(mixmat)
    dataset = Matrix{Float64}(dataset)
    unstained = Matrix{Float64}(unstained)
    _check_inputs(mixmat, dataset, namevec)
    _check_inputs(mixmat, unstained, namevec)
    autofluorescence in eachindex(namevec) || throw(ArgumentError(
        "autofluorescence must index a mixmat column"
    ))

    calibration = mean_unmix(
        mixmat, unstained; percentile_cutoff=percentile_cutoff
    )
    baseline_matrix = copy(mixmat)
    baseline_matrix[:, autofluorescence] .= 0
    adjusted = dataset .- transpose(baseline_matrix * calibration.means)
    fit = TRU_OLS(
        mixmat, adjusted, calibration.cutoffs, namevec;
        autofluorescence=autofluorescence
    )

    result = zeros(Float64, size(dataset, 1), length(namevec))
    name_index = Dict(String(name) => j for (j, name) in enumerate(namevec))
    for i in axes(result, 1)
        for j in eachindex(fit.names[i])
            result[i, name_index[fit.names[i][j]]] = fit.coefficients[i][j]
        end
    end

    if match
        control_abundances = transpose(mixmat \ transpose(unstained))
        for (j, name) in enumerate(String.(namevec))
            events = findall(i -> haskey(fit.removed[i], name), eachindex(fit.removed))
            isempty(events) && continue
            values = [fit.removed[i][name] for i in events]
            _map_distribution!(values, Vector(control_abundances[:, j]))
            result[events, j] = values
        end
    end

    DataFrame([String(name) => result[:, j] for (j, name) in enumerate(namevec)])
end
