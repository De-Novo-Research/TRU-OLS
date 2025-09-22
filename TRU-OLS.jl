using CSV, DataFrames, LinearAlgebra, StatsBase 

function mean_unmix(mixmat::Matrix, obsmat::Matrix, percen::Float64)
    """
    This function takes in a mixing matrix, a matrix of observations (observations as rows),
    and a percentile.  It returns the mean unmixed abundances in a vector and the percen percentile of each unmixed abundance in a vector.
    Run this with a full mixing matrix and unstained control observations in order to produce the nonspecific observation to subtract from multistained
    data and the cutoffs for TRU-OLS.

    Parameters:
    - mixmat: The full mixing matrix. (detectors x dyes)
    - obsmat: The unstained data matrix. (events x detectors)
    - percen: cutoff percentile
    
    Return:
    - mean unmixed obsevation
    - cutoff values vector
    """
    unmixed = []
    m, n = size(obsmat)
    for i in 1:m
        obs = obsmat[i, :]
        recov = mixmat \ Vector(obs)
        append!(unmixed, [recov])
    end
    a =  size(unmixed[1])[1]
    means = []
    vars = []
    for i in 1:a
        temp = 0
        tempv = []
        for j in 1:m
            temp += unmixed[j][i]
            append!(tempv, unmixed[j][i])
        end
        temp = temp ./ m
        append!(means, temp)
        append!(vars, var(tempv))
    end
    per = []
    for i in 1:a
        temp = []
        for j in 1:m
            append!(temp, unmixed[j][i])
        end
        tempnf = percentile(temp, percen)
        append!(per, tempnf)
    end
    return means, per
end

function TRU_OLS(mixmat::Matrix, dataset::Matrix, threshvec::Array{Float64}, namevec::Array{String})
    """
    This function takes in a mixing matrix, a dataset matrix, a vector of cutoffs, and a vector of names.  The vector
    of names are the names of the endmembers in the mixing matrix.  Must be in same order as endmember columns.
    This function runs TRU-OLS regression.  It returns the abundances that are relevant on each cell in the dataset as a list of lists,
    corresponding names as a separate list of lists, and a dictionary of removed columns with their unmixed values.
        
    Parameters:
    - mixmat: The mixing matrix (detectors x dyes)
    - dataset: The dataset matrix (events x detectors)
    - threshvec: Threshold vector for each endmember
    - namevec: Names of endmembers

    Return:
    - unmixed: A list of lists of relevant unmixed values for each cell
    - namel: A list of lists of names of endmembers with relevant unmixed values for each cell
    - removed_cols_dict: A dictionary associating cell indices to irrelevant endmembers and their umixed values
    """
    m, n = size(dataset)
    unmixed = []
    namel = []
    removed_cols_dict = Dict{Int, Dict{String, Float64}}()  # Dictionary to store removed columns by sample index
    
    for i in 1:m
        a = length(threshvec)
        mixmat2 = copy(mixmat)
        threshvec2 = copy(threshvec)
        namevec2 = copy(namevec)
        mbelow = 0
        temp = Float64[]  # Initialize temp
        tempn = String[]  # Initialize tempn
        removed_cols_dict[i] = Dict{String, Float64}()  # Initialize dictionary for this sample
        v = Vector(dataset[i, :])
        while mbelow == 0
            unmix = (mixmat2' * mixmat2)^(-1) * mixmat2' * v
            exclude_list = []
            
            for j in 1:a
                if unmix[j] < threshvec2[j]
                    # Store the name and unmixed value of columns being removed
                    removed_cols_dict[i][namevec2[j]] = unmix[j]
                    append!(exclude_list, j)
                end
            end
            
            if exclude_list == []
                temp = unmix
                tempn = namevec2
                mbelow = 1
            end
            
            mixmat2 = mixmat2[:, Not(exclude_list)]
            threshvec2 = threshvec2[Not(exclude_list)]
            namevec2 = namevec2[Not(exclude_list)]
            a = length(threshvec2)
        end 
        
        append!(unmixed, [temp])
        append!(namel, [tempn])
    end
    
    return unmixed, namel, removed_cols_dict
end

function mapDistribution!(irrelevantData::Vector{Float64}, distributionToMatch::Vector{Float64})
    """
    This function takes in two lists of values.  THe first is a list of all irrelevant abundances for a single endmember over a dataset.
    The second is a list of unmixed control abundances.  The output is a vector with the irrelevant abundances replaced with
    their percentile match from the control

    Parameters:
    - irrelevantData: a list of the irrelevant unmixed values for a single endmember over a whole dataset
    - distributionToMatch: A list of abundances for the same endmember to match the irrelevant data to

    Return:
    - irrelevantData: Percentile matched single stain data to replace irrelevant data
    """
    # Sort the distribution to match percentiles
    sortedDistribution = sort(distributionToMatch)

    # Sort the incoming data and store original indices
    sorted_indices = sortperm(irrelevantData)
    sortedIrrelevantData = irrelevantData[sorted_indices]

    lenM1 = length(sortedDistribution) - 1

    for cntr = 0:length(irrelevantData)-1
        # Calculate percentile in the irrelevant data
        sortedPerc = cntr / length(sortedIrrelevantData)

        # Find index of corresponding percentile in distribution to match
        ssIdx = round(Int, sortedPerc * lenM1) + 1  # +1 for 1-based indexing in Julia

        # Get the original index (possibly shuffled)
        origIdx = sorted_indices[cntr + 1]  # +1 for 1-based indexing

        # Get the new value from the distribution to match
        newValue = sortedDistribution[ssIdx]

        # Set the value at the original index
        irrelevantData[origIdx] = newValue
    end

    return irrelevantData
end

function organize_by_column_name(removed_cols_dict::Dict{Int, Dict{String, Float64}})
    """
    Takes the dictionary of removed columns from TRU_OLS and reorganizes it to return
    a dictionary where keys are column names and values are lists of all unmixed values
    below threshold for each column across all samples.

    Parameters:
    removed_cols_dict: Dictionary output from TRU_OLS

    Return:
    Dictionary associating names to unmixed irrelevant values, without index
    """
    # Initialize result dictionary
    column_values = Dict{String, Vector{Float64}}()
    
    # Iterate through all samples in the dictionary
    for (sample_idx, removed_columns) in removed_cols_dict
        # For each sample, go through the removed columns
        for (col_name, value) in removed_columns
            # If this column name isn't in our result dict yet, initialize it
            if !haskey(column_values, col_name)
                column_values[col_name] = Float64[]
            end
            
            # Add the value to the list for this column name
            push!(column_values[col_name], value)
        end
    end
    
    return column_values
end

function create_complete_dataframe(mixmat::Matrix, namevec::Array, dataset::Matrix, unstained_dataset::Matrix, match::Bool; percen=0.99)
    """
    This function creates a complete dataframe of unmixed values where columns are the original endmembers.
    For endmembers that survived TRU_OLS, the original unmixed values are used.
    For endmembers removed during TRU_OLS, the values are transformed using map_distribution.
    
    Parameters:
    - mixmat: The mixing matrix (detectors x dyes)
    - namevec: List of endmember names
    - dataset: The dataset matrix (events x detectors)
    - unstained_dataset: unstained data (events x detectors)
    - match: boolean to either match unstained control or not
    - percen: cutoff percentile for unstained data.  defaults to 0.99 and unused if match==false
    
    Returns:
    - A dataframe with all original endmembers, containing either the original unmixed values or transformed values
    """

    neg_abunds, cutoff = mean_unmix(mixmat, unstained_dataset, percen)
    zero_baseline_mat = copy(mixmat)
    zero_baseline_mat[:, end] .= 0.0
    baseline = zero_baseline_mat * neg_abunds

    m,n = size(full_tube_data)
    new_tube = copy(dataset)
    for i in 1:m
        new_tube[i, :] = Vector(dataset[i, :]) .- baseline
    end

    # Run TRU_OLS to get unmixed values and removed columns
    unmixed, namel, removed_cols_dict = TRU_OLS(mixmat, new_tube, cutoff, namevec)
    # Create the dataframe structure
    m = size(dataset, 1)  # Number of samples
    n = length(namevec)   # Number of endmembers
    
    # Initialize the result dataframe with appropriate column names
    result = Dict{String, Vector{Float64}}()

    for col_name in namevec
        result[col_name] = zeros(Float64, m)
    end

    # Fill in the dataframe with unmixed values that survived TRU_OLS
    for i in 1:m
        for j in 1:length(namel[i])
            col_name = namel[i][j]
            result[col_name][i] = unmixed[i][j]
        end
    end

    if match == true
        #unmix negative and store values by column 
        x, y = size(unstained_dataset)
        unstained_df = DataFrame([name => Float64[] for name in namevec])
        for i in 1:x
            push!(unstained_df, mixmat \ unstained_dataset)
        end
        # Process removed columns to organize by column name
        column_removed_values = organize_by_column_name(removed_cols_dict)
        
        # Transform the removed values using map_distribution
        transformed_values = Dict{String, Vector{Float64}}()
        for (col_name, values) in column_removed_values
            transformed_values[col_name] = map_distribution(values, unstained_df[!, col_name])
        end
        
        # Fill in the dataframe with transformed values for removed columns
        for i in 1:m
            removed_in_sample = removed_cols_dict[i]
            for col_name in namevec
                # If this column was removed in this sample and we have a transformed value
                if haskey(removed_in_sample, col_name) && result[col_name][i] == 0.0
                    # Get the index of this column in the transformed values array
                    if haskey(transformed_values, col_name)
                        # Find the index of this sample's removed value in the original array
                        idx = 0
                        for (j, (sample_idx, sample_removed)) in enumerate(removed_cols_dict)
                            if sample_idx == i && haskey(sample_removed, col_name)
                                idx = findfirst(x -> x == removed_in_sample[col_name], column_removed_values[col_name])
                                break
                            end
                        end
                        
                        if idx != nothing && idx > 0 && idx <= length(transformed_values[col_name])
                            result[col_name][i] = transformed_values[col_name][idx]
                        end
                    end
                end
            end
        end
    end
    
    df = DataFrame(result)
    return df
end


