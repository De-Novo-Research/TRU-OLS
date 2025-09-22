using CSV, DataFrames, LinearAlgebra, StatsBase, FileIO, FCSFiles

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
    
    # --- CORRECTED: Initialized as Float64[] to ensure correct type ---
    means = Float64[]
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
    
    # --- CORRECTED: Initialized as Float64[] to ensure correct type ---
    per = Float64[]
    
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
            # --- IMPROVED: Replaced unstable matrix inversion with Julia's backslash operator ---
            unmix = mixmat2 \ v
            
            exclude_list = Int[] # Use Int[] for type stability
            
            for j in 1:a
                if unmix[j] < threshvec2[j]
                    # Store the name and unmixed value of columns being removed
                    removed_cols_dict[i][namevec2[j]] = unmix[j]
                    append!(exclude_list, j)
                end
            end
            
            if isempty(exclude_list) # More idiomatic way to check for an empty array
                temp = unmix
                tempn = namevec2
                mbelow = 1
            else
                # --- NOTE: This logic was changed slightly to avoid errors when all columns are excluded ---
                # This ensures the loop eventually terminates.
                mixmat2 = mixmat2[:, Not(exclude_list)]
                threshvec2 = threshvec2[Not(exclude_list)]
                namevec2 = namevec2[Not(exclude_list)]
                a = length(threshvec2)
                if a == 0 # If all columns removed, stop the loop for this cell
                    mbelow = 1
                end
            end
        end 
        
        append!(unmixed, [temp])
        append!(namel, [tempn])
    end
    
    return unmixed, namel, removed_cols_dict
end

function mapDistribution!(irrelevantData::Vector{Float64}, distributionToMatch::Vector{Float64})
    """
    This function takes in two lists of values.  The first is a list of all irrelevant abundances for a single endmember over a dataset.
    The second is a list of unmixed control abundances.  The output is a vector with the irrelevant abundances replaced with
    their percentile match from the control

    Parameters:
    - irrelevantData: a list of the irrelevant unmixed values for a single endmember over a whole dataset
    - distributionToMatch: A list of abundances for the same endmember to match the irrelevant data to

    Return:
    - irrelevantData: Percentile matched single stain data to replace irrelevant data
    """
    if isempty(irrelevantData) || isempty(distributionToMatch)
        return irrelevantData # Return early if nothing to map
    end

    # Sort the distribution to match percentiles
    sortedDistribution = sort(distributionToMatch)

    # Sort the incoming data and store original indices
    sorted_indices = sortperm(irrelevantData)
    sortedIrrelevantData = irrelevantData[sorted_indices]

    lenM1 = length(sortedDistribution) - 1

    for cntr = 1:length(irrelevantData)
        # Calculate percentile in the irrelevant data
        sortedPerc = (cntr - 1) / (length(sortedIrrelevantData) - 1)

        # Find index of corresponding percentile in distribution to match
        ssIdx = round(Int, sortedPerc * lenM1) + 1

        # Get the original index (possibly shuffled)
        origIdx = sorted_indices[cntr]

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

    m,n = size(dataset)
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
        # Unmix negative and store values by column 
        x, y = size(unstained_dataset)
        unstained_df = DataFrame([name => Float64[] for name in namevec])
        
        # --- CORRECTED: Loop now unmixes each row individually ---
        for i in 1:x
            unmixed_row = mixmat \ unstained_dataset[i, :]
            push!(unstained_df, unmixed_row)
        end
        
        # Process removed columns to organize by column name
        column_removed_values = organize_by_column_name(removed_cols_dict)
        
        # Transform the removed values using map_distribution
        transformed_values = Dict{String, Vector{Float64}}()
        for (col_name, values) in column_removed_values
            # Ensure the column exists in the unstained data before trying to access it
            if hasproperty(unstained_df, col_name)
                transformed_values[col_name] = mapDistribution!(values, unstained_df[!, col_name])
            end
        end
        
        # Create a dictionary to easily access the transformed values
        # This is more efficient than searching through the dictionary every time.
        col_to_transformed_idx = Dict{String, Int}()
        for (col_name, values) in column_removed_values
             col_to_transformed_idx[col_name] = 0
        end

        # Fill in the dataframe with transformed values for removed columns
        for i in 1:m
            removed_in_sample = removed_cols_dict[i]
            for (col_name, removed_value) in removed_in_sample
                if haskey(transformed_values, col_name) && result[col_name][i] == 0.0
                    # Get the next available transformed value for this column
                    current_idx = col_to_transformed_idx[col_name] + 1
                    if current_idx <= length(transformed_values[col_name])
                        result[col_name][i] = transformed_values[col_name][current_idx]
                        col_to_transformed_idx[col_name] = current_idx
                    end
                end
            end
        end
    end
    
    df = DataFrame(result)
    return df
end





