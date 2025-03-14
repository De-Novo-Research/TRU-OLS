The file provides a reference implementation for the TRU_OLS algorithm. The official publication is at <insert link when published.>
 
How To Use
 
#Read in a mixing matrix, an unstained dataset, and a nulticolor dataset from CSV files as DataFrames
mm=CSV.read(mixmat_path,DataFrame)
us=CSV.read(unstained_path,DataFrame)
ms=CSV.read(multi_path,DataFrame)

#Get names from mixing matrix.
namel = names(mixmat)

#Cast all DataFrames to Matrices.
mm = Matrix(mm)
us = Matrix(us)
ms = Matrix(ms)

#Invoke function create_complete_dataframe with those as arguments, using negative matching and a default 99 percentile cutoff
result = create_complete_dataframe(mm, namel, ms, us, true) 
