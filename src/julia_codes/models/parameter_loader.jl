using CSV
using DataFrames

"""
    load_parameteers(file_path::String, μ::Float64)

Given parameter values with CSV format, load the dataframe including
parameter values. The parameter named μ means the administration timing
and it is added at the end of load values.

# Arguments
- `file_path::String`: CSV path containing parameter values
- `μ::Float64`: Additional parameter set by user (administration timing)

# Returns
- `Vector{Any}`: Parameter vector (because u1, u2 are vector)
"""
function load_parameter(file_path::String, μ::Float64)
    # load parameters from the csv file
    all_params = CSV.read(file_path, DataFrame)
    
    # to Float 64 (for easy usage)
    all_params.Value = map(p -> tryparse(Float64, p) !== nothing ? parse(Float64, p) : eval(Meta.parse(p)), all_params.Value)

    # set parameter μ as the last value
    all_params[end, :Value] = μ

    # return parameter values
    return all_params.Value
end