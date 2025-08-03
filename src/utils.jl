using DataFrames # , ElasticArrays

export get_untreated_profiles, get_treated_profiles, split_by, expr


Base.getindex(d::Lincs, sym::Symbol) = d.inst[!, sym]
Base.getindex(d::Lincs, sym::Symbol, value::Symbol) = d.inst[!, sym] .== value
Base.getindex(d::Lincs, sym::Symbol, values::Vector{Symbol}) = [v in values for v in d.inst[!, sym]]
Base.getindex(d::Lincs, sym::Symbol, value::String) = d.inst[!, sym] .== value
Base.getindex(d::Lincs, sym::Symbol, values::Vector{String}) = [v in values for v in d.inst[!, sym]]
Base.getindex(d::Lincs, v::BitVector) = d.inst[v, :]

Base.getindex(df::DataFrame, sym::Symbol, value::Symbol) = df[!, sym] .== value
Base.getindex(df::DataFrame, sym::Symbol, values::Vector{Symbol}) = [v in values for v in df[!, sym]]
Base.getindex(df::DataFrame, sym::Symbol, value::String) = df[!, sym] .== value
Base.getindex(df::DataFrame, sym::Symbol, values::Vector{String}) = [v in values for v in df[!, sym]]
Base.getindex(df::DataFrame, v::BitVector) = df[v, :]

split_by(d::Lincs, sym::Symbol) = [String(v) => d[sym, v] for v in unique(d[sym])] |> Dict
expr(d::Lincs, i::Vector{Int64}) = d.expr[:,i]

function create_filter(lm_data::Lincs, criteria::Dict{Symbol, Vector{Symbol}})
    # Returns a bit vector: 1 for experiments that satisfy all the criteria, else 0.
    filters = []
    for (k, v) in criteria
        f = lm_data[k, v]
        push!(filters, f)
    end
    return reduce(.&, filters)
end




function get_untreated_profiles(lm_data::Lincs)
    #= Returns a df of 979 columns: cell line, untreated expression profile (978 genes).
    If several untrt profiles are available for a given cell line, they are all stored in the returned df. =#
    criteria = Dict{Symbol, Vector{Symbol}}(:qc_pass => [Symbol("1")], :pert_type => [:ctl_untrt, :ctl_vehicle])
    f = create_filter(lm_data, criteria)
    untrt_experiments = lm_data[f] 
    untrt_profiles = lm_data.expr[:, f]
    df = DataFrame(transpose(untrt_profiles), :auto)
    insertcols!(df, 1, :cell_line => untrt_experiments.cell_iname)
    # Add header
    col_names = ["cell_line"; lm_data.gene.gene_symbol]
    rename!(df, col_names)
    return df
end


function get_treated_profiles(lm_data::Lincs)
    #= Returns a dataFrame of 982 columns: cell line, dose, exposure time, compound, expression profile (978 genes).
    If a compound is used several times for a given experimental setup (cell line, dose, exposure time), 
    all the associated profiles are in the returned df. =#

    criteria = Dict{Symbol, Vector{Symbol}}(:qc_pass => [Symbol("1")], :pert_type => [:trt_cp])
    f = create_filter(lm_data, criteria) 
    trt_experiments = lm_data[f] 
    trt_profiles = lm_data.expr[:, f]
    df = DataFrame(transpose(trt_profiles), :auto)

    # Create a dict pert_id to SMILES from lm_data.compound:
    pert_id_to_smiles = Dict(zip(lm_data.compound.pert_id, lm_data.compound.canonical_smiles))

    smiles = [pert_id_to_smiles[pert_id] for pert_id in trt_experiments.pert_id]

    insertcols!(df, 1, :cell_line => trt_experiments.cell_iname, 
                        :dose => trt_experiments.pert_idose, 
                        :exposure_time => trt_experiments.pert_itime, 
                        :smiles => smiles)      

    # Add header
    col_names = ["cell_line"; "dose"; "exposure_time"; "smiles"; lm_data.gene.gene_symbol]
    rename!(df, col_names)

    # Remove the rows where the SMILES, the dose or the exposure time is missing:
    filtered_df = filter(row -> row.smiles != "" && row.dose != Symbol("") && row.exposure_time != Symbol(""), eachrow(df)) |> DataFrame

    return filtered_df
end




