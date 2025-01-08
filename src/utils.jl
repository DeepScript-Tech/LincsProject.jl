using DataFrames, ElasticArrays

export get_neutral_profiles, get_treatment_data


Base.getindex(d::Lincs, sym::Symbol, value::Symbol) = d.inst[!, sym] .== value
Base.getindex(d::Lincs, sym::Symbol, values::Vector{Symbol}) = filter(row -> row[sym] in values, d.inst)
Base.getindex(d::Lincs, v::BitVector) = d.inst[v, :]

Base.getindex(df::DataFrame, sym::Symbol, value::Symbol) = df[!, sym] .== value
Base.getindex(df::DataFrame, sym::Symbol, values::Vector{Symbol}) = filter(row -> row[sym] in values, df)
Base.getindex(df::DataFrame, v::BitVector) = df[v, :]


function create_filter(lm_data::Lincs, criteria::Dict{Symbol, Symbol})
    # Returns a bit vector: 1 for experiments that satisfy all the criteria, else 0.
    filters = []
    for (k, v) in criteria
        f = lm_data[k, v]
        push!(filters, f)
    end
    return reduce(.&, filters)
end


function get_neutral_profiles(lm_data::Lincs)
    #= Returns a df of 979 columns: cell line, neutral expression profile (978 genes).
    If several neutral profiles are available for a given cell line, they are all stored in the returned df. =#

    cell_line_arr = ElasticArray{Symbol}(undef, 1, 0)
    profile_arr = ElasticArray{Float32}(undef, size(lm_data.gene)[1], 0) 
    
    criteria = Dict{Symbol, Symbol}(:qc_pass => Symbol("1"), :pert_type => :ctl_untrt)
    f = create_filter(lm_data, criteria)
    neutral_experiments = lm_data[f] 
    
    for cl in unique(neutral_experiments.cell_iname)
        cl_filter = (f .& lm_data[:cell_iname, cl]) 
        cl_indices = findall(x -> x == 1, cl_filter) 
        cl_profiles = lm_data.expr[:, cl_indices] 
        append!(cell_line_arr, fill(cl, length(cl_indices))) 
        append!(profile_arr, cl_profiles)
    end
    
    # Create df
    profiles = transpose(profile_arr) |> Array 
    df = DataFrame(profiles, :auto)
    # Insert cell_line as the first column
    cell_lines = String.(vec(collect(cell_line_arr)))
    insertcols!(df, 1, :cell_line => cell_lines)
    # Add header
    col_names = ["cell_line"; lm_data.gene.gene_symbol]
    rename!(df, col_names)

    return df
end


function get_treatment_data(lm_data::Lincs)
    #= Returns a dataFrame of 982 columns: cell line, dose, exposure time, compound, expression profile (978 genes).
    If a compound is used several times for a given experimental setup (cell line, dose, exposure time), 
    all the associated profiles are in the returned df. =#

    cell_line_arr = ElasticArray{Symbol}(undef, 1, 0)
    dose_arr = ElasticArray{Symbol}(undef, 1, 0)
    time_arr = ElasticArray{Symbol}(undef, 1, 0)
    compound_arr = ElasticArray{Symbol}(undef, 1, 0)
    profile_arr = ElasticArray{Float32}(undef, size(lm_data.gene)[1], 0) 

    criteria = Dict{Symbol, Symbol}(:qc_pass => Symbol("1"), :pert_type => :trt_cp)
    f = create_filter(lm_data, criteria) 
    selected_experiments = lm_data[f] 

    for cl in unique(selected_experiments.cell_iname)
        cl_filter = (f .& lm_data[:cell_iname, cl]) 
        cl_indices = findall(x -> x == 1, cl_filter)
        cl_experiments = lm_data[cl_filter]
        cl_profiles = lm_data.expr[:, cl_indices] 
        append!(cell_line_arr, fill(cl, length(cl_indices))) 
        append!(dose_arr, cl_experiments.pert_idose)
        append!(time_arr, cl_experiments.pert_itime)
        append!(compound_arr, cl_experiments.pert_id)
        append!(profile_arr, cl_profiles)
    end

    # Get SMILES 
    compounds = vec(collect(compound_arr))
    row_indices = Int32[]
    for cp in compounds
        push!(row_indices, (findall(x -> x == cp, lm_data.compound.pert_id))[1])
    end
    smiles = lm_data.compound[row_indices, :canonical_smiles]
    
    # Create df
    profiles = transpose(profile_arr) |> Array 
    df = DataFrame(profiles, :auto)
    # Insert columns cell_line, dose, exposure_time and smiles
    cell_lines = String.(vec(collect(cell_line_arr)))
    doses = String.(vec(collect(dose_arr)))
    times = String.(vec(collect(time_arr)))
    insertcols!(df, 1, :cell_line => cell_lines, :dose => doses, :exposure_time => times, :smiles => smiles)
    # Add header
    col_names = ["cell_line"; "dose"; "exposure_time"; "smiles"; lm_data.gene.gene_symbol]
    rename!(df, col_names)

    return df
end




