using DataFrames, ElasticArrays

export get_neutral_profiles, get_treatment_data


function create_filter(lm_data::Lincs, criteria::Dict{Symbol, String})
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

    cell_line_arr = ElasticArray{String}(undef, 1, 0)
    profile_arr = ElasticArray{Float32}(undef, size(lm_data.gene_df)[1], 0) 
    
    criteria = Dict{Symbol, String}(:qc_pass => "1", :pert_type => "ctl_untrt")
    f = create_filter(lm_data, criteria)
    neutral_experiments = lm_data[f] 
    
    for cl in unique(neutral_experiments.cell_iname)
        cl_str = lm_data.inst_si[cl]
        cl_filter = (f .& lm_data[:cell_iname, cl_str]) 
        cl_indices = findall(x -> x == 1, cl_filter) 
        cl_profiles = lm_data.expr[:, cl_indices] 
        append!(cell_line_arr, fill(cl_str, length(cl_indices))) 
        append!(profile_arr, cl_profiles)
    end
    
    # Create df
    profiles = transpose(profile_arr) |> Array 
    df = DataFrame(profiles, :auto)
    # Insert cell_line as the first column
    cell_lines = vec(collect(cell_line_arr))
    insertcols!(df, 1, :cell_line => cell_lines)
    # Add header
    gene_names = lm_data.gene_df.gene_symbol
    col_names = ["cell_line"; gene_names]
    rename!(df, Symbol.(col_names))

    return df
end


function get_treatment_data(lm_data::Lincs)
    #= Returns a dataFrame of 982 columns: cell line, dose, exposure time, compound, expression profile (978 genes).
    If a compound is used several times for a given experimental setup (cell line, dose, exposure time), 
    all the associated profiles are in the returned df. =#

    cell_line_arr = ElasticArray{String}(undef, 1, 0)
    dose_arr = ElasticArray{String}(undef, 1, 0)
    time_arr = ElasticArray{String}(undef, 1, 0)
    compound_arr = ElasticArray{String}(undef, 1, 0)
    profile_arr = ElasticArray{Float32}(undef, size(lm_data.gene_df)[1], 0) 

    criteria = Dict{Symbol, String}(:qc_pass => "1", :pert_type => "trt_cp")
    f = create_filter(lm_data, criteria) 
    selected_experiments = lm_data[f] 

    for cl in unique(selected_experiments.cell_iname)
        cl_str = lm_data.inst_si[cl]
        cl_filter = (f .& lm_data[:cell_iname, cl_str]) 
        cl_indices = findall(x -> x == 1, cl_filter)
        cl_experiments = lm_data[cl_filter]
        cl_profiles = lm_data.expr[:, cl_indices] 
        append!(cell_line_arr, fill(cl_str, length(cl_indices)) ) 
        append!(dose_arr, lm_data.inst_si[cl_experiments.pert_idose])
        append!(time_arr, lm_data.inst_si[cl_experiments.pert_itime])
        append!(compound_arr, lm_data.inst_si[cl_experiments.pert_id])
        append!(profile_arr, cl_profiles)
    end

    # Get SMILES 
    compounds = vec(collect(compound_arr))
    row_indices = Int32[]
    for cp in compounds
        push!(row_indices, (findall(x -> x == cp, lm_data.compound_df.pert_id))[1])
    end
    smiles = lm_data.compound_df[row_indices, :canonical_smiles]
    
    # Create df
    profiles = transpose(profile_arr) |> Array 
    df = DataFrame(profiles, :auto)
    # Insert columns cell_line, dose, exposure_time and smiles
    cell_lines = vec(collect(cell_line_arr))
    doses = vec(collect(dose_arr))
    times = vec(collect(time_arr))
    insertcols!(df, 1, :cell_line => cell_lines, :dose => doses, :exposure_time => times, :smiles => smiles)
    # Add header
    gene_names = lm_data.gene_df.gene_symbol
    col_names = ["cell_line"; "dose"; "exposure_time"; "smiles"; gene_names]
    rename!(df, Symbol.(col_names))

    return df
end




