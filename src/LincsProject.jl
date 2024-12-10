module LincsProject

    export Lincs

    include("structs.jl")
    using .Structs

    include("utils.jl")

end

#NEXT: fetch the important functionalities below

# working on new branch?

# using HDF5, CSV, DataFrames
# using ProgressBars

# export Data, getindex, length # , StrIndex
# export compound, gene

## Provide a 2-way indexing between string and int

# ## remove to replace by symbols
# struct StrIndex
#     str2id::Dict{String, Int32}
#     id2str::Vector{String}

#     StrIndex(vs::Vector{String}) = new(Dict(vs[i] => i for i = 1:length(vs)), vs) # must be uniqued!
#     StrIndex(ds::HDF5.Dataset) = StrIndex(ds[:])
# end

## Indexing

# Base.getindex(idx::StrIndex, s::String) = idx.str2id[s]
# Base.getindex(idx::StrIndex, i::Integer) = idx.id2str[i]
# Base.getindex(idx::StrIndex, v::AbstractVector{String}) = [idx[s] for s in v]
# Base.getindex(idx::StrIndex, v::AbstractVector{<:Integer}) = [idx[i] for i in v]
# Base.getindex(idx::StrIndex, df::AbstractDataFrame) = mapcols(col -> idx[col], df)    
# Base.length(idx::StrIndex) = length(idx.id2str)

## HDF5 IO

# Base.setindex!(f::HDF5.File, s::StrIndex, k::String) = setindex!(f, s.id2str, k)

# function Base.setindex!(f::HDF5.File, df::AbstractDataFrame, k::String)
#     g = create_group(f, k)
#     for (name, vec) in pairs(eachcol(df))
#         g[String(name)] = vec
#     end
# end

# function DataFrames.DataFrame(g::HDF5.Group)
#     convert(p) = (p.first, p.second[:]) # To pull data from the HDF5 dataset
#     return DataFrame(Dict(map(convert, pairs(g))))
# end

## LINCS data structure

# function Lincs(fn::String)
#     println("Loading from cache...")

#     f = h5open(fn)
#     gene_df = DataFrame(f["gene_df"])
#     gene_si = StrIndex(f["gene_si"])
#     compound_df = DataFrame(f["compound_df"])
#     compound_si = StrIndex(f["compound_si"])
#     inst_df = DataFrame(f["inst_df"])
#     inst_si = StrIndex(f["inst_si"])
#     expr = f["expr"][:,:]
#     return Lincs(expr, gene_df, gene_si, compound_df, compound_si, inst_df, inst_si)
# end

# Base.getindex(d::Lincs, sym::Symbol) = d.inst_df[!, sym]
# Base.getindex(d::Lincs, sym::Symbol, value::Int32) = d.inst_df[!, sym] .== value
# Base.getindex(df::DataFrame, sym::Symbol, value::Int32) = df[!, sym] .== value
# Base.getindex(d::Lincs, sym::Symbol, value::String) = d[sym, d.inst_si[value]]

# Base.getindex(d::Lincs, sym::Symbol, v::AbstractVector{<:Integer}) = reduce(.|, [d[sym, value] for value = v])
# Base.getindex(d::Lincs, sym::Symbol, v::AbstractVector{String}) = d[sym, d.inst_si[v]]
# Base.getindex(d::Lincs, v::BitVector) = d.inst_df[v,:]
# Base.getindex(d::Lincs, v::BitVector, z) = d.inst_df[v,z]

# Base.unique(d::Lincs, sym::Symbol) = unique(d.inst_df[!, sym])
# (d::Lincs)(value::Integer) = d.inst_si[value]
# (d::Lincs)(v::AbstractVector{<:Integer}) = d.inst_si[v]

# compound(d::Lincs) = d.compound_df
# compound(d::Lincs, id::String) = d.compound_df[d.compound_si[id],:]

# gene(d::Lincs) = d.gene_df
# gene(d::Lincs, sym::String) = d.gene_df[d.gene_si[sym], :]
