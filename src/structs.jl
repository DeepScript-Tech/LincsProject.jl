module Structs
export Lincs

using DataFrames

struct Lincs
    expr::Matrix{Float32}
    gene::DataFrame
    # gene_si::StrIndex ## Used as an index to _df
    compound::DataFrame
    # compound_si::StrIndex ## Used as an index to _df
    inst::DataFrame ## Converted to identifiers only, use inst_si to convert back
    # inst_si::StrIndex
end

include("parser.jl") ## constructors

end