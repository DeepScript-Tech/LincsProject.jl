module Structs
export Lincs

using DataFrames

struct Lincs
    expr::Matrix{Float32}
    gene::DataFrame
    compound::DataFrame
    inst::DataFrame ## Converted to identifiers only, use inst_si to convert back
end

include("parser.jl") ## constructors

end