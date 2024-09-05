using LincsProject, HDF5, CSV, DataFrames, ProgressBars

# z = Lincs("data/lincs_all_genes.h5")
include("../src/parallel_sort.jl")

# z = Lincs("/home/muninn/scratch/lincs_beta/", "level3_beta_all_n3026460x12328.gctx", "data/out.jld2")

# prefix, gctx = "/home/muninn/scratch/lincs_beta/", "level3_beta_all_n3026460x12328.gctx"
prefix, gctx = "/scratch/lincs_beta/", "level3_beta_all_n3026460x12328.gctx"

f = h5open(prefix * gctx)
expr = f["0/DATA/0/matrix"]
exprGene_id = Symbol.(f["0/META/ROW/id"][:])

## Compound information

println("Loading compound annotations...")
compound_df = CSV.File(prefix * "compoundinfo_beta.txt", 
                       delim="\t", types=String, missingstring=nothing, pool=false, ntasks=1) |> DataFrame
gdf = groupby(compound_df, [:pert_id, :canonical_smiles, :inchi_key])
# pert_id unique per smiles/key. some have multiple cmap_name, keep only the first

compound_df = combine(gdf, :cmap_name => (x -> first(x)) => :first_name)
compound_df[!,:pert_id] = Symbol.(compound_df.pert_id)

## Gene and sample annotations

println("Loading gene and sample annotations...")
gene_df = CSV.File(prefix * "geneinfo_beta.txt",
                   delim="\t", types=String, missingstring=nothing, pool=false) |> DataFrame
gene_df.gene_id = Symbol.(gene_df.gene_id)
gene_df.gene_type = Symbol.(gene_df.gene_type)
gene_df.src = Symbol.(gene_df.src)
gene_df.feature_space = Symbol.(gene_df.feature_space)

inst_df = CSV.File(prefix * "instinfo_beta.txt", delim = '\t', types=String, missingstring=nothing, pool=false) |> DataFrame

@Threads.threads for col in names(inst_df)
    inst_df[!, col] = Symbol.(inst_df[!, col]) # 
end

expr_id = Symbol.(f["0/META/COL/id"][:]) # Order of samples in the matrix

e2s = psortperm(expr_id)                # Permutation to sort expr_id into standard order
i2s = psortperm(inst_df.sample_id)      # Permutation to sort inst_df.sample_id into standard order
s2e = psortperm(e2s)                    # Permutation to sort standard order into the order of expr_id
i2e = i2s[s2e]                          # Permutation to sort inst_df.sample_id into the order of expr_id

inst_df = inst_df[i2e,:]

println("Subsetting landmark genes")
g = groupby(gene_df, :feature_space)
lm_id = get(g, (feature_space=:landmark,), nothing).gene_id
# gene_df = gene_df[dfGene_idx[lm_id],:]
# lm_sym = StrIndex(gene_df.gene_symbol)
# lm_row = exprGene_idx[lm_id] ## convert to the matrix rows
# lm_df = get(g, (feature_space=:landmark,), nothing)
lm_row = [findfirst(id -> id == sym, exprGene_id) for sym in lm_id]

chunk_size = 32*1000
ngene, ninst = size(expr)
nlm = length(lm_row)

final = Matrix{Float32}(undef, (nlm, ninst)) 

function load!(final, expr, lm_row, r::UnitRange{T}) where T <: Integer
    # println(r)
    # println(expr[:,r])
    slab = expr[:,r]
    final[:,r] = slab[lm_row,:]
end

function test(chunk_size, ninst, final, expr, lm_row)
    ## This actually loads the file (about 15 minutes)
    @Threads.threads for start in ProgressBar(1:chunk_size:ninst)
        # println("start=$start")
        r = start:min(start+chunk_size-1, ninst)
        load!(final, expr, lm_row, r)
        # println("end=$start")
    end
end

test(chunk_size, ninst, final, expr, lm_row)

lincs = Lincs(final, gene_df, compound_df, inst_df)
jldsave("/home/iorek/lemieuxs/tmp/save_test.jld2"; lincs)


# out_fn = "data/save_test.h5"
# f_out = h5open(out_fn, "w")
# f_out["gene_df"] = gene_df
# f_out["compound_df"] = compound_df
# f_out["inst_df"] = inst_df
    
# # About 11GB for the expressions (landmark genes only)
# f_out["expr"] = final
# close(f_out)

res = Lincs(expr, gene_df, compound_df, inst_df)


 

# @code_warntype test(chunk_size, ninst, final, expr, lm_row)
# @code_warntype load!(final, expr, lm_row, 1:chunk_size)

# ###

# chars = ['a':'z'; 'A':'Z']
# l = [Symbol(String(rand(chars, 1000))) for _=1:1000]
# using JLD2
# save("bing.jld2", Dict("l" => l))