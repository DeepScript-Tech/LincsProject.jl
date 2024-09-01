using JLD2

function Lincs(fn::String)
    println("Loading from cache...")
    d = load(fn)
    return d
end

function Lincs(prefix::String, gctx::String, out_fn::String)
    isfile(out_fn) && return Lincs(out_fn)
    println("Parsing from LINCS files...")

    f_out = h5open(out_fn, "w")

    println("Loading from original files...")
    f = h5open(prefix * gctx)
    expr = f["0/DATA/0/matrix"]
    exprGene_idx = StrIndex(f["0/META/ROW/id"][:])
    
    ## Compound information

    println("Loading compound annotations...")
    compound_df = CSV.File(prefix * "compoundinfo_beta.txt", 
                           delim="\t", types=String, missingstring=nothing, pool=false) |> DataFrame
    gdf = groupby(compound_df, [:pert_id, :canonical_smiles, :inchi_key])
    # pert_id unique per smiles/key. some have multiple cmap_name, keep only the first
    
    compound_df = combine(gdf, :cmap_name => (x -> first(x)) => :first_name)
    compound_si = StrIndex(compound_df.pert_id)

    ## Gene and sample annotations
    
    println("Loading gene and sample annotations...")
    gene_df = CSV.File(prefix * "geneinfo_beta.txt",
                       delim="\t", types=String, missingstring=nothing, pool=false) |> DataFrame
    dfGene_idx = StrIndex(gene_df[:,"gene_id"])
    
    inst_df = CSV.File(prefix * "instinfo_beta.txt", delim = '\t', types=String, missingstring=nothing, pool=false) |> DataFrame
    dfInst_idx = StrIndex(inst_df[:,"sample_id"])
    exprInst_idx = StrIndex(f["0/META/COL/id"][:])
    exprInst_o = dfInst_idx[exprInst_idx.id2str]
    inst_df = inst_df[exprInst_o,:] ## Reorder the df to fit the matrix
    
    println("Preparing sample annotation global StrIndex...")
    inst_si = StrIndex(unique(reduce(vcat, [unique(c) for c in eachcol(inst_df)])))
    for i in names(inst_df)
        inst_df[!, i] = inst_si[inst_df[!,i]]
    end

    println("Subsetting landmark genes")
    g = groupby(gene_df, "feature_space")
    lm_id = get(g, (feature_space="landmark",), nothing).gene_id
    gene_df = gene_df[dfGene_idx[lm_id],:]
    lm_sym = StrIndex(gene_df.gene_symbol)
    lm_row = exprGene_idx[lm_id] ## convert to the matrix rows
    
    z = zeros(Float32, (size(expr)[1], length(lm_row)))
    for i=1:length(lm_row)
        z[lm_row[i],i] = 1
    end

    chunk_size = 8 * 4096
    ngene, ninst = size(expr)

    nlm = length(lm_row)

    final = Matrix{Float32}(undef, (nlm, ninst)) 

    ## This actually loads the file (about 15 minutes)
    for start in ProgressBar(1:chunk_size:ninst)
        r = start:min(start+chunk_size-1, ninst)
        final[:, r] = z' * expr[:,r]
    end
    
    println("Saving parsed dataset to $(out_fn)")
    
    # About 1GB for all indices and dataframes
    f_out["gene_df"] = gene_df
    f_out["gene_si"] = lm_sym
    f_out["compound_df"] = compound_df
    f_out["compound_si"] = compound_si
    f_out["inst_df"] = inst_df
    f_out["inst_si"] = inst_si
    
    # About 11GB for the expressions (landmark genes only)
    f_out["expr"] = final
    close(f_out)
    
    return Lincs(final, gene_df, lm_sym, compound_df, compound_si, inst_df, inst_si)
end
