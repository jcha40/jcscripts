using CSV
using DataFrames
using SparseArrays
using DelimitedFiles
using ArgParse

function pileup(scidx_fn::String, beds_fns::String, chrom_sizes::String, out_fn::String; window::Int=201)
    sizes = Dict{String, Int64}()
    open(chrom_sizes, "r") do f
        while !eof(f)
            line = split(readline(f))
            sizes[line[1]] = parse(Int64, line[2])
        end
    end

    scidx = CSV.read(scidx_fn, DataFrame, delim="\t", comment="#")
    contig = scidx[1, "chrom"]
    start = 1
    chrom_dict = Dict{String, SparseMatrixCSC{Int64, Int64}}()
    n = nrow(scidx)
    for i in 1:n
        if scidx[i, "chrom"] != contig
            idx = scidx[start:i - 1, "index"]
            arr = sparse(vcat(idx, idx), vcat(fill(1, i - start), fill(2, i - start)), vcat(scidx[start:i - 1, "forward"], scidx[start:i - 1, "reverse"]), sizes[contig], 2)
            chrom_dict[contig] = arr
            contig = scidx[i, "chrom"]
            start = i
        end
    end
    idx = scidx[start:n, "index"]
    chrom_dict[contig] = sparse(vcat(idx, idx), vcat(fill(1, n - start + 1), fill(2, n - start + 1)), vcat(scidx[start:n, "forward"], scidx[start:n, "reverse"]), sizes[contig], 2)

    f = open(beds_fns, "r")
    beds = readlines(f)
    close(f)

    pileup_mat = zeros(Int64, (length(beds) * 2, window))
    for i in 1:length(beds) # Not sure why there are race conditions when I try to multithread this
        bed = CSV.read(beds[i], DataFrame, delim="\t", header=["Chromosome", "Start", "End", "Name", "Score", "Strand"])
        strand = Vector{Int64}(undef, nrow(bed))
        @simd for j in 1:nrow(bed)
            strand[j] = convert(Int64, bed[j, "Strand"] == "-")
        end

        for j in 1:nrow(bed)
            if bed[j, "Start"] > 0 && bed[j, "End"] <= sizes[bed[j, "Chromosome"]]
                strand_idx = strand[j]
                chrom_sparse = chrom_dict[bed[j, "Chromosome"]]
                start = bed[j, "Start"] - 1
                @simd for k in bed[j, "Start"]:bed[j, "End"]
                    pileup_mat[2 * i + strand_idx - 1, k - start] += chrom_sparse[k, strand_idx + 1]
                end
            end
        end
    end

    writedlm(out_fn, pileup_mat)
end

function main()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--scidx-fn", "-s"
            help = "scidx tab file for sample"
            required = true
        "--beds-fns", "-b"
            help = "file containing list of bed filenames to pileup on"
            required = true
        "--chrom-sizes", "-c"
            help = "chromosome sizes file"
            required = true
        "--out", "-o"
            help = "output file name"
            default = "bulk_pileup.txt"
        "--window", "-w"
            help = "window size in bp"
            arg_type = Int
            default = 201
    end

    args = parse_args(s)

    pileup(args["scidx-fn"], args["beds-fns"], args["chrom-sizes"], args["out"], window=args["window"])
end

main()
