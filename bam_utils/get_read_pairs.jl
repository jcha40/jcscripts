using XAM

function get_read_pairs(bam::String)
    record_dict = Dict{String, BAM.Record}()
    return Channel() do c
        open(BAM.Reader, bam) do reader
            for record in reader
                if !BAM.isprimary(record) || !BAM.ismapped(record)
                    continue
                end
                readname = BAM.tempname(record)
                if !haskey(record_dict, readname)
                    record_dict[readname] = record
                    continue
                end
                mate = pop!(record_dict, readname)
                if BAM.flag(record) & 0x40 == 0x40
                    push!(c, tuple(record, mate))
                else
                    push!(c, tuple(mate, record))
                end
            end
        end
    end
end
