module PowerParser

import JLD2
import Downloads

include("parser.jl")

global SILENCED = false;

function parse_pglib(
    ::Type{T},
    dataset :: String,
    datadir :: String;
    out_type=Data{T}
) :: Union{Data{T}, NamedTuple} where T <: Real
    mkpath(datadir)
    data_path = joinpath(datadir, dataset)
    fname = if isfile(data_path)
        data_path
    else
        Downloads.download(
            "https://raw.githubusercontent.com/power-grid-lib/pglib-opf/master/$dataset",
            data_path,
        )
    end
    parse_file(T, fname; datadir, out_type)
end

function parse_file(
    ::Type{T},
    fname :: String;
    datadir=nothing,
    out_type=Data{T}
) :: Union{Data{T}, NamedTuple} where T <: Real
    if out_type != Data && out_type != NamedTuple
        @error "Argument out_type must have value NamedTuple | PowerParser.Data"
    end
    _, f = splitdir(fname)
    name, _ = splitext(f)
    cached_path = nothing
    named_tuple = out_type != Data{T}
    if isnothing(datadir)
        cached_path = joinpath(datadir, name * "_ " * string(T) * ".jld2")
    elseif !isdir(datadir)
        mkdir(datadir)
    end

    if !isnothing(cached_path) && isfile(cached_path)
        SILENCED || @info "Loading cached JLD2 file at " * cached_path
        data = JLD2.load(cached_path)["data"]
        return named_tuple ? struct_to_nt(data) : data
    else
        SILENCED || @info "Loading MATPOWER file at " * fname
        data = process_ac_power_data(T, fname)
        if !isnothing(cached_path)
            SILENCED || @info "Caching parsed matpower file to " * cached_path
            JLD2.save(cached_path, "data", data)
        end
        return named_tuple ? struct_to_nt(data) : data
    end
end

function silence()
    @info "PowerParser.jl has been silenced for the rest of the session."
    global SILENCED = true
end

export parse_file, parse_pglib, Data, BusData, BranchData, StorageData, GenData

end # module PowerParser
