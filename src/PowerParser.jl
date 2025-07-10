module PowerParser

import JLD2
import Downloads

include("parser.jl")

function parse_matpower_pglib(
    ::Type{T},
    dataset :: String,
    datadir :: String;
    named_tuple=false
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
    parse_matpower_file(T, fname; datadir, named_tuple)
end

function parse_matpower_file(
    ::Type{T},
    fname :: String;
    datadir=nothing,
    named_tuple=false
) :: Union{Data{T}, NamedTuple} where T <: Real
    _, f = splitdir(fname)
    name, _ = splitext(f)
    cached_path = nothing
    isnothing(datadir) || cached_path = joinpath(datadir, name * ".jld2")

    if !isnothing(cached_path) && isfile(cached_path)
        @info "Loading cached JLD2 file at " * cached_path
        data = JLD2.load(cached_path)["data"]
        return named_tuple ? struct_to_nt(data) : data
    else
        @info "Loading MATPOWER file"
        data = process_ac_power_data(T, fname)
        if !isnothing(cached_path)
            @info "Caching parsed matpower file to "
            JLD2.save(cached_path, "data", data)
        end
        return named_tuple ? struct_to_nt(data) : data
    end
end


end # module PowerParser
