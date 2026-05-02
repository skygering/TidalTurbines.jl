export load_all_snapshots, load_snapshot, dimensionalize, export_all_snapshots, snapshot_differences

struct Snapshot
    h::Vector{Float64}
    v1::Vector{Float64}
    v2::Vector{Float64}
    b::Vector{Float64}
    t::Float64
    file::String
end

function load_all_snapshots(folder)
    files = sort(readdir(folder))
    snaps = Snapshot[]

    for f in files
        if startswith(f, "solution_") & endswith(f, ".h5")
            println(f)
            push!(snaps, load_snapshot(joinpath(folder, f)))
        end
    end
    return snaps
end

function load_snapshot(file)
    h5open(file, "r") do f

        t = read(attributes(f)["time"])

        vars = Dict{String,Vector{Float64}}()
        for i in 1:4
            ds = f["variables_$i"]
            name = read(attributes(ds)["name"])
            vars[name] = read(ds)
        end

        return Snapshot(
            vars["H"],
            vars["v1"],
            vars["v2"],
            vars["b"],
            t,
            file
        )
    end
end

function dimensionalize(snapshot; g = 9.81, D = 30)
    h  = snapshot.h .* D
    v1  = (snapshot.v1 ./ snapshot.h) .* sqrt(g * D)
    v2  = (snapshot.v2 ./ snapshot.h) .* sqrt(g * D)
    b = snapshot.b .* D
    t = snapshot.t * sqrt(D / g)

    return Snapshot(h, v1, v2, b, t, snapshot.file)
end

function snapshot_differences(a, b; atol = 1e-6)
    out = Snapshot[]

    for sa in a
        # find closest time in b
        idx = argmin(abs.(getfield.(b, :t) .- sa.t))
        sb = b[idx]

        # only accept if times are close enough
        if abs(sa.t - sb.t) ≤ atol
            diff = Snapshot(
                sa.h .- sb.h,
                sa.v1 .- sb.v1,
                sa.v2 .- sb.v2,
                sa.b,
                sa.t,
                sa.file
            )
            push!(out, diff)
        end
    end

    return out
end

function write_snapshot_like_template(template_file, out_file, snap::Snapshot)
    # copy original file structure
    println(template_file)
    println(out_file)
    cp(template_file, out_file; force=true)

    # overwrite values
    h5open(out_file, "r+") do f
        attrs = attributes(f)
        # overwrite time
        write(attrs["time"], snap.t)
        # overwrite variables
        for i in 1:4
            ds = f["variables_$i"]
            name = read(attributes(ds)["name"])

            if name == "H"
                write(ds, snap.h)
            elseif name == "v1"
                write(ds, snap.v1)
            elseif name == "v2"
                write(ds, snap.v2)
            elseif name == "b"
                write(ds, snap.b)
            else
                error("Unknown variable name: $name")
            end
        end
    end
end

function export_all_snapshots(snaps, outdir)
    mkpath(outdir)
    # copy mesh file to output output directory
    input_dir = dirname(snaps[1].file)
    files = filter(f -> occursin("mesh", f) && endswith(f, ".h5"), readdir(input_dir))
    for f in files
        cp(joinpath(input_dir, f), joinpath(outdir, f); force=true)
    end
    # copy over all state files (with snaps data)
    for snap in snaps
        # keep original filename
        fname = basename(snap.file)
        out_file = joinpath(outdir, fname)

        write_snapshot_like_template(snap.file, out_file, snap)
    end
end
