ENV["GKSwstype"] = "100"
using Plots
using LinearAlgebra
using ProgressMeter
using DelimitedFiles

mutable struct Positions
    pos::Vector{Float64}
end

traj_file = "trj.xyz"
out_file = "rdf.txt"
println("How many frames would you like to analyze?")
n_frames = parse(Int64, chomp(readline()))
println("Box length?")
l_box = parse(Float64, chomp(readline()))

open(traj_file) do file
    global N = parse(Int64, readline(file))
    println("$N particles in the system.")
end

# Initialize variables
box = [l_box, l_box, l_box]
ρ = N / (box[1] * box[2] * box[3])
n_bins = 100
bins = zeros(Float64, n_bins)
d_bins = box[1] / 2
s_bin = d_bins / n_bins

function process_frame(coords)
    for i in 2:N
        for j in 1:i-1
            ∆d = coords[i].pos .- coords[j].pos
            ∆d .-= box .* floor.((∆d .+ box / 2) ./ box)
            r = norm(∆d)
            ii = Int(floor(r / s_bin))
            if ii <= n_bins - 1
                bins[ii + 1] += 2 / (n_frames * N)
            end
        end
    end
end

# Create a progress bar with an update rate of 10%
progress_bar = Progress(n_frames, 1, "Processing frames ", 10)

open(traj_file) do file
    frame_count = 0
    while !eof(file)
        N = parse(Int64, readline(file))
        read_header = readline(file)

        coords = Vector{Positions}(undef, N)
        for i in 1:N
            line = readline(file)
            pos = parse.(Float64, split(line)[2:4])
            coords[i] = Positions(pos)
        end

        process_frame(coords)
        frame_count += 1

        # Update the progress bar
        next!(progress_bar)

        if frame_count >= n_frames
            break
        end
    end
end

for i in 1:n_bins
    v_bin = ((i + 1)^3 - i^3)*s_bin^3
    nid = (4 / 3) * π * v_bin * ρ
    bins[i] = bins[i] / nid
end

open(out_file, "w") do f
    write(f, "r g\n")
    for i in 1:n_bins
        string_out = string(i * s_bin) * " " * string(bins[i]) * "\n"
        write(f, string_out)
    end
end

# Create the RDF plot
r_values = [i * s_bin for i in 1:n_bins]
p = plot(r_values, bins,
    line=:solid,
    marker=:circle,
    markercolor=:blue,
    markeralpha=0.5,
    markersize=4,
    xlabel="Distance",
    ylabel="RDF",
    title="Radial Distribution Function",
    legend=:none,
    grid=:off,
    framestyle=:box)

# Save the RDF plot
savefig(p, "rdf.png")