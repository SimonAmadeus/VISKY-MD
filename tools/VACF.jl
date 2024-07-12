# This code calculates the velocity autocorrelation function of a given trajectory.

ENV["GKSwstype"] = "100"
using Plots
using LinearAlgebra
using ProgressMeter
using Statistics
#println("Enter input file (XYZ format).")
#traj_file = chomp(readline())
#println("Enter name of text out-file (e.g. vacf.txt).")
#out_file = chomp(readline())
#println("Enter name of plot output file (eg:plot_vacf.png).")
#file0 = chomp(readline())
#println("How many frames would you like to analyze?")
#n_frames = parse(Int64, chomp(readline()))

#traj_file = "/run/media/simon/Doc-Stick/PhD/Projects/MPMDP/Testruns/hPF_test/trj.xyz"
traj_file = "trj.xyz"
#out_file = "/run/media/simon/Doc-Stick/PhD/Projects/MPMDP/Testruns/hPF_test/vel_autocorr.txt"
out_file = "VACF.txt"
n_frames = 50
dt = 1

fh = 2 # File header.
n_particles = 0 # Number of particles.
total_frames = 1

open(traj_file) do file
    global N = parse(Int64, readline(file))
    println("$N particles in the system.")
end

vx = zeros(n_frames, N)
vy = zeros(n_frames, N)
vz = zeros(n_frames, N)

open(traj_file) do file
    for (t, line) in enumerate(eachline(file))
        m = mod(t, fh + N)
        frame = (t ÷ (fh + N)) + 1
        if (m != 1 && m != 2 ) && frame <= n_frames
            a = split(line)
            ind = m - fh
            if m == 0
                ind = N
            end
            vx[frame, ind] = parse(Float64, a[5])
            vy[frame, ind] = parse(Float64, a[6])
            vz[frame, ind] = parse(Float64, a[7])
        end
    end
end

# Substract center of mass velocity.
for i in 1:n_frames
    vx[i, :] .-= mean(vx[i, :])
    vy[i, :] .-= mean(vy[i, :])
    vz[i, :] .-= mean(vz[i, :])
end

# Calculate the velocity autocorrelation function (VACF).
vacf = zeros(n_frames)
@showprogress for tau in 1:n_frames
    sum_vdotv = 0.0
    for i in 1:n_frames - tau
        sum_vdotv += dot(vx[i, :], vx[i + tau, :])
        sum_vdotv += dot(vy[i, :], vy[i + tau, :])
        sum_vdotv += dot(vz[i, :], vz[i + tau, :])
    end
    vacf[tau] = sum_vdotv / (3 * (N - tau) * n_frames)
end

out_file = open(out_file, "w")
for tau in 1:n_frames
    global string_out = string(tau * dt) * " " * string(vacf[tau]) * "\n"
    write(out_file, string_out)
end

# plotting:
default(size=(800, 600), linewidth=2, markersize=4, legendfontsize=12, guidefontsize=14, tickfontsize=12)
gr()

time = collect(1:n_frames)
time *= dt

# Create the VACF plot
p = plot(time, vacf,
    line=:solid,
    #marker=:circle,
    #markercolor=:blue,
    #markeralpha=0.5,
    #markersize=4,
    xlabel="Timesteps",
    ylabel="VACF",
    #title="Velocity Autocorrelation Function",
    legend=:none,
    #grid=:off,
    framestyle=:box)

# Save the VACF plot
savefig(p, "VACF.png")