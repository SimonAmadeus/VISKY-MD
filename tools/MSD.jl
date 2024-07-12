using Plots
using LinearAlgebra
using ProgressMeter

#println("Enter input file (XYZ format).")
#traj_file = chomp(readline())
#println("Enter out-file (e.g. msd.txt).")
#out_file = chomp(readline())
#println("How many frames would you like to analyze?")
#n_frames = parse(Int64, chomp(readline()))

traj_file = "trj.xyz"
out_file = "msd.txt"

n_frames = 500
dt = 1

fh = 2 # File header.
N = 0 # Number of particles.

open(traj_file) do file
    global N = parse(Int64, split(readline(file))[1])
    global n_chains = Int(N / len_chn)
    println("$n_chains chains in the system.")
end


x = zeros(n_frames, N)
y = zeros(n_frames, N)
z = zeros(n_frames, N)

open(traj_file) do file
    t = 0
    for line in eachline(file)
        t += 1
        m = mod(t, fh + N)
        frame = (t ÷ (fh + N)) + 1
        if (m != 1 && m != 2 )
            a = split(line)
            ind = m - fh
            if(m == 0)
                ind = N
                frame -= 1
            end
            x[frame, ind] = parse(Float64, a[2])
            y[frame, ind] = parse(Float64, a[3])
            z[frame, ind] = parse(Float64, a[4])                  
        end        
    end
end

l_msd = n_frames ÷ 2
msd = zeros(l_msd)
for i in 1:(n_frames - l_msd)
    l = 0
    for j in i + 1:(i + l_msd)
        l += 1
        for k in 1:N
            msd[l] += (x[j, k] - x[i, k])^2 + (y[j, k] - y[i, k])^2 + (z[j, k] - z[i, k])^2
        end
    end         
end

time = collect(1:l_msd)
time *= dt

msd /= (N * (n_frames - l_msd))

plot(time, msd)
xaxis!("timesteps")
yaxis!("MSD")
title!("MSD")
savefig("msd.png")

# Compute the diffusion coefficient
t_fit_start = 1000  # Starting point for the linear fit
t_fit_end = 2000    # Ending point for the linear fit

#slope, _ = linreg(time[t_fit_start:t_fit_end], msd[t_fit_start:t_fit_end])
slope, _ = linreg(time, msd)

D = slope / 6

println("Diffusion coefficient: $D")

out_file = open(out_file, "w")
string_out = "t msd\n"
write(out_file, string_out)
for i in 1:length(msd)
    global string_out = string(time[i]) * " " * string(msd[i]) * "\n"
    write(out_file, string_out)
end
close(out_file)

out_file_2 = open("diffusion_coefficient.txt", "w")
write(out_file_2, "Diffusion coefficient: $D")
close(out_file_2)