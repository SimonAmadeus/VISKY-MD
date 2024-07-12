using LinearAlgebra
using LaTeXStrings
using Plots
using Measures

filename = "restart.data"

mutable struct Velocity
    velocity::Vector{Float64}
end

open(filename) do file
    global number_particles = 0
    counter = 0
    counter_check = "check"
    for line in eachline(file)
        
        a = split(line) # create array from line.

        # read header
        if length(a) > 1
            if a[2] == "atoms"
                number_particles = parse(Int64, a[1])
                println("$number_particles particles in the system.")
                # define arrays for particles, velocities, bonds, angles
                global velocity_vec = Vector{Velocity}(undef, number_particles)
            end
        end

        if counter >= 2 && counter_check == "Velocities"
            counter += 1
            index = parse(Int64, a[1])
            vel_i = [parse(Float64, a[2]), parse(Float64, a[3]), parse(Float64, a[4])]
            velocity_vec[index] = Velocity(vel_i)
            if counter == (number_particles + 2)
                counter = 0
                counter_check = "check"
            end
        end
        if counter == 1 && counter_check == "Velocities"
            counter = 2
        end
        if isempty(a) == false && a[1] == "Velocities"
            counter_check = "Velocities"
            counter = 1
        end 
    end
end

v_vec = zeros(Float64, length(velocity_vec))

for i in 1:length(velocity_vec)
    v_vec[i] = sqrt(velocity_vec[i].velocity[1]^2 + velocity_vec[i].velocity[2]^2 + velocity_vec[i].velocity[3]^2)
end

# parameters for the Maxwell-Boltzmann distribution
#m = 1.0  # mass of a particle
#k = 1.380649e-23 # Boltzmann constant in J/K
#T = 300.0 # temperature in K

# function to calculate the Maxwell-Boltzmann distribution
#P(v; m=m, T=T) = 4pi * (m / (2pi * k * T))^(3/2) * v^2 * exp(-m * v^2 / (2 * k * T))

P(v) = (4/pi)^(1/2) * v^2 * exp(-v^2 / 2)

function P(v; T=0.7)
    prefactor = (1 / (2π * T))^(3/2)
    return 4π * prefactor * v^2 * exp(-v^2 / (2 * T))
end


# calculate the range of velocities to plot the Maxwell-Boltzmann distribution
v_min, v_max = minimum(v_vec), maximum(v_vec)
v_range = range(v_min, v_max, length=100)

# Calculate the Maxwell-Boltzmann distribution for the range of velocities
P_v = P.(v_range)

#println(P_v, v_range)

function plot_velocity_distribution()

    p = histogram(v_vec, bins=30, norm=:pdf, alpha=0.5,
        label="sim.. data",
        size=(1300, 800), xlabel=L"v", ylabel=L"P(v)",
        title="Velocity Distribution",
        framestyle=:box,
        #legend=:outerright,
        legendfontsize=18, titlefontsize=20, guidefontsize=18, tickfontsize=16,
        xticks=:mirror, yticks=:mirror,
        margin=10mm)
    
    #histogram(v_vec, bins=30)
    plot!(v_range, P_v, line=(:red, 2), label="Maxwell-Boltzmann Distribution")

    savefig(p, "velocity_distribution.png")
end

plot_velocity_distribution()