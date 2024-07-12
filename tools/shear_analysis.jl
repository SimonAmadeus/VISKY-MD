using CSV
using DataFrames
using Statistics
using Plots
using LsqFit
using LaTeXStrings

# This code will calculate the viscosity from a single shear flow simulation,
# and plot the velocity profile

#gr()
pyplot()

n_slabs = 20
z_length = 30
d = z_length / n_slabs

z = zeros(n_slabs)
v_x = zeros(n_slabs)

momentum = 0
viscosity = 0
viscosity_error = 0

# read momentum from file
infile = "momentum.txt"
if isfile(infile)
    open(infile, "r") do file
        for line in eachline(file)
            if occursin("Total momentum exchange j:", line)
                _, val = split(line, ':')
                global momentum = parse(Float64, strip(val))
            end
        end
    end
else
    println(infile, " does not exist.")
end

# read velocity profile from file
infile = "velocity_profile.txt"
if isfile(infile)
    f = open(infile, "r")
    data = CSV.read(f, DataFrame)
    if isempty(data)
        println("Data is empty.")
    else
        z = data[!, "z"]
        v_x = data[!, "v_x"]
    end
    close(f)
else
    println(infile, " does not exist.")
end

# Function to plot shear profile 
function plot_shear_profile(z, v_x)
    # Initialize plot with ACS style settings
    p = plot(size=(1300, 800), 
             xlabel=L"z / \sigma", 
             ylabel=L"$v_x (m / \epsilon)^{1/2}$",
             framestyle=:box,
             margin=10Plots.mm,  # Add margin to avoid label cutoff
             left_margin=12Plots.mm,
             legendfontsize=18, 
             titlefontsize=20, 
             guidefontsize=18, 
             tickfontsize=16,
             xticks=:mirror, 
             yticks=:mirror,
             grid=:on)  # Enable grid

    plot!(p, z[1:11] .- 1.5, v_x[1:11], label="")

    # Save the plot with higher resolution
    savefig(p, "shear_profile.png")
end

# define function for fitting
function func_lin(x, p)
    return p[1] .* x .+ p[2]
end

# call the function
plot_shear_profile(z, v_x)

fit = curve_fit(func_lin, z[2:9], v_x[2:9], [0.0, 0.0])
m = fit.param[1]
dm = estimate_errors(fit)[1]
viscosity_value = abs(momentum / m)
viscosity_error = abs(momentum / m^2) * dm

println("Viscosity: ", viscosity_value)
println("Error: ", viscosity_error)
