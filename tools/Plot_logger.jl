using DelimitedFiles
using Plots

# Read the data from a file (assuming the data is in a file called "md_output.txt")
data = readdlm("logger_MD.txt", skipstart=1)

# Extract the columns of data
t = data[:, 1]
T = data[:, 2]
V = data[:, 3]
K = data[:, 4]
p = data[:, 5]

# Create a plot with multiple subplots
p1 = plot(t, T, xlabel="Time", ylabel="Temperature", legend=false)
p2 = plot(t, V, xlabel="Time", ylabel="Potential Energy", legend=false)
p3 = plot(t, K, xlabel="Time", ylabel="Kinetic Energy", legend=false)
p4 = plot(t, p, yscale=:log10, xlabel="Time", ylabel="Total momentum", legend=false)
plot(p1, p2, p3, p4, layout=(4,1), size=(800,800))

# Save the figure as an image
savefig("logger.png")
