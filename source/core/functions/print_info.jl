function write_final_forces(force_vec)
    open("final_forces.txt", "w") do file
        for (i, force) in enumerate(force_vec)
            x, y, z = force
            println(file, "$i $x $y $z")
        end
    end
end

function write_forces(force_file, force_vec)
    for (i, force) in enumerate(force_vec)
        x, y, z = force
        println(force_file, "$i $x $y $z")
    end
end

function write_drift(momentum_drift)
    open("momentum_drift.txt", "w") do file
        m  = mean(momentum_drift)
        sd = std(momentum_drift)
        mx = maximum(abs.(momentum_drift))
        mn = minimum(abs.(momentum_drift))

        println(file, "Mean: $m")
        println(file, "Std: $sd")
        println(file, "Max (|drift|): $mx")
        println(file, "Min (|drift|): $mn")
    end
end