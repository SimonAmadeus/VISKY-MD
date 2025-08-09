export Logger, logger_init

mutable struct Logger_struct <: Logger 
    Timestep::Bool
    Temp::Bool
    Pressure::Bool
    PotentialEnergy::Bool
    KineticEnergy::Bool
    Momentum::Bool
end

function logger_init()
    Timestep = true
    Temp = true
    Pressure = false
    PotentialEnergy = false
    KineticEnergy = false
    Momentum = false
    return Logger_struct(Timestep, Temp, Pressure, PotentialEnergy, KineticEnergy, Momentum)
end

function apply_logger(fileO, first_step::Bool, step::Int64, velocity_vec, energy::Float64, curr_p::Float64, particle_vec::Vector{Particle}, nl::No_Logger) end

function apply_logger(fileO, first_step::Bool, step::Int64, velocity_vec, energy::Float64, curr_p::Float64, particle_vec::Vector{Particle}, log::Logger_struct) 
    totalkineticenergy = 0
    # First, print header.
    if first_step == true
        string_out=""
        if log.Timestep == true
            string_out*="timestep"*" | "
        end
        if log.Temp == true
            string_out*="temperature"*" | "
        end
        if log.Pressure == true
            string_out*="pressure"*" | "
        end
        if log.PotentialEnergy == true
            string_out*="pot. energy"*" | "
        end
        if log.KineticEnergy == true
            string_out*="kin. energy"*" | "
        end
        if log.Momentum == true
            string_out*="tot. momentum"
        end
        string_out*="\n"
        write(fileO, string_out)

        flush(fileO)
    else
        string_out=""
        if log.Timestep == true
            string_out*=string(step)*" "
        end
        if log.Temp == true
            totalkineticenergy = calc_total_kinetic_energy(particle_vec, velocity_vec)
            inst_temp = calc_inst_temp(particle_vec, totalkineticenergy)
            string_out*=format(inst_temp,precision=6)*" "
        end
        if log.Pressure == true
            string_out*=format(curr_p, precision=6)*" "
        end
        if log.PotentialEnergy == true
            #total_energy = sum(energy_vec)
            string_out*=format(energy,precision=6)*" "
        end
        if log.KineticEnergy == true
            string_out*=format(totalkineticenergy,precision=6)*" "
        end
        if log.Momentum == true
            #momentum = calc_totalmomentum(particle_vec, velocity_vec)
            momentum = calc_momentum_mag(particle_vec, velocity_vec)
            string_out*=format(momentum, precision=16)*" "
        end
        string_out*="\n"
        write(fileO, string_out)
        flush(fileO)
    end
end
