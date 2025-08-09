function print_info(args::System)
    
    println("You are currently using VISKY-MD, written by Simon Alberti (Version July24).")

    if args.non_bonded_interactions != No_Non_Bonded()
        println("Non-bonded interactons: ", args.non_bonded_interactions)
    else
        println("No non bonded interactions are used.")
    end
    
    if args.bonded_interactions != No_Bonds()
        println("Bonded interactions: ", args.bonded_interactions)
    else
        # println("No bonds are used.")
    end
    
    if args.angle_interactions != No_Angles()
        println("Angle interactions: ", args.angle_interactions)
    else
        # println("No angle interactions are used.")
    end
    
    if typeof(args.non_bonded_interactions) == DPD
        if args.thermostat != No_Thermostat()
            println("While using DPD, no additional thermostat is needed!")
            args.thermostat = No_Thermostat()
        end
    elseif args.thermostat != No_Thermostat()
        println("Thermostat: ", args.thermostat)
    else
        println("No thermostat is used.")
    end

    if args.collisions != No_Collisions()
        println("Collisions are activated: ", args.collisions)
    end
end