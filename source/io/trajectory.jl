export XYZTrajectoryDump, LammpsTrajectoryDump

function apply_dump(fileO, step::Int64, args::System, particle_vec::Vector{Particle}, force_vec, velocity_vec, ::No_TrajectoryDump) end

struct XYZTrajectoryDump <: TrajectoryDump end

struct LammpsTrajectoryDump <: TrajectoryDump end

function apply_dump(fileO, step::Int64, args::System, particle_vec::Vector{Particle}, force_vec, velocity_vec, ::XYZTrajectoryDump) 
    
    x = 0 
    y = 0
    z = 0
    
    string_out=""
    string_out=string(length(particle_vec))
    string_out*="\n"
    write(fileO, string_out)
    string_out=""
    string_out*=string(step)
    string_out*="\n"
    write(fileO, string_out)
    for i = 1:length(particle_vec)
        string_out=""
        type_ = particle_vec[i].spec
        x = particle_vec[i].pos[1] + (particle_vec[i].image[1] * args.box[1])
        y = particle_vec[i].pos[2] + (particle_vec[i].image[2] * args.box[2])
        z = particle_vec[i].pos[3] + (particle_vec[i].image[3] * args.box[3])
        vx = velocity_vec[i][1]
        vy = velocity_vec[i][2]
        vz = velocity_vec[i][3]
        #string_out = string(type_)*" "*string(x)*" "*string(y)*" "*string(z)
        string_out = string(type_)*" "*string(x)*" "*string(y)*" "*string(z)*" "*string(vx)*" "*string(vy)*" "*string(vz)
               string_out*="\n"
        write(fileO, string_out)
    end
end

function apply_dump(fileO, step::Int64, args::System, particle_vec::Vector{Particle}, force_vec, velocity_vec, ::LammpsTrajectoryDump) 
    
    x = 0 
    y = 0
    z = 0
    ix = 0
    iy = 0
    iz = 0
    
    # Header.

    string_out = ""
    string_out *= "ITEM: TIMESTEP"
    string_out *= "\n"
    write(fileO, string_out)
    string_out = ""
    string_out *= string(step)
    string_out *= "\n"
    write(fileO, string_out)
    string_out = ""
    string_out *= "ITEM: NUMBER OF ATOMS"
    string_out *= "\n"
    write(fileO, string_out)
    string_out = ""
    string_out *= string(length(particle_vec))
    string_out *= "\n"
    write(fileO, string_out)
    string_out = ""
    string_out *= "ITEM: BOX BOUNDS pp pp pp"
    string_out *= "\n"
    write(fileO, string_out)
    string_out = ""
    string_out *= string(- 0.5 * args.box[1])*" "*string(0.5 * args.box[1])
    string_out *= "\n"
    write(fileO, string_out)
    string_out = ""
    string_out *= string(- 0.5 * args.box[2])*" "*string(0.5 * args.box[2])
    string_out *= "\n"
    write(fileO, string_out)
    string_out = ""
    string_out *= string(- 0.5 * args.box[3])*" "*string(0.5 * args.box[3])
    string_out *= "\n"
    write(fileO, string_out)
    string_out = ""
    string_out *= "ITEM: ATOMS id type x y z ix iy iz"
    string_out *= "\n"
    write(fileO, string_out)

    # Trajectory.
    for i = 1:length(particle_vec)
        string_out=""
        type_ = particle_vec[i].spec
        x = particle_vec[i].pos[1]
        y = particle_vec[i].pos[2]
        z = particle_vec[i].pos[3]
        ix = particle_vec[i].image[1]
        iy = particle_vec[i].image[2]
        iz = particle_vec[i].image[3]
        string_out *= string(i)*" "*string(type_)*" "*string(x)*" "*string(y)*" "*string(z)*" "*string(ix)*" "*string(iy)*" "*string(iz)
        string_out*="\n"
        write(fileO, string_out)
    end
end