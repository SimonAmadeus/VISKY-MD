export Harmonic_Angles

struct Harmonic_Angles <: Angle_Interactions
    κ_a
    θ_0::Float64 # Equilibrium angle.
end

function apply_angle_interactions!(args::System, particle_vec::Vector{Particle}, force_vec, energy::Float64, stress::Vector{Float64}, angle_vec::Vector{Angle}, ::No_Angles)
    return energy, stress
end

function apply_angle_interactions!(args::System, particle_vec::Vector{Particle}, force_vec, energy::Float64, stress::Vector{Float64}, angle_vec::Vector{Angle}, harmonic_angle::Harmonic_Angles)
    θ_sum = 0
    for i in 1:length(angle_vec)

        pos_1 = particle_vec[angle_vec[i].p_1].pos .+ particle_vec[angle_vec[i].p_1].image .* args.box
        pos_2 = particle_vec[angle_vec[i].p_2].pos .+ particle_vec[angle_vec[i].p_2].image .* args.box
        pos_3 = particle_vec[angle_vec[i].p_3].pos .+ particle_vec[angle_vec[i].p_3].image .* args.box

        vec_21 = pos_1 - pos_2 
        vec_23 = pos_3 - pos_2
        vec_32 = pos_2 - pos_3
        # Scalar product.
        scal_p = vec_21[1] * vec_23[1] + vec_21[2] * vec_23[2] + vec_21[3] * vec_23[3]
        θ = scal_p / (norm(vec_21) * norm(vec_23))

        θ = acos(clamp(θ, -1.0, 1.0)) * 180 / π
        θ_sum += θ

        pvec_1 = cross(vec_21, cross(vec_21, vec_23)) / norm(cross(vec_21, cross(vec_21, vec_23)))
        pvec_3 = cross(vec_32, cross(vec_21, vec_23)) / norm(cross(vec_32, cross(vec_21, vec_23))) 

        force_1 = - 2 * (harmonic_angle.κ_a * (θ - harmonic_angle.θ_0) / norm(vec_21)) * pvec_1
        force_3 = - 2 * (harmonic_angle.κ_a * (θ - harmonic_angle.θ_0) / norm(vec_23)) * pvec_3
        force_2 = - force_1 - force_3 

        force_vec[angle_vec[i].p_1].force += force_1
        force_vec[angle_vec[i].p_2].force += force_2
        force_vec[angle_vec[i].p_3].force += force_3

        energy += harmonic_angle.κ_a * (θ - harmonic_angle.θ_0) ^ 2

    end
    return energy, stress
end
