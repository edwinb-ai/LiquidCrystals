using FastGaussQuadrature: gausslegendre
using LinearAlgebra: dot

const npoints = 100

function compute_kernel(x, w, theta_k, theta_j)
    term = @. cos(π * x)
    term = @. cos(theta_k) * cos(theta_j) - sin(theta_k) * sin(theta_j) * term
    result = dot(w, sqrt.(max.(0.0, 1 .- term .^ 2)))

    return π * result
end

function compute_kernel_matrix(x, w, θ::Vector{Float64})
    Nθ = length(θ)
    K = zeros(Nθ, Nθ)

    for k in 1:Nθ
        for j in 1:Nθ
            K[k, j] = compute_kernel(x, w, θ[k], θ[j])
        end
    end

    return K
end

function solve_ψ(
    c::Float64, wθ::Vector{Float64}, K::Matrix{Float64}; tol=1e-3, max_iter=1000
)
    # We initialize with the isotropic distribution
    τ = ones(npoints) ./ (4.0 * π)
    A = zeros(npoints)

    for iter in 1:max_iter
        for k in 1:npoints
            A[k] = 16.0 * c * sum(wθ .* K[k, :] .* τ) / π
        end
        new_ψ = @. exp(-A)
        Z = dot(wθ, new_ψ) * 4.0 * π
        new_ψ ./= Z

        # Convergence check
        δ = maximum(abs.(new_ψ .- τ))
        if δ < tol
            println("Converged in $iter iterations with δ = $δ")
            return new_ψ
        end
        τ .= new_ψ
    end

    println("Warning: did not converge within $max_iter iterations")
    return τ
end

function compute_pressure(c, kernel, orient, theta_w)
    u = orient .* theta_w
    pressure = c + 32.0 * c^2 * dot(u, kernel * u)

    return pressure
end

function main()
    # Dimensionless concentration
    c = 9.5

    # For the azimutal angle, these do not have to be transformed
    (phi_x, phi_w) = gausslegendre(npoints)

    # For the radial angle
    (theta_x, theta_w) = gausslegendre(npoints)
    # Do the correct transformation for the nodes
    θ_x = acos.((theta_x .+ 1.0) ./ 2.0)
    θ_w = theta_w ./ 2.0

    # Precompute the kernel matrix
    kernel = compute_kernel_matrix(phi_x, phi_w, θ_x)

    # Compute the orientational distribution
    orient = solve_ψ(c, θ_w, kernel; tol=1e-10, max_iter=1000)

    # Compute the pressure
    pressure = compute_pressure(c, kernel, orient, θ_w)

    println("Pressure: $pressure")

    return nothing
end

main()
