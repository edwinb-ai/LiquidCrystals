using CairoMakie

# Define the variables that give the description of the system
const L = 5.0
const D = 1.0
const vsph = (π * D^3 / 6.0) + (π * L * D^2 / 4.0)
# The number of points for dicretization
const N = 320
# The number of iterations
const niter = 1000
const tol = 1e-9

compute_concentration(rho) = π * L^2 * D * rho / 4.0

pressure(c) = c + c^2

function g_kl(ϕ, x, y)
    term_1 = cos(x) * cos(y)
    term_2 = sin(x) * sin(y) * cos(ϕ)

    return sqrt(1.0 - (term_1 + term_2)^2)
end

gaussian_tau(θ, c) = (c / π)^2 * exp(-2.0 * (c * θ)^2 / π)

α_correction(eta) = (4.0 * eta - 3.0 * eta^2) / (1.0 - eta)^2

function compute_kernel(phi, theta)
    kernel = zeros(N, N)
    for i in axes(kernel, 1)
        for k in axes(kernel, 2)
            term_1 = 3.0 * g_kl(phi[1], theta[i], theta[k]) / 2.0
            term_2 = 0.0
            for l in 2:(N - 1)
                term_2 += g_kl(phi[l], theta[i], theta[k])
            end
            term_3 = 3.0 * g_kl(phi[end], theta[i], theta[k]) / 2.0
            kernel[i, k] = 2.0 * π * (term_1 + term_2 + term_3) / (N + 1.0)
        end
    end

    return kernel
end

function compute_orientation(concentration, kernel, θk, alfa)
    # Now we need the weights for the trapezoidal rule
    Δk = zeros(N)
    Δk[1] = 1.0 - (cos(θk[1]) + cos(θk[2])) / 2.0
    for i in 2:(N - 1)
        Δk[i] = (cos(θk[i - 1]) - cos(θk[i + 1])) / 2.0
    end
    Δk[N] = (cos(θk[N]) + cos(θk[N - 1])) / 2.0
    sum_Δk = sum(Δk)
    # Check that they are normalized
    @assert sum_Δk ≈ 1.0

    # We generate a trial solution
    τk = gaussian_tau.(θk, Ref(concentration))

    # We now begin the loop
    error = 0.0
    ψk = zeros(N)
    Ak = zeros(N)

    for _ in 1:niter
        for i in 1:N
            Ak[i] = sum(Δk .* kernel[i, :] .* τk) * alfa * 16.0 * concentration / π
        end
        Z = sum(exp.(-Ak) .* Δk) * 4.0 * π
        ψk = exp.(-Ak) ./ Z
        # Compute the error
        error = maximum(abs.(ψk .- τk))
        if error < tol
            break
        else
            τk = copy(ψk)
        end
    end

    return ψk, Δk
end

function plot_orientations()
    concentrations = [4.5, 5.0, 6.0, 7.0, 10.0]

    f = Figure()
    ax = Axis(f[1, 1])

    for c in concentrations
        (orientation_dist, _, _) = compute_orientation(c)
        lines!(ax, orientation_dist; label="c = $c")
    end

    axislegend(; position=:rt)
    save("psi.png", f)

    return nothing
end

function compute_pressure(concentration, orientation_dist, kernel, delta_k)
    pressure = 0.0

    for i in axes(kernel, 1)
        for j in axes(kernel, 2)
            pressure +=
                delta_k[i] *
                delta_k[j] *
                kernel[i, j] *
                orientation_dist[i] *
                orientation_dist[j]
        end
    end
    pressure *= 32.0 * concentration^2

    return concentration + pressure
end

function eta2concentration(eta)
    # Define the density
    rho = eta / vsph
    # We compute the dimensionless concentration
    concentration = compute_concentration(rho)

    return concentration
end

function main()
    # Define the packing fraction
    eta = 0.01:0.001:0.4
    pressure = similar(eta)

    # Now we need to define( the angles
    θk = @. π * range(1, N) / (2.0 * (N + 1.0))
    # We also need azimuthal angles
    ϕj = @. 2.0 * π * range(1, N) / (N + 1.0)

    # Compute the kernel
    kernel = compute_kernel(ϕj, θk)
    exc_volume = π * L^2 / (4.0 * D)

    for i in eachindex(eta)
        alfa = α_correction(eta[i])
        # alfa = 1.0
        c = eta2concentration(eta[i])
        (ψk, Δk) = compute_orientation(c, kernel, θk, alfa)
        pressure[i] = compute_pressure(c, ψk, kernel, Δk) * vsph / exc_volume
    end

    f = Figure()
    ax = Axis(f[1, 1]; xlabel=L"\eta", ylabel=L"P")
    scatter!(ax, collect(eta), pressure)
    save("eos.png", f)

    return nothing
end

main()
