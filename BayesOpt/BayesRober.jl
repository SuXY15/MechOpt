using Random
using StatsPlots
using LinearAlgebra
using DifferentialEquations
using Turing
using PyPlot
Random.seed!(42)

# the ODE function
function robertson(dy, y, p, t)
    k = [0.04, 3e7, 1e4] .* p
    dy[1] = -k[1]*y[1]+k[3]*y[2]*y[3]
    dy[2] =  k[1]*y[1]-k[2]*y[2]^2-k[3]*y[2]*y[3]
    dy[3] =  k[2]*y[2]^2
end

# settings
y0 = [1.0, 0.0, 0.0]
p = [1., 1., 1.]
datasize = 50
tsteps = 10 .^ range(log10(1e-4), log10(1e8), length=datasize);
tspan = (0.0, tsteps[end]+1e-3);
solver = KenCarp4()
prob = ODEProblem(robertson, y0, tspan, p)

# generate data
y_true = solve(prob, solver, p=p, saveat=tsteps)
scale = vec(maximum(y_true, dims=2));
scale = [1, 5e-5, 1]
y_data = Array(y_true) + 0.05 * randn(size(Array(y_true))) .* scale

# Plot simulation and noisy observations
function doplot(tsteps, y_data, y_pred)
    fig, axs = PyPlot.subplots(3,1, figsize=(4,5))
    for i in 1:3
        axs[i].plot(tsteps, y_data[i,:], "k.", label="Data")
        axs[i].plot(tsteps, y_pred[i,:], "-",  label="")
        axs[i].set_xscale("log")
        axs[i].set_ylabel("\$y_{$i}\$")
    end
    axs[2].legend(loc="best", frameon=false)
    axs[3].set_xlabel("Time [s]")
    fig.subplots_adjust(left=0.15, right=0.98, hspace=0.35, top=0.98)
    return fig
end
display(doplot(tsteps, y_data, y_true))

## First attempt: reproduce the webpage

# ================================
# Prepare for optimization
@model function fitlv(data, prob)
    # Prior distributions.
    σ ~ InverseGamma(2, 3)
    k1 ~ truncated(LogNormal( log(1.0), 0.5), 0.2, 4.)
    k2 ~ truncated(LogNormal( log(1.0), 0.5), 0.2, 4.)
    k3 ~ truncated(LogNormal( log(1.0), 0.5), 0.2, 4.)

    # Simulate ODE
    p = [k1, k2, k3]
    predicted = solve(prob, solver; p=p, saveat=tsteps)

    # Observations
    for i in 1:length(predicted)
        data[:, i] ~ MvNormal(predicted[i], diagm(σ^2 * scale.*scale))
    end
    return nothing
end

model = fitlv(y_data, prob)
# Sample 3 independent chains with forward-mode automatic differentiation (the default).
chain = sample(model, NUTS(0.65), MCMCSerial(), 1000, 1; progress=false)
Plots.plot(chain)

idxs = rand(1:1000,100)
posterior_samples = chain[[:k1, :k2, :k3]][idxs,:,1]

# Plot results
fig, axs = PyPlot.subplots(3,1, figsize=(4,5))
for i in 1:3
    axs[i].plot(tsteps, y_data[i,:], "k.", label="Data")
    axs[i].set_xscale("log")
    axs[i].set_ylabel("\$y_{$i}\$")
end
axs[2].legend(loc="best", frameon=false)
axs[3].set_xlabel("Time [s]")
fig.subplots_adjust(left=0.15, right=0.98, hspace=0.35, top=0.98)
for p in eachrow(Array(posterior_samples))
    y_pred = solve(prob, solver; p=p, saveat=tsteps)
    axs[1].plot(tsteps, y_pred[1,:], alpha=0.1, color="#BBBBBB")
    axs[2].plot(tsteps, y_pred[2,:], alpha=0.1, color="#BBBBBB")
    axs[3].plot(tsteps, y_pred[3,:], alpha=0.1, color="#BBBBBB")
end
display(fig)

Plots.plot(chain)
