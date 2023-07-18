using LinearAlgebra
using YAML
using BSON: @save, @load
using Plots, Printf, Random
using StatsBase
using Statistics
using Distributions
using Base.Threads
using LatinHypercubeSampling
using DelimitedFiles
Random.seed!(0x7777777);

## settings
N = 10; # para dims
M = 100; # number of data samples
S = 1000; # number of para samples
n_epoch = 100;
σ_data = 0.1;
σ_para = 0.2;
UF_true = 1 .+ 1.0 * rand(N);
UF_init = 1.5 .+ 2 * rand(N);

A = rand(N, N);
function η(x, θ)
    return sin.(x') * exp.(θ);
    # return x' * θ .+ θ' * θ;
    # return x' * A * θ;
end

## generate training data
X = rand(M, N);
function LogMvNormal(n; μ=0, UF=2)
    N = length(μ);
    σ_UF = log.(sqrt.(UF)) ./ 3;
    r = hcat([rand(truncated(Normal(0, σ_UF[i]), -3*σ_UF[i], 3*σ_UF[i]), n) for i in 1:N]...)
    return (μ .* exp.(r'))';
end
θ_true = 1 .+ rand(N);
η_true = η(X', θ_true);
η_data = η_true .+ randn(M)*σ_data;

θ_init = θ_true .+ randn(N)*σ_para;
η_init = η(X', θ_init);

θs_true = LogMvNormal(S; μ=θ_true, UF=UF_true);
θs_init = LogMvNormal(S; μ=θ_init, UF=UF_init);

function π_θη(θ)
    η_pred = η(X', θ);
    ϵ = abs.(η_data .- η_pred) ./ σ_data / 10;
    return exp(-0.5 * mean(ϵ.^2));
end
@show π_θη(θ_true)
@show π_θη(θ_init)
@show π_θη(θ_pred)

θ_pred = deepcopy(θ_init);
UF_pred = deepcopy(UF_init);
for k in 1:n_epoch
    θs_pred = LogMvNormal(S; μ=θ_pred, UF=UF_pred.+0.1);
    p = [π_θη(θs_pred[i,:]) for i in 1:S];
    p_samples = p ./ maximum(p);
    θ_samples = θs_pred[rand(S) .< p_samples,:];
    θ_pred = vcat(mean(θ_samples, dims=1)...);
    UF_pred = vcat(maximum(θ_samples, dims=1) ./ minimum(θ_samples, dims=1)...);
end
η_pred = η(X', θ_pred);
θs_pred = LogMvNormal(S; μ=θ_pred, UF=UF_pred);

include("visual_bayes.jl")
plot_dots(η_data, η_init, η_pred; ylabel="η")
plot_dots(θ_true, θ_init, θ_pred; ylabel="θ")
plot_pdf(θs_true, θs_init, θs_pred; i=1, j=2)
plot_pdf(θs_true, θs_init, θs_pred; i=3, j=4)
plot_bars(UF_true, UF_init, UF_pred; ylabel="UF")
