using Arrhenius
using Sundials
using LinearAlgebra
using OrdinaryDiffEq
using ForwardDiff
using ForwardDiff: jacobian, jacobian!
using Flux, DiffEqFlux
using Flux.Optimise:update!
using Flux.Losses:mae
using DiffEqSensitivity
using YAML
using BSON: @save, @load
using Plots, Printf, Random
using StatsBase
using Statistics
using Distributions
using ProgressBars
using Base.Threads
using BandedMatrices
using SparseArrays
using LatinHypercubeSampling
using DelimitedFiles
# using Interpolations
using Dierckx # for interpolation

# set for parallel threads number
Threads.nthreads() = 4;

# GKS settings
ENV["GKSwstype"] = "100";

# load input.yaml
runtime = YAML.load_file("./input.yaml");
expr_name = runtime["expr_name"];
is_restart = runtime["is_restart"];

# load $expr_name/config.yaml
conf = YAML.load_file("$expr_name/config.yaml");

# Mention: Julia Dict from YAML do not keep order in YAML
mechs = sort(collect(conf["mechs"]));
name_arr = [n for (n,m) in mechs];
mech_arr = [m for (n,m) in mechs];

# Load settings from config
fuel = conf["fuel"];
oxyd = conf["oxyd"];

P_arr = conf["P_arr"];
T0_arr = conf["T0_arr"];
phi_arr = conf["phi_arr"];

dT = Float64(conf["dTign"]);
dTabort = Float64(conf["dTabort"]);

# print restart info
if is_restart
    println("Continue to run $expr_name ...\n");
else
    println("Runing $expr_name ...\n");
end

# """
#     Create File Path with `is_restart` input
#     if is_restart = true, do not delete old dirs
#     if is_restart = false, delete old dirs and create new dirs
# """
function createPath(path, is_restart)
    if ispath(path)
        if !is_restart
            rm(path, recursive=true);
            sleep(0.1);
            mkdir(path);
        end
    else
        mkdir(path);
    end
end

# path settings
exp_path = string(expr_name);
fig_path = string("$exp_path/figs");
ckpt_path = string("$exp_path/checkpoint");

# create paths
createPath(fig_path, is_restart);
createPath(ckpt_path, is_restart);
createPath("$fig_path/conditions", is_restart);

# global variables, used in simulator.jl
global P = P_arr[1] * one_atm;
global gas = CreateSolution(mech_arr[1]);
global ns = gas.n_species;
global nr = gas.n_reactions;
