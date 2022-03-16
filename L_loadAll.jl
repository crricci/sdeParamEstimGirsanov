using BenchmarkTools
using Distributions
using ForwardDiff
using KernelDensity
using Optim
using Parameters
using PyPlot
using LinearAlgebra

include("L_parameters.jl")
include("L_generateData.jl")
include("L_SDE.jl")
include("L_kde.jl")
include("L_plotLib.jl")