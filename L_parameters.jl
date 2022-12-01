
@with_kw struct SDEmodel{T}

    # MODEL
    ν::T = 1.0
    μ::T = 1.0
    σ::T = 1.0

    d::Int = 2
    T_end::T = 1.0
    X₀::Vector{T} = ones(d)

    # NUMERICAL
    Δt_SDE::T = 1e-3
    h_SDE::T = sqrt(Δt_SDE)
    Ns::Int = 1000
    
    # DATA
    t_data::Vector{T} = collect(0.1:0.1:T_end)

    # AUXILIARY 
    Bt::MvNormal{T, Distributions.PDMats.PDiagMat{T, Vector{T}}, Vector{T}} = 
                                         MultivariateNormal(zeros(T,d),I(d))
    A::Matrix{T} = diagm([T(k)^2 for k in 1:d])
end


function B₀(x)
    """vector field elementwise"""
    return [B_comp(xi) for xi in x]
end

function B_comp(x)
    """scalar nonlinearity"""
    return sin(x)
end

function drift_Linear(x,ν,A)
    """linear vector value drift""" 
    return - ν * A * x 
end

function drift_nonLinear(x,μ)
    """nonlinear vector value drift""" 
    return  μ * B₀(x)
end

function drift(x,ν,μ,A)
    """full vector value drift""" 
    return drift_Linear(x,ν,A) + drift_nonLinear(x,μ)
end

function diffusion(x,σ,d)
    """full vector value diffusion"""
    return σ * I(d)
end

