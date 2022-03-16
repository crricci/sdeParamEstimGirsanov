

@with_kw struct SDEmodel

    # MODEL
    ν = 1.0
    μ = 1.0
    σ = 0.01

    d = 2
    T = 5.0
    X₀ = ones(d)

    # NUMERICAL
    Δt_SDE = 1e-3
    h_SDE = sqrt(Δt_SDE)
    Ns = 100
    
    # DATA
    t_data = 0.1:0.1:T

    # AUXILIARY 
    Bt = MultivariateNormal(zeros(d),I(d))
    A = diagm([k^2 for k in 1:d])
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
    return - p.ν * p.A * x 
end

function drift_nonLinear(x,μ)
    """nonlinear vector value drift""" 
    return  p.μ * B₀(x)
end

function drift(x,ν,μ,A)
    """full vector value drift""" 
    return drift_Linear(x,ν,A) + drift_nonLinear(x,μ)
end

function diffusion(x,σ,d)
    """full vector value diffusion"""
    return p.σ * I(p.d)
end

