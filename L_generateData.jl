

function generateData(p)
    """
    generate data from the equation, saving the results only at t_data
    t_data is a vector of times, and is assumed to contain 0.0 (initial time)
    """
    @unpack ν,μ,σ,d,T,X₀,t_data,Bt,A = p
    Δt, h = p.Δt_SDE, p.h_SDE

    Nt = round(Int,T/Δt) + 1


    X = zeros(Nt,d)
    X[1,:] = X₀

    for t = 1:Nt-1
        X[t+1,:] = X[t,:] + Δt * drift(X[t,:],ν,μ,A) + 
                        h * diffusion(X[t,:],σ,d) * rand(Bt)
    end
    
    i_save = findall(ti -> ti in t_data, 0:Δt:T)
    data = [vec(X[i,:]) for i in i_save]

    return data
end
