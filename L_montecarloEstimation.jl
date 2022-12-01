

function MC_singleInterval(p,data; lower = 0.0, upper = 2.0)
    @unpack_SDEmodel p 

    t_start = 0.0
    t_stop = t_data[end]
    t_save = t_data
    X_start = repeat(X₀,1,Ns)
    
    Nt = round(Int,(t_stop - t_start)/Δt_SDE) + 1
    Bt_samples = rand(Bt,Ns,Nt)
    X = zeros(Nt,d,Ns)
    samples = [zeros(d,Ns) for t in t_data]

    function f(μₑ)
        solveSDE!(p.ν,μₑ,p.σ,t_start,t_stop,t_save,X_start,Bt_samples,X,samples,p)
        return pseudoML(samples,data,p)
    end

    result = optimize(f,lower,upper; show_trace=true)

    return result,samples
end


function MC_multiInterval(p,data; lower = 0.0, upper = 2.0)

    @unpack_SDEmodel p 

    τ = [(t_data[i],t_data[i+1]) for i in 1:length(t_data)-1]
    τ = [(0.0,t_data[1]),τ...]

    Nt = [round(Int,(τᵢ[2]-τᵢ[1])/Δt_SDE + 1) for τᵢ in τ];
    X = [zeros(Ntᵢ,d,Ns) for Ntᵢ in Nt];
    Bt_samples = [rand(Bt,Ns,Ntᵢ) for Ntᵢ in Nt];
    samples = [zeros(d,Ns) for t in t_data];

    function f(μₑ)
        t_start = 0.0
        t_stop = t_data[1]
        X_start = repeat(X₀,1,Ns)        
        solveSDE!(p.ν,μₑ,p.σ,t_start,t_stop,t_stop,X_start,
                                    Bt_samples[1],X[1],samples[1],p)
        for i in 2:length(τ)
            τᵢ = τ[i]
            t_start = τᵢ[1]
            t_stop = τᵢ[2]
            X_start = repeat(data[i-1],1,Ns)
            solveSDE!(p.ν,μₑ,p.σ,t_start,t_stop,t_stop,X_start,
                                    Bt_samples[i],X[i],samples[i],p)
        end
        return pseudoML(samples,data,p)
    end

    result = optimize(f,lower,upper; show_trace=true)
    
    f(Optim.minimizer(result))
    return result,samples
end

