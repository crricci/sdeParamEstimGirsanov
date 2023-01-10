
function solveSDE!(ν,μ,σ,t_start,t_stop,t_save,X_start,Bt_samples,X,samples,p)
    """ overwrite samples -> t_save x Matrix(d x Ns)
        ν,μ,σ = parameters of the SDE,
        (t_start,t_stop) = time range
        t_save = times at which to overwrite samples
        X_start = starting condition -> d x Ns, 
        Bt_samples = increments Brownian Motion -> Ns x Nt 
        X = used as preallocation -> Nt x d x Ns
        p = SDEmodel
    """ 
    Δt, h = p.Δt_SDE, p.h_SDE
    Nt = round(Int,(t_stop - t_start)/Δt) + 1

    X[1,:,:] .= X_start   

    for i = 1:p.Ns
        for t=1:Nt-1
            X[t+1,:,i] .= X[t,:,i] + Δt * drift(X[t,:,i],ν,μ,p.A) + 
                                    h * diffusion(X[t,:,i],σ,p.d) * Bt_samples[i,t]
        end
    end

    i_save = findall(ti -> ti in t_save, t_start:Δt:t_stop)
    if length(i_save) > 1
        samples .= [X[i,:,:] for i in i_save]
    else 
        samples .= X[i_save[1],:,:]
    end
end


function solveODE!(ν,μ,t_start,t_stop,t_save,x_start,x,x_save,p)
    """ solve ODE without noise, x_save (preallocated) -> t_save x d 
    x = (preallocated) solution of the ODE -> Nt x d 
    ν,μ = parameters of the ODE,
    (t_start,t_stop) = time range
    t_save = times at which to overwrite samples
    x_start = starting condition -> d 
    p = SDEmodel
    """ 
    Δt = p.Δt_SDE
    Nt = round(Int,(t_stop - t_start)/Δt) + 1

    x[1,:] .= x_start

    for t=1:Nt-1
        x[t+1,:] .= x[t,:] + Δt * drift(x[t,:],ν,μ,p.A) 
    end

    i_save = findall(ti -> ti in t_save, t_start:Δt:t_stop)
    if length(i_save) > 1
        x_save .= [x[i,:] for i in i_save]
    else 
        x_save .= x[i_save[1],:]
    end
    
end

