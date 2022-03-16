include("L_loadAll.jl")


function main()

end


function solveSDE!Example(p)
    @unpack_SDEmodel p 
    Nt = round(Int,T/Δt_SDE) + 1

    t_start = 0.0
    t_stop = T
    t_save = t_data
    X_start = repeat(X₀,1,Ns)
    Bt_samples = [rand(Bt) for i in zeros(Ns,Nt)]
    X = zeros(Nt,d,Ns)
    samples = [zeros(d,Ns) for t in t_data]

    solveSDE!(p.ν,p.μ,p.σ,t_start,t_stop,t_save,X_start,Bt_samples,X,samples,p)
    return samples 
end