include("L_loadAll.jl")


function main(p)
    data = generateData(p)
    result, samples = MC_singleInterval(p,data)
    result, samples = MC_multiInterval(p,data)
end


function solveSDE!Example(p)
    @unpack_SDEmodel p 

    t_start = 0.0
    t_stop = t_data[end]
    t_save = t_data
    Nt = round(Int,(t_stop - t_start)/Δt_SDE) + 1
    X_start = repeat(X₀,1,Ns);
    Bt_samples = [rand(Bt) for i in zeros(Ns,Nt)];
    X = zeros(Nt,d,Ns);
    samples = [zeros(d,Ns) for t in t_data];

    solveSDE!(p.ν,p.μ,p.σ,t_start,t_stop,t_save,X_start,Bt_samples,X,samples,p)
    return samples 
end
