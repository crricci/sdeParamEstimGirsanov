

function GirsanovTransformation(μ,σ,t_start,t_stop,X_start,Bt_samples,X,GirsValues,p)
    """ overwrite samples -> t_save x Matrix(d x Ns)
    μ,σ = parameters of the transformation,
    (t_start,t_stop) = time range
    t_save = times at which to overwrite samples
    X_start = starting condition -> d x Ns, 
    Bt_samples = increments Brownian Motion -> Ns x Nt 
    X = used as preallocation -> Nt x d x Ns
    p = SDEmodel
""" 


    return ρ
end

function Girsanov_singleInterval(p,data; lower = 0.0, upper = 2.0)

end

function Girsanov_multiInterval(p,data; lower = 0.0, upper = 2.0)

end