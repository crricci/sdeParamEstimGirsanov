
function plotSampleData(samples,data,p)
    @assert p.d == 1
    @assert length(data[1]) == 1
    @assert size(samples[1],1) == 1

    scatter(p.t_data,data,color="red")
    grid(true)

    for (i,t) in enumerate(p.t_data)
        sample = samples[i]
        scatter(t*ones(p.Ns),sample,0.5*ones(p.Ns),color="blue")
    end
end
