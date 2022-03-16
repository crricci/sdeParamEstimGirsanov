function pseudoML(samples,data,p)
    """ samples -> t_save x Matrix(d x Ns)
        data -> t_save x Vector(d)
        return the psudoLogLikelyhood of  the samples computed at data
    """
    NData = length(data)
    return sum([-log(KDEstimator(samples[t],data[t],p)) for t in 1:NData])
end

function KDEstimator(samples::Matrix{Float64},point,p)
    """samples -> d x Ns, point -> d"""
    SD = vec(std(samples,dims=2))   # standard deviations for each component 
    invHsq = diagm([1 / S_ruleOfThumb(sd,p.Ns,p.d) for sd in SD])   # 1/√Hᵢᵢ
    
    return 1/p.Ns * sum([ K(point.-samples[:,i],invHsq,p) for i in 1:p.Ns])
end

function S_ruleOfThumb(sd,Ns,d)
    return (4/(d+2))^(1/(d+4)) * Ns^(-1/(d+4)) * sd
end

function K(x,p)
    """gaussian kernel"""
    return pdf(p.Bt,x)
end

function K(x,invHsq,p)
    """kernel with rescaling
    invHsq is the square root of the inverse of diagonal matrix with 
    sd on diagonal"""
    return det(invHsq) * K(invHsq * x,p)
end