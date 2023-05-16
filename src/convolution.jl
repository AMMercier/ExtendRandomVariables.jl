struct convolutionRV <: RandomVariable
    distr
    id::Int64
end

abstract type Convolution{T} <: UnivariateDistribution{T} end

struct ConvolutionDiscrete <: Convolution{Discrete}
    pdf
    l::Float64
    u::Float64
end

struct ConvolutionContinuous <: Convolution{Continuous}
    pdf
    l::Float64
    u::Float64
end

function continuous_mult_convolution(dx::ContinuousDistribution, dy::ContinuousDistribution)
    #Order has to be quite high for some convolutions
    function pdf(z)
        if z == 0
            val = 0
        else
            val = quadgk(x -> (Distributions.pdf(dx,x) * Distributions.pdf(dy,z/x))*(1/abs(x)), -Inf, 1e-10, order=200)[1]
            val += quadgk(x -> (Distributions.pdf(dx,x) * Distributions.pdf(dy,z/x))*(1/abs(x)), 1e-10, Inf, order=200)[1]
        end
    end
    return pdf
end

function discrete_mult_convolution(dx::DiscreteDistribution, dy::DiscreteDistribution)
    l = *(quantile(dx, 1e-10), quantile(dy, 1e-10))
    u = *(quantile(dx, 1-1e-10), quantile(dy, 1-1e-10))
    sup = collect(l:u)
    pmf(z) = sum(Distributions.pdf.(dx, sup).*Distributions.pdf.(dy, z ./ sup))
    return pmf, l, u
end

function continuous_add_convolution(dx::ContinuousDistribution, dy::ContinuousDistribution)
    #Order has to be quite high for some convolutions
    pdf(z) = quadgk(x -> (Distributions.pdf(dx,z-x) * Distributions.pdf(dy,x)), -Inf, Inf, order=200)[1] 
    return pdf
end

function discrete_add_convolution(dx::DiscreteDistribution, dy::DiscreteDistribution)
    l = sum([quantile(dx, 1e-10), quantile(dy, 1e-10)])
    u = sum([quantile(dx, 1-1e-10), quantile(dy, 1-1e-10)])
    sup = collect(l:u)
    pmf(z) = sum(Distributions.pdf.(dx, z .- sup).*Distributions.pdf.(dy, sup))
    return pmf, l, u
end

function continuous_sub_convolution(dx::ContinuousDistribution, dy::ContinuousDistribution)
    #Order has to be quite high for some convolutions
    pdf(z) = quadgk(x -> (Distributions.pdf(dx,z+x) * Distributions.pdf(dy,x)), -Inf, Inf, order=200)[1] 
    return pdf
end

function discrete_sub_convolution(dx::DiscreteDistribution, dy::DiscreteDistribution)
    l = -([quantile(dx, 1e-10), quantile(dy, 1e-10)])
    u = -([quantile(dx, 1-1e-10), quantile(dy, 1-1e-10)])
    sup = collect(l:u)
    pmf(z) = sum(Distributions.pdf.(dx, z .+ sup).*Distributions.pdf.(dy, sup))
    return pmf, l, u
end

(Base.:*)(X::RV, Y::RV) = begin
    if X.id == Y.id
        error("Dependent variables not yet supported.")
    elseif X.distr isa ContinuousDistribution && Y.distr isa DiscreteDistribution
        error("Mixed distributions (Continuous & Discrete) not yet supproted.")
    elseif X.distr isa DiscreteDistribution && Y.distr isa ContinuousDistribution
        error("Mixed distributions (Continuous & Discrete) not yet supported.")
    end
    min = minimum(X.distr) * minimum(Y.distr)
    max = maximum(X.distr) * maximum(Y.distr)
    if X.distr isa ContinuousDistribution && Y.distr isa ContinuousDistribution
        distr = continuous_mult_convolution(X.distr, Y.distr)
        conv = ConvolutionContinuous(distr, min, max)
    elseif X.distr isa DiscreteDistribution && Y.distr isa DiscreteDistribution
        distr = discrete_mult_convolution(X.distr, Y.distr)
        conv = ConvolutionDiscrete(distr[1], distr[2], distr[3])
    end
    return RV(conv)
end

(Base.:+)(X::RV, Y::RV) = begin
    if X.id == Y.id
        error("Dependent variables not yet supported.")
    elseif X.distr isa ContinuousDistribution && Y.distr isa DiscreteDistribution
        error("Mixed distributions (Continuous & Discrete) not yet supproted.")
    elseif X.distr isa DiscreteDistribution && Y.distr isa ContinuousDistribution
        error("Mixed distributions (Continuous & Discrete) not yet supported.")
    end
    min = minimum(X.distr) + minimum(Y.distr)
    max = maximum(X.distr) + maximum(Y.distr)
    if X.distr isa ContinuousDistribution && Y.distr isa ContinuousDistribution
        distr = continuous_add_convolution(X.distr, Y.distr)
        conv = ConvolutionContinuous(distr, min, max)
    elseif X.distr isa DiscreteDistribution && Y.distr isa DiscreteDistribution
        distr = discrete_add_convolution(X.distr, Y.distr)
        conv = ConvolutionDiscrete(distr[1], distr[2], distr[3])
    end
    return RV(conv)
end

(Base.:-)(X::RV, Y::RV) = begin
    if X.id == Y.id
        error("Dependent variables not yet supported.")
    elseif X.distr isa ContinuousDistribution && Y.distr isa DiscreteDistribution
        error("Mixed distributions (Continuous & Discrete) not yet supproted.")
    elseif X.distr isa DiscreteDistribution && Y.distr isa ContinuousDistribution
        error("Mixed distributions (Continuous & Discrete) not yet supported.")
    end
    min = minimum([minimum(X.distr), minimum(Y.distr)])
    max = maximum([maximum(X.distr), maximum(Y.distr)])
    if X.distr isa ContinuousDistribution && Y.distr isa ContinuousDistribution
        distr = continuous_sub_convolution(X.distr, Y.distr)
        conv = ConvolutionContinuous(distr, min, max)
    elseif X.distr isa DiscreteDistribution && Y.distr isa DiscreteDistribution
        distr = discrete_sub_convolution(X.distr, Y.distr)
        conv = ConvolutionDiscrete(distr[1], distr[2], distr[3])
    end
    return RV(conv)
end

function rv_max(dx::UnivariateDistribution, dy::UnivariateDistribution)
    function pdf(z)
        temp(x) = cdf(dx, x)*cdf(dy, x)
        return central_fdm(12, 1)(temp, z)
    end
    return pdf
end

(Base.max)(X::RV, Y::RV) = begin
    if X.id == Y.id
        error("Dependent variables not yet supported.")
    elseif X.distr isa ContinuousDistribution && Y.distr isa DiscreteDistribution
        error("Mixed distributions (Continuous & Discrete) not yet supproted.")
    elseif X.distr isa DiscreteDistribution && Y.distr isa ContinuousDistribution
        error("Mixed distributions (Continuous & Discrete) not yet supported.")
    end
    min = minimum([minimum(X.distr), minimum(Y.distr)])
    max = maximum([maximum(X.distr), maximum(Y.distr)])
    if X.distr isa ContinuousDistribution && Y.distr isa ContinuousDistribution
        distr = rv_max(X.distr, Y.distr)
        conv = ConvolutionContinuous(distr, min, max)
    elseif X.distr isa DiscreteDistribution && Y.distr isa DiscreteDistribution
        distr = rv_max(X.distr, Y.distr)
        conv = ConvolutionDiscrete(distr, min, max)
    end
    return RV(conv)
end

function rv_min(dx::UnivariateDistribution, dy::UnivariateDistribution)
    function pdf(z)
        temp(x) = 1 - (1 - cdf(X.distr, x))*(1 - cdf(Y.distr, x))
        return central_fdm(12, 1)(temp, z)
    end
    return pdf
end

(Base.min)(X::RV, Y::RV) = begin
    if X.id == Y.id
        error("Dependent variables not yet supported.")
    elseif X.distr isa ContinuousDistribution && Y.distr isa DiscreteDistribution
        error("Mixed distributions (Continuous & Discrete) not yet supproted.")
    elseif X.distr isa DiscreteDistribution && Y.distr isa ContinuousDistribution
        error("Mixed distributions (Continuous & Discrete) not yet supported.")
    end
    min = minimum([minimum(X.distr), minimum(Y.distr)])
    max = maximum([maximum(X.distr), maximum(Y.distr)])
    if X.distr isa ContinuousDistribution && Y.distr isa ContinuousDistribution
        distr = rv_min(X.distr, Y.distr)
        conv = ConvolutionContinuous(distr, min, max)
    elseif X.distr isa DiscreteDistribution && Y.distr isa DiscreteDistribution
        distr = rv_min(X.distr, Y.distr)
        conv = ConvolutionDiscrete(distr, min, max)
    end
    return RV(conv)
end

#Extend cdf and pdf for convolution case
function Distributions.cdf(d::ConvolutionContinuous, x::Real)::Real
    #Check properties of CDF
    if x == Inf
        ret = 1
    elseif x == -Inf
        ret = 0
    else
        #Possible optimization here
        ret = quadgk(z -> d.pdf(z), -Inf, x, order=100)[1]
    end
    ret
end

function Distributions.cdf(d::T1, x::T2)::Real where {T1 <: ConvolutionDiscrete, T2 <: Union{Float64, Int}}
    #Check properties of CDF
    if x == Inf
        ret = 1
    elseif x == -Inf
        ret = 0
    else
        sup = collect(d.l:x)
        ret = sum(d.pdf.(sup))
    end
    ret
end

function Distributions.pdf(d::T, x::Real)::Real where {T <: Convolution}
    #Check PMF/PDF properties (otherwise not a mass/density function)
    if x == Inf || x == -Inf
        ret = 0
    else
        ret = d.pdf(x)
    end
    ret
end

function Distributions.maximum(d::T)::Real where {T <: Convolution}
    d.u
end

function Distributions.minimum(d::T)::Real where {T <: Convolution}
    d.l
end

function Distributions.mean(d::T)::Real where {T <: ConvolutionContinuous}
    ret = quadgk(x -> x*d.pdf(x), minimum(d), maximum(d))[1]
    ret
end

function Distributions.mean(d::T)::Real where {T <: ConvolutionDiscrete}
    sup = collect((quantile(d, 1e-10)):(quantile(d, 1 - 1e-10)))
    return sum(pdf.(d, sup).*(sup))
end

function Distributions.quantile(d::ConvolutionContinuous, p::Real)::Real 
    ret = find_zero(x -> cdf(d,x) - p, 0)
    ret
end

function Distributions.quantile(d::ConvolutionDiscrete, p::Real)::Real 
    #quantile for p is smallest value of x for which CDF(x) ≥ p
    sup = collect(d.l:d.u)
    temp = vcat(cdf(d, sup), 1)
    ret = findfirst(temp .≥ p) 
    ret - 1
end