function plotPDF(X::RV, x::Real, y::Real; kwargs...)
    if typeof(X.distr) <: DiscreteDistribution
        sup = collect(x:y)
        ys = pdf(X.distr, sup)
        plot(sup, ys, seriestype=:scatter; kwargs...)
    else
        plot(x -> pdf(X.distr, x), x, y; kwargs...)
    end
end    

function plotPDF(X::RV; kwargs...) 
    if typeof(X.distr) <: DiscreteDistribution
        sup = collect((quantile(X.distr, 1e-10)):(quantile(X.distr, 1 - 1e-10)))
        ys = pdf(X.distr, sup)
        plot(sup, ys, seriestype=:scatter; kwargs...)
    else
        plot(x -> pdf(X.distr, x), quantile(X.distr, 0.01), quantile(X.distr, 0.99); kwargs...)
    end
end

function plotCDF(X::RV, x::Real, y::Real; kwargs...)
    plot(x -> cdf(X.distr, x), x, y; kwargs...)
end

plotCDF(X::RV; kwargs...) = plot(x -> cdf(X.distr, x); kwargs...)


function plotPDF!(X::RV, x::Real, y::Real; kwargs...)
    if typeof(X.distr) <: DiscreteDistribution
        sup = collect(x:y)
        ys = pdf(X.distr, sup)
        plot!(sup, ys, seriestype=:scatter; kwargs...)
    else
        plot!(x -> pdf(X.distr, x), x, y; kwargs...)
    end
end    

function plotPDF!(X::RV; kwargs...) 
    if typeof(X.distr) <: DiscreteDistribution
        sup = collect((quantile(X.distr, 1e-10)):(quantile(X.distr, 1 - 1e-10)))
        ys = pdf!(X.distr, sup)
        plot(sup, ys, seriestype=:scatter; kwargs...)
    else
        plot!(x -> pdf(X.distr, x), quantile(X.distr, 0.01), quantile(X.distr, 0.99); kwargs...)
    end
end

function plotCDF!(X::RV, x::Real, y::Real; kwargs...)
    plot!(x -> cdf(X.distr, x), x, y; kwargs...)
end

plotCDF!(X::RV; kwargs...) = plot!(x -> cdf(X.distr, x); kwargs...)
