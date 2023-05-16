function plotPDF(X::RV, x::Real, y::Real; kwargs...)
    plot(x -> pdf(X.distr, x), x, y; kwargs...)
end    

plotPDF(X::RV; kwargs...) = plot(x-> pdf(X.distr, x); kwargs...)

function plotCDF(X::RV, x::Real, y::Real; kwargs...)
    plot(x -> cdf(X.distr, x), x, y; kwargs...)
end

plotCDF(X::RV; kwargs...) = plot(x -> cdf(X.distr, x); kwargs...)

function plotPDF(X::Convolution, x::Real, y::Real; kwargs...)
    plot(x -> pdf(X.distr, x), x, y; kwargs...)
end    

plotPDF(X::Convolution; kwargs...) = plot(x-> pdf(X.distr, x); kwargs...)

function plotCDF(X::Convolution, x::Real, y::Real; kwargs...)
    plot(x -> cdf(X.distr, x), x, y; kwargs...)
end

plotCDF(X::Convolution; kwargs...) = plot(x -> cdf(X.distr, x); kwargs...)

function plotPDF!(X::RV, x::Real, y::Real; kwargs...)
    plot!(x -> pdf(X.distr, x), x, y; kwargs...)
end    

plotPDF!(X::RV; kwargs...) = plot!(x-> pdf(X.distr, x); kwargs...)

function plotCDF!(X::RV, x::Real, y::Real; kwargs...)
    plot!(x -> cdf(X.distr, x), x, y; kwargs...)
end

plotCDF!(X::RV; kwargs...) = plot!(x -> cdf(X.distr, x); kwargs...)

function plotPDF!(X::Convolution, x::Real, y::Real; kwargs...)
    plot!(x -> pdf(X.distr, x), x, y; kwargs...)
end    

plotPDF(X::Convolution; kwargs...) = plot(x-> pdf(X.distr, x); kwargs...)

function plotCDF!(X::Convolution, x::Real, y::Real; kwargs...)
    plot!(x -> cdf(X.distr, x), x, y; kwargs...)
end

plotCDF!(X::Convolution; kwargs...) = plot!(x -> cdf(X.distr, x); kwargs...)
