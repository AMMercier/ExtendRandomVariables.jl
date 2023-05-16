module ExtendRandomVariables

using Distributions, RandomVariables, QuadGK, Plots, Roots

import Base.\, Base.diff, Base.intersect, Base.length, Base.iterate
import Base.copy, Base.ndims, Base.reduce, Base.issubset, Base.union, Base.xor
import Base.<, Base.>, Base.<=, Base.>=, Base.==, Base.!=, Base.in
import Base.!, Base.|, Base.&, Base.⊻
import Base.+, Base.-, Base.*, Base./, Base.inv, Base.exp, Base.log, Base.sqrt
import Base.abs, Base.^, Base.adjoint, Base.max, Base.min
import Distributions.mean, Distributions.var, Distributions.skewness, Distribution.quantile
import QuadGK.quadgk
import Roots.find_zero
import Plots.plot, Plots.plot!

include("condexp.jl")
include("convolution.jl")
include("plotRV.jl")

export 𝔼, |
export *, +, -, min, max, cdf, pdf, quantile
export plotPDF, plotCDF, plotPDF!, plotCDF!

end