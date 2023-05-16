# ExtendRandomVariables

A Julia package to extend the [RandomVariables.jl](https://github.com/ManuelStapper/RandomVariables.jl) package, as an add-on to the Distributions.jl package. The scope of this package is:

- Allow for conditional expectations `X=RV(Normal()); 𝔼(X|(X<3))`
- Allow for addition of two independent random variables to create a new random variable `Z=X+Y`
- Allow for the multiplcation of two independent random variables to create a new random variable `Z=X*Y`
- Allow for the maximum of two independent random variables to create a new random variable `Z=max(X,Y)`
- A suite of plotting functions

## Installation

Both Distributions.jl and RandomVariables.jl must be installed and loaded to use ExtendRandomVariables.jl

```
julia> import Pkg
julia> Pkg.add("Distributions")
julia> Pkg.add("RandomVariables")
julia> Pkg.add https://github.com/AMMercier/ExtendRandomVariables.jl.git
julia> using RandomVariables, Distributions, ExtendRandomVariables
```
## Introduction



## Random Variables
Random variables can be defined by the `RV()` from the RandomVariables.jl package function from any univariate distribution taken from the Distriutions.jl package

```
X = RV(Normal()) #Random variable X~N(0,1)
Y = RV(Poisson(4)) #Random variable Y~Poisson(4)
```

All random variables have two components: 1) an underlying distribution and 2) an ID to keep track of dependent vs. independent variables.

## Events
Probabilites are only relevant when considering events. That is, the probability function is from the event space, $\mathcal{F}$, to the reals, $\mathbb{R}$ (for the scope of this package), reflected by closed, open, or "clopen" (both closed and open) intervals. Consequently, for each event, the RandomVariables.jl package only keeps track of the random variables in the event and the intervals defined by the event, denoted cc, oo, and co or oc for closed, open, or clopen intervals, respectivly. The collection of intervals defining the event are called a "box". For $n$ independent random variables, we have an $n$ dimensional box. If we have access to the cumulative distribution function (CDF), $F$, then

$P(a< X\leq b) = F_X(b) - F_X (a)$

We define an event through various operators such as `X>3` and also combine events. For example, we can write

```
#Events
A = X > 1
B = X ≤ 3
#Combinations
A∩B
```

## Transformed Random Variables
Given a random variable $X$, we can transform that random variable with a function $f$, $f(X)$. Correspondingly, we can track what occurs to the invertal of an event with $f(X)$ such as `f(X)<x`. Therefore, to find `P(f(X)<x)` we have $P(f(X)\leq x)$

$P(f(X)\leq x) = P(X \leq f^{-1}(x)) = F_X (f^{-1}(x))$ 

Therefore, in addition to the distribution and ID, we must keep track of the function $f$ and the inverse function for the interval, `fInv`. In sum, we can write expressions such as `P(exp(X)<3)`.

The following functions are supported by RandomVariables.jl, for $a\in\mathbb{R}$:
- `X+a` 
- `X*a` 
- `inv(X)`, i.e. $Z=X^{-1}$
- `X/a` or `a/X`
- `exp(X)`
- `sqrt(X)`
- `abs(X)`
- `X^a` 

## Conditional Events
We can recall that the probability *of* event *A* conditional on the *given* event $B$ is defined as

$P(A\mid B) = \frac{P(A\cap B)}{P(B)}$

That is, the probability of event $A$ *and* $B$ normalized that event $B$ has already occured. As, the RandomVariables.jl type of events allows for `P(A∩B)`, we can easily write `P(A|B)`. 

## Conditional Expectation
The conditional expectation differs from the conditional event in that the first argument, $X$ in `𝔼(X|A)` , is a random variable and not an event. Because we usually require the conditional PDF (or PMF) to compute the conditional expectation and we have direct access to the CDF of any odd combination of random variables through the expression `P(X<x)`, we compute the conditional expectation by 

$\mathbb{E}(X|A) = \int_{-\infty}^0 1-F_{X|A}(x) dx - \int_0^{\infty} F_{X|A}(x) dx$

Thus, we can compute the conditional expectation for either continuous or discrete random variables using the `𝔼` expression (`\bbE` followed by `tab` in Julia):
```
𝔼(X)
𝔼(X|(X<3))
𝔼(X|(exp(X)<2)∩(abs(X)<1.5))
```
It should be noted that if $X$ is independent of $A$, then $\mathbb{E}(X\mid A) = \mathbb{E}(X)$.

Lastly, mathematically, for $\mathbb{E}(X\mid Y)$ where both $X$ and $Y$ are random variables, we have that $\mathbb{E}(X\mid Y) = g(Y)$. That is, we obtain a function of the given random variable $Y$. This last application of the conditional expectation has been left for future work, as it would require an explicit dependence of $X$ and $Y$ where so far we have assumed $X$ and $Y$ are independent. 

## Addition, Multiplication, and Max/Min
We can find the addition or multiplication of two random variables through a "convolutional" framework. For addition of two independent random variables, we have that 

$P(X+Y \leq z) = \int_{-\infty}^{\infty} f_X(z-x)f_Y(x) dx$
$P(X-Y \leq z) = \int_{-\infty}^{\infty} f_X(z+x)f_Y(x) dx$

which allows `X+Y` and `X-Y`. Similarly, for the multiplication of two random variables $X$ and $Y$

$P(XY \leq z) = \int_{-\infty}^{\infty} f_X(x)f_Y(z/x)\frac{1}{|x|} dx$

This is implemented using Gaussian quadrature from the QuadGK.jl package. It is worth noting that for some random variables, the order of the corresponding polynomial is quite large and computation time relativly slow. 

For `max` or `min` of two independent random variables $X$ and $Y$, we have that we can most easily identify the CDF. Letting $Z_1=\max(X,Y)$ and $Z_2 =\min(X,Y)$, we have that

$F_{Z_1}(z) = F_X(z)F_Y(z) $
$F_{Z_2}(z) = 1 - (1-F_X(z))(1-F_Y(z))$

Therefore, if we wish to have the PDF, we must compute $f_Z(z) = \frac{d}{dz} F_Z(z)$. As the RandomVariable.jl types do not support automatic differentiation, we can use finite differencing methods provided by the FiniteDifferences.jl package. 

Relatedly, if we wish to determine the quantile function for our random variables, we can find the inverse of the corresponding CDF, $F^{-1}(x)$. For continuous CDF, we can do this through a root finding method using the Roots.jl package. That is, for the p-th quantile, we solve for the root $F(x) - p$. This is accomplished through by the code `find_zero(x -> cdf(d,x) - p, 0)`. For a discrete CDF, we have that for the p-th quantile, we find the smallest value of $x$ for which $F(x) \geq p$. We can accomplish this computationally by running over the support and determining the first element of the support that fulfills the condition $F(x) \geq p$. If the lower bound or the upper bound of the support is infinite, we can run from or to the `1e-10` quantile or `1-1e-10` quantile, respectively.   

## Plots
When considering a random variable, there are two primary descriptions we would wish to plot: the cummulative distribution function (CDF) and the probability density (or mass) function (PDF or PMF). Consequently, we can write:
```
X = RV(Normal())
Y = RV(Exponential(3))
Z = X+Y
plotCDF(Z, -6, 10, xlabel="test")
plotPDF(Z, -6, 10, xlabel="test")
```

This provides the following plots:

![Alt text](https://i.ibb.co/PG9KgmK/CDFPlot.png "CDF")
![Alt text](https://i.ibb.co/rdRQYvk/PDFPlot.png "PDF")

## Future Work
The largest component of future work is to implement dependent random variables, which would allow more interesting conditional moments. This could be accomplished through a variety of methods, including covariance structure or including interpolated, emperical distributions. Additionally, more arithmetic operations on random variables (such as `X/Y`) would also be appreciated. Moreover, more general optimizations could be made, especially with respect to the addition and multiplcation of random variables. Lastly, reworking the base types from RandomVariables.jl to be compatiable with automatic differentiation, including differentiation of random variables (i.e., Malliavin calculus), would be a more ambiuous extension of the RandomVariables.jl.

## References:
[RandomVariables.jl package](https://github.com/ManuelStapper/RandomVariables.jl)
[Distributions.jl package](https://github.com/JuliaStats/Distributions.jl)
[QuadGK.jl package](https://github.com/JuliaMath/QuadGK.jl)
[Roots.jl package](https://github.com/JuliaMath/Roots.jl)
[Plots.jl package](https://github.com/JuliaPlots/Plots.jl)
