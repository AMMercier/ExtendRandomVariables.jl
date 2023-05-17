# ExtendRandomVariables

A Julia package to extend the [RandomVariables.jl](https://github.com/ManuelStapper/RandomVariables.jl) package, as an add-on to the Distributions.jl package. The scope of this package is:

- Allow for conditional expectations `X=RV(Normal()); 𝔼(X|(X<3))`
- Allow for addition of two independent random variables to create a new random variable `Z=X+Y`
- Allow for the multiplication of two independent random variables to create a new random variable `Z=X*Y`
- Allow for the maximum and minimum of two independent random variables to create a new random variable `Z=max(X,Y)`, `Z=min(X,Y)`
- A plotting functions for the CDF and PDF of random variables

## Installation
Both Distributions.jl and RandomVariables.jl must be installed and loaded to use ExtendRandomVariables.jl
```
julia> import Pkg
julia> Pkg.add("Distributions")
julia> Pkg.add("RandomVariables")
julia> Pkg.add("QuadGK")
julia> Pkg.add("Roots")
julia> Pkg.add("Plots")
julia> Pkg.add(url = "https://github.com/AMMercier/ExtendRandomVariables.jl.git")
julia> using RandomVariables, Distributions, ExtendRandomVariables, Plots
```
## Introduction
### *"In order to make a perfect and beautiful machine, it is not requisite to know how to make it."* 
##### - Robert Beverley MacKenzie

Mathematics is often a "beautiful machine", where we set definite, abstract objects, and allow them to interact in strange, and emergent ways. Consequently, it is often not necessary to enumerate all the pathways by which these constructs interact, but rather simply provide the ruleset by which they interrelate and allow the user to mix and match as they please. This "toy box" approach to mathematics is easily applicable to - and perhaps even highly beneficial for - computer programming. In this way, the goal is not to create a single program, but rather an ecosystem of interacting types and operations to allow emergent behavior. This synergy between mathematics and programming presents an opportunity in a critical area to modern pure and applied mathematics: random variables.

Random variables, a key concept in probability theory, hold paramount importance in various domains, ranging from finance to biology. By their very nature, random variables introduce uncertainty into mathematical models, allowing us to capture the inherent unpredictability of real-world phenomena. Random variables are a function which map from a sample space, $\Omega$, to a measurable space, $E$, which for the remainder of this paper will be the measurable space induced by the reals, $\mathbb{R}$. This links the "real world" (the sample space) to a mathematically rich, and analyzable structure (the measurable space). This means that Julia, with its handling of functions, is a natural setting to attempt to implement random variables as their own abstract type. Consequently, through statistical analysis and probabilistic reasoning, mathematicians can harness the power of random variables to make informed decisions, quantify risk, and gain deeper insights into complex systems.

Therefore, the goal of this project is to extend the usefulness of random variables as types form the RandomVariables.jl in Julia to allow operations on random variables in an analogous manner as floats, with easily approachable and familiar syntax, and conditional moments. Additionally, we introduce plotting of the CDF and PDF or PMF of the distribution given by the random variable. In this way, the rules for a "beautiful machine" can be set such that any given user can utilize rules to craft their own emergent machine. 

While the majority of the following discussion applies to continuous random variables, the same rules analogously  apply to discrete random variables, as well.

## Random Variables
Random variables can be defined by the `RV()` from the RandomVariables.jl package function from any univariate distribution taken from the Distriutions.jl package

```
X = RV(Normal()) #Random variable X~N(0,1)
Y = RV(Poisson(4)) #Random variable Y~Poisson(4)
```

All random variables have two components: 1) an underlying distribution and 2) an ID to keep track of dependent vs. independent variables.

## Events
Probabilites are only relevant when considering events. That is, the probability function is from the event space, $\mathcal{F}$, to the reals, $\mathbb{R}$ (for the scope of this package), reflected by closed, open, or "clopen" (both closed and open) intervals. Consequently, for each event, the RandomVariables.jl package only keeps track of the random variables in the event and the intervals defined by the event, denoted `cc`, `oo`, and `co` or `oc` for closed, open, or clopen intervals, respectively. The collection of intervals defining the event is called a "box". For $n$ independent random variables, we have an $n$ dimensional box. If we have access to the cumulative distribution function (CDF), $F$, then

$P(a< X\leq b) = F_X(b) - F_X (a)$

We define an event through various operators such as `X>3` and combine events. For example, we can write

```
#Events
A = X > 1
B = X ≤ 3
#Combinations
A∩B
```

## Transformed Random Variables
Given a random variable $X$, we can transform that random variable with a function $f$, $f(X)$. Correspondingly, we can track what occurs to the interval of an event with $f(X)$ such as `f(X)<x`. Therefore, to find `P(f(X)<x)` we have $P(f(X)\leq x)$

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

That is, the probability of event $A$ *and* $B$ normalized that event $B$ has already occurred. As, the RandomVariables.jl type of events allows for `P(A∩B)`, we can easily write `P(A|B)`. 

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

Lastly, mathematically, for $\mathbb{E}(X\mid Y)$ where both $X$ and $Y$ are random variables, we have that $\mathbb{E}(X\mid Y) = g(Y)$. That is, we obtain a function of the given random variable $Y$. This last application of the conditional expectation has been left for future work, as it would require an explicit dependence of $X$ and $Y$, where so far we have assumed $X$ and $Y$ are independent. 

## Addition, Multiplication, and Max/Min
We can find the addition or multiplication of two random variables through a "convolutional" framework. For addition of two independent random variables, we have that 

$P(X+Y \leq z) = \int_{-\infty}^{\infty} f_X(z-x)f_Y(x) dx$
$P(X-Y \leq z) = \int_{-\infty}^{\infty} f_X(z+x)f_Y(x) dx$

which allows `X+Y` and `X-Y`. Similarly, for the multiplication of two random variables $X$ and $Y$

$P(XY \leq z) = \int_{-\infty}^{\infty} f_X(x)f_Y(z/x)\frac{1}{|x|} dx$

This is implemented using Gaussian quadrature from the QuadGK.jl package. It is worth noting that for some random variables, the order of the corresponding polynomial is quite large and computation time relatively slow. Other more computationally efficient options to complete these convolutions could be explored in the future, such as FTT.

For `max` or `min` of two independent random variables $X$ and $Y$, we have that we can most easily identify the CDF. Letting $Z_1=\max(X,Y)$ and $Z_2 =\min(X,Y)$, we have that

$F_{Z_1}(z) = F_X(z)F_Y(z) $
$F_{Z_2}(z) = 1 - (1-F_X(z))(1-F_Y(z))$

Therefore, if we wish to have the PDF, we must compute $f_Z(z) = \frac{d}{dz} F_Z(z)$. As the RandomVariable.jl types do not support automatic differentiation, we can use finite differencing methods provided by the FiniteDifferences.jl package. A future reworking of the RandomVariables.jl types

Relatedly, if we wish to determine the quantile function for our random variables, we can find the inverse of the corresponding CDF, $F^{-1}(x)$. For continuous CDF, we can do this through a root finding method using the Roots.jl package. That is, for the p-th quantile, we solve for the root of $F(x) - p$. This is accomplished through by the code `find_zero(x -> cdf(d,x) - p, 0)` where `d` is our distribution. For a discrete CDF, we have that for the p-th quantile, we find the smallest value of $x$ for which $F(x) \geq p$. We can accomplish this computationally by running over the support and determining the first element of the support that fulfills the condition $F(x) \geq p$. If the lower bound or the upper bound of the support is infinite, we can run from or to the `1e-10` quantile or `1-1e-10` quantile, respectively.   

Taken together, we can do some preliminary checks. Consider the following examples:

```
#Continuous random variables
julia> X1 = RV(Normal()); Y1 = RV(Normal(3))
RV(Normal{Float64}(μ=3.0, σ=1.0), 11)

julia> Z1=X1+Y1; 𝔼(Z1)
3.0000000000000226
julia> W1=X1*Y2; 𝔼(W1)
-8.719078556333654e-16

#Discrete random variables
julia> X2 = RV(Poisson(3)); Y2 = RV(Poisson(4))
RV(Poisson{Float64}(λ=4.0), 17)

julia> Z2 = X2 + Y2; 𝔼(Z2)
6.99999999697663
julia> W2 = X2*Y2; 𝔼(W2)
11.999971100191376

#Misc. Examples:
julia> 𝔼(W1|(exp(W1) < 1))
-2.3942635174048403
julia> 𝔼(W2|(exp(W2) < 3))
0.38603161430513455
```

## Plots
When considering a random variable, there are two primary descriptions we would wish to plot: the cumulative distribution function (CDF) and the probability density (or mass) function (PDF or PMF). Consequently, we can write:
```
X = RV(Normal())
Y = RV(Exponential(3))
Z = X+Y
plotCDF(Z, -6, 10, xlabel="test")
plotPDF(Z, -6, 10, xlabel="test")
```

This provides the following plots of the CDF and PDF, respectively.

![CDF](https://github.com/AMMercier/ExtendRandomVariables.jl/blob/main/images/CDFPlot.png "CDF")
![PDF](https://github.com/AMMercier/ExtendRandomVariables.jl/blob/main/images/PDFPlot.png "PDF")

We can also use the "bang" operator `!` via Julia convention for plotting: `plotCDF!` and `plotPDF!`.

## Future Work
The largest component of future work is to implement dependent random variables, which would allow more interesting conditional moments. This could be accomplished through a variety of methods, including covariance structure or including interpolated, empirical distributions. Moreover, more arithmetic operations on random variables (such as `X/Y`) would also be appreciated. THis extends to transformations of combinations of random variables, such as `exp(X+Y)` or `exp(X*Y)`. Broadly speaking, for more complex distributions, Monte Carlo methods could be used at the cost of either precision, accuracy, or performance. Additionally, more general optimizations could be made, especially with respect to the addition and multiplication of random variables. Furthermore, random variables from mixed distributions (both continuous and discrete) would also be greatly appreciated and would allow for even greater ease of use for the user. Lastly, reworking the base types from RandomVariables.jl to be compatible with automatic differentiation, including differentiation of random variables (i.e., Malliavin calculus), would be a more ambitious extension of the RandomVariables.jl. It is a hope that at some point in the future random variables may be implemented, up to their limitations, similar to integers, floats, or other types used daily.

## Author:
Alexander M. Mercier

## References:
MacKenzie, Robert Beverley. "The Darwinian theory of the transmutation of species examined." London: Nisbet & Co (1868): 318.

[RandomVariables.jl package](https://github.com/ManuelStapper/RandomVariables.jl)

[Distributions.jl package](https://github.com/JuliaStats/Distributions.jl)

[QuadGK.jl package](https://github.com/JuliaMath/QuadGK.jl)

[Roots.jl package](https://github.com/JuliaMath/Roots.jl)

[Plots.jl package](https://github.com/JuliaPlots/Plots.jl)
