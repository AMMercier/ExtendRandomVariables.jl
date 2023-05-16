struct condRVevent
    of::RV
    given::event
    id::Int64
end

(Base.:|)(X::T1, A::T2) where {T1 <: RandomVariable, T2 <: event} = begin
    condRVevent(X, A, X.id)
end

𝔼(X::RV) = E(X)

function 𝔼(X::condRVevent)
    Y = X.of
    A = X.given
    if typeof(X.of.distr) <: ContinuousDistribution
    #Continuous expectation: 𝔼[X] = ∫(0,∞) 1-P(X≤x) dx - ∫(-∞,0) P(X≤x) dx 
        val = quadgk(y -> P((X.of≥y)|X.given), 0, Inf, order=200)[1]
        val -= quadgk(y -> P((X.of≤y)|X.given), -Inf, 0, order=200)[1] 
    elseif typeof(X.of.distr) <: DiscreteDistribution
        val = condexp(X.of, X.given)
    end
    val
end

function condexp(X::RV, A::event) #For discrete
    sup = collect(quantile(X.distr, 1e-10):quantile(X.distr, 1-1e-10))
    val = 0
    for x in sup
        val += x * P((X == x) ∩ A)
    end
    val * 1/P(A)
end

function condexp(X::RVtransformed, A::event) #For discrete
    sup = collect(quantile(X.distr, 1e-10):quantile(X.distr, 1-1e-10))
    val = 0
    for x in sup
        val += X.f(x) * P((X == x) ∩ A)
    end
    val * 1/P(A)
end 
