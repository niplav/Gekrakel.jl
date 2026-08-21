# Summary statistics. Each question is answered from the closed form when one
# survived, and from the particles otherwise.

# A closed form need not answer every question. `Truncated` has exact
# quantiles, cdf and pdf but, outside the normal case, no mean or variance in
# Distributions.jl. Those questions fall back to the particles — which for a
# truncation are quantile-spaced from the exact distribution, making this
# midpoint quadrature on the quantile function rather than a Monte Carlo
# estimate.
hasmoments(::Nothing) = false
hasmoments(::Distribution) = true
hasmoments(::Truncated) = false
hasmoments(::Truncated{<:Normal}) = true

Distributions.mean(n::Node) = hasmoments(n.dist) ? mean(n.dist) : pmean(parts(n))
Distributions.var(n::Node) = hasmoments(n.dist) ? var(n.dist) : pvar(parts(n))
Distributions.std(n::Node) = hasmoments(n.dist) ? std(n.dist) : pstd(parts(n))
Distributions.median(n::Node) = exact(n) ? median(n.dist) : pquantile(parts(n), 0.5)

Distributions.quantile(n::Node, p::Real) = exact(n) ? quantile(n.dist, p) : pquantile(parts(n), p)
Distributions.quantile(n::Node, ps::AbstractVector) = [quantile(n, p) for p in ps]

"""
	stdev(x)

Squiggle spells `std` this way.
"""
stdev(x) = std(x)

Distributions.cdf(n::Node, x::Real) = exact(n) ? cdf(n.dist, x) : count(<=(x), samples(n)) / length(samples(n))

function Distributions.pdf(n::Node, x::Real)
	exact(n) && return pdf(n.dist, x)
	s = samples(n)
	h = 0.9 * min(std(s), (pquantile(parts(n), 0.75) - pquantile(parts(n), 0.25)) / 1.34) * length(s)^(-0.2)
	iszero(h) && return x == s[1] ? Inf : 0.0
	sum(exp(-0.5 * ((x - v) / h)^2) for v in s) / (length(s) * h * sqrt(2pi))
end

"""
	inv_cdf(x, p)

The `p`-quantile of `x`. Squiggle calls this `inv`.
"""
inv_cdf(n::Node, p::Real) = quantile(n, p)

"""
	rand(x)
	rand(x, k)

Draw from `x`'s materialised particles. (Spelled `rand` rather than Squiggle's
`sample` because `StatsBase.sample` already owns that name.)
"""
Base.rand(n::Node) = rand(samples(n))
Base.rand(n::Node, k::Integer) = rand(samples(n), k)

"""
	credible(x, mass=0.9)

The central credible interval of `x`, as a tuple.
"""
function credible(n::Node, mass::Real=0.9)
	t = (1 - mass) / 2
	(quantile(n, t), quantile(n, 1 - t))
end
