# Distribution constructors, named after their Squiggle counterparts.

const Z95 = 1.6448536269514722  # quantile(Normal(), 0.95)

normal(mu::Real, sigma::Real) = leaf(Normal(mu, sigma))
lognormal(mu::Real, sigma::Real) = leaf(LogNormal(mu, sigma))
uniform(lo::Real, hi::Real) = leaf(Uniform(lo, hi))
cauchy(loc::Real, scale::Real) = leaf(Cauchy(loc, scale))
gamma(shape::Real, scale::Real) = leaf(Gamma(shape, scale))
exponential(rate::Real) = leaf(Exponential(1 / rate))
bernoulli(p::Real) = leaf(Bernoulli(p))
poisson(lambda::Real) = leaf(Poisson(lambda))
triangular(lo::Real, mid::Real, hi::Real) = leaf(TriangularDist(lo, mid, hi))
beta(a::Real, b::Real) = leaf(Beta(a, b))
pointmass(x::Real) = leaf(Dirac(float(x)))

# 90% credible interval forms, e.g. normal(p5=1, p95=10).
function normal(; p5::Real, p95::Real)
	p5 < p95 || error("need p5 < p95")
	normal((p5 + p95) / 2, (p95 - p5) / (2Z95))
end

function lognormal(; p5::Real, p95::Real)
	0 < p5 < p95 || error("need 0 < p5 < p95")
	lognormal((log(p5) + log(p95)) / 2, (log(p95) - log(p5)) / (2Z95))
end

# Squiggle's beta({mean, stdev}), same closed form.
function beta(; mean::Real, stdev::Real)
	0 < mean < 1 || error("beta mean must lie in (0,1)")
	v = stdev^2
	v < mean * (1 - mean) || error("stdev too large for that mean")
	k = mean * (1 - mean) / v - 1
	beta(mean * k, (1 - mean) * k)
end

"""
	to(lo, hi)
	lo … hi

The distribution whose 5th and 95th percentiles are `lo` and `hi`: lognormal
when `lo > 0`, normal otherwise. Matches Squiggle's `to`.

!!! warning "Precedence"
	`…` binds *looser* than `+` and `*`, so `10…100 * 2` parses as `10…200`,
	not `(10…100) * 2`. Squiggle's `to` behaves the same way. Parenthesise
	whenever you combine.
"""
function to(lo::Real, hi::Real)
	lo < hi || error("need lo < hi, got $lo … $hi")
	lo > 0 ? lognormal(p5=lo, p95=hi) : normal(p5=lo, p95=hi)
end

…(lo::Real, hi::Real) = to(lo, hi)
