# Adapter onto StochasticOrders.jl.
#
# The orders themselves are not about Gekrakel: they need two index-aligned
# sample vectors and nothing else, so they live in their own package and work
# equally well on an MCMC chain or a CSV column. What Gekrakel adds is the three
# things that *do* need a `Node`:
#
#   1. unwrapping — `samples(x)`, so the orders take distributions directly;
#   2. the integrability guard, which needs the closed form to detect;
#   3. `certify`, which decides dominance from parameters instead of estimating.
#
# Splitting it this way forces the honest statement of what depends on
# Gekrakel's provenance machinery and what is plain marginal statistics. The
# answer turned out to be: almost nothing does.

# Import every name Gekrakel re-exports or adds a method to, explicitly. A
# `using` would bind them too, but then `Gekrakel.jl`'s `export` list would
# depend on what StochasticOrders happens to export, and an upstream rename
# would change Gekrakel's public surface without saying so.
import StochasticOrders: relate, fosd, sosd, icx, statewise, pgeq, wass,
	margin, DomMargin, lub, glb, ⊔, ⊓, envelope, regret, domgraph, hasse,
	domlevels, undominated, dominating, margins, robust, score

# The second-order and convex orders integrate the cdf, and those integrals
# diverge for a non-integrable tail: population SOSD is undefined for a Cauchy,
# while the empirical partial sums are always finite and would return a
# confident verdict about a quantity that does not exist. A category error, not
# noise — and `algebra.jl` keeps Cauchy closed under `+` and scalar `*`, so such
# nodes are easy to build by accident.
#
# Only the closed-form case is detectable. A *sampled* node descended from a
# Cauchy looks like any other vector, and `Node` stores a materialisation thunk
# rather than an operation, so there is no expression tree to walk. Documented
# rather than guessed at.
function _integrable(n::Node)
	exact(n) && !isfinite(mean(n.dist)) && error(
		"$(nameof(typeof(n.dist))) has no finite mean, so ∫F diverges and " *
		"second-order dominance is undefined for it. Use :first, which is " *
		"mean-free.")
	nothing
end

_guard(n::Node, order::Symbol) = (order in (:second, :convex) && _integrable(n); n)
score(n::Node, order::Symbol) = score(samples(_guard(n, order)), order)

relate(x::Node, y::Node; order::Symbol=:first) =
	relate(samples(_guard(x, order)), samples(_guard(y, order)); order)

for f in (:fosd, :sosd, :icx, :statewise)
	ord = f === :fosd ? :(:first) : f === :sosd ? :(:second) :
		f === :icx ? :(:convex) : :(:statewise)
	@eval $f(x::Node, y::Node) = relate(x, y; order=$ord) in (:dominates, :equal)
end

pgeq(x::Node, y::Node) = pgeq(samples(x), samples(y))
wass(x::Node, y::Node) = wass(samples(x), samples(y))
margin(x::Node, y::Node; order::Symbol=:first) =
	margin(samples(_guard(x, order)), samples(_guard(y, order)); order)

# The lattice operations come back as sorted vectors; wrap them into nodes.
# Shuffled, because a sorted particle vector would correlate perversely with
# anything it later met — the same reason `trunc.jl` shuffles. The result is
# nobody's function, so it gets a fresh leaf and no closed form.
_node(v::Vector{Float64}) = gek(shuffle!(v))

lub(x::Node, y::Node) = _node(lub(samples(x), samples(y)))
glb(x::Node, y::Node) = _node(glb(samples(x), samples(y)))
envelope(xs::AbstractVector{Node}) = _node(envelope([samples(x) for x in xs]))
regret(xs::AbstractVector{Node}) = regret([samples(x) for x in xs])

for f in (:domgraph, :hasse, :domlevels, :undominated, :dominating, :margins)
	@eval $f(xs::AbstractVector{Node}; order::Symbol=:first) =
		$f([samples(_guard(x, order)) for x in xs]; order)
end

robust(xs::AbstractVector{Node}; gammas=0:0.5:4) =
	robust([samples(x) for x in xs]; gammas)

# --------------------------------------------------------- certified verdicts

"""
	certify(x, y)

Decide dominance from closed forms when possible: `true`/`false` is a *proof*,
`nothing` means undecided. Never estimates — that is what [`fosd`](@ref) is for.

Only pairs where both nodes are `exact` can be decided, and only within a
family. A recursive certifier over the expression DAG would decide far more
(monotone maps preserve dominance, independent sums preserve it, and
`isdisjoint(a.deps, b.deps)` is already the independence test) — but `Node`
stores a materialisation thunk rather than an operation and its children, so
there is no tree to walk. That would need a change to `node.jl`.
"""
function certify(x::Node, y::Node)
	(exact(x) && exact(y)) || return nothing
	_cert(x.dist, y.dist)
end

_cert(a, b) = nothing
_cert(a::Dirac, b::Dirac) = a.value >= b.value

# Location-scale on unbounded support: (t-μ₁)/σ₁ ≤ (t-μ₂)/σ₂ for all t needs
# equal slopes, so unequal scale forces a crossing — and `false` there is a
# genuine proof of *non*-dominance, not an undecided verdict.
#
# `≈` rather than `==` on the scale: it comes out of `(log(hi) - log(lo)) / 2Z`
# for a credible-interval constructor, so `1…10` and `2…20` agree
# mathematically but differ in the last bits.
_scaleeq(x, y) = isapprox(x, y; rtol = 1e-12)
_cert(a::Normal, b::Normal) = _scaleeq(a.σ, b.σ) ? a.μ >= b.μ : false
_cert(a::LogNormal, b::LogNormal) = _scaleeq(a.σ, b.σ) ? a.μ >= b.μ : false
_cert(a::Cauchy, b::Cauchy) = _scaleeq(a.σ, b.σ) ? a.μ >= b.μ : false

# Anchored at zero, pure scale, so the cdfs cannot cross.
_cert(a::Exponential, b::Exponential) = a.θ >= b.θ
_cert(a::Gamma, b::Gamma) = a.α >= b.α && a.θ >= b.θ
_cert(a::Uniform, b::Uniform) = a.a >= b.a && a.b >= b.b
_cert(a::Beta, b::Beta) = a.α >= b.α && a.β <= b.β
_cert(a::Bernoulli, b::Bernoulli) = a.p >= b.p
_cert(a::Poisson, b::Poisson) = a.λ >= b.λ

# `n₁ ≥ n₂ ∧ p₁ ≥ p₂` is sufficient but not necessary — Binomial(100, 0.5)
# dominates Binomial(1, 0.9). The support is finite, so decide it exactly.
function _cert(a::Binomial, b::Binomial)
	for k in 0:max(a.n, b.n)
		cdf(a, k) > cdf(b, k) && return false
	end
	true
end

# Same base distribution: the quantile function is F⁻¹ of an affine
# interpolation between F(lo) and F(hi), so ordering both bounds suffices.
function _cert(a::Truncated, b::Truncated)
	a.untruncated == b.untruncated || return nothing
	minimum(a) >= minimum(b) && maximum(a) >= maximum(b)
end
