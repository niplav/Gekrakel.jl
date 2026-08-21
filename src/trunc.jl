# Truncation and clamping. These are different operations and Squiggle keeps
# them apart too, but the distinction is easy to lose:
#
#   truncate  conditions on lo <= x <= hi, renormalising. Mass outside the
#             bounds is *removed*, so the result is a new source of randomness
#             and gets a fresh leaf id — it is no longer correlated with its
#             parent, because it is not a function of it.
#   clamp     is pointwise min/max. Mass outside the bounds piles up on the
#             bounds. Correlation with the parent is preserved.
#
# If you want "at least 0.05 years" you want clamp; if you want "given that it
# was under 1" you want truncate.

"""
	truncate(x, lo, hi)
	truncateLeft(x, lo)
	truncateRight(x, hi)

Condition `x` on lying in the given range. Exact when a closed form survives
(via `Distributions.truncated`), otherwise by resampling the particles that
fall inside.

Conditioning severs the link to the parent: `truncate(x, 0, 1) - x` is not
zero. Use [`clamp`](@ref) for the pointwise operation.
"""
function Base.truncate(n::Node, lo::Real, hi::Real)
	lo < hi || error("need lo < hi, got $lo … $hi")
	if exact(n)
		t = Node(truncated(n.dist, lo, hi), BitSet([newid()]), () -> nothing, nothing)
		t.mk = () -> _qparts(t.dist, NSAMP[])
		return t
	end
	s = filter(x -> lo <= x <= hi, samples(n))
	isempty(s) && error("truncating to [$lo, $hi] left no samples")
	d = Node(nothing, BitSet([newid()]), () -> nothing, nothing)
	d.mk = () -> Particles([s[i] for i in _spread(length(s), NSAMP[])])
	d
end

truncateLeft(n::Node, lo::Real) = truncate(n, lo, Inf)
truncateRight(n::Node, hi::Real) = truncate(n, -Inf, hi)

# Stretch m surviving samples back out to k, evenly rather than randomly: the
# survivors are already quantile-spaced, so walking them in order preserves
# that and avoids adding Monte Carlo noise on top of the conditioning.
_spread(m::Integer, k::Integer) = [1 + ((i - 1) * m) ÷ k for i in 1:k]

# Quantile-spaced sampling, then shuffled — MonteCarloMeasurements' own
# systematic sampler does the same, but it asks the distribution for its mean
# first, which `Truncated` does not implement. The shuffle is what makes two
# independently built nodes uncorrelated.
_qparts(d::Distribution, n::Integer) =
	Particles(shuffle!([quantile(d, (i - 0.5) / n) for i in 1:n]))

"""
	clamp(x, lo, hi)

Pointwise `min(max(x, lo), hi)`, sample by sample. Correlation with `x` is
preserved, and mass outside the bounds piles up on them rather than being
renormalised away — the difference from [`truncate`](@ref).
"""
Base.clamp(n::Node, lo::Real, hi::Real) = over(x -> clamp(x, lo, hi), n)
Base.max(n::Node, lo::Real) = over(x -> max(x, lo), n)
Base.min(n::Node, hi::Real) = over(x -> min(x, hi), n)
Base.max(a::Real, n::Node) = max(n, a)
Base.min(a::Real, n::Node) = min(n, a)
Base.max(a::Node, b::Node) = over(max, a, b)
Base.min(a::Node, b::Node) = over(min, a, b)
