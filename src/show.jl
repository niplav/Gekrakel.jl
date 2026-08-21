# Display. MonteCarloMeasurements shows `mean ± std`, which badly misdescribes
# the heavy-tailed lognormals that dominate Fermi estimates; Squiggle's
# percentile triple is the right summary, so that is what we print.

const BLOCKS = ('▁', '▂', '▃', '▄', '▅', '▆', '▇', '█')

function num(x::Real)
	!isfinite(x) && return string(x)
	a = abs(x)
	(a != 0 && (a < 1e-3 || a >= 1e5)) ? string(round(x, sigdigits=3)) :
	string(round(x, sigdigits=4))
end

"""
	spark(x; bins=28)

A unicode-block histogram of `x`, clipped to its 1st..99th percentile so a
single outlier cannot flatten the whole picture. Exact continuous nodes are
drawn from their density directly, which is both smoother and cheaper than
materialising particles for the sake of a picture.
"""
function spark(n::Node; bins::Int=28)
	lo, hi = quantile(n, 0.01), quantile(n, 0.99)
	hi > lo || return "\u2581"^bins
	mid(i) = lo + (hi - lo) * (i - 0.5) / bins

	# Exact discrete distributions get one block per outcome; binning them
	# on a continuous grid produces a comb of empty bins.
	if exact(n) && n.dist isa DiscreteUnivariateDistribution
		ks = floor(Int, lo):ceil(Int, hi)
		h = if length(ks) <= 2bins
			[pdf(n.dist, k) for k in ks]
		else
			c = zeros(Float64, bins)
			for k in ks
				c[clamp(1 + floor(Int, bins * (k - lo) / (hi - lo + 1)), 1, bins)] += pdf(n.dist, k)
			end
			c
		end
		m = maximum(h)
		return m > 0 ? String([BLOCKS[clamp(1 + round(Int, 7v / m), 1, 8)] for v in h]) : "\u2581"^length(h)
	end

	h = if exact(n) && n.dist isa ContinuousUnivariateDistribution
		[pdf(n.dist, mid(i)) for i in 1:bins]
	else
		c = zeros(Float64, bins)
		for v in samples(n)
			lo <= v <= hi || continue
			c[clamp(1 + floor(Int, bins * (v - lo) / (hi - lo)), 1, bins)] += 1
		end
		c
	end

	m = maximum(h)
	m > 0 || return "\u2581"^bins
	String([BLOCKS[clamp(1 + round(Int, 7 * v / m), 1, 8)] for v in h])
end

kind(n::Node) = exact(n) ? string(nameof(typeof(n.dist))) : "sampled"

function Base.show(io::IO, n::Node)
	p5, p50, p95 = quantile(n, [0.05, 0.5, 0.95])
	print(io, num(p5), " … ", num(p50), " … ", num(p95))
end

function Base.show(io::IO, ::MIME"text/plain", n::Node)
	p5, p50, p95 = quantile(n, [0.05, 0.5, 0.95])
	println(io, kind(n), exact(n) ? " (exact)" : " ($(length(samples(n))) samples)")
	println(io, "  ", spark(n))
	print(io, "  5% ", num(p5), "   50% ", num(p50), "   95% ", num(p95),
		"   mean ", num(mean(n)))
end

# The `MIME"text/html"` method lives in svg.jl, next to the drawing code.
