"""
	mixture(parts, weights=uniform)
	mx(parts, weights=uniform)

Mix distributions (or plain numbers, which become point masses). Weights are
normalised. The mixture stays exact when every component is exact and no two
components share a leaf.
"""
function mixture(cs::AbstractVector, ws::AbstractVector=fill(1.0, length(cs)))
	length(cs) == length(ws) || error("need one weight per component")
	length(cs) > 0 || error("mixture of nothing")
	all(>=(0), ws) || error("weights must be non-negative")
	ns = Node[c isa Node ? c : pointmass(c) for c in cs]
	w = ws ./ sum(ws)

	deps = BitSet([newid()])          # the selector is itself a source of randomness
	for n in ns
		union!(deps, n.deps)
	end

	shared = sum(length(n.deps) for n in ns) != length(union((n.deps for n in ns)...))
	d = (!shared && all(exact, ns) && mixable(ns)) ?
		MixtureModel([n.dist for n in ns], w) : nothing

	derive(d, deps, function ()
		ps = [parts(n) for n in ns]
		N = length(ps[1].particles)
		cut = cumsum(w)
		out = Vector{Float64}(undef, N)
		for i in 1:N
			r = rand()
			k = findfirst(>=(r), cut)
			out[i] = ps[k === nothing ? length(ps) : k].particles[i]
		end
		Particles(out)
	end)
end

const mx = mixture

# `MixtureModel` cannot combine continuous with discrete components, and
# `pointmass` is a `Dirac`, which Distributions.jl classifies as discrete. So
# `mx([normal(0,1), 5])` — a distribution mixed with a plain number, the case
# the docs advertise — has no closed form even though every component has one.
# Without this check it threw. The particle path handles it fine.
mixable(ns) = all(n -> n.dist isa ContinuousUnivariateDistribution, ns) ||
	all(n -> n.dist isa DiscreteUnivariateDistribution, ns)
