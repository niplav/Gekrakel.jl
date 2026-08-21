# The core value type.
#
# A Node is a lazily-materialised random quantity. It carries three things:
#
#   dist   a closed-form Distribution when one is known, else nothing.
#          Used *only* to answer summary questions (mean, quantile, cdf...)
#          exactly. Never sampled from.
#   deps   the set of leaf ids this value depends on. Two values may be
#          combined analytically only when their dep sets are disjoint;
#          otherwise they are correlated and must be combined numerically.
#   mk     a thunk rebuilding this value's particles from its *children's*
#          cached particles. Materialising through the tree (rather than
#          resampling `dist`) is what keeps correlations honest.

const NSAMP = Ref(2000)
const NEXTID = Ref(0)

newid() = (NEXTID[] += 1)

abstract type Gek <: Real end

mutable struct Node <: Gek
	dist::Union{Distribution,Nothing}
	deps::BitSet
	mk::Function
	cache::Any
end

# A fresh independent source of randomness.
function leaf(d::Distribution)
	n = Node(d, BitSet([newid()]), () -> nothing, nothing)
	n.mk = () -> Particles(NSAMP[], d)
	n
end

derive(d, deps, mk) = Node(d, deps, mk, nothing)

"""
	parts(x)

The particle representation of `x`, materialised and memoised on first use.
"""
function parts(n::Node)
	if n.cache === nothing
		n.cache = n.mk()
	end
	n.cache
end
parts(x::Real) = x
parts(x::AbstractParticles) = x

"""
	samples(x)

The raw sample vector backing `x`, as a plain `Vector{Float64}`.
"""
samples(n::Node) = parts(n).particles
samples(p::AbstractParticles) = p.particles

"""
	nsamp()
	setnsamp!(n)

Global particle count. Nodes already materialised keep the count they were
built with (`Particles{T,N}` puts N in the type, so mixing counts is an
error); call this before building anything.
"""
nsamp() = NSAMP[]
setnsamp!(n::Integer) = (NSAMP[] = Int(n))

# Does `x` still have an exact answer available?
exact(n::Node) = n.dist !== nothing

"""
	gek(x)

Wrap something as a `Node`: a Distributions.jl distribution, raw
`Particles`, a plain sample vector, or a number. This is the door the rest of
the ecosystem comes in through — Turing chains, fetched forecasts, bootstrap
resamples, anything you already have samples of.

Sample vectors and particles arrive without provenance, so they get a fresh
leaf id and no closed form.
"""
gek(d::Distribution) = leaf(d)
gek(n::Node) = n
function gek(p::AbstractParticles)
	length(p.particles) == NSAMP[] || error(
		"got $(length(p.particles)) samples but the global count is $(NSAMP[]); " *
		"`Particles{T,N}` puts N in the type, so these could not be combined. " *
		"Call setnsamp!($(length(p.particles))) before building anything else.")
	Node(nothing, BitSet([newid()]), () -> p, p)
end
gek(v::AbstractVector{<:Real}) = gek(Particles(collect(float.(v))))

"""
	over(f, x...)

Apply `f` sample-by-sample. Needed for anything that has to see one number at
a time — branches, lookups, `if` — since a distribution cannot answer `x > 0`.

	over(x -> x > 1e5 ? 0.45x : 0.25x, income)

Correlations are preserved: the i-th sample of every argument is passed
together. Plain numbers may be mixed in freely; they are passed through
unchanged to every call.
"""
function over(f, xs::Union{Node,Real}...)
	ns = filter(x -> x isa Node, xs)
	isempty(ns) && error("over needs at least one Node argument")
	derive(nothing, union((n.deps for n in ns)...), function ()
		# a Real becomes a constant stream; a Node its sample vector
		ss = [x isa Node ? samples(x) : nothing for x in xs]
		k = length(samples(ns[1]))
		Particles([f((s === nothing ? xs[j] : s[i] for (j, s) in enumerate(ss))...)
			for i in 1:k])
	end)
end
