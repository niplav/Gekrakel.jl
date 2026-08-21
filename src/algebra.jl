# Analytic closure table.
#
# `closure(op, a, b)` returns the exact distribution of `op(a, b)` when the two
# operands are *independent* and a closed form is known, and `nothing`
# otherwise. Callers must check independence themselves (see `binop`) — this
# table assumes it.

closure(op, a, b) = nothing

# Point masses behave like the number they sit on.
closure(op, a::Dirac, b) = closure(op, a.value, b)
closure(op, a, b::Dirac) = closure(op, a, b.value)
closure(op, a::Dirac, b::Dirac) = Dirac(op(a.value, b.value))

closure(::typeof(+), a::Normal, b::Normal) = Normal(a.μ + b.μ, hypot(a.σ, b.σ))
closure(::typeof(-), a::Normal, b::Normal) = Normal(a.μ - b.μ, hypot(a.σ, b.σ))
closure(::typeof(+), a::Normal, c::Real) = Normal(a.μ + c, a.σ)
closure(::typeof(+), c::Real, a::Normal) = Normal(c + a.μ, a.σ)
closure(::typeof(-), a::Normal, c::Real) = Normal(a.μ - c, a.σ)
closure(::typeof(-), c::Real, a::Normal) = Normal(c - a.μ, a.σ)
closure(::typeof(*), a::Normal, c::Real) = Normal(a.μ * c, abs(c) * a.σ)
closure(::typeof(*), c::Real, a::Normal) = Normal(c * a.μ, abs(c) * a.σ)
closure(::typeof(/), a::Normal, c::Real) = iszero(c) ? nothing : Normal(a.μ / c, a.σ / abs(c))

closure(::typeof(*), a::LogNormal, b::LogNormal) = LogNormal(a.μ + b.μ, hypot(a.σ, b.σ))
closure(::typeof(/), a::LogNormal, b::LogNormal) = LogNormal(a.μ - b.μ, hypot(a.σ, b.σ))
closure(::typeof(*), a::LogNormal, c::Real) = c > 0 ? LogNormal(a.μ + log(c), a.σ) : nothing
closure(::typeof(*), c::Real, a::LogNormal) = c > 0 ? LogNormal(log(c) + a.μ, a.σ) : nothing
closure(::typeof(/), a::LogNormal, c::Real) = c > 0 ? LogNormal(a.μ - log(c), a.σ) : nothing
closure(::typeof(/), c::Real, a::LogNormal) = c > 0 ? LogNormal(log(c) - a.μ, a.σ) : nothing
closure(::typeof(^), a::LogNormal, c::Real) = iszero(c) ? nothing : LogNormal(a.μ * c, abs(c) * a.σ)

closure(::typeof(+), a::Cauchy, b::Cauchy) = Cauchy(a.μ + b.μ, a.σ + b.σ)
closure(::typeof(+), a::Cauchy, c::Real) = Cauchy(a.μ + c, a.σ)
closure(::typeof(+), c::Real, a::Cauchy) = Cauchy(c + a.μ, a.σ)
closure(::typeof(*), a::Cauchy, c::Real) = iszero(c) ? nothing : Cauchy(a.μ * c, abs(c) * a.σ)
closure(::typeof(*), c::Real, a::Cauchy) = iszero(c) ? nothing : Cauchy(c * a.μ, abs(c) * a.σ)

closure(::typeof(+), a::Gamma, b::Gamma) = a.θ ≈ b.θ ? Gamma(a.α + b.α, a.θ) : nothing
closure(::typeof(*), a::Gamma, c::Real) = c > 0 ? Gamma(a.α, a.θ * c) : nothing
closure(::typeof(*), c::Real, a::Gamma) = c > 0 ? Gamma(a.α, c * a.θ) : nothing
closure(::typeof(+), a::Poisson, b::Poisson) = Poisson(a.λ + b.λ)
closure(::typeof(+), a::Binomial, b::Binomial) = a.p ≈ b.p ? Binomial(a.n + b.n, a.p) : nothing

# Unary closures.
uclosure(f, d) = nothing
uclosure(::typeof(-), d::Normal) = Normal(-d.μ, d.σ)
uclosure(::typeof(-), d::Cauchy) = Cauchy(-d.μ, d.σ)
uclosure(::typeof(-), d::Dirac) = Dirac(-d.value)
uclosure(::typeof(exp), d::Normal) = LogNormal(d.μ, d.σ)
uclosure(::typeof(log), d::LogNormal) = Normal(d.μ, d.σ)
uclosure(::typeof(inv), d::LogNormal) = LogNormal(-d.μ, d.σ)
uclosure(::typeof(sqrt), d::LogNormal) = LogNormal(d.μ / 2, d.σ / 2)

# The provenance gate: take the closed form only when the operands share no
# leaf. Either way the thunk recurses through the tree, so the particle
# representation stays correctly correlated.
function binop(op, a::Node, b::Node)
	d = isdisjoint(a.deps, b.deps) ? closure(op, a.dist, b.dist) : nothing
	derive(d, union(a.deps, b.deps), () -> op(parts(a), parts(b)))
end

binop(op, a::Node, c::Real) = derive(closure(op, a.dist, c), copy(a.deps), () -> op(parts(a), c))
binop(op, c::Real, a::Node) = derive(closure(op, c, a.dist), copy(a.deps), () -> op(c, parts(a)))

# Raw particles have unknown provenance, so no closed form survives contact
# with them; the thunk still keeps the numerics right.
binop(op, a::Node, p::AbstractParticles) = derive(nothing, copy(a.deps), () -> op(parts(a), p))
binop(op, p::AbstractParticles, a::Node) = derive(nothing, copy(a.deps), () -> op(p, parts(a)))

unop(f, a::Node) = derive(uclosure(f, a.dist), copy(a.deps), () -> f(parts(a)))

for op in (:+, :-, :*, :/, :^)
	@eval begin
		Base.$op(a::Node, b::Node) = binop($op, a, b)
		Base.$op(a::Node, b::Real) = binop($op, a, b)
		Base.$op(a::Real, b::Node) = binop($op, a, b)
		Base.$op(a::Node, b::AbstractParticles) = binop($op, a, b)
		Base.$op(a::AbstractParticles, b::Node) = binop($op, a, b)
	end
end

# Julia lowers `x^2` to Base.literal_pow; without this it would bypass the
# closure table entirely.
Base.literal_pow(::typeof(^), a::Node, ::Val{n}) where {n} = binop(^, a, n)

Base.:-(a::Node) = unop(-, a)
Base.:+(a::Node) = a

const UNARY = (:exp, :log, :log2, :log10, :log1p, :exp2, :exp10, :expm1,
	:sqrt, :cbrt, :abs, :inv, :sin, :cos, :tan, :asin, :acos, :atan,
	:sinh, :cosh, :tanh, :sign, :floor, :ceil, :round)

for f in UNARY
	@eval Base.$f(a::Node) = unop($f, a)
end

# Comparisons are left to MonteCarloMeasurements, whose error message explains
# why comparing uncertain values is ill-defined and what to do instead.
for op in (:<, :<=, :>, :>=, :(==))
	@eval begin
		Base.$op(a::Node, b::Node) = $op(parts(a), parts(b))
		Base.$op(a::Node, b::Real) = $op(parts(a), b)
		Base.$op(a::Real, b::Node) = $op(a, parts(b))
	end
end

Base.promote_rule(::Type{Node}, ::Type{<:Real}) = Node
Base.convert(::Type{Node}, x::Real) = pointmass(x)
Base.convert(::Type{Node}, n::Node) = n
