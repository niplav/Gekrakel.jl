using Gekrakel
using MonteCarloMeasurements
using Distributions
using Statistics
using Test

@testset "Gekrakel" begin

@testset "credible-interval constructors" begin
	x = 10…100
	@test quantile(x, [0.05, 0.95]) ≈ [10, 100] rtol = 1e-9
	@test x.dist isa LogNormal
	@test pquantile(parts(x), [0.05, 0.95]) ≈ [10, 100] rtol = 0.01
	y = to(-5, 5)
	@test y.dist isa Normal
	@test quantile(y, [0.05, 0.95]) ≈ [-5, 5] rtol = 1e-9
	@test_throws ErrorException to(10, 1)
	@test quantile(normal(p5=-1, p95=1), 0.5) ≈ 0 atol = 1e-12
	@test mean(beta(mean=0.3, stdev=0.1)) ≈ 0.3
end

@testset "correlation is preserved" begin
	x = 10…100
	@test all(iszero, samples(x - x))
	@test all(≈(1), samples(x / x))
	@test cor(samples(x), samples(2x)) ≈ 1
	@test pstd(parts(x + x - 2x)) ≈ 0 atol = 1e-9
end

@testset "provenance gates the closed form" begin
	a = normal(0, 1)
	b = normal(0, 1)
	# Independent: closed form applies.
	@test exact(a + b)
	@test (a + b).dist.σ ≈ sqrt(2)
	# Correlated: it must not. Squiggle gets this case wrong.
	@test !exact(a + a)
	@test pstd(parts(a + a)) ≈ 2 rtol = 0.01
	@test pstd(parts(a - a)) ≈ 0 atol = 1e-12
end

@testset "materialisation walks the tree" begin
	a = normal(0, 1)
	b = normal(0, 1)
	z = a + b
	# If z resampled its closed form instead of recursing, this would be √3.
	@test pstd(parts(z - a)) ≈ 1 rtol = 0.01
	@test pstd(parts(z)) ≈ sqrt(2) rtol = 0.05
	@test samples(z) ≈ samples(a) .+ samples(b)
end

@testset "analytic closures" begin
	# Bit-exact means prove no sampling happened.
	@test mean(normal(0, 1) + normal(0, 1)) === 0.0
	@test ((10…100) * (2…5)).dist isa LogNormal
	@test ((10…100) / (2…5)).dist isa LogNormal
	@test (2 * normal(1, 3)).dist == Normal(2, 6)
	@test (normal(1, 3) - 1).dist == Normal(0, 3)
	@test exp(normal(0, 1)).dist == LogNormal(0, 1)
	@test log(lognormal(0, 1)).dist == Normal(0, 1)
	@test (cauchy(0, 1) + cauchy(0, 2)).dist == Cauchy(0, 3)
	@test (poisson(2) + poisson(3)).dist == Poisson(5)
	@test (gamma(2, 3) + gamma(4, 3)).dist == Gamma(6, 3)
	@test !exact(gamma(2, 3) + gamma(4, 7))   # scales differ, no closed form
	@test ((10…100)^2).dist isa LogNormal
	@test mean(pointmass(3) + pointmass(4)) == 7
end

@testset "particle fallback" begin
	s = sin(10…100)
	@test !exact(s)
	@test length(samples(s)) == nsamp()
	@test all(v -> -1 <= v <= 1, samples(s))
	@test !exact(normal(0, 1) * normal(0, 1))     # no closed form for this
	@test pmean(parts(normal(5, 1) * normal(3, 1))) ≈ 15 rtol = 0.05
end

@testset "mixtures" begin
	m = mixture([normal(0, 1), normal(10, 1)], [0.3, 0.7])
	@test mean(m) ≈ 7
	@test exact(m)
	@test pmean(parts(m)) ≈ 7 rtol = 0.05
	@test mean(mixture([0, 1], [0.5, 0.5])) ≈ 0.5
	# A component sharing a leaf with another kills the closed form.
	x = normal(0, 1)
	@test !exact(mixture([x, x + 1]))
	@test_throws ErrorException mixture([normal(0, 1)], [1.0, 2.0])
end

@testset "gek and over" begin
	@test gek(Normal(0, 1)).dist == Normal(0, 1)
	@test !exact(gek(randn(nsamp())))
	@test length(samples(gek(randn(nsamp())))) == nsamp()
	@test_throws ErrorException gek(randn(nsamp() + 1))
	x = 10…100
	@test all(iszero, samples(x - over(identity, x)))
	@test pmean(parts(over(+, normal(3, 1), normal(4, 1)))) ≈ 7 rtol = 0.05
	@test all(v -> v in (0.0, 1.0), samples(over(v -> v > 0 ? 1.0 : 0.0, normal(0, 1))))
	@test_throws ErrorException over(identity)
end

@testset "mixture of a distribution with a plain number" begin
	# `pointmass` is a Dirac, which Distributions.jl calls discrete, so
	# MixtureModel refuses to combine it with a continuous component. Used to
	# throw; must fall through to particles instead.
	m = mx([normal(0, 1), 5.0], [0.9, 0.1])
	@test !exact(m)
	@test pmean(parts(m)) ≈ 0.5 rtol = 0.25
	@test count(==(5.0), samples(m)) / nsamp() ≈ 0.1 atol = 0.03
	@test exact(mx([poisson(3), 5.0], [0.9, 0.1]))          # both discrete
	@test exact(mx([normal(0, 1), normal(5, 1)], [0.9, 0.1]))  # both continuous
end

@testset "over mixes nodes and numbers" begin
	@test pmean(parts(over((a, b) -> a * b, normal(3, 1), 2.0))) ≈ 6 rtol = 0.05
	@test pmean(parts(over((a, b) -> a * b, 2.0, normal(3, 1)))) ≈ 6 rtol = 0.05
	@test_throws ErrorException over(identity, 1.0)
end

@testset "truncate conditions, clamp does not" begin
	x = normal(0, 1)
	t = truncateRight(x, 1)
	@test exact(t)
	@test quantile(t, 0.999) <= 1
	@test mean(t) < 0                       # the upper tail is gone
	@test maximum(samples(truncateLeft(0…100, 1))) > 1
	@test minimum(samples(truncateLeft(0…100, 1))) >= 1
	@test length(samples(truncate(sin(10…100), -0.5, 0.5))) == nsamp()
	@test_throws ErrorException truncate(x, 1, 0)

	# clamp keeps the correlation truncate destroys
	c = clamp(x, -0.5, 0.5)
	@test all(v -> -0.5 <= v <= 0.5, samples(c))
	@test all(iszero, samples(clamp(x, -Inf, Inf) - x))
	@test all(iszero, samples(max(x, -Inf) - x))
	@test pmean(parts(max(normal(0, 1), normal(0, 1)))) > 0.4
end

@testset "summaries" begin
	x = 10…100
	@test cdf(x, median(x)) ≈ 0.5 rtol = 1e-9
	@test inv_cdf(x, 0.5) ≈ median(x)
	@test credible(x) == (quantile(x, 0.05), quantile(x, 0.95))
	@test std(x) ≈ std(LogNormal(x.dist.μ, x.dist.σ))
	@test stdev(x) == std(x)
	@test cdf(sin(10…100), 1.0) == 1.0
	@test pdf(normal(0, 1), 0.0) ≈ 1 / sqrt(2pi)
	@test pdf(sin(10…100), 0.0) > 0
	@test length(rand(x, 7)) == 7
end

@testset "precedence gotcha is real" begin
	# `…` binds looser than `*`, exactly as Squiggle's `to` does. Pinned by a
	# test so the documentation cannot silently rot.
	@test (10…100 * 2).dist == (10…200).dist
	@test ((10…100) * 2).dist != (10…200).dist
end

@testset "display" begin
	@test occursin("…", sprint(show, 10…100))
	@test occursin("5%", sprint(show, MIME"text/plain"(), 10…100))
	@test occursin("<div", sprint(show, MIME"text/html"(), 10…100))
	@test length(spark(normal(0, 1))) == 28
	@test length(spark(sin(10…100))) == 28
	@test spark(pointmass(1)) isa String
end

@testset "svg" begin
	h = sprint(show, MIME"text/html"(), 10…100)
	@test occursin("<svg", h) && occursin("polyline", h)
	@test occursin("polygon", h)                       # shaded 90% band
	@test occursin("<line", svgdens(poisson(4)))       # discrete gets sticks
	@test occursin("<svg", svgdens(pointmass(3)))      # degenerate support
	@test occursin("<svg", svgdens(sin(10…100)))       # sampled, via the KDE
	# labels must be inside the viewBox, not relying on overflow:visible
	m = match(r"viewBox='0 0 (\d+) (\d+)'", svgdens(normal(0, 1); w = 340, h = 78))
	@test parse(Int, m[2]) > 78
end

@testset "distfn" begin
	p = distfn(a -> normal(a, 1), 0.0, 10.0; n = 12, xlab = "a")
	@test length(p.xs) == 12
	@test issorted(p.med)                              # median tracks the mean
	@test all(p.lo[1] .<= p.lo[2]) && all(p.hi[2] .<= p.hi[1])  # 90% ⊇ 50%
	@test p.mass == [0.9, 0.5]                         # widest first
	@test occursin("<svg", sprint((io, v) -> show(io, MIME"text/html"(), v), p))
	@test occursin("a ∈", sprint(show, MIME"text/plain"(), p))
	@test distfn(a -> 2a, 0.0, 1.0; n = 4).med ≈ [0.0, 2/3, 4/3, 2.0]
	# the zero line is drawn only when zero is in range
	crossing = distfn(a -> normal(a - 5, 1), 0.0, 10.0; n = 8)
	@test occursin("stroke-dasharray='3,3'",
		sprint((io, v) -> show(io, MIME"text/html"(), v), crossing))
	@test_throws ErrorException distfn(identity, 1.0, 0.0)
end

@testset "dominance adapter" begin
	# The orders themselves are StochasticDominance.jl's business and are
	# tested there. What belongs here is that `Node`s reach them correctly and
	# that the three Gekrakel-specific parts work: unwrapping, the
	# integrability guard, and returning `Node`s from the lattice ops.
	@test fosd(normal(1, 1), normal(0, 1))
	@test !fosd(normal(0, 2), normal(0, 1))
	@test relate(normal(1, 1), normal(0, 1)) === :dominates
	@test margin(normal(0, 2), normal(0, 1)).v1 > 0
	@test wass(normal(1, 1), normal(0, 1)) ≈ 1 rtol = 1e-6
	@test 0 <= pgeq(normal(0, 1), normal(0, 1)) <= 1

	# the coupling, which is Gekrakel's own contribution: shared leaves make
	# the test enormously sharper than independent leaves with equal marginals
	z = normal(0, 1) + normal(0, 1) + normal(0, 1)
	@test fosd(z + 0.01, z)
	@test !fosd(normal(0, 1) + normal(0, 1) + normal(0, 1) + 0.01, z)
	@test statewise(z + 1, z) && pgeq(z + 1, z) == 1.0

	# the guard needs a closed form to detect, so it lives on this side
	@test_throws ErrorException sosd(cauchy(0, 1), cauchy(1, 1))
	@test fosd(cauchy(1, 1), cauchy(0, 1))       # first order is mean-free
	@test sosd(normal(0, 1), normal(0, 2))
	@test icx(normal(0, 2), normal(0, 1))

	# the lattice ops must hand back `Node`s, not sorted vectors
	a, b = 1…10, 2…20
	@test (a ⊔ b) isa Node && (a ⊓ b) isa Node && envelope([a, b]) isa Node
	@test length(samples(a ⊔ b)) == nsamp()
	@test fosd(a ⊔ b, a) && fosd(a ⊔ b, b)
	@test fosd(a, a ⊓ b) && fosd(b, a ⊓ b)

	# and the collection functions must accept a Vector{Node}
	ch = [normal(i, 1) for i in 0:4]
	@test undominated(ch) == [5] && dominating(ch) == 5
	@test domlevels(ch) == [4, 3, 2, 1, 0]
	@test robust(ch) == [5]
	@test length(regret(ch)) == 5 && size(margins(ch)) == (5, 5)
end

@testset "certified dominance" begin
	@test certify(normal(1, 1), normal(0, 1)) === true
	@test certify(normal(0, 2), normal(0, 1)) === false
	# 1…10 and 2…20 are scale multiples, so their log-scales are *equal* and
	# dominance is decidable either way
	@test certify(2…20, 1…10) === true
	@test certify(1…10, 2…20) === false
	# unequal scale forces a crossing, so `false` is a *proof* of
	# non-dominance rather than an undecided verdict
	@test certify(1…10, 1…100) === false
	@test certify(normal(0, 1), 1…10) === nothing   # different families
	@test certify(gek(LogNormal(1, 2)), gek(LogNormal(0, 2))) === true
	@test certify(poisson(4), poisson(3)) === true
	@test certify(pointmass(3), pointmass(2)) === true
	# n>=n and p>=p is sufficient but not necessary; the support is finite
	@test certify(gek(Binomial(100, 0.5)), gek(Binomial(1, 0.9))) === true
	@test certify(sin(10…100), normal(0, 1)) === nothing    # not exact
	# where it decides, it agrees with the particle test
	for (x, y) in ((normal(1, 1), normal(0, 1)), (normal(0, 2), normal(0, 1)),
			(poisson(4), poisson(3)))
		c = certify(x, y)
		c === nothing || @test c == fosd(x, y)
	end
end

@testset "sweep" begin
	base = normal(0, 1)                       # one leaf, shared by every cell
	g = sweep((x, y) -> base * x + y, (1.0, 4.0), (0.0, 3.0); n=(5, 4),
		xlab="a", ylab="b", slab="med")
	@test size(g.nodes) == (5, 4) && size(g.vals) == (5, 4)
	@test length(g.xs) == 5 && length(g.ys) == 4
	@test !isempty(g.front)
	@test g.vals[1, 4] ≈ 3 rtol = 0.2         # median of base*1 + 3

	h = sprint((io, v) -> show(io, MIME"text/html"(), v), g)
	@test occursin("<svg", h) && occursin("<style>", h)
	@test occursin("prefers-color-scheme:dark", h)   # both modes declared
	@test occursin("var(--gk", h)
	@test occursin("<title>", h)                     # per-cell hover
	@test occursin("rotate(-90)", h)                 # y-axis label
	@test occursin("b ∈", sprint(show, MIME"text/plain"(), g))
	@test_throws ErrorException sweep((x, y) -> base, (4.0, 1.0), (0.0, 3.0))

	# cells built from their own fresh leaves are mutually independent, which
	# makes the coupling and the frontier meaningless — must warn
	@test_logs (:warn,) match_mode=:any sweep((x, y) -> normal(x, 1) + y,
		(1.0, 2.0), (0.0, 1.0); n=(2, 2))
end

@testset "comparisons refuse to guess" begin
	@test_throws Exception normal(0, 1) > 0
end

end
