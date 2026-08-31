"""
Gekrakel — Squiggle's estimation vocabulary as a Julia library.

	using Gekrakel
	revenue = 10…100          # 90% credible interval, lognormal
	cost    = normal(5, 2)
	profit  = revenue - cost

Distributions are ordinary numbers here, so they flow through any Julia code.
Shared variables stay correlated (`x - x` is exactly zero), and closed forms
are used wherever they exist and the operands are independent.
"""
module Gekrakel

using MonteCarloMeasurements
using Distributions
using RecipesBase
using Random: shuffle!

import Distributions: mean, var, std, median, quantile, cdf, pdf

# Gekrakel's own vocabulary.
export Node, Gek,
	normal, lognormal, uniform, cauchy, gamma, exponential, bernoulli,
	poisson, triangular, beta, pointmass, to, …,
	mixture, mx, gek, over, truncateLeft, truncateRight,
	mean, var, std, stdev, median, quantile, cdf, pdf, inv_cdf,
	samples, parts, credible, spark, svgdens, distfn, DistFn, exact,
	nsamp, setnsamp!,
	certify, sweep, Sweep

# Re-exported from StochasticOrders, which owns the orders themselves. Gekrakel
# only adds `Node` methods to them (see `dominance.jl`), so a user who has
# `using Gekrakel` needs these names to reach those methods at all. Every one
# is explicitly imported there rather than picked up from `using`, so a rename
# upstream fails loudly here instead of quietly changing what this exports.
export relate, fosd, sosd, icx, statewise, pgeq, wass, margin, DomMargin,
	lub, glb, ⊔, ⊓, envelope, regret,
	domgraph, hasse, domlevels, undominated, dominating, margins, robust

include("node.jl")
include("ctor.jl")
include("algebra.jl")
include("trunc.jl")
include("summary.jl")
include("mixture.jl")
include("show.jl")
include("svg.jl")
include("dominance.jl")
include("sweep.jl")
include("plots.jl")

end # module
