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

export Node, Gek,
	normal, lognormal, uniform, cauchy, gamma, exponential, bernoulli,
	poisson, triangular, beta, pointmass, to, …,
	mixture, mx, gek, over, truncateLeft, truncateRight,
	mean, var, std, stdev, median, quantile, cdf, pdf, inv_cdf,
	samples, parts, credible, spark, svgdens, distfn, DistFn, exact,
	nsamp, setnsamp!

include("node.jl")
include("ctor.jl")
include("algebra.jl")
include("trunc.jl")
include("summary.jl")
include("mixture.jl")
include("show.jl")
include("svg.jl")
include("plots.jl")

end # module
