# Gekrakel.jl

Squiggle's estimation vocabulary, as a Julia library rather than a language.

```julia
using Gekrakel

pop       = 2.6e6
per_house = 2…4             # 90% credible interval
pianos    = 0.002…0.02
tunings   = 0.5…2
per_day   = 3…6
workdays  = 200…250

tuners = pop / per_house * pianos * tunings / (per_day * workdays)
```

```
LogNormal (exact)
  ▇█▇▆▄▄▃▃▂▂▂▂▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁
  5% 1.46   50% 6.128   95% 25.73   mean 8.965
```

That chain is entirely closed-form; nothing was sampled.

Squiggle's limitation is being a language: a Turing.jl posterior, a CSV, a
Metaculus fetch, or a search over inputs all live outside it. Here
distributions are ordinary `Real`s, so they flow through any Julia code.
[MonteCarloMeasurements.jl](https://github.com/baggepinnen/MonteCarloMeasurements.jl)
does the sampling; Gekrakel adds credible-interval constructors, closed forms,
mixtures, plotting, and Squiggle's display conventions.

## Correlation and closed forms

Shared variables stay shared:

```julia
x = 10…100
x - x        # exactly 0
x / x        # exactly 1
```

Closed forms are used only when the operands share no source of randomness;
otherwise Gekrakel falls back to correlated samples. (Squiggle gets this wrong
and reports `x - x` for a symbolic normal as `normal(0, √2)`.)

```julia
a, b = normal(0, 1), normal(0, 1)
exact(a + b)    # true:  Normal(0, √2)
exact(a + a)    # false: sampled, and std(a + a) == 2
```

## Writing models

**Precedence.** `…` binds *looser* than `+` and `*`, just like Squiggle's `to`:

```julia
10…100 * 2      # means 10…200
(10…100) * 2    # what you probably wanted
to(10, 100)     # plain function, no precedence trouble
```

**Branches.** A distribution can't answer `x > 0`, so comparisons raise. Use
`over`, which applies a function sample by sample and keeps correlations:

```julia
income = 50_000…300_000
tax = over(x -> x > 1e5 ? 0.45x : 0.25x, income)
```

**Truncate vs. clamp.** `truncate(x, lo, hi)` conditions on the range and
renormalises, producing a new, uncorrelated distribution. `clamp(x, lo, hi)`
is pointwise, piles mass up at the bounds, and stays correlated with `x`.

**Outside samples.** `gek` wraps anything:

```julia
gek(Normal(0, 1))    # a Distributions.jl distribution
gek(chain[:θ])       # a Turing posterior, or any sample vector
```

Sample vectors must have `nsamp()` elements (default 2000); call
`setnsamp!(n)` before building anything else to change that.

## Plotting

In Pluto, VS Code or IJulia every distribution renders itself as a density
plot, with no `plot` call. For functions of one or two variables:

```julia
distfn(h -> value(h, 90), 0.5, 20; xlab="hr/day")    # ℝ → distribution: quantile bands
sweep((h, d) -> value(h, d), (0.5, 20), (5, 90))     # ℝ² → distribution: heatmap
```

`sweep` also outlines the inputs whose outputs nobody stochastically
dominates. Build the distributions *outside* the function you pass, so every
cell shares the same randomness; `sweep` warns if they don't. A Plots.jl
recipe is included as well.

## Comparing options

The dominance orders live in
[StochasticOrders.jl](../StochasticOrders.jl), which works on any sample
vectors; Gekrakel re-exports them for distributions.

```julia
fosd(x, y)          # does x first-order dominate y?
pgeq(x, y)          # P(x ≥ y), using the shared randomness
undominated(xs)     # options nothing dominates
margin(x, y)        # why dominance failed, and by how much
robust(xs)          # options that survive a sweep of risk aversions
certify(x, y)       # decide from closed forms, or `nothing`
```

Expect dominance to fail almost everywhere: options that differ in risk
generally can't dominate one another. `margin`, `regret` and `robust` are the
tools that discriminate in practice.

## Caveats

**Prefer medians for heavy tails.** Samples sit on a quantile grid, so
quantiles are very accurate but the *mean* of a tail-heavy distribution is
biased low, and consistently so: `lognormal(0, 4)` has mean 2981, but its
2000 samples average 1259. `exact(x)` tells you whether a mean was computed
analytically or estimated from samples.

**Name collisions.** `beta` and `gamma` collide with SpecialFunctions; use
`import Gekrakel as G` if you need both.

## Reference

| | |
|---|---|
| `lo…hi`, `to(lo, hi)` | 90% CI; lognormal if `lo > 0`, else normal |
| `normal`, `lognormal` | `(μ, σ)` or `(p5=, p95=)`; lognormal takes log-space parameters |
| `uniform`, `cauchy`, `gamma`, `exponential`, `beta` | `exponential` takes a rate; `beta(mean=, stdev=)` works too |
| `bernoulli`, `poisson`, `triangular`, `pointmass` | |
| `mixture(cs, ws)`, `mx` | numbers become point masses |
| `mean`, `var`, `std`, `median`, `quantile`, `cdf`, `pdf` | exact when a closed form survives |
| `credible(x, mass)`, `inv_cdf`, `stdev` | |
| `over(f, xs...)` | sample-wise application |
| `truncate`, `truncateLeft`, `truncateRight`, `clamp` | |
| `gek(x)`, `samples(x)`, `rand(x, k)`, `exact(x)` | |
| `spark(x)`, `svgdens(x)`, `distfn`, `sweep` | display |
| `nsamp()`, `setnsamp!(n)` | particle count |
| `fosd`, `sosd`, `statewise`, `pgeq`, `margin`, `undominated`, `robust`, … | see StochasticOrders.jl |
| `certify(x, y)` | dominance from closed forms |

## Installation

Unregistered, and depends on the sibling package StochasticOrders.jl:

```julia
using Pkg
Pkg.develop(path="../StochasticOrders.jl")
Pkg.develop(path="../Gekrakel.jl")
Pkg.test("Gekrakel")
```
