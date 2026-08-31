# Gekrakel.jl

Squiggle's estimation vocabulary, as a Julia library rather than a language.

```julia
using Gekrakel

pop        = 2.6e6
per_house  = 2…4               # 90% credible interval
pianos     = 0.002…0.02
tunings    = 0.5…2
per_day    = 3…6
workdays   = 200…250

households = pop / per_house
tuners     = households * pianos * tunings / (per_day * workdays)
```

```
LogNormal (exact)
  ▇█▇▆▄▄▃▃▂▂▂▂▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁
  5% 1.46   50% 6.128   95% 25.73   mean 8.965
```

That whole chain is closed-form; nothing was sampled.

## Why a library

Squiggle's defect is being a language: a Metaculus fetch, a Turing.jl
posterior, a CSV, or a search for stochastically dominating inputs all live
outside it. Here distributions are ordinary `Real`s, so they flow through any
Julia code you already have. `MonteCarloMeasurements.jl` does the heavy
lifting; Gekrakel adds credible-interval constructors, closed forms, mixtures,
and Squiggle's display conventions.

## Two things it gets right

**Correlation.** Shared variables stay shared.

```julia
x = 10…100
x - x        # exactly 0, not a distribution around 0
x / x        # exactly 1
```

**Provenance.** Closed forms are applied only between operands that share no
source of randomness. Squiggle gets this wrong: given symbolic
`x = normal(0,1)`, it reports `x - x` as `normal(0, √2)`.

```julia
a, b = normal(0, 1), normal(0, 1)
(a + b).dist          # Normal(0, 1.414...) — exact, independent
exact(a + a)          # false — correlated, so it falls back to sampling
std(a + a)            # 2.0, not 1.414
```

Materialisation walks the expression tree rather than resampling the closed
form, so `z = a + b; z - a` still recovers `a`'s partner exactly.

## Fat tails: trust the median, and check `exact`

MonteCarloMeasurements samples on a quantile grid rather than iid, which makes
marginals far more accurate than iid at the same N — `quantile` is good to ~6
digits. But the grid puts exactly one point in the outermost stratum, so when
a distribution's *mean* is tail-dominated, the particle mean is biased low.
For `lognormal(0, 4)`, true mean 2981:

```julia
mean(lognormal(0, 4))                 # 2981   — closed form
mean(over(identity, lognormal(0, 4))) # 1259   — particles, same distribution
```

Worse, the grid is deterministic, so the error is *quiet*: repeat it six times
and you get 1259 every time. Iid sampling would at least be loudly noisy.
Convergence is slow — at σ=3 the ratio to the true mean goes 0.68 → 0.95 as N
goes 500 → 50000.

A `mixture` adds an iid selector on top, which does not remove the bias but
does add variance, and the variance can be large when the tail lives in a
low-weight arm: `mx([2e3…4e4, 3e4…1e6], [0.95, 0.05])` has an exact mean of
28147, and its particle mean moves over a 20% range between rebuilds.

Two consequences, and the package hands you the tool for both:

- Prefer the median, or `1 - cdf(x, 0)`, over the mean as a decision statistic
  for anything heavy-tailed. Quantiles are what the grid is good at.
- `exact(x)` tells you whether the mean you just read was computed or
  estimated. If a chain of closed forms survived, `mean` is analytic and none
  of this applies. Once anything breaks the chain — a branch through `over`, a
  correlated operand, `sin` — the mean quietly becomes an estimate.

## Stochastic dominance

The thing a language can't do: hand a function to the library and ask which
inputs give outputs that dominate.

The orders live in their own package,
[StochasticOrders.jl](../StochasticOrders.jl), because they need nothing
but two index-aligned sample vectors — they work equally well on an MCMC chain
or a CSV column. Gekrakel adds the three parts that do need a `Node`:
unwrapping, the integrability guard (detectable only from a closed form), and
`certify`.

```julia
fosd(x, y)          # does x first-order dominate y (weakly)?
undominated(xs)     # indices nothing dominates — the maximal antichain
dominating(xs)      # the greatest element, or `nothing`
margin(x, y)        # why it failed, and by how much
```

For equal-N empirical measures this needs no grid and no KDE. `F̂_X ≤ F̂_Y`
pointwise iff the quantile functions reverse, and both are step functions with
jumps at the same `i/N`, so the whole thing collapses to `sort(xs) .>= sort(ys)`
— O(N log N), and exact on the empirical measures. Everything else is the same
vector under a different reduction: `mean(abs.(d))` is Wasserstein-1,
`max.(sort, sort)` is the lattice join (`⊔`), `cumsum` gives second-order
dominance.

**Three orders, and they answer different questions.**

| | |
|---|---|
| `fosd` | better for every increasing utility. Marginal. |
| `sosd` | better for every increasing *concave* utility — every risk-averse preference |
| `icx` | the risk-seeking counterpart |
| `statewise`, `pgeq` | better in every *world*, using the coupling |

`statewise` is the one this library can ask and Squiggle cannot. `fosd`
compares marginals, so sorting throws the index alignment away; `statewise`
keeps it, and the index *is* the model's joint structure — same draw of every
shared input. Prefer `pgeq`, which returns `P(x ≥ y)`: a probability-one claim
isn't estimable from N draws, so a true 0.9995 reads as `false` and tells you
nothing.

The coupling matters more than it looks. Shared leaves make the two sorted
vectors positively dependent, so `x = y + δ` passes at *one* particle; built
from independent leaves with the same marginals it fails unless δ clears the
noise in the extreme order statistics. Which is why `sweep` insists you build
the leaves **once, outside** the function.

**Expect nothing to dominate.** Dominance requires the extreme order statistics
to be ordered, and between two location-scale distributions it requires *equal
scale* — it is essentially a same-shape-shifted relation. Any two options
differing in spread are incomparable after two comparisons, so in a model whose
point is that options differ in risk, the frontier is nearly everything. That
is information. `undominated` is never empty and needs no sentinel: "nothing
dominates" and "everything is undominated" are the same sentence.

So the deliverables that actually discriminate are `margin`, `regret`, and
`robust` — the last sweeping risk aversion over a CARA family and keeping what
survives every setting.

`certify` decides from closed forms where it can, and returns `nothing` rather
than guessing. A `true` there is a proof, not an estimate.

## Functions of two variables

```julia
sweep(f, (0.5, 20), (5, 90); n=(16,16), stat=median, xlab="hr/day", ylab="days")
```

A heatmap of one summary per cell with the undominated cells outlined — two
input dimensions plus a distribution is four-ish dimensions and no single glyph
survives that. Diverging blue↔red on a neutral grey midpoint, symmetric about
zero (a diverging scale has to be, or the midpoint stops meaning "nothing"), the
red arm derived by rotating the blue steps in OKLCH so the arms match in
lightness. Both themes are selected rather than flipped, and the poles clear
CVD separation in both.

## Precedence gotcha

`…` sits in Julia's precedence group 10, **looser** than `+` (11) and `*` (12):

```julia
10…100 * 2       # parses as 10…200
(10…100) * 2     # what you meant
```

Squiggle's `to` behaves identically. Parenthesise whenever you combine.
`to(10, 100)` is a plain function and has no such problem.

## Branches

A distribution cannot answer `x > 0`, so comparisons raise. Use `over`, which
applies a function sample-by-sample while keeping correlations intact:

```julia
income = 50_000…300_000
tax = over(x -> x > 1e5 ? 0.45x : 0.25x, income)
```

## Truncating vs. clamping

Easy to conflate, and the difference is the whole answer:

```julia
truncate(x, 0.05, Inf)   # condition on x ≥ 0.05; mass below is removed and
                         # the rest renormalised. A new source of randomness,
                         # so it is no longer correlated with x.
clamp(x, 0.05, Inf)      # pointwise max; mass below piles up at 0.05, and
                         # correlation with x survives.
```

`truncate` is exact when a closed form survives, via `Distributions.truncated`.

## Bringing in outside samples

```julia
gek(Normal(0, 1))         # any Distributions.jl distribution
gek(chain[:θ])            # a Turing posterior
gek(fetch_metaculus(...)) # anything you have samples of
```

Sample vectors arrive without provenance, so they get a fresh leaf and no
closed form. They must have `nsamp()` elements — `Particles{T,N}` puts the
count in the type, so mixing counts is a type error; call `setnsamp!(n)` first.

## Notebooks

Pluto, VS Code, IJulia and Documenter all display a value by asking it for
`MIME"text/html"`, and `Node` answers with an inline SVG density — shaded
across the 90% interval, dashed at the median. So every named distribution in
a notebook draws itself, with no `plot` call. That is the one thing Squiggle's
playground does that a REPL does not.

For a function from a number to a distribution — Squiggle's `Plot.distFn` —
there is `distfn`, which draws nested quantile bands around the median:

```julia
distfn(h -> value(h, 90).net, 0.5, 20; xlab="hr/day")
```

A zero line appears whenever zero is in range, since that is the decision
boundary whenever the y-axis is money. `f` is evaluated once per point, so
keep `n` small while iterating. Band edges are jagged where the model is
sample-based: neighbouring x are independent draws, and that jaggedness is an
honest picture of the Monte Carlo error rather than something to smooth away.

Everything is drawn in `currentColor`, so it inherits the host page's
foreground and needs no light/dark handling.

Pluto gotcha, unrelated to this package: a name may only be assigned in one
cell, which fights the `x = ...; x = truncate(x, ...)` style. Use `let` blocks.

## Reference

| | |
|---|---|
| `lo…hi`, `to(lo,hi)` | 90% CI; lognormal if `lo > 0`, else normal |
| `normal(μ,σ)`, `normal(p5=,p95=)` | |
| `lognormal(μ,σ)`, `lognormal(p5=,p95=)` | log-space parameters, as in Squiggle |
| `uniform`, `cauchy`, `gamma`, `exponential` | `exponential` takes a *rate* |
| `bernoulli`, `poisson`, `triangular` | |
| `beta(a,b)`, `beta(mean=,stdev=)` | |
| `pointmass(x)` | |
| `mixture(cs, ws)`, `mx` | weights normalised; numbers become point masses |
| `mean`, `var`, `std`/`stdev`, `median` | exact when a closed form survives |
| `quantile(x,p)`, `inv_cdf(x,p)`, `credible(x,mass)` | |
| `cdf(x,v)`, `pdf(x,v)` | empirical / Gaussian-KDE when sampled |
| `rand(x)`, `rand(x,k)`, `samples(x)`, `parts(x)` | `rand`, not `sample` — `StatsBase` owns that name |
| `over(f, xs...)` | sample-wise application, for branchy code; numbers may be mixed in |
| `truncate(x,lo,hi)`, `truncateLeft`, `truncateRight` | condition on a range |
| `clamp(x,lo,hi)`, `min`, `max` | pointwise, unlike `truncate` |
| `gek(x)` | wrap a distribution, particles, or a sample vector |
| `exact(x)` | is a closed form still available? |
| `spark(x)` | the unicode histogram |
| `fosd`, `sosd`, `icx`, `statewise`, `pgeq` | the dominance orders |
| `relate(x,y)`, `margin(x,y)` | the verdict, and why it failed |
| `x ⊔ y`, `x ⊓ y`, `envelope(xs)` | lattice join, meet, and the join of a family |
| `undominated`, `dominating`, `domlevels`, `hasse`, `domgraph` | the order on a family |
| `margins`, `regret`, `robust` | the rankings for when nothing dominates |
| `certify(x,y)` | decide from closed forms; `nothing` if undecided |
| `wass(x,y)` | Wasserstein-1 |
| `sweep(f, xr, yr)` | ℝ²→distribution heatmap with the frontier outlined |
| `svgdens(x)` | inline SVG density; what `MIME"text/html"` shows |
| `distfn(f, lo, hi; n, mass, xlab)` | plot ℝ → distribution, as quantile bands |
| `nsamp()`, `setnsamp!(n)` | global particle count, default 2000 |

Closed forms cover: sums/differences of normals, products/quotients/powers of
lognormals, affine scaling of normals and lognormals, sums of Cauchys,
Poissons, equal-scale Gammas and equal-`p` Binomials, `exp` of a normal, `log`
of a lognormal, and mixtures of independent exact components. Everything else
falls through to 2000 particles, which is not a fallback so much as the normal
case.

## Name collisions

`beta` and `gamma` collide with `SpecialFunctions`; `import Gekrakel as G` if
you need both.

## Tests

```fish
julia --project=. -e 'using Pkg; Pkg.test()'
```
