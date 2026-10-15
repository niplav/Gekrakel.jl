# When does the next Millennium Prize problem fall?
# Navier–Stokes fell 2026-09-08; five remain. Competing risks: ordinary
# mathematics, versus an AI horizon crossing the residual difficulty.

using Gekrakel, Distributions
setnsamp!(20_000)

# AI Impacts, Kaplan–Meier over Wikipedia conjectures: time-to-proof is
# exponential, half-life ~117y ⇒ mean 169y. Memoryless, so Riemann's 167 years
# of survival say nothing — residual difficulty is distributed as on day one.
# Millennium problems draw more effort than the average conjecture (2 of 7 in
# 26y ⇒ λ̂ ≈ 2× base, on n=2), so scale the clock.
attn = 1…4
years = over(a -> rand(Exponential(169 / (5a))), attn)   # min over the five
human = 12 * years                          # months to resolution, if nothing changes

# The same difficulty in expert-hours: effort arrives at `effort` hours a year.
effort = 3e3…3e4
difficulty = years * effort * attn  # mean ~5e5 h, exponential tail

# Frontier reach now, same units, anchored on what fell routinely in 2026:
# unit-distance (May), Grothendieck & Jacobian (Jul). N–S is not an anchor —
# one hit is a draw, not the frontier.
frontier = 1e3…3e4
doubling = (3…5) / (1.2…3)               # METR months ÷ verifiability uplift

# The trend can break either way. Downward: it stalls outright. Upward: each
# successive doubling takes `rho` times as long as the last, the superexponential
# branch AI 2027 argues for. rho = 1 recovers a plain exponential.
rho = mx([pointmass(1.0), 0.8…0.95], [0.6, 0.4])
# Time for `count` doublings when the first takes `first` and each is `ratio`
# times the last: a geometric sum, and plain `first * count` when ratio is 1.
elapsed(first, ratio, count) =
	first * (isone(ratio) ? count : (1 - ratio^count) / (1 - ratio))

machine = over(doubling, rho, difficulty, frontier) do first, ratio, difficulty, frontier
	elapsed(first, ratio, log2(max(difficulty / frontier, 1)))
end
machine = mx([machine, machine + (18…120)], [0.85, 0.15])     # 15%: it stalls

months = over(min, machine, human)
println("months: ", months)
for (label, cutoff) in ("2026-12"=>3, "2027-12"=>15, "2029-12"=>39, "2036-12"=>123)
	println(label, "  P = ", round(cdf(months, cutoff), digits=2))
end
println("human-only 2029-12: P = ", round(cdf(human, 39), digits=3))
