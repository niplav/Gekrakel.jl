# Inline SVG rendering.
#
# This is what makes a reactive notebook into a Squiggle playground. Pluto,
# VS Code, IJulia and Documenter all display a value by asking it for
# `MIME"text/html"`, so a good method here means every named distribution in a
# notebook draws itself with no `plot` call — which is the one thing Squiggle's
# viewer does that a plain REPL does not.
#
# Everything is drawn in `currentColor` at varying opacity, so it inherits the
# host page's foreground and needs no light/dark handling of its own. No
# dependencies: the shapes are strings.

_r(x) = string(round(float(x), digits = 2))
_pts(ps) = join((_r(x) * "," * _r(y) for (x, y) in ps), " ")
_svg(w, h, body) = string("<svg viewBox='0 0 ", w, " ", h, "' width='100%' ",
	"style='max-width:", w, "px;display:block;overflow:visible'>", body, "</svg>")
# Axis labels have a fixed gutter to fit in, so they get a more compact
# rendering than `num`: "-1.7e5" rather than "-167000.0".
function _axnum(x::Real)
	a = abs(x)
	(a >= 1e4 || (a > 0 && a < 1e-2)) || return num(x)
	e = floor(Int, log10(a))
	string(round(x / 10.0^e, digits = 1), "e", e)
end

_lbl(x, y, s; anchor = "start", op = ".55") = string(
	"<text x='", _r(x), "' y='", _r(y), "' text-anchor='", anchor,
	"' font-size='9' font-family='monospace' fill='currentColor' ",
	"fill-opacity='", op, "'>", s, "</text>")

# ---------------------------------------------------------------- densities

function _curve(n::Node, lo, hi, q5, q95, w, h)
	k = 121
	xs = range(lo, hi, length = k)
	ys = [pdf(n, x) for x in xs]
	m = maximum(ys)
	m > 0 || (m = 1.0)
	X(x) = (x - lo) / (hi - lo) * w
	Y(y) = h - 1 - (h - 3) * y / m
	curve = [(X(x), Y(y)) for (x, y) in zip(xs, ys)]
	area = vcat([(0.0, float(h))], curve, [(float(w), float(h))])
	inner = [(X(x), Y(y)) for (x, y) in zip(xs, ys) if q5 <= x <= q95]
	band = isempty(inner) ? "" : string("<polygon points='",
		_pts(vcat([(inner[1][1], float(h))], inner, [(inner[end][1], float(h))])),
		"' fill='currentColor' fill-opacity='.22'/>")
	string("<polygon points='", _pts(area), "' fill='currentColor' fill-opacity='.09'/>",
		band, "<polyline points='", _pts(curve),
		"' fill='none' stroke='currentColor' stroke-opacity='.75' stroke-width='1.5'/>")
end

# One stick per outcome. Binning a discrete support on a continuous grid
# produces a comb of empty bins, the same reason `spark` special-cases it.
function _sticks(n::Node, lo, hi, w, h)
	ks = floor(Int, lo):ceil(Int, hi)
	length(ks) > 60 && return _curve(n, lo, hi, lo, hi, w, h)
	ps = [pdf(n.dist, k) for k in ks]
	m = maximum(ps)
	m > 0 || (m = 1.0)
	X(k) = length(ks) == 1 ? w / 2 : (k - lo) / (hi - lo) * w
	join(string("<line x1='", _r(X(k)), "' y1='", h - 1, "' x2='", _r(X(k)),
		"' y2='", _r(h - 1 - (h - 3) * p / m),
		"' stroke='currentColor' stroke-opacity='.7' stroke-width='2'/>")
		for (k, p) in zip(ks, ps))
end

"""
	svgdens(x; w=340, h=78)

`x` as an inline SVG density, shaded across its 90% credible interval, with a
dashed median. Clipped to the 0.5th..99.5th percentile so one outlier cannot
flatten the picture. Exact distributions are drawn from their density; sampled
ones from the KDE in `pdf`.
"""
function svgdens(n::Node; w::Int = 340, h::Int = 78)
	lo, hi = quantile(n, 0.005), quantile(n, 0.995)
	q5, q50, q95 = quantile(n, [0.05, 0.5, 0.95])
	hi > lo || return _svg(w, h + 12, string("<line x1='", w ÷ 2, "' y1='", h - 1,
		"' x2='", w ÷ 2, "' y2='2' stroke='currentColor' stroke-width='2'/>",
		_lbl(w ÷ 2, h + 9, num(lo); anchor = "middle")))
	body = (exact(n) && n.dist isa DiscreteUnivariateDistribution) ?
		_sticks(n, lo, hi, w, h) : _curve(n, lo, hi, q5, q95, w, h)
	mx = (q50 - lo) / (hi - lo) * w
	# the viewBox is 12 taller than the plot so the axis labels are inside it;
	# `overflow:visible` is not honoured everywhere HTML is rendered
	_svg(w, h + 12, string(body,
		"<line x1='", _r(mx), "' y1='", h - 1, "' x2='", _r(mx), "' y2='2' ",
		"stroke='currentColor' stroke-opacity='.5' stroke-dasharray='2,2'/>",
		"<line x1='0' y1='", h - 1, "' x2='", w, "' y2='", h - 1,
		"' stroke='currentColor' stroke-opacity='.25'/>",
		_lbl(0, h + 9, _axnum(lo)), _lbl(w, h + 9, _axnum(hi); anchor = "end")))
end

function Base.show(io::IO, ::MIME"text/html", n::Node)
	p5, p50, p95 = quantile(n, [0.05, 0.5, 0.95])
	print(io, "<div style='font-family:monospace;font-size:.85em;padding-bottom:12px'>",
		svgdens(n),
		"<div style='padding-top:12px'><b>5%</b> ", num(p5),
		" &nbsp; <b>50%</b> ", num(p50), " &nbsp; <b>95%</b> ", num(p95),
		" &nbsp; <b>mean</b> ", num(mean(n)),
		"</div><div style='opacity:.55'>", kind(n),
		exact(n) ? " (exact)" : " ($(length(samples(n))) samples)",
		"</div></div>")
end

# ------------------------------------------------------------ ℝ → Dist plots

"""
	distfn(f, lo, hi; n=48, mass=(0.9, 0.5), xlab="x")

Plot a function from a number to a distribution — Squiggle's `Plot.distFn`.
`f` is evaluated on `n` points across `[lo, hi]` and drawn as nested quantile
bands around the median.

	distfn(h -> value(M, P, h, 90).net, 0.5, 20; xlab="hr/day")

Each point is an independent evaluation of `f`, so a model that materialises
particles costs `n` times as much as evaluating it once; keep `n` small while
iterating. Correlations *within* one evaluation are preserved as usual, but
neighbouring x are independent draws, so a band edge is jagged where the model
is sample-based. That jaggedness is an honest picture of the Monte Carlo error
and is worth reading rather than smoothing away.
"""
struct DistFn
	xs::Vector{Float64}
	mass::Vector{Float64}
	lo::Vector{Vector{Float64}}
	hi::Vector{Vector{Float64}}
	med::Vector{Float64}
	xlab::String
end

function distfn(f, a::Real, b::Real; n::Int = 48, mass = (0.9, 0.5), xlab = "x")
	a < b || error("need lo < hi, got $a … $b")
	xs = collect(range(float(a), float(b), length = n))
	ms = sort(collect(float.(mass)), rev = true)     # widest first, so it is
							 # drawn underneath
	ds = map(xs) do x
		v = f(x)
		v isa Node ? v : pointmass(v)
	end
	los = [[quantile(d, (1 - m) / 2) for d in ds] for m in ms]
	his = [[quantile(d, 1 - (1 - m) / 2) for d in ds] for m in ms]
	DistFn(xs, ms, los, his, [median(d) for d in ds], xlab)
end

const PAD = 46      # left gutter for the y labels, inside the viewBox

function Base.show(io::IO, ::MIME"text/html", p::DistFn; w::Int = 380, h::Int = 150)
	ylo = minimum(minimum, p.lo)
	yhi = maximum(maximum, p.hi)
	yhi > ylo || (yhi = ylo + 1)
	X(x) = PAD + (x - p.xs[1]) / (p.xs[end] - p.xs[1]) * w
	Y(y) = (h - 1) - (h - 3) * (y - ylo) / (yhi - ylo)

	body = IOBuffer()
	for (i, m) in enumerate(p.mass)
		up = [(X(x), Y(y)) for (x, y) in zip(p.xs, p.hi[i])]
		dn = [(X(x), Y(y)) for (x, y) in zip(reverse(p.xs), reverse(p.lo[i]))]
		print(body, "<polygon points='", _pts(vcat(up, dn)),
			"' fill='currentColor' fill-opacity='", _r(0.10 + 0.09i), "'/>")
	end
	# zero is the decision boundary whenever the y-axis is money, so draw it
	# when it is in range
	if ylo < 0 < yhi
		print(body, "<line x1='", PAD, "' y1='", _r(Y(0)), "' x2='", PAD + w,
			"' y2='", _r(Y(0)),
			"' stroke='currentColor' stroke-opacity='.45' stroke-dasharray='3,3'/>",
			_lbl(PAD - 4, Y(0) + 3, "0"; anchor = "end", op = ".7"))
	end
	print(body, "<polyline points='", _pts([(X(x), Y(y)) for (x, y) in zip(p.xs, p.med)]),
		"' fill='none' stroke='currentColor' stroke-opacity='.9' stroke-width='1.5'/>")

	# drop an endpoint label that would sit on top of the zero label
	ends = join(_lbl(PAD - 4, Y(y) + o, _axnum(y); anchor = "end")
		for (y, o) in ((yhi, 7), (ylo, 0))
		if !(ylo < 0 < yhi) || abs(Y(y) - Y(0)) > 11)

	print(io, "<div style='font-family:monospace;font-size:.85em;padding:0 0 14px 0'>",
		_svg(PAD + w + 2, h + 13, string(String(take!(body)), ends,
			_lbl(PAD, h + 11, _axnum(p.xs[1])),
			_lbl(PAD + w / 2, h + 11, p.xlab; anchor = "middle", op = ".8"),
			_lbl(PAD + w, h + 11, _axnum(p.xs[end]); anchor = "end"))),
		"<div style='padding-top:14px;opacity:.55'>median, ",
		join((string(round(Int, 100m), "%") for m in p.mass), " and "),
		" bands over ", length(p.xs), " points</div></div>")
end

function Base.show(io::IO, ::MIME"text/plain", p::DistFn)
	println(io, "DistFn over ", p.xlab, " ∈ [", num(p.xs[1]), ", ",
		num(p.xs[end]), "], ", length(p.xs), " points")
	m, M = extrema(p.med)
	b = M > m ?
		String([BLOCKS[clamp(1 + round(Int, 7(v - m) / (M - m)), 1, 8)] for v in p.med]) :
		"▄"^length(p.med)
	println(io, "  ", b)
	print(io, "  median ", num(m), " … ", num(M))
end
