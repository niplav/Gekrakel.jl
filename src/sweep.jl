# Functions from a 2-D box to a distribution.
#
# Two input dimensions plus a distribution is four-ish dimensions, and no
# single glyph survives that. So: a heatmap of one scalar summary per cell,
# with the dominance frontier outlined on it. The frontier is what turns a
# picture of "how big is the answer here" into "which inputs are even
# candidates".
#
# Diverging blue↔red with a neutral grey midpoint, four steps per arm, matched
# in lightness and chroma across the arms (the red arm is the blue steps rotated
# in OKLCH, which is what keeps a diverging scale honest). Both modes are
# *selected*, not flipped: the dark arms keep their chroma up rather than
# running to near-white pastels, where the two poles would collide. The poles
# clear CVD separation and the normal-vision floor in both modes; adjacent steps
# within an arm are deliberately close, as in any sequential ramp, and the
# legend plus the cell tooltips carry the relief.

const POSL = ("#cde2fb", "#86b6ef", "#2a78d6", "#104281")   # light, positive
const NEGL = ("#fad5d7", "#e9999f", "#c64656", "#76202b")   # light, negative
const POSD = ("#184f95", "#256abf", "#3987e5", "#6da7ec")   # dark, positive
const NEGD = ("#892a35", "#b13d4b", "#d65665", "#e4838c")   # dark, negative
const MIDL, MIDD = "#f0efec", "#383835"

"""
	sweep(f, (xlo,xhi), (ylo,yhi); n=(16,16), stat=median, xlab, ylab)

Evaluate `f(x, y)` across the box and render it as a heatmap of `stat`, with
the first-order-undominated cells outlined.

`stat` defaults to `median`. For a decision, `x -> 1 - cdf(x, 0)` is usually
more to the point — the probability the outcome beats zero.

`f` is called `n₁ · n₂` times, so a 16×16 grid is 256 model evaluations.

**Build the leaves once, outside `f`.** If `f` constructs its own
distributions, every cell mints fresh leaf ids, the cells become mutually
independent, and both the coupling and the dominance comparison lose their
meaning. `sweep` checks this and warns.
"""
struct Sweep
	xs::Vector{Float64}
	ys::Vector{Float64}
	nodes::Matrix{Node}
	vals::Matrix{Float64}
	front::Vector{Int}       # linear indices of the undominated cells
	xlab::String
	ylab::String
	slab::String
end

function sweep(f, xr::Tuple{Real,Real}, yr::Tuple{Real,Real};
		n::Tuple{Int,Int}=(16, 16), stat=median,
		xlab::AbstractString="x", ylab::AbstractString="y",
		slab::AbstractString="median", frontier::Bool=true)
	xr[1] < xr[2] && yr[1] < yr[2] || error("need lo < hi on both axes")
	xs = collect(range(float(xr[1]), float(xr[2]), length=n[1]))
	ys = collect(range(float(yr[1]), float(yr[2]), length=n[2]))
	ns = [gek(f(x, y)) for x in xs, y in ys]
	_coupled(ns)
	vals = [float(stat(v)) for v in ns]
	front = frontier ? undominated(vec(ns)) : Int[]
	Sweep(xs, ys, ns, vals, front, xlab, ylab, slab)
end

# Every cell should be a function of the same leaves. If the global
# intersection of dep sets is empty, `f` is building its own randomness and the
# whole comparison is between unrelated draws.
function _coupled(ns::AbstractArray{Node})
	shared = reduce(intersect, (v.deps for v in ns))
	isempty(shared) && @warn """
		sweep cells share no leaf, so they are mutually independent draws.
		Build the distributions once outside `f` and let `f` only combine
		them, or the dominance frontier and every coupled quantity below
		are comparing unrelated randomness."""
	nothing
end

# Bin a value onto the diverging scale: index 0 is the neutral midpoint,
# 1..4 the positive arm, -1..-4 the negative arm.
function _bin(v, lo, hi)
	sp = max(abs(lo), abs(hi))
	sp > 0 || return 0
	k = clamp(ceil(Int, 4 * abs(v) / sp), 0, 4)
	k == 0 && return 0
	v < 0 ? -k : k
end

_binfill(b) = b == 0 ? "var(--gkm)" : string("var(--gk", b < 0 ? "n" : "p", abs(b), ")")

function _gridcss(id)
	decl(pos, neg, mid) = string(
		join(("--gkp$i:$(pos[i]);" for i in 1:4)),
		join(("--gkn$i:$(neg[i]);" for i in 1:4)), "--gkm:$mid;")
	string("<style>",
		".$id{", decl(POSL, NEGL, MIDL), "}",
		"@media (prefers-color-scheme:dark){:root:not([data-theme=\"light\"]) .$id{",
		decl(POSD, NEGD, MIDD), "}}",
		":root[data-theme=\"dark\"] .$id{", decl(POSD, NEGD, MIDD), "}",
		"</style>")
end

function Base.show(io::IO, ::MIME"text/html", g::Sweep)
	nx, ny = length(g.xs), length(g.ys)
	lo, hi = minimum(g.vals), maximum(g.vals)
	cw, ch = 22, 18
	pad = 52
	w, h = nx * cw, ny * ch
	id = "gk" * string(hash((g.xs, g.ys, g.vals)); base=36)

	b = IOBuffer()
	# cells, y ascending upward. A 2px surface gap between fills, per the mark
	# spec — adjacent cells must not bleed into one another.
	for j in 1:ny, i in 1:nx
		v = g.vals[i, j]
		x = pad + (i - 1) * cw
		y = h - j * ch
		lin = (j - 1) * nx + i          # column-major linear index of nodes[i,j]
		print(b, "<rect x='", x + 1, "' y='", y + 1, "' width='", cw - 2,
			"' height='", ch - 2, "' rx='2' fill='", _binfill(_bin(v, lo, hi)),
			"'><title>", g.xlab, " ", num(g.xs[i]), ", ", g.ylab, " ",
			num(g.ys[j]), "\n", g.slab, " ", num(v),
			lin in g.front ? "\nundominated" : "", "</title></rect>")
	end
	# the frontier, drawn as a ring on the cells nothing dominates — a second
	# encoding, so membership is never carried by colour alone
	for lin in g.front
		i = mod1(lin, nx)
		j = 1 + (lin - i) ÷ nx
		print(b, "<rect x='", pad + (i - 1) * cw + 1, "' y='", h - j * ch + 1,
			"' width='", cw - 2, "' height='", ch - 2,
			"' rx='2' fill='none' stroke='currentColor' stroke-width='1.5' ",
			"stroke-opacity='.85'/>")
	end
	# axes
	print(b, _lbl(pad, h + 12, _axnum(g.xs[1])),
		_lbl(pad + w / 2, h + 12, g.xlab; anchor="middle", op=".8"),
		_lbl(pad + w, h + 12, _axnum(g.xs[end]); anchor="end"),
		_lbl(pad - 6, h - 2, _axnum(g.ys[1]); anchor="end"),
		_lbl(pad - 6, 9, _axnum(g.ys[end]); anchor="end"),
		"<g transform='translate(11,", h / 2, ") rotate(-90)'>",
		_lbl(0, 0, g.ylab; anchor="middle", op=".8"), "</g>")

	# Legend: nine swatches. The scale is symmetric about zero — it has to be,
	# or the neutral midpoint would not mean "nothing" — so the ends are ±sp,
	# not the data's own min and max. When the data is lopsided one arm goes
	# largely unused, and saying so is more honest than rescaling the arms
	# independently and quietly moving the zero.
	sp = max(abs(lo), abs(hi))
	ly = h + 26
	for (k, bn) in enumerate((-4, -3, -2, -1, 0, 1, 2, 3, 4))
		print(b, "<rect x='", pad + (k - 1) * 16, "' y='", ly,
			"' width='14' height='10' rx='2' fill='", _binfill(bn), "'/>")
	end
	print(b, _lbl(pad, ly + 22, _axnum(-sp)),
		_lbl(pad + 72, ly + 22, "0"; anchor="middle"),
		_lbl(pad + 144, ly + 22, _axnum(sp); anchor="end"),
		_lbl(pad + 156, ly + 9, g.slab; op=".8"),
		_lbl(pad + 156, ly + 22, string("data ", _axnum(lo), " … ", _axnum(hi))))

	print(io, "<div class='", id, "' style='font-family:monospace;font-size:.85em'>",
		_gridcss(id),
		_svg(pad + w + 120, h + 56, String(take!(b))),
		"<div style='opacity:.55;padding-top:6px'>", nx * ny, " cells, ",
		length(g.front), " undominated", "</div></div>")
end

function Base.show(io::IO, ::MIME"text/plain", g::Sweep)
	println(io, "Sweep ", length(g.xs), "×", length(g.ys), " over ",
		g.xlab, " ∈ [", num(g.xs[1]), ", ", num(g.xs[end]), "], ",
		g.ylab, " ∈ [", num(g.ys[1]), ", ", num(g.ys[end]), "]")
	print(io, "  ", g.slab, " ", num(minimum(g.vals)), " … ",
		num(maximum(g.vals)), ";  ", length(g.front), " of ",
		length(g.vals), " undominated")
end
