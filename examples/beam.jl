# Is it worth beaming humanity's archive at the stars before we're disempowered?
#
# Every option fires 1 kJ, 1 ns pulses at 1.064 μm. They differ in what
# directs the beam (a 10 m mirror on the ground, the same mirror in orbit, or
# the Sun's gravitational lens from 550 AU), in repetition rate, and in where
# they point (a random patch of sky, or down a filament into a rich cluster).
#
# The recipients are Hanson-style grabby aliens. Our own unendorsed successor
# is assumed grabby too, so their front and its front meet at a border; the
# point of a beam is that a photon outruns both fronts and lands deep inside
# their territory, beyond anything our successor can ever reach or veto.

# Run from anywhere, with or without --project: stack the repo's environment.
pushfirst!(LOAD_PATH, dirname(@__DIR__))

using Gekrakel
setnsamp!(20_000)

const ly    = 9.4607e15                 # m
const lam   = 1.064e-6                  # m, Nd:YAG
const hnu   = 6.626e-34 * 2.998e8 / lam # J per photon
const ndet  = 10                        # photons in the 1 ns gate to call it a detection
const pulse = 1e3                       # J
const bits  = 10                        # per detected pulse, with 1024-slot PPM

# Antenna gain of a diffraction-limited aperture, and the Sun's on-axis
# magnification (Turyshev–Toth, with r_s = 2GM/c² = 2953 m).
gain(d) = (pi * d / lam)^2
sgl = 4pi^2 * 2953 / lam                # ≈ 1.2e11

# ---------------------------------------------------------------------------
# Geometry (Hanson, Martin, McCarter & Paulson 2021)
#
# Nearest grabby origin 1–4 Gly; everyone expands at v. The two fronts close
# at 2v, so the border forms at distance v·t_border = D/2 from us. A photon
# leaving now meets their inbound front at D/(1+v): further out than the
# border, and long before it forms.
#
# The paper also quotes contact in 0.2–2 Gyr, which doesn't quite cohere with
# 1–4 Gly at closing speed c. Distance is the sturdier number, so drive from it;
# driving from time instead halves the distance and quadruples the photons.

v    = 0.3…0.8                  # expansion speed, units of c
dist = (1…4) * 1e9              # ly to their origin
t_border = over((v, d) -> d / 2v, v, dist)          # yr
d_border = v * t_border                             # ly
d_meet   = over((v, d) -> d / (1 + v), v, dist)     # ly
lead     = (1 - v) * d_meet     # how far the photon is ahead of our successor then

# ---------------------------------------------------------------------------
# Recipients: a percolating front, not a filled ball
#
# Expansion over a sparse graph of stars and galaxies settles only some of the
# sites it sweeps (Landis 1998), and only some of those carry a big detector
# pointed somewhere useful.

stars   = 2e-11…2e-10           # per ly³, cosmic mean (galaxies × stars/galaxy)
grabbed = 0.3…0.6               # fraction of distant volume already claimed
settled = 0.1…1.0               # of swept sites, fraction actually occupied
watched = 0.003…0.3             # of occupied sites, fraction with a detector
# How far out the cone counts. Beyond ~3 Gly this model's flat-space 1/d² and
# d³ go wrong (flux should use luminosity distance, volume comoving distance),
# so treat the far end as a nuisance parameter. It decides how much aiming
# matters: a long cone dwarfs any single cluster.
reach   = 2e9…1e10              # ly

# Their per-site collecting area is the parameter that decides everything.
# 1e6 m² ≈ a km-scale array, 1e14 ≈ Earth's cross-section, 7e22 ≈ a 1 AU swarm.
area = 1e6…1e23                 # m²

# ...and someone has to be running all-sky nanosecond photometry at all.
looking = 0.2…0.9

# ---------------------------------------------------------------------------
# The options' building blocks

# Strehl ratio: fraction of the ideal gain surviving the atmosphere under
# laser-grade adaptive optics at 1 μm. Without AO it'd be ~(r₀/D)² ≈ 4e-4.
ao = clamp(0.3…0.7, 0, 1)

# Aiming down a filament into a rich cluster core adds an overdense stretch
# to the path: (overdensity ~200–1000) × (core path ~2 Mpc), in ly, plus a
# filament's ~8× over ~30 Mpc. Cluster centres drift ~1 Mpc over the 2.6 Gyr
# between the light we see and our photon's arrival, inside the core, so the
# aim holds; an individual galaxy drifts ~20× its own size and would not.
cluster = 1e9…7e9

tele_ground = 0.3…1             # $B, 10 m class with laser-grade AO, on a mountain
tele_space  = 10…40             # $B, 10 m in orbit plus ~10 MW of solar array
laser_kw    = 0.02…0.1          # $B, 1 kJ @ 1 Hz: a DiPOLE100-class laser with ~10× the pulse energy
laser_mw    = 0.3…3             # $B, 1 kJ @ 1 kHz: ~1000 such beamlines, pump diodes dominate
sgl_craft   = 40…150            # $B, 10 MW fission craft, 30 yr cruise to 550 AU

options = [
	(name = "ground, 1 Hz, random sky",     g = gain(10),       s = ao,  rate = 1,    aim = 0.0,     cost = tele_ground + laser_kw),
	(name = "ground, 1 Hz, at a cluster",   g = gain(10),       s = ao,  rate = 1,    aim = cluster, cost = tele_ground + laser_kw),
	(name = "ground, 1 kHz, at a cluster",  g = gain(10),       s = ao,  rate = 1000, aim = cluster, cost = tele_ground + laser_mw),
	(name = "orbit, 1 kHz, at a cluster",   g = gain(10),       s = 1.0, rate = 1000, aim = cluster, cost = tele_space + laser_mw),
	(name = "Sun lens, 1 kHz, at a cluster", g = gain(1) * sgl, s = 1.0, rate = 1000, aim = cluster, cost = sgl_craft + laser_mw),
]

# ---------------------------------------------------------------------------
# Model

photons(o, d) = pulse * o.g * o.s * area / (4pi * (d * ly)^2 * hnu)

# Monitoring sites the beam passes while it's still bright enough to detect.
# Photons fall as 1/d², so the beam stays detectable out to d_meet·√(N/ndet).
function sites(o; coverage = 1)
	p    = photons(o, d_meet)
	far  = over((p, d, r) -> p >= ndet ? clamp(d * sqrt(p / ndet), d, max(d, r)) : d, p, d_meet, reach)
	spot = 4pi / o.g                        # beam solid angle, sr
	cone = stars * spot * (far^3 - d_meet^3) / 3
	aimed = over((f, d, a) -> f > d ? a : 0.0, far, d_meet, o.aim) * stars * spot * d_meet^2
	(cone + aimed) * grabbed * settled * watched * coverage
end

received(o; coverage = 1) = over((n, l) -> (1 - exp(-n)) * l, sites(o; coverage), looking)

# ---------------------------------------------------------------------------
# Output

gly(x) = x / 1e9
fmt(x) = string(round(x, sigdigits = 2))

println("GEOMETRY  (5% … 50% … 95%)")
println("  their origin, Gly          ", gly(dist))
println("  border forms at, Gly       ", gly(d_border))
println("  photon meets front, Gly    ", gly(d_meet))
println("  photon's lead over our successor then, Gly  ", gly(lead))
println("  → arrives ", round(median(d_meet / d_border), digits = 2), "× deeper than the border, ",
	round(Int, median(t_border - d_meet) / 1e6), " Myr before it forms")
println()

cols = ("option", "photons/site", "P(closes)", "sites|closes", "E[P]", "E[P] patchy", "MB/yr", "\$B", "E[P]/\$B")
w    = (30, 12, 9, 12, 7, 11, 7, 5, 8)
line(xs) = println(join((i == 1 ? rpad(x, w[i]) : lpad(x, w[i]) for (i, x) in enumerate(xs)), " "))

println("OPTIONS")
println("  photons/site: median, at the front.  P(closes): photons ≥ $ndet there.")
println("  sites|closes: median monitoring sites in the beam, in worlds where it closes.")
println("  patchy: coverage ×1e-4.  E[P]/\$B: expected receptions per billion dollars.")
line(cols)
line(map(n -> "─"^n, w))
for o in options
	r = received(o)
	p = photons(o, d_meet)
	n = filter(>(0), samples(sites(o)))
	line((o.name,
		fmt(median(p)),
		fmt(mean(over(x -> x >= ndet, p))),
		isempty(n) ? "–" : fmt(median(n)),
		fmt(mean(r)),
		fmt(mean(received(o; coverage = 1e-4))),
		fmt(o.rate * bits * 3.156e7 / 8e6),
		fmt(median(o.cost)),
		fmt(mean(r / o.cost))))
end
println()

needed = ndet * 4pi * (median(d_meet) * ly)^2 * hnu / (pulse * gain(10) * median(ao))
println("Collecting area per site that the ground options need at the median distance:")
println("  ", fmt(needed), " m²  (Earth's cross-section is 1.3e14); prior mass above it: ",
	fmt(mean(over(a -> a > needed, area))))
