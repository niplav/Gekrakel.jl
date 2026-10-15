# What does the fast-track transmitter actually cost, all in?
#
# The rig: a 1 m commercial telescope, no adaptive optics (without AO the
# atmosphere caps the gain at ~(πr₀/λ)² whatever the aperture, so 1 m is as
# good as 4 m), a catalogue 10 J @ 10 Hz nanosecond Nd:YAG, own site, one year
# of build and one year of transmitting. Everything in $k, 90% intervals.
#
# Prices checked where they're public (Oct 2026): PlaneWave PW1000 $575k list;
# Ash-Dome 20'6" $101–120k. Ekspla sells a 10 J @ 10 Hz NanoFlux but doesn't
# publish a price, so that line is a guess, and the biggest uncertain item.
#
# Own site, not a hosted pier: an uplink beam's Rayleigh/Mie scatter is visible
# to telescopes kilometres away, which is why Mauna Kea and La Palma run laser
# traffic control. A shared site means negotiating that; a quiet one is simpler.

pushfirst!(LOAD_PATH, dirname(@__DIR__))
using Gekrakel
redirect_stdout(devnull) do        # borrow the reception model, mute its table
	include(joinpath(@__DIR__, "beam.jl"))
end

# ---------------------------------------------------------------------------
# Hardware

telescope = 550…700         # PW1000 at list, ± options and discount
dome      = 120…250         # 20'6" Ash-Dome plus installation
laser     = 700…2500        # 10 J @ 10 Hz flashlamp Nd:YAG, unpublished price
chiller   = 30…100          # cooling and power conditioning the laser doesn't ship with
bench     = 100…400         # beam expander, high-damage-threshold mirrors, pupil injection
safety    = 50…250          # aircraft detection (cameras, ADS-B), interlocks, shutters
controls  = 20…80           # GPS timing, PPM modulator, control computers
guider    = 10…40           # acquisition camera, pointing model
spares    = 50…200          # flashlamps (~3e8 shots/yr means several changes), optics

hardware = telescope + dome + laser + chiller + bench + safety + controls + guider + spares

# ---------------------------------------------------------------------------
# Site

civil = 200…900             # foundation, pier, laser room, control room
power = 30…200              # ~20–30 kW hookup, UPS, backup generator
land  = 10…80               # one year's lease

site = civil + power + land

# ---------------------------------------------------------------------------
# People: full-time-equivalent years times loaded cost per FTE-year.
# Build: lead, 1–2 laser engineers, opto-mechanical, controls/software, part-time
# laser safety officer and regulatory work. Operations: mostly automated nights,
# with someone on call and periodic laser maintenance.

loaded  = 120…250           # $k per FTE-year, Europe/US, salary plus overhead
build   = 4…7               # FTE-years
operate = 1.5…3             # FTE-years

people = (build + operate) * loaded

# ---------------------------------------------------------------------------
# Moving things and permission

travel  = 40…150            # factory acceptance (Vilnius, Michigan), site trips, integration stays
freight = 30…150            # ~3 t telescope crate, laser, chillers; insured
duties  = uniform(0, 0.15) * hardware   # zero to ~15% depending on country
legal   = 30…200            # Laser Clearinghouse / aviation filings, permits, insurance
service = (0.03…0.1) * laser            # laser service contract, year one
energy  = 10…50             # ~20–30 kW around the clock, at most

logistics = travel + freight + duties + legal + service + energy

# Integration surprises the line items don't see: things arrive late, benches
# misalign, the first months go to debugging rather than transmitting.
overrun = 1.05…1.6

total = (hardware + site + people + logistics) * overrun

# ---------------------------------------------------------------------------
# What it buys: one year of 10 Hz at a 30–50% duty cycle (nights, weather,
# satellite-avoidance windows), stacked on one 1 MB message.

r0     = 0.15…0.3           # m, turbulence scale at 1 μm
duty   = 0.3…0.5
copies = 10 * bits * 3.156e7 * duty / 8e6
cheap  = (g = gain(1), s = over(r -> min(1, r^2), r0) * (10 / pulse) * copies, aim = cluster)
ep     = received(cheap)

# ---------------------------------------------------------------------------
# Output

usd(k) = k >= 1000 ? string(round(k / 1000, digits = 2), "M") : string(round(Int, k), "k")
ci(x)  = string(usd(quantile(x, 0.05)), " … ", usd(quantile(x, 0.95)))

items = [("laser", laser), ("people", people), ("telescope", telescope),
	("site works", civil), ("beam bench", bench), ("dome", dome),
	("safety", safety), ("duties", duties), ("spares", spares),
	("power", power), ("legal & insurance", legal), ("laser service", service),
	("travel", travel), ("freight", freight), ("chiller", chiller),
	("controls", controls), ("land", land), ("guider", guider), ("energy", energy)]
sort!(items, by = x -> -mean(x[2]))

sumall = sum(mean(x) for (_, x) in items)
println(rpad("item", 20), lpad("median", 9), "   ", rpad("90% interval", 20), lpad("share", 6))
println("─"^60)
for (name, x) in items
	println(rpad(name, 20), lpad(usd(median(x)), 9), "   ", rpad(ci(x), 20),
		lpad(string(round(Int, 100 * mean(x) / sumall), "%"), 6))
end
println("─"^60)
println(rpad("before overrun", 20), lpad(usd(median(total / overrun)), 9), "   ", ci(total / overrun))
println(rpad("TOTAL, with overrun", 20), lpad(usd(median(total)), 9), "   ", ci(total))
println()
println("copies of a 1 MB message sent in the year: ", round(Int, median(copies)))
println("E[P(reception)]:             ", round(mean(ep), digits = 3))
println("expected receptions per \$M:  ", round(mean(ep) / mean(total / 1000), digits = 4))
println("\$M per expected reception:   ", round(mean(total / 1000) / mean(ep), digits = 1))
