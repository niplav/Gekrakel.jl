# How long from "go" to the first bit leaving the fast-track transmitter?
#
# Same rig as fasttrack.jl: 1 m commercial telescope, no AO, catalogue
# 10 J @ 10 Hz Nd:YAG, professional staff hired from scratch, money in hand.
# A small critical-path network; every duration in months, 90% intervals.
#
# Neither PlaneWave nor Ekspla publishes lead times (both are made to order),
# so those are guesses. So is everything about permits, which vary wildly by
# country and turn out to be the thing that matters.
#
# Two ways to get a building:
#   greenfield: lease bare land, get permits, pour a foundation, raise a dome
#   brownfield: an existing pad or vacant observatory building, plus a prefab
#               laser room and a roll-off roof instead of a custom dome

pushfirst!(LOAD_PATH, dirname(@__DIR__))
using Gekrakel
setnsamp!(20_000)

# ---------------------------------------------------------------------------
# Common to both

hire      = 2…5             # professional team on board; the lead sooner, engineers later
spec      = 1…3             # laser spec frozen, quote in, order placed
tel_lead  = 4…12            # PW1000, made to order
las_lead  = 4…9             # 10 J @ 10 Hz NanoFlux-class, made to order
ship      = 0.5…1.5         # freight plus customs, per shipment
tel_inst  = 1…2             # mount, collimate, pointing model
las_inst  = 0.5…1.5         # into the laser room, alignment, acceptance
design    = 2…4             # beam bench, interlocks, aircraft detection
fab       = 2…4             # parts in, bench assembled
integrate = 1…3             # bench meets telescope meets laser
freeze    = 1…3             # design stable enough to file for clearance
review    = 1…6             # Laser Clearinghouse / aviation authority
commish   = 1…4             # low-power tests, pointing, first debugging

# ---------------------------------------------------------------------------
# Site scenarios: (find & lease, permits, civil works, enclosure lead, enclosure install)

greenfield = (site = 2…6,   permit = 2…9,   civil = 2…5,   enc = 3…6, enc_inst = 0.5…1.5)
brownfield = (site = 0.5…2, permit = 0.3…2, civil = 0.5…2, enc = 1…3, enc_inst = 0.3…1)

const paths = ("building", "telescope delivery", "laser", "beam bench", "clearance")

# Returns (months to first bit, index into `paths` of what was critical).
function firstbit(site, permit, civil, enc, enc_inst, hire, spec, tel_lead, las_lead,
		ship, tel_inst, las_inst, design, fab, integrate, freeze, review, commish)
	building  = site + permit + civil
	enclosure = max(building, enc) + enc_inst
	telescope = max(enclosure, tel_lead + ship) + tel_inst
	laser     = max(spec + las_lead + ship, building) + las_inst
	bench     = hire + design + fab
	hardware  = max(telescope, laser, bench) + integrate
	clearance = hire + freeze + review

	crit = if clearance > hardware
		5
	elseif telescope >= max(laser, bench)
		enclosure >= tel_lead + ship ? 1 : 2
	elseif laser >= bench
		3
	else
		4
	end
	(max(hardware, clearance) + commish, crit)
end

function run(s)
	leaves = (s.site, s.permit, s.civil, s.enc, s.enc_inst, hire, spec, tel_lead, las_lead,
		ship, tel_inst, las_inst, design, fab, integrate, freeze, review, commish)
	months = over((xs...) -> firstbit(xs...)[1], leaves...)
	crit   = over((xs...) -> firstbit(xs...)[2], leaves...)
	months, crit
end

# ---------------------------------------------------------------------------

for (name, s) in (("GREENFIELD", greenfield), ("BROWNFIELD", brownfield))
	months, crit = run(s)
	println(name)
	println("  months to first bit   ", months)
	for cutoff in (6, 9, 12, 18, 24)
		println("  P(≤ ", lpad(cutoff, 2), " months)       ", round(cdf(months, cutoff), digits = 2))
	end
	println("  critical path:")
	c = samples(crit)
	for (i, p) in enumerate(paths)
		share = count(==(i), c) / length(c)
		share > 0.005 && println("    ", rpad(p, 20), round(Int, 100 * share), "%")
	end
	println()
end
