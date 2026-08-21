# Density curve drawn from `pdf`, which is analytic for exact nodes and a
# Gaussian KDE otherwise. Doing it here rather than deferring to
# `seriestype=:density` keeps StatsPlots out of the dependency list.

@recipe function f(n::Node; npoints=200)
	lo, hi = quantile(n, 0.005), quantile(n, 0.995)
	legend --> false
	if exact(n) && n.dist isa DiscreteUnivariateDistribution
		ks = floor(Int, lo):ceil(Int, hi)
		seriestype --> :sticks
		markershape --> :circle
		collect(ks), [pdf(n, k) for k in ks]
	else
		hi > lo || return [lo], [1.0]
		xs = range(lo, hi, length=npoints)
		seriestype --> :path
		fillrange --> 0
		fillalpha --> 0.25
		collect(xs), [pdf(n, x) for x in xs]
	end
end
