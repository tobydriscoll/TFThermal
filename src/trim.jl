"Find where the time series is strictly decreasing"
function findends(x)
	d = diff(x)
	ib = findfirst(<(0),d)
	ie = findlast(<(0),d)
    return (isnothing(ib) | isnothing(ie)) ? [] : ib:ie+1
end

function findspikes(x,window,Z=2)
	elim = []
	for i in 1:window:length(x)
		z = zscore( x[i:min(end,i+window)] )
		append!(elim, i-1 .+ findall(@. abs(z)>Z) )
    end
    return elim
end

function smooth(x)
	n = length(x)
	x̄ = [ mean(x[max(1,i-1):min(end,i+1)]) for i in 1:n ]
	# Fit cubic to lower half to smooth aggressively
	m,M = extrema(x̄)
	k = findfirst( <((m+M)/2), x̄)
	p = Polynomials.fit(collect(k:n),x̄[k:n],3)
	return [x̄[1:k]; p.(k+1:n)]
end

"Rayanne's method for finding the best region of a time series"
function rayanne(x, frame_rate=20)
	n = length(x)
	rate(d, i, j) = d < -0.6*(maximum(x) - minimum(x)) ? d/(j-i) : 0
	drop_rate = zeros(n,n)
	for i in 1:n, j in i+frame_rate-1:n
		drop_rate[i,j] = rate(x[j]-x[i], i, j)
	end
	gdbegin,gdend = Tuple(argmin(drop_rate))

	# Walk up the hill to the left
	while x[gdbegin] < x[max(gdbegin-2*frame_rate,1)]
		(gdbegin==1) && break
		gdbegin = gdbegin - 1
	end

	# Walk down the hill to the right
    while x[min(gdend+2*frame_rate,end)] < 0.98x[gdend]
		(gdend==n) && break
		gdend = gdend + 1
	end

	# Maybe there was an earlier bump?
	rise = drop_rate[gdbegin:gdend, gdbegin:gdend] .> 5
	if any(rise)
		pos = findfirst(rise)
		gdend = gdbegin + pos[1]
	end

	return gdbegin:gdend
end

function initialtrim(x, window=5; kw...)
	idx = collect(findends(x))
	xx = Hampel.filter(x[idx], window, boundary=:repeat)
	keep = falses(length(xx))
	# push ends outward
	keep[rayanne(smooth(xx); kw...)] .= true
	idx = idx[keep]

	# go from the max to the min
	isempty(idx) && (return idx)
	M = argmax(x[idx])
	m = argmin(x[idx])
	return idx[M:m]
end

function adjustfit(
	meas::MeasuredValues, t, I, θ;
	restol=0.025,
	frame_rate=20,
	verbose=false
	)
	idx = initialtrim(I; frame_rate)
	n = length(idx)
	if (n < 11) || (maximum(I[idx]) > 2I[idx[1]])
		@warn "Insufficient usable data. Giving up."
		return [], false, idx
	end
	done(f) = f.residual < restol

	bestmodel = FittedModel(TFModelExp, meas, t[idx], I[idx], θ[idx])
	success = done(bestmodel)
	while !success
		if verbose
			@info "n = $(length(idx)), residual: $(bestmodel.residual)"
		end
		idx = idx[1:end-frame_rate]
		if length(idx) < 0.75n
			@warn "Unable to shorten data further. Giving up."
			break
		end
		bestmodel = FittedModel(TFModelExp, meas, t[idx], I[idx], θ[idx])
		success = done(bestmodel)
	end
	return bestmodel, success, idx
end
