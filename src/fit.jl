# create callback to halt solver if dh/dt > 0 or I increases at a step
function solvercb(Ifun)
	CBFOO = [0.,0.]
	function increaseI(u,t,integrator)
		return (u[1] > integrator.uprev[1]) || (Ifun(u[1], u[2]) > Ifun(integrator.uprev...))
	end
	function posdhdt(u,t,integrator)
		integrator.f(CBFOO, u, integrator.p, t)
		return CBFOO[1]-1e-3
	end
	affect!(integrator) = terminate!(integrator)
	solvercb = CallbackSet(DiscreteCallback(increaseI, affect!),ContinuousCallback(posdhdt, affect!))
end

struct FittedModel
	model::AbstractModel
	t::AbstractVector{<:Quantity}
	I::AbstractVector{<:Real}
	Tₕ::AbstractVector{<:Quantity}
	residual::AbstractFloat
end

# create a fitted model from dimensional data
function FittedModel(
	MT::Type{<:AbstractModel},
	measured::MeasuredValues,
	t::AbstractVector{<:Quantity},
	I::AbstractVector{<:Real},
	Tₕ::AbstractVector{<:Quantity};
	kw...
	)
	@assert measured.ts ≈ t[end] - t[1] "Inconsistent time scale."

	# Nondimensionalize data.
	tt = (t .- t[1]) / measured.ts
	II = I / I[1]
	TT = (Tₕ .- uconvert(u"K",Ts)) / (Tb - Ts)

	min_res, min_param = fit(MT, measured, tt, II, TT; kw...)
	model = solve(MT(min_param, measured))
	return FittedModel(model, t, I, Tₕ, min_res)
end

constants(M::FittedModel) = constants(model(M))
parameters(f::FittedModel) = parameters(f.m)
names(f::FittedModel) = names(f.m)
model(f::FittedModel) = f.m

solution(f::FittedModel,t::Real;kwargs...) = solution(model(f),t;kwargs...)
function solution(f::FittedModel,t::Unitful.Time;kwargs...)
	solution(model(f),(t-f.t₀)/f.ts;kwargs...)
end
solution(f::FittedModel,t::AbstractVector;kwargs...) = [ solution(f,t;kwargs...) for t in t ]

intensity(f::FittedModel,t::Real) = f.I₀*intensity(model(f),t)
intensity(f::FittedModel,t::Unitful.Time) = f.I₀*intensity(model(f),(t-f.t₀)/f.ts)
intensity(f::FittedModel,t::AbstractVector) = [ intensity(f,t) for t in t ]

function show(io::IO, M::FittedModel)
	res = round(M.residual, sigdigits=4)
	print(io,"Fitted ", M.model)
	print(";  with residual $res")
end

function Dict(f::FittedModel)
	Dict( ("model"=>Dict(model(f)), "t"=>f.t, "I"=>f.I, "residual"=>f.residual) )
end

# Fit a single model to nondimensional data, giving initial values.
function fit(
	MT::Type{<:AbstractModel},
	measured::MeasuredValues,
	t::AbstractVector{<:Real},
	I::AbstractVector{<:Real},
	Tₕ::AbstractVector{<:Real};
	initpar=missing,
	# method=:LN_NELDERMEAD
	method=:LN_PRAXIS
	)

	@assert all(@. 0 ≤ t ≤ 1 + 2eps(Float32)) "Invalid time vector, $(extrema(t))."
	derived = DerivedParameters(measured)

	function misfit(fun, t, val)
		Q = 0.0
		# trapezoid rule
		y = (fun(t[1]) - val[1]) ^ 2
		x = t[1]
		for (t, I) in zip(t[2:end], val[2:end])
			ynew = (fun(t) - I) ^ 2
			Q += (t - x) * (y + ynew)
			x, y = t, ynew
		end
		return Q / 2
	end

	r = similar(t, length(t)-1)
	r̂ = similar(r)
	log_v = similar(I)
	log_v̂ = similar(I)
	function misfit_ratio(fun, t, val)
		for n in eachindex(t)
			# @show n, Ifun(h(t[n]), f(t[n]))
			log_v̂[n] = log(max(eps(), fun(t[n])))
			log_v[n] = log(val[n])
			if n > 1
				idt = 1 / (t[n] - t[n-1])
				r[n-1] = (log_v[n] - log_v[n-1]) * idt
				r̂[n-1] = (log_v̂[n] - log_v̂[n-1]) * idt
			end
		end
		# trapezoid rule
		Q = 0.0
		y = (r̂[1] - r[1]) ^ 2
		x = t[2]
		for (t, r, r̂) in zip(t[3:end], r[2:end], r̂[2:end])
			ynew = (r̂ - r) ^ 2
			Q += (t - x) * (y + ynew)
			x, y = t, ynew
		end
		return Q / 2
	end

	PT = parameter_type(MT)
	# to be minimized to find model parameters
	function objective(x, grad)
		p̂ = ExpModelParameters(x, measured)
		@debug p̂
		M = solve(MT(p̂, measured))
		sol = M.ode_solution
		# Impose a penalty for early termination.
		if successful_retcode(sol.retcode)
			Î = intensity(M)
			@debug Î(0.1)
			Q1 = misfit(Î, t, I)
			T̂ₕ = solution(M, :Th, dim=false)
			Q2 = misfit(T̂ₕ, t, Tₕ)
			return Q1 + 5Q2
		else
			return 1 / sol.t[end]
		end
	end

	# Set up optimization.
	n = length(units(PT))
	opt = Opt(method, n)
	lower, upper = bounds(PT)
	opt.lower_bounds = nondimensional(lower, measured)
	opt.upper_bounds = nondimensional(upper, measured)
	initpar = something(initpar, [lower + (m/4)*(upper - lower) for m in 1:3])
	opt.xtol_rel = 1e-3
	opt.xtol_abs = 1e-5
	opt.maxeval = 800
	opt.min_objective = objective
	# @show opt.upper_bounds

	# Optimize over all initializations.
	bestmin, bestpar = Inf, nondimensional(initpar[1], measured)
	bestret = []
	for p̂ in initpar
		p = nondimensional(p̂, measured)
		# @show objective(p, [])
		minval, minp, ret = NLopt.optimize(opt, p)
		@debug minval, minp, ret
		if minval < bestmin
			bestmin, bestpar, bestret = minval, minp, ret
		end
	end

	if bestret == :FAILURE
		@warn "optimization failed for $(typeof(M))"
	end

	return bestmin, bestpar
end
