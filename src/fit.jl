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
	t::AbstractVector
	I::AbstractVector
	residual::AbstractFloat
end

function FittedModel(m::AbstractModel, t::AbstractVector{<:Real}, I::AbstractVector)
	inten = intensity(m)
	y = [(I-inten(t))^2 for (t,I) in zip(t,I)]
	misfit = sqrt(trap(t,y))
	return FittedModel(m,1,0,1,t,I,misfit)
end

# create a fitted model from dimensional data
function FittedModel(m::AbstractModel, t::AbstractVector{<:Quantity}, I::AbstractVector{<:Real})
	t₀ = t[1]
	ts = t[end]-t₀
	I₀ = I[1]
	c = constants(m)
	c = ModelConstants(c.h₀,ts,c.f̂₀)
	m = typeof(m)(c,m.parameters,m.solution)
	# get residual
	tn = (t.-t₀)/ts
	In = I/I₀
	model_n = FittedModel(m,tn,In)
	return FittedModel(m,ts,t₀,I₀,t,I,model_n.residual)
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

function show(io::IO,M::FittedModel)
	res = round(M.residual,sigdigits=4)
	print(io,"Fitted ",M.m)
	print(";  with residual $res")
end

function Dict(f::FittedModel)
	Dict( ("model"=>Dict(model(f)), "t"=>f.t, "I"=>f.I, "residual"=>f.residual) )
end

function make_optim_ivp(MT::Type{<:AbstractModel}, measured)
	y = Chebyshev.points(CHEBN, [ỹ_c, 0])
	Dy = Chebyshev.diffmat(CHEBN, [ỹ_c, 0])
	@unpack Pc, Pr, ℒ, k̃, d̃, T̃inf, T̃₀, f₀, Bi = DerivedParameters(measured)
	function ode!(du, u, p, t)
		@unpack h, c, fl, Tc = u
		dimen = dimensional(MT, p, measured)
		trial = TrialParameters(dimen, measured)
		@unpack Pc, Pr, Bi, ℒ, K̃, k̃, d̃, T̃inf, T̃₀, f₀, strain = DerivedParameters(measured, trial)
		T₀ = Tc[end]  # temp at corneal surface
		Th = ( T₀ + Bi * h * T̃inf ) / ( 1 + ( Bi + ℒ / K̃ ) * h )  # TF surface temp
		Je = ( T₀ + T̃inf * Bi * h) / ( K̃ * (1 + Bi * h) + ℒ * h ) # evap rate

		tmp = Je - Pc * (c-1)
		du.h = -strain(t) * h - tmp
		du.c = tmp * c / h
		du.fl = tmp * fl / h

		Tcʹ = Dy * u.Tc
		du.Tc .= Pr * (Dy * Tcʹ)
		du.Tc[1] = Tc[1] - 1
		du.Tc[end] = Tcʹ[end] + (k̃ / d̃) * ( ℒ * Je + Bi * (Th - T̃inf) )
		return du
	end
	Mass = Diagonal([1; 1; 1; 0; ones(CHEBN - 2); 0])
	T₀₀ = T̃₀ * (1 + Bi + ℒ / K̃) - Bi * T̃inf  # corneal surface temp in theory
	T_init = @. T₀₀ + (1 - T₀₀) * (y / ỹ_c)
	u0 = ComponentArray(h = 1., c = 1., fl = f₀, Tc = T_init)
	# p = nondimensional(MT, parameters(MT,measured))
	return ODEProblem(ODEFunction(ode!, mass_matrix=Mass), u0, (0., 1.))
end

# Fit a single model to nondimensional data, giving initial values.
function fit(
	MT::Type{<:AbstractModel},
	measured::MeasuredValues,
	t::AbstractVector{<:Real},
	I::AbstractVector{<:Real},
	initpar;
	method=:LN_NELDERMEAD
	)

	@assert all(@. 0 ≤ t ≤ 1) "Invalid time vector."
	derived = DerivedParameters(measured)
	inten = intensity(derived)

	ivp = make_optim_ivp(MT, measured)
	function misfit(Ifun, h, f, t, I)
		Q = 0.0
		# trapezoid rule
		y = (Ifun(h(t[1]), f(t[1])) - I[1]) ^ 2
		x = t[1]
		for (t, I) in zip(t[2:end], I[2:end])
			ynew = (Ifun(h(t), f(t)) - I) ^ 2
			Q += (t - x) * (y + ynew)
			x, y = t, ynew
		end
		return Q / 2
	end

	# Impose a penalty for early termination.
	function misfit_fail(sol)
		if sol.retcode==:Success
			h = t -> sol(t, idxs=:h)
			f = t -> sol(t, idxs=:f)
			return misfit(inten, h, f, t, I)
		else
			return 1/sol.t[end]
		end
	end

	lossfun = build_loss_objective(
		ivp,
		Tsit5(),
		misfit_fail,
		maxiters = 2000,
		verbose=false,
		verbose_opt=false,
		verbose_steps=10,
		reltol=1e-8,
		abstol=1e-9,
		callback=solvercb(inten)
		)

		# Set up optimization.
		n = length(names(M))
		opt = Opt(method, n)
		PT = parameter_type(MT)
		lower, upper = bounds(PT)
		opt.lower_bounds = nondimensional(lower, measured)
		opt.upper_bounds = nondimensional(upper, measured)
		opt.xtol_rel = 1e-5
		opt.xtol_abs = 1e-7
		opt.maxeval = 10000
		opt.min_objective = lossfun

		# Optimize over all initializations.
		bestmin,bestpar = Inf, nondimensional(initpar[1], measured)
		bestret = []
		for p̂ in initpar
			p = nondimensional(p̂, measured)
			minval, minp, ret = NLopt.optimize(opt, p)
			#(typeof(M)==ModelM) && (@show minval,p,minp,ret)
			if minval < bestmin
				bestmin, bestpar, bestret = minval, minp, ret
			end
	end

	if bestret == :FAILURE
		@warn "optimization failed for $(typeof(M))"
	end

	#@show typeof(M),bestmin,bestpar,dimensionalize(M,bestpar)

	return bestmin, bestpar #FittedModel(solve(M,bestpar),t,I)
end

# Fit all the models, using simpler ones to help initialize the more complex ones.
# function fit(con::ModelConstants,t::AbstractVector{<:Real},I)

# 	function safe(m::AbstractModel)
# 		# get parameters that are safely pushed away from the bounds
# 		p̂ = parameters(m)
# 		p̂ = min.(p̂,0.98*upper(m))
# 		return max.(p̂,0.98*lower(m))
# 	end

# 	# Start with model O
# 	init = [ [v*u"μm/minute"] for v in [1;5:5:35] ]
# 	modO = fit(ModelO(con),t,I,init)

# 	# Model F
# 	init = [ [2u"μm/minute",0.06u"s^-1"],
# 		[2u"μm/minute",-0.06u"s^-1"],
# 		[10u"μm/minute",0.06u"s^-1"],
# 		[10u"μm/minute",-0.06u"s^-1"],
# 		[0u"μm/minute",-0.06u"s^-1"],
# 		[0u"μm/minute",0.06u"s^-1"],
# 		]
# 	# Use O result to initialize model F
# 	p̂ = safe(modO.m)
# 	append!(init,[
# 		[p̂[1],0u"s^-1"], [p̂[1],0.06u"s^-1"], [p̂[1],-0.06u"s^-1"],
# 	 ] )
# 	modF = fit(ModelF(con),t,I,init)

# 	# Model D
# 	init = [
# 		[1u"μm/minute",0.2u"s^-1",0.2u"s^-1"],
# 		[6u"μm/minute",0.2u"s^-1",0.2u"s^-1"],
# 		[0u"μm/minute",0.1u"s^-1",0.2u"s^-1"],
# 		[15u"μm/minute",0.1u"s^-1",0.2u"s^-1"],
# 		[0u"μm/minute",-0.1u"s^-1",0.2u"s^-1"],
# 		]
# 	p̂ = safe(modO.m)
# 	append!(init,[
# 		[p̂[1],0u"s^-1",0u"s^-1"], [p̂[1],0.1u"s^-1",0.1u"s^-1"], [p̂[1],-0.1u"s^-1",0.1u"s^-1"],
# 	 ] )
# 	p̂ = safe(modF.m)
# 	append!(init,[
# 		 [p̂[1],p̂[2],0u"s^-1"], [p̂[1],p̂[2],0.5u"s^-1"],
# 	  ] )
# 	modD = fit(ModelD(con),t,I,init)

# 	# Model M
# 	init = [
# 		[1u"μm/minute",0.1u"s^-1",0.2u"s^-1",0.5u"s^-1"],
# 		[6u"μm/minute",-0.1u"s^-1",0.2u"s^-1",0.5u"s^-1"],
# 		[15u"μm/minute",0.1u"s^-1",0.2u"s^-1",0.5u"s^-1"],
# 		[24u"μm/minute",0.1u"s^-1",-0.2u"s^-1",0.8u"s^-1"],
# 		[20u"μm/minute",0.1u"s^-1",-0.2u"s^-1",0.8u"s^-1"],
# 		[0u"μm/minute",0.2u"s^-1",0.2u"s^-1",0.5u"s^-1"],
# 		]
# 	p̂ = safe(modO.m)
# 	append!(init,[
# 		[p̂[1],0u"s^-1",0u"s^-1",0u"s^-1"], [p̂[1],0.2u"s^-1",-0.1u"s^-1",0.5u"s^-1"],
# 	] )
# 	p̂ = safe(modF.m)
# 	append!(init,[
# 		 [p̂[1],p̂[2],0u"s^-1",0u"s^-1"], [p̂[1],p̂[2],-0.1u"s^-1",0.5u"s^-1"],
# 	] )
# 	p̂ = safe(modD.m)
# 	append!(init,[
# 		[p̂[1],0u"s^-1",p̂[2],p̂[3]],
# 	] )
# 	modM = fit(ModelM(con),t,I,init)

# 	# Re-try model F with what D and M found
# 	init = [ safe(modF.m) ]
# 	p̂ = safe(ModelF(con,modD.m.parameters[1:2]))
# 	push!(init,p̂)
# 	pM = modM.m.parameters[1:3]
# 	p̂ = safe(ModelF(con,[pM[1],pM[2]+pM[3]]))
# 	push!(init,p̂)
# 	try
# 		modF = fit(ModelF(con),t,I,init)
# 	catch
# 		@warn "There was a problem with final model F fit"
# 	end

# 	return (O=modO,F=modF,D=modD,M=modM)
# end

# function fit(con::ModelConstants,t::AbstractVector{<:Quantity},I)
# 	# Rescale time and intensity.
# 	ts = t[end]-t[1]
# 	tt = (t.-t[1])/ts
# 	II = I/I[1]
# 	cc = ModelConstants(con.h₀,ts,con.f̂₀)
# 	return fit(cc,tt,II)
# end