# h(t): film thickness (relative to initial)
# c(t): film osmolarity (relative to initial)

# <:Quantity types have dimensional units, <:Real types are dimensionless


# Create a function that computes FL intensity as a function of h and c.
# Depends on initial FL and Naperian constant
function intensity(C::DerivedParameters)
	ϕ = C.ϕ
	II(h, fl) = (1 - exp(-ϕ * h * fl)) / (1 + fl^2)
	I₀ = II(1, C.f₀)
	return (h, fl) -> II(h, fl) / I₀
end

#######################################################################
# Parent type for all the ODE models
#######################################################################

abstract type AbstractModel end

const CHEBN = 64

# Create an IVP for a given model instance.
function makeivp(M::AbstractModel)
	@unpack Pc, Pr, Bi, ℒ, K̃, k̃, d̃, T̃inf, T̃₀, ỹ_c, f₀, strain = M.derived
	y = Chebyshev.points(CHEBN, [ỹ_c, 0])
	Dy = Chebyshev.diffmat(CHEBN, [ỹ_c, 0])
	function ode!(du, u, p, t)
		@unpack h, c, fl, Tc = u
		@unpack Pc, Pr, Bi, ℒ, K̃, k̃, d̃, T̃inf, T̃₀, ỹ_c, strain, y, Dy = p
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
	ode_param = (; strain, Pc, Pr, Bi, ℒ, K̃, k̃, d̃, T̃inf, T̃₀, ỹ_c, Dy, y)
	Mass = Diagonal([1; 1; 1; 0; ones(CHEBN - 2); 0])
	T₀₀ = T̃₀ * (1 + Bi + ℒ / K̃) - Bi * T̃inf  # corneal surface temp in theory
	T_init = @. T₀₀ + (1 - T₀₀) * (y / ỹ_c)
	u0 = ComponentArray(h = 1., c = 1., fl = f₀, Tc = T_init)
	return ODEProblem(ODEFunction(ode!, mass_matrix=Mass), u0, (0., 1.), ode_param)
end

# Convenience access functions for all models
# constants(M::AbstractModel) = M.derived
# parameters(M::AbstractModel) = uconvert.(units(M), M.parameters)
# nondimensionalize(M::AbstractModel) = nondimensionalize(M, parameters(M))
# isknown(M::AbstractModel) = !isnothing(M.solution)
# strain(M::AbstractModel) = strain(M, nondimensionalize(M, M.parameters))

# Solve a particular model type at given dimensional parameters
import DifferentialEquations.solve
solve(M::AbstractModel, p̂::AbstractVector{<:Quantity}) = solve(typeof(M)(M, p̂))
function solve(M::AbstractModel)
	ivp = makeivp(M)
	solution = DifferentialEquations.solve(ivp, reltol=1e-7, abstol=1e-8)
	if solution.retcode != ReturnCode.Success
		@warn "Solution was terminated without success, final time = $(solution.t[end])"
	end
	return typeof(M)(M, solution)
end

# Convert nondimensional parameters and solve
solve(M::AbstractModel, p::AbstractVector{<:Real}) = solve(M, dimensionalize(M, p))

# Access solution and intensity functions for either dimensional or dimensionless time
solution(M::AbstractModel; kw...) = [solution(M, comp; kw...) for comp in (:h, :c, :fl, :Tc)]
function solution(M::AbstractModel, component::Symbol; dim=false)
	if component == :Tc
		yc = M.derived.ỹ_c
		return function(t, y)
			Tc = M.ode_solution(t, idxs=:Tc)
			T = Chebyshev.interp(Tc, 2*(yc - y) / yc - 1)
			if dim
				T = Ts + (Tb - Ts) * T
				return uconvert(u"°C", T)
			else
				return T
			end
		end
	elseif component == :Th
		@unpack Bi, ℒ, K̃, T̃inf = M.derived
		return function(t)
			h = M.ode_solution(t, idxs=:h)
			T₀ = M.ode_solution(t, idxs=:Tc)[end]
			Th = ( T₀ + Bi * h * T̃inf ) / ( 1 + ( Bi + ℒ / K̃ ) * h )
			if dim
				Th = Ts + (Tb - Ts) * Th
			end
			return uconvert(u"°C", Th)
		end
	else
		if dim
			if component == :h
				factor = M.measured.h₀
			elseif component == :fl
				factor = M.measured.f̂₀
			elseif component == :c
				factor = c₀
			end
		else
			factor = 1
		end
		return t -> factor * M.ode_solution(t, idxs=component)
	end
end
solution(M::AbstractModel, t::Real, component::Symbol; kw...) = solution(M, component; kw...)(t)

function intensity(M::AbstractModel)
	h = solution(M, :h)
	fl = solution(M, :fl)
	I = intensity(M.derived)
	return t -> I(h(t), fl(t))
end

function evap_rate(M::AbstractModel; dim=false)
	h = solution(M, :h, dim=false)
	Tc = solution(M, :Tc, dim=false)
	@unpack Pc, Pr, Bi, ℒ, K̃, k̃, d̃, T̃inf, T̃₀ = M.derived
	@unpack h₀, ts = M.measured
	if dim
		return function(t)
			Je = ( Tc(t, 0) + T̃inf * Bi * h(t)) / ( K̃ * (1 + Bi * h(t)) + ℒ * h(t) )
			return uconvert(u"μm/s", Je * h₀ / ts)
		end
	else
		return t -> ( Tc(t, 0) + T̃inf * Bi * h(t)) / ( K̃ * (1 + Bi * h(t)) + ℒ * h(t) )
	end
end

# solution(M::AbstractModel,t::Real;idxs=nothing) = M.solution(t;idxs)
# solution(M::AbstractModel) = (t;kwargs...) -> solution(M,t;kwargs...)
#
# intensity(M::AbstractModel,t) = intensity(M.derived)(solution(M,t)...)

# Compact display
# function show(io::IO,M::AbstractModel)
# 	if isempty(M.trial)
# 		print("$(typeof(M)) with unknown parameters")
# 	else
# 		print(io,"$(typeof(M)) with parameters")
# 		for (n,p,u) in zip(names(M),M.parameters,units(M))
# 			print(io," $n = $(round(u,p,digits=4)),")
# 		end
# 		print(io,"\b")
# 	end
# end

#######################################################################
# Specific model types
#######################################################################

## Model M
# Strain can relax from a nonzero value to zero
mutable struct ExpModelParameters
	v₀::typeof(1.0u"μm/s")
	b₁::typeof(1.0u"1/s")
	b₂::typeof(1.0u"1/s")
	Bi::Float64
end

function nondimensionalize(p̂::ExpModelParameters, measured::MeasuredValues)
	δ = measured.h₀
	τ = measured.ts
	return uconvert.(Unitful.NoUnits, [p̂.v₀ * τ / δ, p̂.b₁ * τ, p̂.b₂ * τ, p̂.Bi])
end

# create strain function from parameters
function TrialParameters(p̂::ExpModelParameters, meas::MeasuredValues)
	# time is always nondimensional
	ĝ(t) = p̂.b₁ * exp(-p̂.b₂ * t * meas.ts)
	return TrialParameters(v₀=p̂.v₀, Bi=p̂.Bi, ĝ=ĝ)
end

mutable struct TFModelExp <: AbstractModel
	measured::MeasuredValues
	derived::DerivedParameters
	parameters::Union{Missing, ExpModelParameters}
	ode_solution::Union{Missing, ODESolution}
end

names(::TFModelExp) = ["v₀", "b₁", "b₂", "Bi"]
units(::TFModelExp) = [u"μm/s", u"1/s", u"1/s", NoUnits]
parameters(M::TFModelExp) = uconvert.(units(M), M.parameters)

# Constructors
TFModelExp(tri::MeasuredValues) = TFModelExp(tri, DerivedParameters(tri), missing, missing)

function TFModelExp(p̂::ExpModelParameters, meas::MeasuredValues)
	con = DerivedParameters(meas, TrialParameters(p̂, meas))
	return TFModelExp(meas, con, p̂, missing)
end

function TFModelExp(M::TFModelExp, p̂::ExpModelParameters)
	TFModelExp(p̂, M.measured)
end

function TFModelExp(M::TFModelExp, sol::ODESolution)
	M.ode_solution = sol
	return M
end

# Convert to dimensional parameters
function dimensionalize(::Type{TFModelExp}, c::MeasuredValues, p::AbstractVector{<:Real})
	return p / c.ts
end

dimensionalize(M::TFModelExp, p::AbstractVector{<:Real}) = dimensionalize(TFModelExp, M.measured, p)

# Convert to nondimensional parameters
function nondimensionalize(::Type{TFModelExp}, c::MeasuredValues, p̂::AbstractVector{<:Quantity})
	return uconvert.(Unitful.NoUnits, p̂ * c.ts)
end



# function nondimensionalize(M::ModelM,p̂::AbstractVector{<:Real})
# 	c = M.derived
# 	a = ustrip(uconvert(unit(1/units(M)[1]),c.ts/c.h₀))
# 	b = ustrip(uconvert(unit(1/units(M)[2]),c.ts))
# 	return [ p̂[1]*a, p̂[2]*b, p̂[3]*b, p̂[4]*b ]
# end
