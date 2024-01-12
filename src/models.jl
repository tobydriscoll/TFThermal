# h(t): film thickness (relative to initial)
# c(t): film osmolarity (relative to initial)

# <:Quantity types have dimensional units, <:Real types are dimensionless


# Create a function that computes FL intensity as a function of h and c.
# Depends on initial FL and Naperian constant
function intensity(C::ModelConstants)
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
	@unpack Pc, Pr, Bi, ℒ, K̃, k̃, d̃, T̃inf, T̃₀, ỹ_c, f₀ = M.constants
	y = Chebyshev.points(CHEBN, [ỹ_c, 0])
	Dy = Chebyshev.diffmat(CHEBN, [ỹ_c, 0])
	function ode!(du, u, p, t)
		@unpack h, c, fl, Tc = u
		@unpack Pc, Pr, Bi, ℒ, K̃, k̃, d̃, T̃inf, T̃₀, ỹ_c, strain = p
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
	ode_param = (strain=strain(M), Pc, Pr, Bi, ℒ, K̃, k̃, d̃, T̃inf, T̃₀, ỹ_c, Dy, y)
	Mass = Diagonal([1; 1; 1; 0; ones(CHEBN - 2); 0])
	T_00 = T̃₀ * (1 + Bi + ℒ / K̃) - Bi * T̃inf  # corneal surface temp in theory
	T_init = @. T_00 + (1 - T_00) * (y / ỹ_c)
	u0 = ComponentArray(h = 1., c = 1., fl = f₀, Tc = T_init)
	return ODEProblem(ODEFunction(ode!, mass_matrix=Mass), u0, (0., 1.), ode_param)
end

# Convenience access functions for all models
# constants(M::AbstractModel) = M.constants
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
		yc = M.constants.ỹ_c
		return function(t, y)
			Tc = M.ode_solution(t, idxs=:Tc)
			T = Chebyshev.interp(Tc, 2*(yc - y) / yc - 1)
			if dim
				T = Ts + (Tb - Ts) * T
			end
			return uconvert(u"°C", T)
		end
	elseif component == :Th
		@unpack Bi, ℒ, K̃, T̃inf = M.constants
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
				factor = M.trial.h₀
			elseif component == :fl
				factor = M.trial.f̂₀
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
	I = intensity(M.constants)
	return t -> I(h(t), fl(t))
end

# solution(M::AbstractModel,t::Real;idxs=nothing) = M.solution(t;idxs)
# solution(M::AbstractModel) = (t;kwargs...) -> solution(M,t;kwargs...)
#
# intensity(M::AbstractModel,t) = intensity(M.constants)(solution(M,t)...)

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
mutable struct TFModel3 <: AbstractModel
	trial::TrialParameters
	constants::ModelConstants
	parameters::Union{Missing, AbstractVector{<:Quantity}}
	strain::Function
	ode_solution::Union{Missing, ODESolution}
end

names(::TFModel3) = ["b₁","b₂"]
units(::TFModel3) = [u"1/s", u"1/s"]
parameters(M::TFModel3) = uconvert.(units(M), M.parameters)
strain(M::TFModel3) = M.strain

# Constructors
TFModel3(tri::TrialParameters) = TFModel3(tri, ModelConstants(tri), missing, t -> NaN, missing)
function TFModel3(M::TFModel3, p̂::AbstractVector{<:Quantity})
	TFModel3(p̂, M.trial)
end
function TFModel3(M::TFModel3, sol::ODESolution)
	M.ode_solution = sol
	return M
end
function TFModel3(p̂::AbstractVector{<:Quantity}, tri::TrialParameters)
	con = ModelConstants(tri)
	p = nondimensionalize(TFModel3, tri, p̂)
	strain(t) = p[1] * exp(-p[2] * t)
	return TFModel3(tri, con, p̂, strain, missing)
end

# Convert to dimensional parameters
function dimensionalize(::Type{TFModel3}, c::TrialParameters, p::AbstractVector{<:Real})
	return p / c.ts
end

dimensionalize(M::TFModel3, p::AbstractVector{<:Real}) = dimensionalize(TFModel3, M.trial, p)

# Convert to nondimensional parameters
function nondimensionalize(::Type{TFModel3}, c::TrialParameters, p̂::AbstractVector{<:Quantity})
	return uconvert.(Unitful.NoUnits, p̂ * c.ts)
end

nondimensionalize(M::TFModel3, p̂::AbstractVector{<:Quantity}) = nondimensionalize(TFModel3, M.trial, p̂)


# function nondimensionalize(M::ModelM,p̂::AbstractVector{<:Real})
# 	c = M.constants
# 	a = ustrip(uconvert(unit(1/units(M)[1]),c.ts/c.h₀))
# 	b = ustrip(uconvert(unit(1/units(M)[2]),c.ts))
# 	return [ p̂[1]*a, p̂[2]*b, p̂[3]*b, p̂[4]*b ]
# end
