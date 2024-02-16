const ρ = 1e3u"g/L"           # density of water, gm/liter
const P̂c = 12.1u"μm/s"
const Vw = 18e-6u"m^3/mol"    # molar volume of water, m^3/mol
const c₀ = 300u"mol/m^3"      # serum osmolarity, mOsM or Osm/m^3, same thing; 300 ok too
const σ₀ = 0.045u"N/m"
const μ = 1.3e-3u"Pa*s"

# const Pcdim = 12.1e-6u"m/s"
const dc = 5e-4u"m"           # corneal thickness in m
const kc = 0.58u"W/(m*K)"     # corneal conductivity, W/m/(deg K)
const Lm = 2.3e6u"J/kg"       # latent heat of vap, J/kg
const Ts = 20u"°C"            # saturation temp (guess); deg C
const Tb = 37u"°C"            # body temp
const Tinf = 25u"°C"          # room temp
const κc = 1.9e-7u"m^2/s"     # corneal thermal diffusivity, m^2/s
const ỹ_c = -5.      # 5 corneal thicknesses per Li et al

## solute parameters
# fcr is critical concentation for max value in intensity at fixed d
const fcr = 0.002; # material property: critical FL concentration in mole frac
const Mw_FL = 376u"g/mol" # molecular weight of Na_2FL, g/mol
# Naperian extinction coefficient
const eps_f = 1.75e7u"1/(m*M)"
# transmittance parameter
const fcr2 = ρ*fcr/Mw_FL

@with_kw mutable struct MeasuredValues
	h₀::typeof(1.0u"μm") = 4.0u"μm"
	ts::typeof(1.0u"s") = 20.0u"s"
	f̂₀::Float64
	T₀::typeof(1.0u"°C")
end

# function show(io::IO, m::MeasuredValues)
# 	print(io,"Measured trial values:")
# 	print(io,"    h₀ = $(m.h₀), ts = $(m.ts), f̂₀ = $(100m.f̂₀) %, T₀ = $(m.T₀)")
# end

# these are needed to determine the behavior of a model
@with_kw mutable struct TrialParameters
	v₀::typeof(1.0u"μm/s") = 0.1u"μm/s"
	Bi::Float64 = 9e-4
	ĝ::Function = t -> NaN
end

# these are what appear in the equations
@with_kw struct DerivedParameters
	Pc::Float64
	Pr::Float64
	f₀::Float64
	ϕ::Float64
	ℒ::Float64
	T̃inf::Float64
	K̃::Float64
	k̃::Float64
	d̃::Float64
	T̃₀::Float64
	Bi::Float64
	strain::Function
end

function DerivedParameters(p::MeasuredValues, trial::TrialParameters)
	@unpack h₀, ts, f̂₀, T₀ = p
	@unpack v₀, Bi, ĝ = trial
	k = 0.68u"W/(m*K)"
	h₀ = uconvert(u"m", h₀)

	# dimensional parameters computed for trials
	ℓ = uconvert(u"μm", (ts * σ₀ * h₀^3 / μ) ^ 0.25)
	U = ℓ/ts
	ϵ = h₀/ℓ    # aspect ratio
	V = ϵ*U
	capK = (Tb - Ts) / (ρ * v₀)     # gets initial v0 from temperatures

	## dimensionless parameters
	# dimensionless extinction coefficient
	ϕ = uconvert(Unitful.NoUnits, eps_f*fcr2*h₀)
	# corneal permeability (non dimnl)
	Pc = P̂c * Vw * c₀ / V
	f₀ = f̂₀ / fcr

	# cornea-related
	# depth into eye for thermal problem
	k̃ = k / kc       # conductivity ratio
	d̃ = h₀ / dc      # depth ratio
	Pr = κc * ts / dc^2   # Prandtl number

	# evap-related
	# v = uconvert(Unitful.NoUnits, v₀ / (h₀ / ts))
	ℒ = uconvert(Unitful.NoUnits, Lm * (ρ * h₀ / ts) / (k * (Tb - Ts) / h₀))   # ratio of latent to conducted
	K̃ = uconvert(Unitful.NoUnits, capK * (ρ * h₀ / ts) / (Tb - Ts))
    T̃inf = uconvert(Unitful.NoUnits, (Tinf - Ts) / ( Tb - Ts ))

	# nondimensionalize experimental temperatures
	T̃₀ = (T₀ - Ts) / (Tb - Ts)

	# nondimensionalize strain function
	strain = t -> ustrip(ĝ(t) * ts)

	return DerivedParameters(Pc, Pr, f₀, ϕ, ℒ, T̃inf, K̃, k̃, d̃, T̃₀, Bi, strain)
end

# function show(io::IO, c::DerivedParameters)
# 	println(io,"DerivedParameters for")
# 	println(io,"    h₀ = $(c.h₀), ts = $(c.ts), f̂₀ = $(c.f̂₀)")
# 	# println(io,"    Pc = $(c.Pc), f₀ = $(c.f₀), ϕ = $(c.ϕ)")
# end
