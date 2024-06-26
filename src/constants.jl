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

# function Base.show(io::IO, p::MeasuredValues)
	# print(io, "Trial measurements:    ", struct_string(p))
# end

# these are needed to determine the behavior of a model
@with_kw mutable struct TrialParameters
	v₀::typeof(1.0u"μm/s") = 0.1u"μm/s"
	Bi::Float64 = 9e-4
	ĝ::Function = t -> NaN
end

# These are known from MeasuredValues only
@with_kw struct DerivedParameters
	Pc::Float64
	Pr::Float64
	f₀::Float64
	ϕ::Float64
	ℒ::Float64
	T̃inf::Float64
	k̃::Float64
	K̃::Float64
	d̃::Float64
	T̃₀::Float64
	Bi::Float64
	strain::Function
end

function DerivedParameters(p::MeasuredValues)
	# These only require measured values. Needed for defining the intensity function without knowing model parameters.
	@unpack h₀, ts, f̂₀, T₀ = p
	h₀ = uconvert(u"m", h₀)
	k = 0.68u"W/(m*K)"

	ℓ = uconvert(u"μm", (ts * σ₀ * h₀^3 / μ) ^ 0.25)
	U = ℓ / ts
	ϵ = h₀ / ℓ    # aspect ratio
	V = ϵ*U

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
    T̃inf = uconvert(Unitful.NoUnits, (Tinf - Ts) / ( Tb - Ts ))
	ℒ = uconvert(Unitful.NoUnits, Lm * (ρ * h₀ / ts) / (k * (Tb - Ts) / h₀))   # ratio of latent to conducted

	# nondimensionalize experimental temperatures
	T̃₀ = (T₀ - Ts) / (Tb - Ts)

	# a few can't be found without also knowing trial parameters
	return DerivedParameters(Pc, Pr, f₀, ϕ, ℒ, T̃inf, k̃, NaN, d̃, T̃₀, 9e-4, identity)
end

function DerivedParameters(p::MeasuredValues, trial::TrialParameters)
	@unpack h₀, ts, f̂₀, T₀ = p
	@unpack v₀, Bi, ĝ = trial
	h₀ = uconvert(u"m", h₀)

	dv = DerivedParameters(p)
	@unpack Pc, Pr, f₀, ϕ, ℒ, T̃inf, k̃, d̃, T̃₀ = dv

	# evap-related
	capK = (Tb - Ts) / (ρ * v₀)     # gets initial v0 from temperatures
	K̃ = uconvert(Unitful.NoUnits, capK * (ρ * h₀ / ts) / (Tb - Ts))

	# nondimensionalize strain function
	strain = t -> ustrip(ĝ(t) * ts)

	return DerivedParameters(Pc, Pr, f₀, ϕ, ℒ, T̃inf, k̃, K̃, d̃, T̃₀, Bi, strain)
end

function initial_thickness(f̂₀, ΔI, Imax=226.22/255)
	r = ΔI * (1 + (f̂₀ / fcr)^2) / Imax
	# convert to Molar
	f₀ʹ = f̂₀ * 0.53 / 0.2
	return uconvert(u"μm", -log1p(-r) / (eps_f * f₀ʹ * u"M"))
end
