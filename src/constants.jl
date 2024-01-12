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

## solute parameters
# fcr is critical concentation for max value in intensity at fixed d
const fcr = 0.002; # material property: critical FL concentration in mole frac
const Mw_FL = 376u"g/mol" # molecular weight of Na_2FL, g/mol
# Naperian extinction coefficient
const eps_f = 1.75e7u"1/(m*M)"
# transmittance parameter
const fcr2 = ρ*fcr/Mw_FL


# these determine the behavior of every model
@with_kw struct TrialParameters
	h₀::Unitful.Length
	ts::Unitful.Time
	f̂₀::Float64
	T₀::Unitful.Temperature
	v₀::Unitful.Velocity
end

# these are what appear in the equations
@with_kw struct ModelConstants
	Pc::Float64
	Pr::Float64
	f₀::Float64
	ϕ::Float64
	ℒ::Float64
	Bi::Float64
	T̃inf::Float64
	K̃::Float64
	k̃::Float64
	d̃::Float64
	ỹ_c::Float64
	T̃₀::Float64
end

function ModelConstants(p::TrialParameters)
	@unpack h₀, ts, f̂₀, T₀, v₀ = p
	# trial-specific?
	k = 0.68u"W/(m*K)"
	h₀ = uconvert(u"m", h₀)
	hconv = 0.9e-3k / h₀    # convective heat transfer coefficient (9e-4)

	# dimensional parameters computed for trials
	ℓ = uconvert(u"μm", (ts * σ₀ * h₀^3 / μ) ^ 0.25)
	U = ℓ/ts
	ϵ = h₀/ℓ    # aspect ratio
	V = ϵ*U
	f0p2 = f̂₀ * ρ / Mw_FL
	capK = (Tb - Ts) / (ρ * v₀)     # gets initial v0 from temperatures

	## dimensionless parameters
	# dimensionless extinction coefficient
	ϕ = uconvert(Unitful.NoUnits, eps_f*fcr2*h₀)
	# corneal permeability (non dimnl)
	Pc = P̂c * Vw * c₀ / V
	f₀ = f̂₀ / fcr

	# cornea-related
	# depth into eye for thermal problem
	ỹ_c = -5.    # 5 corneal thicknesses per Li et al
	k̃ = k / kc       # conductivity ratio
	d̃ = h₀ / dc      # depth ratio
	Pr = κc * ts / dc^2   # Prandtl number

	# evap-related
	# v = uconvert(Unitful.NoUnits, v₀ / (h₀ / ts))
	ℒ = uconvert(Unitful.NoUnits, Lm * (ρ * h₀ / ts) / (k * (Tb - Ts) / h₀))   # ratio of latent to conducted
	K̃ = uconvert(Unitful.NoUnits, capK * (ρ * h₀ / ts) / (Tb - Ts))
    T̃inf = uconvert(Unitful.NoUnits, (Tinf - Ts) / ( Tb - Ts ))

	# flow-related
	# b₁ = b₁p * ts
	# b₂ = b₂p * ts
	Bi = hconv * h₀ / k       # Biot number

	# nondimensionalize experimental temperatures
	T̃₀ = (T₀ - Ts) / (Tb - Ts)

	return ModelConstants(Pc, Pr, f₀, ϕ, ℒ, Bi, T̃inf, K̃, k̃, d̃, ỹ_c, T̃₀)
end

# function show(io::IO, c::ModelConstants)
# 	println(io,"ModelConstants for")
# 	println(io,"    h₀ = $(c.h₀), ts = $(c.ts), f̂₀ = $(c.f̂₀)")
# 	# println(io,"    Pc = $(c.Pc), f₀ = $(c.f₀), ϕ = $(c.ϕ)")
# end
