using TFThermal, Plots, Unitful
default(linewidth=2, titlefontsize=12)

# Trial measurements
h₀ = 5 * u"μm"            # initial film thickness
ts = 20 * u"s"            # time period of observation
f̂₀ = 0.12 * 0.01          # iniital FL concentration
Tₕ = 34.0 * u"°C"         # initial air/film temperature

tp = MeasuredValues(h₀, ts, f̂₀, Tₕ)

# Tuneable parameters controlling the model
v₀ = 2 * u"1μm/minute"    # nominal evaporation rate (try between 1 and 20)
b₁ = 0.5 * u"1/s"           # initial shear rate  (try between -2 and 2)
b₂ = 0 * u"1/s"         # shear decay rate (try between 0 and 2)
Bi = 9 * 1e-4             # Biot number (probably should be << 1)

param = ExpModelParameters(v₀, b₁, b₂, Bi)

model = TFModelExp(param, tp);
model = solve(model);
h, c, fl, Tc = solution(model, dim=true)  # use dim=false for nondimensional
I = intensity(model)
Th = solution(model, :Th, dim=true)           # film surface temperature
Je = evap_rate(model, dim=true)    # evaporation rate

println("At the final time:")
println("  film thickness = ", h(1))
println("  osmolarity = ", c(1))
println("  intensity = ", I(1))
println("  film surface temperature = ", Th(1))
println("  cornea surface temperature = ", Tc(1, 0))

fig = plot(t -> ustrip(h(t)), 0, 1, label="h", size=(550,700),
    title = "film thickness", ylabel="μm", legend=false, layout=(4,1))
plot!(t -> ustrip(c(t)), 0, 1, label="c", subplot=2,
    title = "osmolarity", ylabel="mol / L", legend=false)
plot!(t -> ustrip(I(t)), 0, 1, label="I", subplot=3,
    title = "relative intensity", ylabel="", legend=false)
plot!(t -> ustrip(Th(t)), 0, 1, label="air/film", subplot=4,
    title = "temperature", ylabel="° C")
plot!(t->ustrip(Tc(t, 0)), 0, 1, label="cornea/film", subplot=4, ylabel="° C")

fig
