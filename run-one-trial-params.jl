using TFThermal, Unitful, JLD2, HaltonSequences

# Trial measurements
h₀ = 5 * u"μm"            # initial film thickness
ts = 20 * u"s"            # time period of observation
f̂₀ = 0.12 * 0.01          # iniital FL concentration
Tₕ = 34.0 * u"°C"         # initial air/film temperature

tp = MeasuredValues(h₀, ts, f̂₀, Tₕ)
bases = [5, 11, 17, 29]
bounds = [(0.05, 40), (0, 20), (0, 20), (1e-4, 1e-3)]
t = (0:100) / 100

params = Dict{Int64, ExpModelParameters}()
intensities = Dict{Int64, Vector{Float32}}()
filmtemps = Dict{Int64, Vector{Float32}}()

# @load "trials1.jld2"

trials = 1:50000
halt = HaltonPoint(bases)

for trial in trials
    vals = Any[]
    for (x, bd, unit) in zip(halt[trial], bounds, units(ExpModelParameters))
        push!(vals, unit * (x * (bd[2] - bd[1]) + bd[1]))
    end
    p̂ = ExpModelParameters(vals...)

    model = TFModelExp(p̂, tp);
    model = solve(model);
    I = intensity(model)
    Th = solution(model, :Th, dim=true)           # film surface temperature

    params[trial] = p̂
    intensities[trial] = I.(t)
    filmtemps[trial] = ustrip.(Th.(t))
end

##
I = hcat([intensities[i] for i in trials]...)
T = hcat([filmtemps[i] for i in trials]...)
P = [params[i] for i in trials]

@save "inverse/trials5.jld2" I T P tp bases bounds t
