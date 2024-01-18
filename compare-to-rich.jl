using Unitful, TFThermal

tp = MeasuredValues(
    h₀=6.6063u"μm",
    ts=0.05*420*u"s",
    f̂₀=0.10516*0.01,
    T₀=33.5u"°C",
    )

phat = ExpModelParameters((1/6)*u"μm/s", 0.04u"1/s", 0.04u"1/s", 7.5e-4)
model = TFModelExp(phat, tp)
s = solve(model);
h, c, fl, Tc = solution(s, dim=true)
I = intensity(s)
Th = solution(s, :Th, dim=true)
