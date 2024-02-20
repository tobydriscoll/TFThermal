using Unitful, TFThermal

tp = MeasuredValues(
    h₀=6.6063u"μm",
    ts=0.05*420*u"s",
    f̂₀=0.10516*0.01,
    T₀=33.5u"°C",
    )

phat = ExpModelParameters((1/6)*u"μm/s", 0.04u"1/s", 0.04u"1/s")
model = TFModelExp(phat, tp)
s = solve(model);
h, c, fl, Tc = solution(s, dim=true)
I = intensity(s)
Th = solution(s, :Th, dim=true)

##

t_ = collect(0:0.025:1)
I_ = I.(t_)
initpar = ExpModelParameters(2.0*u"μm/s", 1.0u"1/s", 0.0u"1/s")
fit(TFModelExp, tp, t_, I_, initpar)
