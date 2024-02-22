using Unitful, TFThermal

tp = MeasuredValues(
    h₀=6.6063u"μm",
    ts=0.05*420*u"s",
    f̂₀=0.100516*0.01,
    T₀=33.5u"°C",
    )

phat = ExpModelParameters((20)*u"μm/minute", -0.2u"1/s", 0.8u"1/s", 9e-4)
model = TFModelExp(phat, tp)
s = solve(model);
h, c, fl, Tc = solution(s, dim=true)
I = intensity(s)
Th = solution(s, :Th, dim=false)

##

t_ = collect(0:0.025:1)
I_ = I.(t_)
Th_ = Th.(t_)
initpar = ExpModelParameters(2*u"μm/minute", 0u"1/s", 0u"1/s", 6e-4)
bestmin, bestparam = fit(TFModelExp, tp, t_, I_, Th_, [initpar], method=:LN_NELDERMEAD)
bestmod = solve(TFModelExp(ExpModelParameters(bestparam, tp), tp))
