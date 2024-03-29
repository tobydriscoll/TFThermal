module TFThermal

# using Polynomials
using Unitful, Parameters
using FFTW, LinearAlgebra, ComponentArrays
using NLopt
using DifferentialEquations
using SciMLBase: successful_retcode
# using DiffEqParamEstim
# using Statistics
# using StatsBase
# using Distributions: Uniform,Normal,Exponential,InverseGamma,truncated

import Base: show
# import Polynomials: fit
import DifferentialEquations: solve
export MeasuredValues, TrialParameters, DerivedParameters, ExpModelParameters, units, bounds, TFModelExp, solution, evap_rate, dimensional, nondimensional, solve, intensity
export fit, FittedModel

export Chebyshev
include("chebyshev.jl")
include("constants.jl")
include("models.jl")
include("fit.jl")

end # module TFThermal
