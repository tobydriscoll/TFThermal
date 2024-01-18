module TFThermal

# using Polynomials
using Unitful, Parameters
using FFTW, LinearAlgebra, ComponentArrays
# using NLopt
using DifferentialEquations
# using DiffEqParamEstim
# using Statistics
# using StatsBase
# using Distributions: Uniform,Normal,Exponential,InverseGamma,truncated

import Base: show, names
# import Polynomials: fit
import DifferentialEquations: solve
export MeasuredValues, TrialParameters, DerivedParameters, ExpModelParameters, TFModelExp, solution, evap_rate, dimensionalize, nondimensionalize, solve, intensity

export Chebyshev
include("chebyshev.jl")
include("constants.jl")
include("models.jl")

end # module TFThermal
