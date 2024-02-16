### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 63e6a76f-cc50-438a-bc6f-fd74f166d6bc
import Pkg

# ╔═╡ ed76c79e-d5c5-4106-ace9-b958dd521674
Pkg.activate("/Users/driscoll/Dropbox/research/tearfilm/thermal/TFThermal")

# ╔═╡ 6d283ad2-021a-4899-b015-847fc15f3582
using TFThermal, Plots, Unitful, PlutoUI, LaTeXStrings

# ╔═╡ d0b6473e-c4a2-4bbf-831f-a74c13a4f6d0
default(linewidth=2, titlefontsize=12)

# ╔═╡ b35cca52-eeca-4668-9475-2da756c14f70
md"""
initial film thickness, h₀ $(@bind h₀_ Slider(2:0.5:10, default=5))\
time period of observation, ts $(@bind ts_ Slider(5:5:50, default=20))\
iniital FL concentration, f̂₀ $(@bind f̂₀_ Slider(0.1:0.02:0.2, default=0.1))\
initial air/film temperature, Tₕ $(@bind Tₕ_ Slider(32:0.2:36, default=34))
"""

# ╔═╡ 2c753b5a-f0ab-47ba-8484-e566947b449f
md"""
nominal evaporation rate, v₀ $(@bind v₀_ Slider(1:1:20, default=2))\
Biot number, Bi $(@bind Bi_ Slider(1:1:20, default=9))\
initial shear rate, b₁ $(@bind b₁_ Slider(-2:0.2:2, default=0))\
shear decay rate, b₂ $(@bind b₂_ Slider(0:.1:2, default=1))
"""

# ╔═╡ 29bb0eb9-3507-4e6a-927b-d7e672487756
begin
	v₀ = v₀_ * u"1μm/minute"
	b₁ = b₁_ * u"1/s"
	b₂ = b₂_ * u"1/s"
	Bi = Bi_ * 1e-4
# 	md"""v₀ = $v₀\
# Bi = $(round(1000*Bi, digits=2)) × 10⁻³\
# b₁ = $b₁\
# b₂ = $b₂
# """
	nothing
end

# ╔═╡ be194a89-0310-4a29-9ef8-8009a28faeda
begin
	h₀ = h₀_ * u"μm"
	ts = ts_ * u"s"
	f̂₀ = f̂₀_/100
	Tₕ = Tₕ_ * u"°C"
	# md"""
	# h₀ = $h₀\
	# ts = $ts\
	# f̂₀ = $f̂₀ %\
	# Tₕ = $Tₕ
	# """
	nothing
end

# ╔═╡ 7ea5f3a7-38f0-4b95-88cb-50da5370cee7
md"""
|  |  |  |  |
|---|---|---|---|
| h₀ = $h₀ | ts = $ts | f̂₀ = $f̂₀% | Tₕ = $Tₕ |
| v₀ = $v₀ | Bi = $(round(1000*Bi, digits=2)) × 10⁻³ | b₁ = $b₁ | b₂ = $b₂ |
"""

# ╔═╡ b6654860-ae4a-449b-af70-820c217ba21a
begin
	tp = MeasuredValues(h₀,ts,f̂₀,Tₕ)
	phat = ExpModelParameters(v₀, b₁, b₂, Bi)
	model = TFModelExp(phat, tp)
	s = solve(model);
	h, c, fl, Tc = solution(s, dim=true)
	I = intensity(s)
    Th = solution(s, :Th, dim=true)
end

# ╔═╡ cced9f4e-7924-420c-b8b1-33e0f0c33ccb
begin
	plot(t->ustrip(h(t)), 0, 1, label="h", size=(500,650), 
		title = "film thickness", ylabel="μm", legend=false, layout=(4,1))
	plot!(t->ustrip(c(t)), 0, 1, label="c", subplot=2, 
		title = "osmolarity", ylabel="mol / L", legend=false)
	plot!(t->ustrip(I(t)), 0, 1, label="I", subplot=3, 
		title = "intensity", ylabel="", legend=false)
	plot!(t->ustrip(Th(t)), 0, 1, label="air/film", subplot=4, 
		title = "temperature", ylabel="° C")
	plot!(t->ustrip(Tc(t, 0)), 0, 1, label="cornea/film", subplot=4, 
		ylabel="° C")
end

# ╔═╡ 9ba69f4e-d8f2-4287-9b00-9e5b82d9f6e8
function evap_rate(t) 
	p = model.derived
	h = solution(s, :h)(t)
	T₀ = solution(s, :Tc)(t, 0)
	return ( T₀ + p.T̃inf * p.Bi * h) / ( p.K̃ * (1 + p.Bi * h) + p.ℒ * h ) 
end

# ╔═╡ a0c7fc03-ef14-4bd0-8505-a8f2ba1a97d8
evap_rate(0)

# ╔═╡ 64a030fa-bb63-4731-9de0-fee7bbe186f2
solution(s, :Tc, dim=true)(0, 0)

# ╔═╡ a448274f-6ee5-4d19-8369-bbfecfad7881
Bi

# ╔═╡ 21b91c95-6ff5-4092-919d-34d2489dd71f
Bi_

# ╔═╡ Cell order:
# ╠═63e6a76f-cc50-438a-bc6f-fd74f166d6bc
# ╠═ed76c79e-d5c5-4106-ace9-b958dd521674
# ╠═6d283ad2-021a-4899-b015-847fc15f3582
# ╠═d0b6473e-c4a2-4bbf-831f-a74c13a4f6d0
# ╟─b35cca52-eeca-4668-9475-2da756c14f70
# ╟─2c753b5a-f0ab-47ba-8484-e566947b449f
# ╟─7ea5f3a7-38f0-4b95-88cb-50da5370cee7
# ╟─cced9f4e-7924-420c-b8b1-33e0f0c33ccb
# ╟─29bb0eb9-3507-4e6a-927b-d7e672487756
# ╟─be194a89-0310-4a29-9ef8-8009a28faeda
# ╠═b6654860-ae4a-449b-af70-820c217ba21a
# ╠═9ba69f4e-d8f2-4287-9b00-9e5b82d9f6e8
# ╠═a0c7fc03-ef14-4bd0-8505-a8f2ba1a97d8
# ╠═64a030fa-bb63-4731-9de0-fee7bbe186f2
# ╠═a448274f-6ee5-4d19-8369-bbfecfad7881
# ╠═21b91c95-6ff5-4092-919d-34d2489dd71f
