### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ a70cef7d-2a2f-4155-bdf3-fec9df94c63f
begin
    using Pkg
    Pkg.activate(joinpath(@__DIR__,".."))
    using PlutoUI, HypertextLiteral, UUIDs
    using LinearAlgebra
    using Interpolations
    using VoronoiFVM, GridVisualize, ExtendableGrids
    using LaTeXStrings
    using LessUnitful, Unitful
    using PreallocationTools
    using CairoMakie
	using LaTeXStrings
	using DoubleFloats
	using ForwardDiff 
    default_plotter!(CairoMakie)
    CairoMakie.activate!(type = "svg")
  #  TableOfContents()
end

# ╔═╡ 1e7979cc-1946-4070-9df2-9aa0924efe32
md"""
## Constants and units

See LessUnitful.jl  `ph` and `ufac` and Unitful.jl for `u`


"""

# ╔═╡ 927fac38-2d84-4ae5-984b-732f3b035420
begin
    const χ_S = 78.49-1
    const F = ph"N_A" * ph"e"
    const K = ufac"K"
    const nm = ufac"nm"
	const m=ufac"m"
    const dm = ufac"dm"
    const V = ufac"V"
    const mol = ufac"mol"
	const T=(273.15+25)*ufac"K"
	const RT = ph"R" * T
end

# ╔═╡ ce01181a-81c0-4507-9896-4e5f4aa9b8a8
md"""
## Modified Poisson - Boltzmann equation

"""

# ╔═╡ fd602033-23d2-4d0b-961d-f37c9869a89a
md"""

Let ``\bar c`` be the summary concentration, and ``y_i`` be the solute molar fractions, such that for the species molar concentrations one has ``c_i=\bar c y_i``. Then the Poisson equation states:
```math
\begin{aligned}
  -\nabla ⋅ ((1+χ)ε_0 \nabla \phi) &= F\sum_{i=1}^n z_i \bar c y_i\\
\end{aligned}
```

The mole fractions can be calculated as

```math
\begin{aligned}
   y_i&= \frac{\frac{c_{i,bulk}}{c_{0,bulk}}\exp\left(z_i\phi \frac{F}{RT}\right)}{1+ f_{model}\sum_{j=1}^n \frac{c_{j,bulk}}{c_{0,bulk}}\exp\left(z_i\phi \frac{F}{RT}\right)}
\end{aligned}
```
where ``f_{model}=1`` indicates the Bikerman model and  ``f_{model}=0`` leads to the unmodified Poisson-Boltzmann equation.
"""


# ╔═╡ 0b7f1ab4-7efa-4975-9e31-478a400ef329
md"""
Boundary conditions - fixed surface charges:
"""

# ╔═╡ 1e161882-8327-40ee-949d-c245fb63c841
md"""
``ε\partial_n \phi = \pm q``
"""

# ╔═╡ 9552ea4b-7a78-4dd4-8cd5-435e7945792b
md"""
### WIP: overview on modifications
"""

# ╔═╡ c0b5535a-2ba4-43b8-b20b-6ed017d7e519
md"""
#### Steric interaction
 - (*) Bikerman/Freise, Borukhov/Andelman/Orland
 - Dreyer/Guhlke/Landstorfer/Müller with additional account for pressure

#### Dielectric decrement
 - (*) Field strength dependency
 - Concentration dependency
 
#### Overscreening 
 - (*) Bazant/Storey/Kornyshev - add 4th-Order term
 - Hyon/Eisenberg,  Schammer et al: add repulsive potential term
 - Schammer et al: correlation with AFM experiment
 - De Souza+Bazant: modify BSK
 - Baskin+Prendergast
 - Gupta et al

"""

# ╔═╡ ada81bca-124d-40a2-9964-92fafebdcbc6
md"""
### Field dependent dielectric decrement

```math
χ=χ(E)= \frac{χ_S}{\sqrt{1+ a|E|^2}}
```
"""

# ╔═╡ f1e8c7c3-98d8-45fa-85d3-c4ebafecdd1d
md"""
Modification by Bazant/Storey/Kornyshev 2011:

```math
\begin{aligned}
ε (l_c^2 \nabla^2-1 )\nabla^2 \phi &= \rho\\
-ε\nabla^2 \phi + ε l_c^2 \nabla^4\phi &= \rho\\
- \nabla\cdot (ε \nabla \phi) + \nabla \cdot(ε l_c^2\nabla \Delta \phi) &= \rho\\
- \nabla\cdot (ε \nabla \phi) + \nabla \cdot(ε l_c^2\nabla \psi) &= \rho\\
\end{aligned}
```
with 
```math
-\Delta \phi + \psi = 0
```
$\Rightarrow$ 
```math
\begin{aligned}
- \nabla\cdot (ε \nabla \phi - εl_c^2\nabla \psi) &= \rho\\
-\Delta \phi + \psi &= 0
\end{aligned}
```

"""

# ╔═╡ d4c9d16a-98e8-44e0-bd15-f81dfc71a75d
md"""
## Parameters

"""

# ╔═╡ c2a991c2-c09e-4cbc-8595-de31b11f258c
const f_model =1

# ╔═╡ 9e5c58ed-af67-4ea4-bae8-e4b75ecd08cd
begin
	const floattype=Float64
    const L = 10.0nm # computational domain size
	const nel =20 # number of electrons/nm^2
	const M_bulk = 1 # (bulk) molarity at center of domain
	const E0=floattype(10V/nm) # decrement parameter
	const a=5.0/E0^2 # decrement parameter in χ(E)
    const q = nel * ph"e" / ufac"nm^2" # surface charge
    const c_bulk = [M_bulk, M_bulk] * mol / dm^3 # bulk  concentrations
    const z = floattype[-1, 1] # species charge numbers
    const c̄ = 55.508mol / dm^3 # summary molar concentration
    const ε_0 = floattype(ph"ε_0")

end;

# ╔═╡ a5f44f72-608c-4503-a90a-6fcbf72c2b71
const xscale=identity;  # x axis scale for graphics

# ╔═╡ b3c378db-3d96-4d67-bc81-f585160bb6b2
md"""
Check for bulk electroneutrality:
"""

# ╔═╡ 67562aff-bc9e-4d21-b569-15d6a0c58e0f
@assert dot(c_bulk, z) == 0

# ╔═╡ 6993be42-526d-47bf-8633-1dee7e0a0eab
begin
    const c0_bulk = c̄ - sum(c_bulk) # solvent bulk molar concentration
    const l_debye = sqrt((1+χ_S) * ε_0  * RT / (F^2 * c_bulk[1])) # Debye length
    const dlcap0 = sqrt(2 * (1+χ_S) * ε_0 * F^2 * c_bulk[1] / RT) # Double layer capacitance at point of zero charge (0V)
end;

# ╔═╡ bf9ccc9e-bbdb-4511-a28f-4f1802da9f35
md"""
- Debye length= $(l_debye |> x->round(x,sigdigits=5) |>u"nm")
- Double layer capacitance at zero voltage for symmetric binary electrolyte = $(dlcap0 |> x->round(x,sigdigits=5) |>u"μF/cm^2")
"""

# ╔═╡ 6ae24348-1d9b-4ec9-905d-b3087e93fd9e
md"""
### Discretization grid
"""

# ╔═╡ 53ec6dd6-7554-4d4f-bdb6-e574df37b9b1
const nref = 4 # grid refinement level

# ╔═╡ f6be79eb-626e-4a47-b82a-85fd78e0498f
begin
	const hmin = 1.0e-1 * nm * 2.0^(-nref) # grid size at working electrode
    const hmax = 1.0 * nm * 2.0^(-nref) # grid size at bulk
 
	δx=1.0e-3*nm*0 # X offset for logarithmic in x plots
    X0 = geomspace(δx, L / 2, hmin, hmax)
    X1 = geomspace(L / 2, L-δx, hmax, hmin)
    X = glue(X0, X1)
end

# ╔═╡ c5bb8f4a-d146-4f50-b9e7-06611dc825ee
begin
	grid = simplexgrid(X)
	bfacemask!(grid, [L/2], [L/2],3, tol=1.0e-10*nm)

end

# ╔═╡ 2b7f9ae7-5580-4c06-ac22-d3c0f8325193
gridplot(grid, size = (600, 200))

# ╔═╡ efdf11c6-b75c-4663-86d5-82ead82c397f
md"""
### Problem implementation
"""

# ╔═╡ 6897f195-eb28-48b4-b946-52120796fde5
const Y = DiffCache(ones(floattype,length(z))) # place for temporary data in callbacks

# ╔═╡ 6d26ffef-2456-4f63-b6d7-fc391b760cb7
function molfractions!(y, ϕ)
    N = length(z)
    for i in 1:N
        y[i] = exp(-z[i] * ϕ * F / RT) * c_bulk[i] / c0_bulk
    end
    denom = 1.0/(one(ϕ) + f_model * sum(y))
    for i in 1:N
        y[i] = y[i] * denom
    end
    return nothing
end

# ╔═╡ ffbf3c15-9397-4a07-a48f-458c963fe613
function flux!(y, u, edge, data)
	eins=one(eltype(u))
	h=floattype(edgelength(edge))
	E=(u[1, 1] - u[1, 2])/h
	χ=χ_S/sqrt(eins+a*E^2)
	ε=(eins+χ) * ε_0  
	y[1] = ε* ( (u[1, 1] - u[1, 2]) )
	return nothing
end

# ╔═╡ d0ea64ed-d351-401e-bfb7-bdf1fbc7227d
function spacecharge!(y,ϕ)
    N = length(z)
    molfractions!(y, ϕ)
    sumyz = zero(ϕ)
    for i in 1:N
        sumyz += z[i] * y[i]
    end
    return F * c̄ * sumyz
end

# ╔═╡ 4b497d98-e874-4218-b908-19b792eb9bfa
function reaction!(y, u, node, data)
	tmp=get_tmp(Y, u)
    y[1] = -spacecharge!(tmp,u[1])
    return nothing
end

# ╔═╡ b6f0a027-c4ed-4a9d-81a2-033421a2beea
function bcondition!(y, u, bnode, data)
    boundary_neumann!(y, u, bnode, species = 1, region = 2, value = -q)
    boundary_neumann!(y, u, bnode, species = 1, region = 1, value = q)
	return nothing
end

# ╔═╡ c92d5e01-6988-4cbf-bd2f-bdef7bafa86d
pbsystem = VoronoiFVM.System(
    grid;
    reaction = reaction!,
    flux = flux!,
    bcondition = bcondition!,
    species = [1],
	valuetype=floattype
)

# ╔═╡ 5ffcb2d2-26c8-4e56-8ca8-247f415d6f82
md"""
## Solution
"""

# ╔═╡ f81e55a2-a4fd-4f58-a071-b7f012368dba
sol=solve(pbsystem, inival=0.1, verbose="n", damp_initial=0.1, maxiters=1000)

# ╔═╡ ebf55449-0e3e-4a95-8010-ee6c3862cfdf
md"""
## Postprocessing
"""

# ╔═╡ 33c57d77-17cf-4dfa-8dc8-2346623dbae9
function concentrations(sol)
	n=size(sol,2)
	N=length(z)
	c=zeros(N,n)
	y=zeros(N)
	for i=1:n
		molfractions!(y,sol[1,i])
		for j=1:N
			c[j,i]=c̄*y[j]
		end
	end
	return c
end

# ╔═╡ 00e6a252-2896-40b8-a34f-55fb2780c30c
function bee!(y,ϕ)
	N=length(z)
	for i=1:N
		y[i]= RT* log(c_bulk[i] / c0_bulk)/F - z[i]*ϕ
	end
    return nothing
end

# ╔═╡ bbf14b15-c668-40a9-86b2-8275f055f05f
function bee(sol)
	n=size(sol,2)
	N=length(z)
	e=zeros(N,n)
	y=zeros(N)
	for i=1:n
	  bee!(y,sol[1,i])
		for j=1:N
			e[j,i]=y[j]
		end
	end
	return e

end

# ╔═╡ bd6b0867-c3c8-4d3c-aee6-213cd87dbfc2
 nv=nodevolumes(pbsystem)

# ╔═╡ 950a9a27-1326-40fb-9f83-1b19a9c91828
function spacecharges(sol)
	c=concentrations(sol)
	n=size(sol,2)
	nv0=copy(nv)
	nv0[n÷2+2:end].=0
	nv0[n÷2+1]/=0.5
	nvl=copy(nv)
	tmp=zeros(length(z))
	nvl[1:n÷2].=0
	nvl[n÷2+1]/=0.5
	cdens=[ spacecharge!(tmp,sol[1,i]) for i=1:n]	
	cdens⋅nv0, cdens⋅nvl
end

# ╔═╡ 43facc9c-dee8-4c59-9ff9-85588a7a1614
md"""
Compare calculated space charges with electrode charges
"""

# ╔═╡ e6ccd7bd-ade2-452c-9142-ff7f9a5bd858
spacecharges(sol)

# ╔═╡ 9fab97f7-cd40-41b7-bf7b-e9f7fa8ef8ad
(q,-q)

# ╔═╡ 3c8ef4db-b32e-43df-bb7e-2b698446c9ba
c=concentrations(sol)

# ╔═╡ d74a2d3a-1f00-44f4-8f9d-afa4ef5561b1
ε_r=χ_S./(a*((sol[1,2:end]-sol[1,1:end-1])./(X[2:end]-X[1:end-1])).^2 .+ 1).+1

# ╔═╡ 8cf198aa-d4a5-46bc-af89-fc4a3256ceb6
md"""
## Plotting
"""

# ╔═╡ 4426f9d5-b580-42e3-9a19-f1f03b33e0b1
function plotsol(sol; size=(600,400))
	c=concentrations(sol)
	cm=c[1,:]⋅nv/(mol/dm^3)/L
	cp=c[2,:]⋅nv/(mol/dm^3)/L
	c0=-(sum(c, dims=1).-c̄)
	e=bee(sol)
	
	fig=Figure(;size)
	ax1=Axis(fig[1,1]; xlabel="z/nm",ylabel="ϕ/V", 
			 title="ϕ∈$(round.(Float64.(extrema(sol[1,:])),sigdigits=3)), ε_r ∈$(round.(Float64.(extrema(ε_r)),sigdigits=3))",
			xscale)
	ylims!(ax1,(-10,10))

	ax1r=Axis(fig[1,1];
		yaxisposition=:right,
		ylabel=L"ε_r",
	    xscale
		)
	ylims!(ax1r,(0,100))
	data1=[
	lines!(ax1, X/nm,sol[1,:], color=:green,     linewidth=2),#
	lines!(ax1, X/nm,e[1,:], color=:blue,     linewidth=2, linestyle=:dot),
    lines!(ax1, X/nm,e[2,:], color=:red,     linewidth=2, linestyle=:dot),
		lines!(ax1r, X[1:end-1]/nm,ε_r, color=:pink, linewidth=3)]
	#axislegend(ax1, data1, ["ϕ", "ε_r"], position=:cb, backgroundcolor=:transparent)

	
	ax2=Axis(fig[2,1]; xlabel="z/nm",ylabel="c/(mol/L)",
			title="M_avg=$(round.((cm,cp),sigdigits=3))",
			 xscale
			)
	ylims!(ax2,(0,60))
	data2=[
	lines!(ax2, X/nm,c[1,:]/(mol/dm^3), color=:blue, linewidth=2)
	lines!(ax2, X/nm,c[2,:]/(mol/dm^3), color=:red, linewidth=2)
	lines!(ax2, X/nm,c0[1,:]/(mol/dm^3), color=:green, linewidth=2)]
	axislegend(ax2, data2, [L"c^-",L"c^+", L"c_{solvent}"], backgroundcolor=:transparent, position=:rb)
	fig
end

# ╔═╡ 33de30bb-1195-4b92-a8cf-70dfe4878755
plotsol(sol)

# ╔═╡ 483da60f-8695-45cd-8ac7-6da49405030d
md"""
## TODO
- Additional terms (Yochelis & co)
- Bulk concentration is the concentration at ϕ=0. This is not the average concentration. With the current setting, the model "generates" as much charge as necessary in order to compensate the surface charges. This compensation may be not possible with an a priori given amount of species. 
- It would be interesting to use the amount of species as a constraint.
- The limiting case is the complete depletion of the bulk in order to have
  enough charge in the boundary layer. This would result in some minimal
   initial molarity
- What is the activity coefficient here ?
- How does this translate to semiconductor theory (band gaps etc.)
- Expressions for ε, lennard - jones, see Schammer et al
- Landstorfer + Müller paper
"""

# ╔═╡ 8af12f1c-d35b-4cc9-8185-1bb5adbb69e8
html"""<hr>"""

# ╔═╡ 784b4c3e-bb2a-4940-a83a-ed5e5898dfd4
html"""<style>.dont-panic{ display: none }</style>"""

# ╔═╡ afe4745f-f9f1-4e23-8735-cbec6fb79c41
begin
    function floataside(text::Markdown.MD; top = 1)
        uuid = uuid1()
        return @htl(
            """
            		<style>


            		@media (min-width: calc(700px + 30px + 300px)) {
            			aside.plutoui-aside-wrapper-$(uuid) {

            	color: var(--pluto-output-color);
            	position:fixed;
            	right: 1rem;
            	top: $(top)px;
            	width: 400px;
            	padding: 10px;
            	border: 3px solid rgba(0, 0, 0, 0.15);
            	border-radius: 10px;
            	box-shadow: 0 0 11px 0px #00000010;
            	/* That is, viewport minus top minus Live Docs */
            	max-height: calc(100vh - 5rem - 56px);
            	overflow: auto;
            	z-index: 40;
            	background-color: var(--main-bg-color);
            	transition: transform 300ms cubic-bezier(0.18, 0.89, 0.45, 1.12);

            			}
            			aside.plutoui-aside-wrapper > div {
            #				width: 300px;
            			}
            		}
            		</style>

            		<aside class="plutoui-aside-wrapper-$(uuid)">
            		<div>
            		$(text)
            		</div>
            		</aside>

            		"""
        )
    end
    floataside(stuff; kwargs...) = floataside(md"""$(stuff)"""; kwargs...)
end;


# ╔═╡ 31013a97-fc82-4ad2-a14e-81bc93e36b8d
let
	mb=M_bulk
	floataside(md"""
#### Parameters			   
- Surface charge: q=$(q |> x->round(x,sigdigits=4)|> u"C/cm^2") 
- Model: $(f_model == 1 ? "Bikerman" : "Boltzmann/GuoyChapman")
- Dielectric decrement parameter: a=$(round(Float64(a), sigdigits=5))
- Bulk molarity: M_bulk = $(mb)
""")
end

# ╔═╡ 124fec69-a922-48f4-8449-7960fdee42be
floataside(plotsol(sol; size=(400,400)), top=200)

# ╔═╡ b8fd36a7-d8d1-45f7-b66e-df9132168bfc
# https://discourse.julialang.org/t/adding-a-restart-process-button-in-pluto/76812/5
restart_button() = html"""
<script>
	const button = document.createElement("button")

	button.addEventListener("click", () => {
		editor_state_set(old_state => ({
			notebook: {
				...old_state.notebook,
				process_status: "no_process",
			},
		})).then(() => {
			window.requestAnimationFrame(() => {
				document.querySelector("#process_status a").click()
			})
		})
	})
	button.innerText = "Restart notebook"

	return button
</script>
""";


# ╔═╡ Cell order:
# ╠═a70cef7d-2a2f-4155-bdf3-fec9df94c63f
# ╟─1e7979cc-1946-4070-9df2-9aa0924efe32
# ╠═927fac38-2d84-4ae5-984b-732f3b035420
# ╟─ce01181a-81c0-4507-9896-4e5f4aa9b8a8
# ╟─fd602033-23d2-4d0b-961d-f37c9869a89a
# ╟─0b7f1ab4-7efa-4975-9e31-478a400ef329
# ╟─1e161882-8327-40ee-949d-c245fb63c841
# ╟─9552ea4b-7a78-4dd4-8cd5-435e7945792b
# ╟─c0b5535a-2ba4-43b8-b20b-6ed017d7e519
# ╟─ada81bca-124d-40a2-9964-92fafebdcbc6
# ╟─f1e8c7c3-98d8-45fa-85d3-c4ebafecdd1d
# ╟─d4c9d16a-98e8-44e0-bd15-f81dfc71a75d
# ╠═c2a991c2-c09e-4cbc-8595-de31b11f258c
# ╠═9e5c58ed-af67-4ea4-bae8-e4b75ecd08cd
# ╟─31013a97-fc82-4ad2-a14e-81bc93e36b8d
# ╠═a5f44f72-608c-4503-a90a-6fcbf72c2b71
# ╟─b3c378db-3d96-4d67-bc81-f585160bb6b2
# ╠═67562aff-bc9e-4d21-b569-15d6a0c58e0f
# ╠═6993be42-526d-47bf-8633-1dee7e0a0eab
# ╟─bf9ccc9e-bbdb-4511-a28f-4f1802da9f35
# ╟─6ae24348-1d9b-4ec9-905d-b3087e93fd9e
# ╠═53ec6dd6-7554-4d4f-bdb6-e574df37b9b1
# ╠═f6be79eb-626e-4a47-b82a-85fd78e0498f
# ╠═c5bb8f4a-d146-4f50-b9e7-06611dc825ee
# ╠═2b7f9ae7-5580-4c06-ac22-d3c0f8325193
# ╟─efdf11c6-b75c-4663-86d5-82ead82c397f
# ╠═6897f195-eb28-48b4-b946-52120796fde5
# ╠═6d26ffef-2456-4f63-b6d7-fc391b760cb7
# ╠═ffbf3c15-9397-4a07-a48f-458c963fe613
# ╠═d0ea64ed-d351-401e-bfb7-bdf1fbc7227d
# ╠═4b497d98-e874-4218-b908-19b792eb9bfa
# ╠═b6f0a027-c4ed-4a9d-81a2-033421a2beea
# ╠═c92d5e01-6988-4cbf-bd2f-bdef7bafa86d
# ╟─5ffcb2d2-26c8-4e56-8ca8-247f415d6f82
# ╠═f81e55a2-a4fd-4f58-a071-b7f012368dba
# ╟─ebf55449-0e3e-4a95-8010-ee6c3862cfdf
# ╠═33c57d77-17cf-4dfa-8dc8-2346623dbae9
# ╠═00e6a252-2896-40b8-a34f-55fb2780c30c
# ╠═bbf14b15-c668-40a9-86b2-8275f055f05f
# ╠═bd6b0867-c3c8-4d3c-aee6-213cd87dbfc2
# ╠═950a9a27-1326-40fb-9f83-1b19a9c91828
# ╟─43facc9c-dee8-4c59-9ff9-85588a7a1614
# ╠═e6ccd7bd-ade2-452c-9142-ff7f9a5bd858
# ╠═9fab97f7-cd40-41b7-bf7b-e9f7fa8ef8ad
# ╠═3c8ef4db-b32e-43df-bb7e-2b698446c9ba
# ╠═d74a2d3a-1f00-44f4-8f9d-afa4ef5561b1
# ╟─8cf198aa-d4a5-46bc-af89-fc4a3256ceb6
# ╠═33de30bb-1195-4b92-a8cf-70dfe4878755
# ╠═124fec69-a922-48f4-8449-7960fdee42be
# ╠═4426f9d5-b580-42e3-9a19-f1f03b33e0b1
# ╠═483da60f-8695-45cd-8ac7-6da49405030d
# ╟─8af12f1c-d35b-4cc9-8185-1bb5adbb69e8
# ╟─784b4c3e-bb2a-4940-a83a-ed5e5898dfd4
# ╟─afe4745f-f9f1-4e23-8735-cbec6fb79c41
# ╟─b8fd36a7-d8d1-45f7-b66e-df9132168bfc
