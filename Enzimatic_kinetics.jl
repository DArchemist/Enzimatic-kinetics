### A Pluto.jl notebook ###
# v0.16.0

using Markdown
using InteractiveUtils

# ╔═╡ 404f1675-9b0a-49cd-82cd-16d9ee906e51
md"""
# Enzimatic kinetics homework
"""

# ╔═╡ f7a77b64-0490-4879-b48d-49c9ecf652f0
md"""
A new enzyme (E) has been isolated and purified which converts a single substrate (S) into a single product (P). A Chemical Engineer Student determined Mr (molar mass) by gel filtration as ~ 46,400. However, in SDS gel electrophoresis, a molecular mass of ~ 23 kDa was indicated for the single protein band observed. A solution of the enzyme was analyzed in the following way. The absorbance at 280 nm was found to be 0.512. A 1.00 ml portion of the same solution was subjected to amino acid analysis and was found to contain 71.3 nmol of tryptophan. N-terminal analysis on the same volume of enzyme revealed 23.8 nmol of N-terminal alanine. The following information can be useful to calculate the enzyme concentration and minimum molar mass taking account into each amino acid by separated.
"""

# ╔═╡ c21f8f81-7357-4f8c-9f81-6408368b9e5f
md"""
$[E]_{A = 280 nm} = 1.55 A^{1cm}_{280nm} [mg/ml]$
$Mr_{minimum} = \frac{m_{E}Mr_{AA}}{m_{AA}}$
"""

# ╔═╡ f98db414-36d9-480b-aef2-3e9fcd97ed05
md"""
Where $m_E$ = Enzyme mass; $m_{AA}$ = Aminoacid mass and $Mr_{AA}$ = Aminoacid molar mass.
"""

# ╔═╡ 66aeb19b-e263-4c1d-93ec-6cffcd434d8e


# ╔═╡ a732527f-2c34-4e67-9ba9-9a9681c6adb6
md"""
>1. What is the approximate molecular mass of the enzyme? Discuss this answer. Be sure to use an appropriate number of significant figures in this and other calculations.
"""

# ╔═╡ 9719e542-99c6-472f-96a1-42c5bc1dc4c0
md"""
**Answer:** Both analysis methods proposed (namely, gel filtration and SDS gel electrophoresis) present several pitfalls; including the need that the analyzed protein be of a similar shape and size from those of the standards used to calibrate the method. One should not deem the molar mass estimated by these methods as accurate. To get a better estimate of our protein’s molar mass, we proceed to compute the minimal molecular mass from both, the tryptophan content and the N-terminal alanine content.

First, we need to find the actual protein mass used during the test (we can assume that the protein is pure due to the fact that during gel electrophoresis only a single band was observed).

**Enzime concentration $[mg/ml]$** = $1.55 * 0.512$ = $0.794$ mg/ml 

**Enzime concentration $[g/ml]$** = $0.794e-3$ g/mL

**mass of enzime in 1 mL $[g]$** = $0.794e-3$ g
"""

# ╔═╡ 1904e04a-cbc7-4003-928e-2d8b36a8e037


# ╔═╡ d8e4e024-0b4e-4e6c-ab6e-353a49feb5b7
md"""
Now, there is no clue to assert that one molecule of protein contains only one single molecule of tryptophan, so the minimal molar mass of protein obtained from the tryptophan content would just be an integer fraction of the actual protein’s molar mass.

$Mr [Try] = \frac{0.794e-3 g}{ 71.3e-9 mol} \approx 11100 g/mol$
"""

# ╔═╡ 11f0ae2f-4138-4b3d-9306-e43ca002a20e


# ╔═╡ 1f3bf542-59dd-4e5a-8546-9e9df8881378
md"""
On the other hand, we know that there is only one N-terminal alanine residual per unit of protein. This fact indicates that the minimal molar mass obtained from the content of N-terminal alanine corresponds to the actual protein’s molar mass. It also needs to be an integer multiple of the minimal molar mass obtained from the tryptophan content.

$Mr[Ala] = \frac{0.794e-3 g}{23.8e-9 mol} \approx 33300 g/mol$
"""

# ╔═╡ 8ea53f56-8fe1-4b6d-ac16-b99d1f44eb3e


# ╔═╡ 973eceed-c377-4144-ae7a-1cf4ea79b383
md"""
One can then conclude that the approximate molar mass for our protein is 33300 $Da$, and also that our protein has three residuals of tryptophan in it.
"""

# ╔═╡ 65bfc707-7157-4011-a92a-b6a1d45eb6ff


# ╔═╡ b46db3e0-6ced-4be1-bff7-52bbb00ab7d6


# ╔═╡ c9278605-66a3-4088-b9e1-ff85a20599cc
md"""
>3. What is the molar extinction coefficient ε at 280 $nm$ where $A = εcl$; $A = absorbance$, $c [=] mol / liter$, and l = cell width in $cm$. Assume that all spectrophotometric measurements are made in 1.00 cm cuvettes.
"""

# ╔═╡ 8ba69bfd-3755-4725-93ec-aaa627ec741e
md"""
**Answer:** Assuming that our protein follows the Beer-Lambert law for radiation absorption, we can use the previously computed value for the protein’s molar mass to obtain the molar concentration of protein, like this:

**Molar concentration of protein** = $\frac{0.794 g/L}{(33300 g/mol)} = 2.38e-5 M$

**$\epsilon$** = $\frac{A}{cl}$$ = \frac{0.512}{(2.38e-5 M * 1 cm)} \approx 21500 L/(mol * cm)$
"""

# ╔═╡ Cell order:
# ╟─404f1675-9b0a-49cd-82cd-16d9ee906e51
# ╟─f7a77b64-0490-4879-b48d-49c9ecf652f0
# ╟─c21f8f81-7357-4f8c-9f81-6408368b9e5f
# ╟─f98db414-36d9-480b-aef2-3e9fcd97ed05
# ╟─66aeb19b-e263-4c1d-93ec-6cffcd434d8e
# ╟─a732527f-2c34-4e67-9ba9-9a9681c6adb6
# ╟─9719e542-99c6-472f-96a1-42c5bc1dc4c0
# ╟─1904e04a-cbc7-4003-928e-2d8b36a8e037
# ╟─d8e4e024-0b4e-4e6c-ab6e-353a49feb5b7
# ╟─11f0ae2f-4138-4b3d-9306-e43ca002a20e
# ╟─1f3bf542-59dd-4e5a-8546-9e9df8881378
# ╟─8ea53f56-8fe1-4b6d-ac16-b99d1f44eb3e
# ╟─973eceed-c377-4144-ae7a-1cf4ea79b383
# ╟─65bfc707-7157-4011-a92a-b6a1d45eb6ff
# ╟─b46db3e0-6ced-4be1-bff7-52bbb00ab7d6
# ╟─c9278605-66a3-4088-b9e1-ff85a20599cc
# ╟─8ba69bfd-3755-4725-93ec-aaa627ec741e
