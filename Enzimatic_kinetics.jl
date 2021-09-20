### A Pluto.jl notebook ###
# v0.16.0

using Markdown
using InteractiveUtils

# ╔═╡ 02061597-15e8-431d-8994-a54ae576d0f4
using Plots

# ╔═╡ b74eaa04-7538-442c-9683-97b6e44321d1
using CurveFit

# ╔═╡ fa932dea-bd6f-4e17-91d9-0075229ba65d
using FiniteDifferences

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

# ╔═╡ a5047aa5-3381-49a8-8f74-ae367efd58e6


# ╔═╡ e4962d17-742d-407d-bf97-0ff1282243ea
md"""
A second preparation of the enzyme had an absorbance at 280 nm of 0.485. This
enzyme was diluted very carefully: 1.00 ml into 250 ml and this diluted enzyme
was used for the following experiments (I to III).
"""

# ╔═╡ 3832254e-2350-4844-a37d-34fd81c06517
md"""
**Enzime concentration after dilution**

Initial concentration [M] = $\frac{0.485}{21500 M^{-1} cm^{-1} * 1 cm} = 2.25e-5 M$

Total number of moles in 1 mL [mol] = $2.25e-8 mol$

**Concentration after dilution [M]** = $\frac{2.25e-8 mol}{0.250 L} = 9.02e-8 M$

"""

# ╔═╡ 50353ea5-e4a2-479a-b10f-ce8eaab41e56


# ╔═╡ f078dd41-f689-463d-bcd0-e362f075fd28
md"""
**Experiment I**
>A 1.00 ml portion of the diluted enzyme was added to 250 ml of buffered substrate at pH 7.0 and was mixed rapidly. The resulting initial substrate concentration [S] o was 1.000 mM. This reaction mixture was held at 25.0°C and portions were removed periodically at time t for analysis of the product P formed. The results follow. Plot [P] vs. time. What is $k_{m}$, $v_{m}$ and $k_{cat}$?
"""

# ╔═╡ 5e07759e-2c20-49f5-9b0a-46dbf6af7d43
md"""
**Answer:**

**Enzime concentration [M]** = $\frac{9.02e-8 M * 1e-3 L}{0.25 L} = 3.61e-10 M$

The following is a plot of $[P]$ vs. $t$:

"""

# ╔═╡ 24192914-fe6b-4ab8-b49b-b70458b29280
c_P = [0.104,
0.208,
0.392,
0.554,
0.695,
0.8,
0.881,
0.93,
0.96,
0.978,
0.994,
0.99952]


# ╔═╡ 332356b7-fde4-45fe-8a1d-62a5517450e6
c_S = [0.896,
0.792,
0.608,
0.446,
0.305,
0.2,
0.119,
0.07,
0.04,
0.022,
0.006,
0.00048]

	

# ╔═╡ bb92048e-ffdd-4dbb-91f8-dc4799b812b1
t = [200,
400,
800,
1200,
1600,
2000,
2400,
2800,
3200,
3600,
4400,
6000]


# ╔═╡ ed7b38c9-84f9-430c-a4e0-376f68536aeb
p = plot(t, c_P, xlabel = "t, [s]", ylabel = "[P], [mM]", legend= false)

# ╔═╡ ce5c745f-8459-442d-95f9-674a3c73d242
md"""
Now we adjust $[P]$ vs $t$ to a third-degree polynomial so we can differenciate this function with respect to $t$.
"""

# ╔═╡ 63a39f4e-6597-49bb-9156-429aa040792d
coef = poly_fit(t, c_P, 3)

# ╔═╡ daeb8a07-4847-474a-b458-722c54e1f003
md"""
We then create a function with the computed coefficients for a third-degree polynomial:
"""

# ╔═╡ 08aa1c75-8180-41ab-91b6-0b7c3757657b
c_P_fitted(t) = coef[4] * t^3 + coef[3] * t^2 + coef[2] * t + coef[1]

# ╔═╡ 2b46a058-e494-499b-87e0-745f81cd2b66
md"""
We calculate [P] from the fitted polynomial:
"""

# ╔═╡ 287d691f-4e8e-4707-b8cd-9a489b0b449a
c_P_fitted_data = c_P_fitted.(t)

# ╔═╡ 95dfe925-9338-41fc-a8c4-2b902ac5958f


# ╔═╡ ac22ff09-fad3-433b-b7ca-161b8198daf7
md"""
Now we plot both, experimental data and fitted data to confirm that $[P]$ vs $t$ has a good fitting with a third-degree polynomial. Like this:

"""

# ╔═╡ 0ec0bdaa-9cfc-4afc-918c-2c6316e045a3
plot!(p, t, c_P_fitted_data)

# ╔═╡ ba6a3a35-57ff-49e2-a0dd-480c1f9bfd16
md"""
Now, we differentiate our fitted polynomial in order to compute $\frac{d[P]}{dt}$, like this:
"""

# ╔═╡ a4f5ecfe-5dc6-4330-b5d3-3b7d1ef7eea2
dP_dt = central_fdm(5, 1).(c_P_fitted, t)

# ╔═╡ 73260c32-08ec-4876-9278-07de82d60678
md"""
Inverting the Michaelis-Menten equation produces the simpler Lineweaver-Burk linear equation:

$\frac{1}{r} = \frac{K_m}{V_{max}} * \frac{1}{[S]} + \frac{1}{Vmax}$
"""

# ╔═╡ 703e4b73-52b7-427e-8d5f-318f97158f30
md"""
We then compute $\frac{1}{r}$ and $\frac{1}{[S]}$, like this:
"""

# ╔═╡ 23e7be1b-4a86-46dc-b7c1-cd433d895b1c
__r = [1706.48,
1862.9,
2248.29,
2767.88,
3492.71,
4548.49,
6176.07,
8889.13,
13966.8,
25498.0]

# ╔═╡ 228c7a9c-4f91-4110-9f2b-54e107835183
__S = [1.11607,
1.26263,
1.64474,
2.24215,
3.27869,
5.0,
8.40336,
14.2857,
25.0,
45.4545]

# ╔═╡ 47707f72-e980-4b2a-9c2c-f226bc9d22b0
md"""
Ploting $\frac{1}{r}$ vs $\frac{1}{[S]}$ results in the following:
"""

# ╔═╡ fbdd90d6-b56b-4f78-aaeb-06fea32340ce
p2 = plot(__S, __r, seriestype = :scatter, xlabel = "1/[S], [1/mM]", ylabel = "1/r, [s/mM]", legend= false)

# ╔═╡ c172337c-be49-4f09-90b8-2f58fd32ac93


# ╔═╡ d30980b1-862e-4b1e-b3f4-137127b983a4
md"""
We then fit $\frac{1}{r}$ vs. $\frac{1}{[S]}$ to a linear curve, like this:
"""

# ╔═╡ 9b00fb90-d930-4f36-8c0a-b9b0514386d7
coef2 = linear_fit(__S, __r)

# ╔═╡ b6d6c4a5-db74-4d58-a311-b0be087841a8
md"""
The resulting linear function is:

$\frac{1}{r} = 523.162 * \frac{1}{[S]} + 1481.86$
"""

# ╔═╡ b1fca0a2-1ae4-4651-bcdd-d87964be781d
linear_fitted(__s) = coef2[2] * __s + coef2[1]

# ╔═╡ 6cc6f205-c8cb-4780-a5dd-c4a24e7ebf3f
linear_fitted_data = linear_fitted.(__S)

# ╔═╡ 9c800be4-a018-4371-a66c-33737b038dea
plot!(p2, __S, linear_fitted_data)

# ╔═╡ deea50a9-e6f7-49ad-93a7-383f0263a172
md"""
After fitting the data to a linear curve, we get the kinetic parameters $K_m$ and $V_{max}$:
"""

# ╔═╡ c0e40b15-aae5-4bc6-b296-05e43425d9e0
v_max = 1 / coef2[1]

# ╔═╡ b9c8c2ba-2f8c-40ad-9675-08ef6d93f530
k_m = coef2[2] * v_max 

# ╔═╡ 75361c13-92c3-4e1a-83eb-9c0d7fe481dd
md"""
$v_{max} = 0.000675 \frac{mM}{s}$

$k_m = 0.353 mM$
"""

# ╔═╡ 05510af9-1d2e-42af-bccd-9d683469c4e5


# ╔═╡ a8875485-7f02-4813-b5b2-85620f57c6d6
md"""
Now we have gathered enough information to compute $K_{cat}$:

$K_{cat} = \frac{V_{max}}{[E]} = \frac{0.000000675 M/s}{ 3.61e-10 M} \approx 1870 s^{-1}$
"""

# ╔═╡ 19862db5-8c9f-40ae-99e3-71996d12ae6e


# ╔═╡ f5ecae38-c8b6-420a-a05c-0813de54e742
md"""
**Experiment 2**
>In a second experiment, a series of test tubes were set up, each containing a different amount of buffered substrate at pH 7 but each in a volume of exactly 4.00 ml. The same enzyme solution used above (absorbance at 280 nm = 0.485) was diluted 2.00 ml in 250 ml as in I, then again 2.00 ml in 200 ml. Portions of 1.00 ml of this twice diluted enzyme were added at t = 0 to each of the test tubes of buffered substrate. The reaction was stopped in just 10.0 minutes by adding perchloric acid; a suitable reagent was added to provide for a colorimetric determination of the product. The results were as follows:
"""

# ╔═╡ 5e7f008a-7606-4626-9509-e4f1e1bbac24
md"""
| $[S]$ (mM)    | Amount of product ($\mu mol$/tube) |
| :-----:     | :---------:                        |
|10.0         | 2.29                               |
|5.00         | 2.18                               |
|2.50         | 2.00                               |
|1.20         | 1.69                               |
|0.80         | 1.48                               |
|0.60         | 1.31                               |
|0.40         | 1.07                               |
|0.20         | 0.686                              |
|0.10         | 0.400                              |
"""

# ╔═╡ a294ff66-bd1a-4d99-8963-e0971253cbad
md"""
>1. Plot 1/v vs 1/ [S] where v is in units of µmol per tube and [S] in millimoles/liter. Evaluate $k_{m}$ , $v_{m}$ , and $k_{cat}$ from this plot.
"""

# ╔═╡ 7d910fcb-3f3b-46d3-8680-0fbddac7ba0c
md"""
**Answer:**
"""

# ╔═╡ 71461c4c-21e6-49c2-b075-9ba5e5f32350
md"""
**Enzime concentration after dilution**

Initial concentration [M] = $\frac{0.485}{21500 M^{-1} cm^{-1} * 1 cm} = 2.25e-5 M$

Total number of moles in 2 mL [mol] = $4.5e-8 mol$

**Concentration after dilution [M]** = $\frac{4.50e-8 mol}{0.200 L} = 2.25e-7 M$

"""

# ╔═╡ 7e29a6e2-93b9-4366-a8f7-c43aff7e8029
c_S_2 = [10.0,
5.00,	
2.50,
1.20,	
0.80,	
0.60,	
0.40,	
0.20,	
0.10]

# ╔═╡ f8ed4e5c-8900-48dc-9aac-fffdc38af58c
v = [2.29,
2.18,
2.00,
1.69,
1.48,
1.31,
1.07,
0.686,
0.400]

# ╔═╡ 6b4e692f-87f9-47e0-8ba2-17c5ac307d3c
__v = 1 ./ v

# ╔═╡ 4e02c512-9603-40ba-a8ba-18c1bd934084
__S_2 = 1 ./ c_S_2

# ╔═╡ a04dd81d-6efa-409a-b36a-a150815d37ef
p3 = plot(__S_2, __v, seriestype = :scatter, legend=false, ylabel="1/v, [tube/microMol]", xlabel="1/[S], [L/mMol]")

# ╔═╡ 78b50b4c-83eb-41bc-9b53-4759ee5bf7a1
md"""
Now we fit the data to a linear curve, like this:
"""

# ╔═╡ edf80058-be22-4489-afb8-b52ae2ad4fd5
coef3 = linear_fit(__S_2, __v)

# ╔═╡ 0110e998-39f9-4ec2-ba2a-997110db4791
md"""
The resulting linear function is:

$\frac{1}{v} = 0.208 * \frac{1}{[S]} + 0.416$
"""

# ╔═╡ 7c396b13-ca52-47ea-b457-b8d64edd56fe
linear_fitted_2(__s) = coef3[2] * __s + coef3[1]

# ╔═╡ d093413f-9cf4-4314-988a-b13b51069854
linear_fitted_data_2 = linear_fitted_2.(__S_2)

# ╔═╡ 8ed99266-6d34-4497-99da-0d41806cabd7
plot!(p3, __S_2, linear_fitted_data_2)

# ╔═╡ 281daf6b-7132-42cf-8be5-2b11e96bada9
md"""
After fitting the data to a linear curve, we get the kinetic parameters $K_m$ and $v_{max}$:

"""

# ╔═╡ 01efbbd7-ffef-427d-984c-825758e9e616
v_max_2 = 1 / coef3[1]

# ╔═╡ 76504409-9674-4684-9209-18d2fe8ca7b5
k_m_2 = coef3[2] * v_max_2 

# ╔═╡ 673c3d08-e202-4045-b59c-b44c930851e8
md"""
$v_{max} = 2.40 \frac{\mu M}{tube} = 0.24 \frac{\mu M}{min}$

$k_m = 0.500 \mu M$
"""

# ╔═╡ ac6cf72b-c1ed-42d5-bcbd-3334be59b6da
md"""
Now we have gathered enough information to compute $K_{cat}$:

$K_{cat} = \frac{V_{max}}{[E]} = \frac{0.000000240 M/min}{ 2.25e-7 M} \approx 1.067 min^{-1}$
"""

# ╔═╡ a0290f9e-162f-4eb8-941c-127ddca49749


# ╔═╡ 8983e7f9-8222-41c4-a9b0-c651cb491308
md"""
>2. Plot the same data as v/[S] vs. v. Again evaluate $k_m$ and $v_m$.
"""

# ╔═╡ c2ac22b4-422d-4845-a57f-f3ea033db460
md"""
Multiplying the Lineweaver-Burk by $r$, and carrying out some arithmetical transformations, we arrive to the Eadie-Mofstee equation:

$v = v_{max} - k_m * \frac{v}{[S]}$
"""

# ╔═╡ a67adc45-81f6-402c-985b-06ddfb0a7886
v__S = v ./ c_S_2

# ╔═╡ 75e1a539-1dc6-4fdb-a829-d6554a5abac0
md"""
Then we plot $v$ vs. $\frac{v}{[S]}$, like this:
"""

# ╔═╡ 1d49c42a-6e6f-4e76-bcb9-3979b90fc3c6
p4 = plot(v__S, v, seriestype = :scatter, legend=false, ylabel="v, [microMol/tube]", xlabel="v/[S], [microMol/(mMol * tube)]")

# ╔═╡ cf6a1af5-7b29-4995-805e-b18c80877bb1
md"""
Now we fit the data to a linear curve, like this:
"""

# ╔═╡ c964ad11-733b-42cb-98fd-ee203aa9f2e1
coef4 = linear_fit(v__S, v)

# ╔═╡ d17ce2a7-8571-40ad-b517-15dd191ce6e8
md"""
The resulting linear function is:

$v = 2.40 - 0.499 * \frac{v}{[S]}$
"""

# ╔═╡ 1c8b6607-186b-4aa3-bd68-ddff821a64a3
linear_fitted_3(v__S) = coef4[1] + coef4[2] * v__S

# ╔═╡ 5f5d6b94-5a0c-42bb-a17e-8a20ed61bf7d
linear_fitted_data_3 = linear_fitted_3.(v__S)

# ╔═╡ b43edf17-204c-4fdb-82df-ac3668c016a6
plot!(p4, v__S, linear_fitted_data_3)

# ╔═╡ 8ca98a28-6efd-400d-a900-2922418d9822
md"""
After fitting the data to a linear curve, we get the kinetic parameters $K_m$ and $v_{max}$:

$v_{max} = 2.40 \frac{\mu M}{tube} = 0.240 \frac{mM}{min}$

$k_m = 0.499 \mu M$

"""

# ╔═╡ 1b0f4f50-2a89-470b-8034-ed8405b871d8


# ╔═╡ 22658bba-060e-41ea-b4c2-a86647fb56c2
md"""
**Experiment 3**
>The preceding experiment was repeated but an inhibitor was present in each tube in a concentration equal to 5.00 mM, 10.0 mM, or 25.0 mM. Two different inhibitors were used, A and B. The following results were obtained.
"""

# ╔═╡ eb72068d-8a09-4817-a254-15901abda65b
md"""
|[S], [mM] | Inhibitor A [I] = 5.00 mM | Inhibitor B [I] = 5.00 mM | Inhibitor A [I] = 25.0 mM | Inhibitor B [I] = 10 mM|
| :----: | :----: | :----: | :----: | :----: |
| 5.00   | 2.00   | 1.09   | 1.50   | 0.727  |
| 2.50   | 0.71   | 1.00   | 1.09   | 0.667  |
| 1.20   | 1.31   | 0.848  | 0.686  | 0.565  |
| 0.80   | 1.07   | 0.739  | 0.505  | 0.492  |
| 0.60   | 0.900  | 0.654  | 0.400  | 0.436  |
| 0.40   | 0.686  | 0.533  | 0.282  | 0.356  |
| 0.20   | 0.400  | 0.343  | 0.150  | 0.229  |
| 0.10   | 0.218  | 0.200  | 0.077  | 0.133  |

"""

# ╔═╡ 05bf4f24-3265-46c2-9772-8829d0d2278e
md"""
>1. Plot l/v vs 1/[S] for each of these sets of data. For each case evaluate $v_m$, apparent $k_m$, and inhibitor constants $K_I$.
"""

# ╔═╡ 4519686a-cdde-4a3e-b2d7-d84828c05527
md"""
**Answer:**
For this task we are going to create a function that carries out the needed computations on a given dataset.
"""

# ╔═╡ 55ce8dc1-1c58-4c14-91e7-4beb5da94ddd
function kinetics(v, s) 
	__v = 1 ./ v;
	__s = 1 ./ s;
	px = plot(__s, __v, seriestype = :scatter, legend=false, ylabel="1/v [tube/microMol]", xlabel="1/[S], [L/mMol]")
	coef_x = linear_fit(__s, __v);
	fit_data = coef_x[2] .* __s .+ coef_x[1];
	v_max_x = (1 / coef_x[1]) / 10;
	k_m_x = coef_x[2] * v_max_x * 10;
	plot!(px, __s, fit_data)
	return px, v_max_x, k_m_x
end
	
	
	
	

# ╔═╡ f8ba971d-6736-4b94-88c7-2998a38ce894
md"""
Now we apply this function to every dataset, like this:
"""

# ╔═╡ 31a47595-7db2-4d2c-a77d-77d1b863af78
md"""
**Inhibitor A [I] = 5.00 mM**
"""

# ╔═╡ a7c652dc-018e-4b12-9ea2-6a71c92e5f1a
plot1, v_max1, k_m1 = kinetics([2.00, 0.71, 1.31, 1.07, 0.900, 0.686, 0.400, 0.218], [5.00, 2.50, 1.20, 0.800, 0.600, 0.400, 0.200, 0.100])

# ╔═╡ 5e6414c1-0e5b-4715-be50-e84c02478437
plot1

# ╔═╡ dbab97ca-fabc-4bf2-b84b-416a89ac1172
"V_max Aparente = $v_max1 mM/min"

# ╔═╡ 51ac8189-a162-45f3-9815-d76fe7702bf5
"k_m Aparente = $k_m1 microM"

# ╔═╡ 089a0161-bdb8-431b-b5b6-1ba0a95101b5
md"""
Both, $v_{max}$ and $k_m$ changed, which indicates **Uncompetitive inhibition**. 

We can compute $K_I$ like this:

$K_I = \frac{[I]}{\frac{k_d}{k_{dapp}} - 1} = \frac{5.0 mM}{\frac{0.50}{0.67} - 1}$
"""

# ╔═╡ 67d53def-87cd-47c7-8ace-e8249dbabda0
5 / (0.5/0.67 - 1)

# ╔═╡ ece4dec4-2631-403b-ba01-f506d0c03ce1


# ╔═╡ 697f8523-25f6-4413-86c2-5ebeaacf02af
md"""
**Inhibitor B [I] = 5.00 mM**
"""

# ╔═╡ d9976a5d-06f2-4d04-b711-e64651f04ba2
plot2, v_max2, k_m2 = kinetics([1.09, 1.00, 0.848, 0.739, 0.654, 0.533, 0.343, 0.200], [5.00, 2.50, 1.20, 0.800, 0.600, 0.400, 0.200, 0.100])

# ╔═╡ 971b4a12-60f4-4901-8f06-3ea86c2700f7
plot2

# ╔═╡ 297d81f1-f798-4bc1-9cb2-625960325ec8
"V_max Aparente = $v_max2 mM/min"

# ╔═╡ 84b0982a-096f-4d13-ab75-f98341ecc25c
"k_m Aparente = $k_m2 microM"

# ╔═╡ 0f8be6a6-8997-4289-b439-ee165e94ae2a


# ╔═╡ fdf10260-362c-44cb-8729-8e6e51dd6d92
md"""
**Inhibitor A [I] = 25.0 mM**
"""

# ╔═╡ 9b12ec6b-7ab8-4138-bb94-95286bf29d1d
plot3, v_max3, k_m3 = kinetics([1.50, 1.09, 0.686, 0.505, 0.400, 0.282, 0.150, 0.077], [5.00, 2.50, 1.20, 0.800, 0.600, 0.400, 0.200, 0.100])

# ╔═╡ 5483b0b1-b46b-4d1d-a32b-cd81b9c29252
plot3

# ╔═╡ 06da0f32-2b5e-4235-b5da-9e9769cac275
"V_max Aparente = $v_max3 mM/min"

# ╔═╡ 17f16e09-6b7c-475d-ae93-cf66ee0b04d5
"k_m Aparente = $k_m3 microM"

# ╔═╡ 466e3ac3-cd11-447f-9f91-2ce4c999e9ef


# ╔═╡ d45561e3-7b58-472d-bbf7-fe157973da9d
md"""
**Inhibitor B [I] = 10 mM**
"""

# ╔═╡ 90165c8e-2e4f-42ce-bb62-f600f60e2360
plot4, v_max4, k_m4 = kinetics([0.727, 0.667, 0.565, 0.492, 0.436, 0.356, 0.229, 0.133], [5.00, 2.50, 1.20, 0.800, 0.600, 0.400, 0.200, 0.100])

# ╔═╡ 6f10c9a0-5b62-4fdb-9b9a-20b99df78c34
plot4

# ╔═╡ a72e519d-07bd-47b1-b990-009fc7c693f3
"V_max Aparente = $v_max4 mM/min"

# ╔═╡ 8b2beb4d-d265-4dc0-8371-c0179c47fe9a
"k_m Aparente = $k_m4 microM"

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CurveFit = "5a033b19-8c74-5913-a970-47c3779ef25c"
FiniteDifferences = "26cc04aa-876d-5657-8c51-4c34ba976000"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"

[compat]
CurveFit = "~0.3.5"
FiniteDifferences = "~0.12.18"
Plots = "~1.22.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "84918055d15b3114ede17ac6a7182f68870c16f7"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.1"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "f2202b55d816427cd385a9a4f3ffb226bee80f99"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+0"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "4ce9393e871aca86cc457d9f66976c3da6902ea7"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.4.0"

[[ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "9995eb3977fbf67b86d0a0a0508e83017ded03f2"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.14.0"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "4866e381721b30fac8dda4c8cb1d9db45c8d2994"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.37.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[CurveFit]]
deps = ["LinearAlgebra", "Polynomials"]
git-tree-sha1 = "eb0e0c7f3c49611223c7ed2d75ad98cfe6a9b18c"
uuid = "5a033b19-8c74-5913-a970-47c3779ef25c"
version = "0.3.5"

[[DataAPI]]
git-tree-sha1 = "bec2532f8adb82005476c141ec23e921fc20971b"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.8.0"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "7d9d316f04214f7efdbb6398d545446e246eff02"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.10"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b3bfd02e98aedfa5cf885665493c5598c350cd2f"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.2.10+0"

[[ExprTools]]
git-tree-sha1 = "b7e3d17636b348f005f11040025ae8c6f645fe92"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.6"

[[FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[FiniteDifferences]]
deps = ["ChainRulesCore", "LinearAlgebra", "Printf", "Random", "Richardson", "StaticArrays"]
git-tree-sha1 = "9a586f04a21e6945f4cbee0d0fb6aebd7b86aa8f"
uuid = "26cc04aa-876d-5657-8c51-4c34ba976000"
version = "0.12.18"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "dba1e8614e98949abfa60480b13653813d8f0157"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.5+0"

[[GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "c2178cfbc0a5a552e16d097fae508f2024de61a3"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.59.0"

[[GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "ef49a187604f865f4708c90e3f431890724e9012"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.59.0+0"

[[GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "58bcdf5ebc057b085e58d95c138725628dd7453c"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.1"

[[Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "7bf67e9a481712b3dbe9cb3dac852dc4b1162e02"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+0"

[[Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "60ed5f1643927479f845b0135bb369b031b541fa"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.14"

[[HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "8a954fed8ac097d5be04921d595f741115c1b2ad"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+0"

[[IniFile]]
deps = ["Test"]
git-tree-sha1 = "098e4d2c533924c921f9f9847274f2ad89e018b8"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.0"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[Intervals]]
deps = ["Dates", "Printf", "RecipesBase", "Serialization", "TimeZones"]
git-tree-sha1 = "323a38ed1952d30586d0fe03412cde9399d3618b"
uuid = "d8418881-c3e1-53bb-8760-2df7ec849ed5"
version = "1.5.0"

[[IterTools]]
git-tree-sha1 = "05110a2ab1fc5f932622ffea2a003221f4782c18"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.3.0"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d735490ac75c5cb9f1b00d8b5509c11984dc6943"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.0+0"

[[LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[LaTeXStrings]]
git-tree-sha1 = "c7f1c695e06c01b95a67f0cd1d34994f3e7db104"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.2.1"

[[Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a4b12a1bd2ebade87891ab7e36fdbce582301a92"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.6"

[[LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "761a393aeccd6aa92ec3515e428c26bf99575b3b"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+0"

[[Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "340e257aada13f95f98ee352d316c3bed37c8ab9"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+0"

[[Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "5a5bc6bf062f0f95e62d0fe0a2d99699fed82dd9"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.8"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[Mocking]]
deps = ["ExprTools"]
git-tree-sha1 = "748f6e1e4de814b101911e64cc12d83a6af66782"
uuid = "78c3b35d-d492-501b-9361-3d52fe80e533"
version = "0.7.2"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "3927848ccebcc165952dc0d9ac9aa274a87bfe01"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "0.2.20"

[[NaNMath]]
git-tree-sha1 = "bfe47e760d60b82b66b61d2d44128b62e3a369fb"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.5"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7937eda4681660b4d6aeeecc2f7e1c81c8ee4e2f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+0"

[[OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "15003dcb7d8db3c6c857fda14891a539a8f2705a"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.10+0"

[[Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "438d35d2d95ae2c5e8780b330592b6de8494e779"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.0.3"

[[Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "2537ed3c0ed5e03896927187f5f2ee6a4ab342db"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.0.14"

[[Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs"]
git-tree-sha1 = "4c2637482176b1c2fb99af4d83cb2ff0328fc33c"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.22.1"

[[Polynomials]]
deps = ["Intervals", "LinearAlgebra", "MutableArithmetics", "RecipesBase"]
git-tree-sha1 = "0bbfdcd8cda81b8144de4be8a67f5717e959a005"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "2.0.14"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "ad368663a5e20dbb8d6dc2fddeefe4dae0781ae8"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+0"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[RecipesBase]]
git-tree-sha1 = "44a75aa7a527910ee3d1751d1f0e4148698add9e"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.1.2"

[[RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "7ad0dfa8d03b7bcf8c597f59f5292801730c55b8"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.4.1"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[Richardson]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "e03ca566bec93f8a3aeb059c8ef102f268a38949"
uuid = "708f8203-808e-40c0-ba2d-98a6953ed40d"
version = "1.4.0"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "3240808c6d463ac46f1c1cd7638375cd22abbccb"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.12"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
git-tree-sha1 = "1958272568dc176a1d881acb797beb909c785510"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.0.0"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "8cbbc098554648c84f79a463c9ff0fd277144b6c"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.10"

[[StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "2ce41e0d042c60ecd131e9fb7154a3bfadbf50d3"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.3"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "1162ce4a6c4b7e31e0e6b14486a6986951c73be9"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.5.2"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[TimeZones]]
deps = ["Dates", "Future", "LazyArtifacts", "Mocking", "Pkg", "Printf", "RecipesBase", "Serialization", "Unicode"]
git-tree-sha1 = "6c9040665b2da00d30143261aea22c7427aada1c"
uuid = "f269a46b-ccf7-5d73-abea-4c690281aa53"
version = "1.5.7"

[[URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll"]
git-tree-sha1 = "2839f1c1296940218e35df0bbb220f2a79686670"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.18.0+4"

[[XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "cc4bf3fdde8b7e3e9fa0351bdeedba1cf3b7f6e6"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.0+0"

[[libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "c45f4e40e7aafe9d086379e5578947ec8b95a8fb"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+0"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╠═02061597-15e8-431d-8994-a54ae576d0f4
# ╠═b74eaa04-7538-442c-9683-97b6e44321d1
# ╠═fa932dea-bd6f-4e17-91d9-0075229ba65d
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
# ╟─a5047aa5-3381-49a8-8f74-ae367efd58e6
# ╟─e4962d17-742d-407d-bf97-0ff1282243ea
# ╟─3832254e-2350-4844-a37d-34fd81c06517
# ╟─50353ea5-e4a2-479a-b10f-ce8eaab41e56
# ╟─f078dd41-f689-463d-bcd0-e362f075fd28
# ╟─5e07759e-2c20-49f5-9b0a-46dbf6af7d43
# ╟─24192914-fe6b-4ab8-b49b-b70458b29280
# ╟─332356b7-fde4-45fe-8a1d-62a5517450e6
# ╟─bb92048e-ffdd-4dbb-91f8-dc4799b812b1
# ╠═ed7b38c9-84f9-430c-a4e0-376f68536aeb
# ╟─ce5c745f-8459-442d-95f9-674a3c73d242
# ╠═63a39f4e-6597-49bb-9156-429aa040792d
# ╟─daeb8a07-4847-474a-b458-722c54e1f003
# ╠═08aa1c75-8180-41ab-91b6-0b7c3757657b
# ╟─2b46a058-e494-499b-87e0-745f81cd2b66
# ╠═287d691f-4e8e-4707-b8cd-9a489b0b449a
# ╟─95dfe925-9338-41fc-a8c4-2b902ac5958f
# ╟─ac22ff09-fad3-433b-b7ca-161b8198daf7
# ╠═0ec0bdaa-9cfc-4afc-918c-2c6316e045a3
# ╟─ba6a3a35-57ff-49e2-a0dd-480c1f9bfd16
# ╠═a4f5ecfe-5dc6-4330-b5d3-3b7d1ef7eea2
# ╟─73260c32-08ec-4876-9278-07de82d60678
# ╟─703e4b73-52b7-427e-8d5f-318f97158f30
# ╟─23e7be1b-4a86-46dc-b7c1-cd433d895b1c
# ╟─228c7a9c-4f91-4110-9f2b-54e107835183
# ╟─47707f72-e980-4b2a-9c2c-f226bc9d22b0
# ╠═fbdd90d6-b56b-4f78-aaeb-06fea32340ce
# ╟─c172337c-be49-4f09-90b8-2f58fd32ac93
# ╟─d30980b1-862e-4b1e-b3f4-137127b983a4
# ╠═9b00fb90-d930-4f36-8c0a-b9b0514386d7
# ╟─b6d6c4a5-db74-4d58-a311-b0be087841a8
# ╠═b1fca0a2-1ae4-4651-bcdd-d87964be781d
# ╟─6cc6f205-c8cb-4780-a5dd-c4a24e7ebf3f
# ╠═9c800be4-a018-4371-a66c-33737b038dea
# ╟─deea50a9-e6f7-49ad-93a7-383f0263a172
# ╟─c0e40b15-aae5-4bc6-b296-05e43425d9e0
# ╟─b9c8c2ba-2f8c-40ad-9675-08ef6d93f530
# ╟─75361c13-92c3-4e1a-83eb-9c0d7fe481dd
# ╟─05510af9-1d2e-42af-bccd-9d683469c4e5
# ╟─a8875485-7f02-4813-b5b2-85620f57c6d6
# ╟─19862db5-8c9f-40ae-99e3-71996d12ae6e
# ╟─f5ecae38-c8b6-420a-a05c-0813de54e742
# ╟─5e7f008a-7606-4626-9509-e4f1e1bbac24
# ╟─a294ff66-bd1a-4d99-8963-e0971253cbad
# ╟─7d910fcb-3f3b-46d3-8680-0fbddac7ba0c
# ╟─71461c4c-21e6-49c2-b075-9ba5e5f32350
# ╟─7e29a6e2-93b9-4366-a8f7-c43aff7e8029
# ╟─f8ed4e5c-8900-48dc-9aac-fffdc38af58c
# ╠═6b4e692f-87f9-47e0-8ba2-17c5ac307d3c
# ╠═4e02c512-9603-40ba-a8ba-18c1bd934084
# ╠═a04dd81d-6efa-409a-b36a-a150815d37ef
# ╟─78b50b4c-83eb-41bc-9b53-4759ee5bf7a1
# ╠═edf80058-be22-4489-afb8-b52ae2ad4fd5
# ╟─0110e998-39f9-4ec2-ba2a-997110db4791
# ╠═7c396b13-ca52-47ea-b457-b8d64edd56fe
# ╠═d093413f-9cf4-4314-988a-b13b51069854
# ╠═8ed99266-6d34-4497-99da-0d41806cabd7
# ╟─281daf6b-7132-42cf-8be5-2b11e96bada9
# ╠═01efbbd7-ffef-427d-984c-825758e9e616
# ╠═76504409-9674-4684-9209-18d2fe8ca7b5
# ╟─673c3d08-e202-4045-b59c-b44c930851e8
# ╟─ac6cf72b-c1ed-42d5-bcbd-3334be59b6da
# ╟─a0290f9e-162f-4eb8-941c-127ddca49749
# ╟─8983e7f9-8222-41c4-a9b0-c651cb491308
# ╟─c2ac22b4-422d-4845-a57f-f3ea033db460
# ╠═a67adc45-81f6-402c-985b-06ddfb0a7886
# ╟─75e1a539-1dc6-4fdb-a829-d6554a5abac0
# ╟─1d49c42a-6e6f-4e76-bcb9-3979b90fc3c6
# ╟─cf6a1af5-7b29-4995-805e-b18c80877bb1
# ╠═c964ad11-733b-42cb-98fd-ee203aa9f2e1
# ╟─d17ce2a7-8571-40ad-b517-15dd191ce6e8
# ╠═1c8b6607-186b-4aa3-bd68-ddff821a64a3
# ╠═5f5d6b94-5a0c-42bb-a17e-8a20ed61bf7d
# ╠═b43edf17-204c-4fdb-82df-ac3668c016a6
# ╟─8ca98a28-6efd-400d-a900-2922418d9822
# ╟─1b0f4f50-2a89-470b-8034-ed8405b871d8
# ╟─22658bba-060e-41ea-b4c2-a86647fb56c2
# ╟─eb72068d-8a09-4817-a254-15901abda65b
# ╟─05bf4f24-3265-46c2-9772-8829d0d2278e
# ╟─4519686a-cdde-4a3e-b2d7-d84828c05527
# ╠═55ce8dc1-1c58-4c14-91e7-4beb5da94ddd
# ╟─f8ba971d-6736-4b94-88c7-2998a38ce894
# ╟─31a47595-7db2-4d2c-a77d-77d1b863af78
# ╠═a7c652dc-018e-4b12-9ea2-6a71c92e5f1a
# ╠═5e6414c1-0e5b-4715-be50-e84c02478437
# ╠═dbab97ca-fabc-4bf2-b84b-416a89ac1172
# ╠═51ac8189-a162-45f3-9815-d76fe7702bf5
# ╠═089a0161-bdb8-431b-b5b6-1ba0a95101b5
# ╠═67d53def-87cd-47c7-8ace-e8249dbabda0
# ╟─ece4dec4-2631-403b-ba01-f506d0c03ce1
# ╟─697f8523-25f6-4413-86c2-5ebeaacf02af
# ╠═d9976a5d-06f2-4d04-b711-e64651f04ba2
# ╠═971b4a12-60f4-4901-8f06-3ea86c2700f7
# ╠═297d81f1-f798-4bc1-9cb2-625960325ec8
# ╠═84b0982a-096f-4d13-ab75-f98341ecc25c
# ╟─0f8be6a6-8997-4289-b439-ee165e94ae2a
# ╟─fdf10260-362c-44cb-8729-8e6e51dd6d92
# ╠═9b12ec6b-7ab8-4138-bb94-95286bf29d1d
# ╠═5483b0b1-b46b-4d1d-a32b-cd81b9c29252
# ╠═06da0f32-2b5e-4235-b5da-9e9769cac275
# ╠═17f16e09-6b7c-475d-ae93-cf66ee0b04d5
# ╟─466e3ac3-cd11-447f-9f91-2ce4c999e9ef
# ╟─d45561e3-7b58-472d-bbf7-fe157973da9d
# ╠═90165c8e-2e4f-42ce-bb62-f600f60e2360
# ╠═6f10c9a0-5b62-4fdb-9b9a-20b99df78c34
# ╠═a72e519d-07bd-47b1-b990-009fc7c693f3
# ╠═8b2beb4d-d265-4dc0-8371-c0179c47fe9a
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
