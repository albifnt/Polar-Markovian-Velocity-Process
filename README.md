# Polar-Markovian-Velocity-Process
<p align="justify">
Groundwater-storing subsurface formations result from sedimentation processes,
which produce soil patterns characterized by high spatial variability. The expensive
and difficult measurement procedure makes measurements of soil properties available
at few locations, making the deterministic characterization of aquifers an elusive
goal. To assess the resulting uncertainty affecting the transport of tracers in subsurface flows, 
Monte Carlo (MC) methods are typically applied. As an alternative to the
MC simulation, a recently developed stochastic modeling approach based on Polar
Markovian Velocity Processes (PMVP) allows to obtain accurate transport statistics
at a computational expense that is roughly 1000 times lower than the conventional
MC. The computationally expensive MC step of solving the Darcy flow problem for
a large number of realizations is omitted in PMVP method, by attempting to model
tracer transport subjecting the particle path in the mean field to a random walk.
Previous studies have shown that the overall spreading of the expected value of
tracer concentration is well-predicted by the PMVP method, both in homogeneous
and inhomogeneous settings. However, the conditioning of transmissivity fields on
low-value conductivity measurements leads to a distinct local misestimation of the
concentration, affecting the statistics of the flow. The aim of this project is to analyse
the velocity statistics on highly heterogeneous and statistically inhomogeneous domains by 
means of a light-weight parametrization of the PMVP model. Alternative
computationally cheap approaches like Low Order Approximation (LOA) and interpolating functions 
are investigated for the estimation of velocity statistical moments.
Based on their performance, each model is affected by singular weaknesses in specific
conditioned log-conductivity scenarios, but the option given by interpolation represents a viable solution.
</p>

## Darcy Flow in the Porous Subsurface
<p align="justify">
The physics governing groundwater flows in aquifers has been known since the
19th century, however unfortunately it was mostly studied through a deterministic
approach. The heterogeneity of hydrogeologic parameters is so complex and difficult
to describe quantitatively that even when the physical process is well understood, as
the Darcian process of single-phase incompressible flows, the deterministic prediction
of the advective transport of a solute plume is impossible to define.

The basic equation governing the fluid flow in porous media is Darcy’s Law:
</p>

<p align="center">
$` q(x, t) = −K(x)∇h(x, t) `$
</p>

<p align="justify">
Darcy found a linear relation between the gradient of the hydraulic head <b>h(x, t)</b> and
the specific discharge <b>q(x, t)</b> based on the hydraulic conductivity <b>K(x)</b>. This linear
relation holds as long as the flow in the porous media is laminar; this condition
is satisfied for most subsurface flows [1]. If the hydraulic head and the conductivity are 
known at each point in space, Darcy's Law allows to determine
the flow field <b>U(x, t) = q(x, t)/n(x)</b> by means of the spatially dependent porosity
<b>n(x)</b>. Subjecting the specific discharge to the continuity equation for incompressible
single-phase flow:
</p>

<p align="center">
$` ∇ · q(x, t) = s(x) `$
</p>

<p align="justify">
and thus:
</p>

<p align="center">
$` ∇ · ( K(x)∇h(x, t) ) = s(x). `$
</p>

<p align="justify">
This is an inhomogeneous elliptic Partial Differential Equation (PDE) for the hydraulic head. 
The right-hand side takes into account the possible re- and discharge of wells. In the absence 
of sources, the divergence of the discharge field reduces to zero.

The general advection-dispersion equation in an Eulerian reference frame reads
</p>

<p align="center">
$` \frac{∂C}{∂t} + \frac{∂}{∂xi}(uiC) − \frac{∂}{∂xj}(Dij \frac{∂Ci}{∂xi}) = 0 `$
</p>

<p align="justify">
written using the Einstein notation, where <b>i, j = 1, . . . , d</b> represent the dimensions.
Here <b>C(x, t)</b> is the concentration of a passive scalar, while the <b>Dij</b> represents the
dispersive tensor accounting for the Pore-Scale-Dispersion contribution (the combined effect of 
mechanical dispersion and molecular diffusion).

If soil properties were known at each point in space, the theory outlined would be enough to 
characterize subsurface flows. However, measurements of hydraulic conductivity values in groundwater 
aquifers have revealed that K-values typically vary over several orders of magnitude at different 
spatial locations [2]. For this reason, common interpolation algorithms cannot be applied and the
resulting field is uncertain. This uncertainty is translated through Darcy's Law equation into 
uncertainty of the velocity field <b>U(x, t)</b> and in turn into uncertainty of the concentration 
<b>C(x, t)</b>. As a consequence, a probabilistic approach is required.
</p>

## The Polar Markovian Velocity Process Model Derived from Perturbation Theory
<p align="justify">
The starting point of the polar reformulation of the Markovian Process is the
temporal RW model derived from perturbation theory in [3]. The model
is defined in terms of linear stochastic differential equations (SDEs) for the Cartesian 
velocity components
</p>

<p align="center">
$` dU_1 = - \frac{8U}{15l_Y}(U_1 - U)dt + \sqrt{ 2 \frac{U^3 σ^2_Y}{5l_Y} }(dW) `$
</p>
<p align="center">
$` dU_2 = - (\frac{16U}{15l_Y}U_2 + ω^2_0 y)dt + \sqrt{ 2 \frac{2 U^3 σ^2_Y}{15l_Y} }(dW) `$
</p>

<p align="justify">
where <b>dy = U2(t)dt</b>, <b>ω0 = 0.338 · 16U/(15lY )</b>, <b>U</b> is the mean flow velocity in the
x1-direction and <b>dW</b> is a Wiener Process [4].  Despite these equations showed an accurate agreement 
in both longitudinal and transverse macrodispersion, further inspection revealed that the temporal evolution of
both velocity components along particle path lines is strongly correlated, especially
for large $` σ^2_Y `$. If we switch to a polar coordinate system, where the particle 
displacement is defined by the velocity magnitude and the directional angle, instead
of the classical Cartesian velocity components, the problem of statistical dependence
vanishes, regardless of the conductivity field heterogeneity. Manipulating and re-writing the SDEs
for the angle <b>θ</b> and the non-dimensional velocity magnitude <b>U'</b>, we obtain
</p>

<p align="center">
$` dU' = - \frac{8}{15}(U' - 1)dt' + \sqrt{ \frac{2}{5} }σ_Y dW' `$
</p>
<p align="center">
$` θ = - (\frac{16}{15}θ + (ω_0 \frac{l_Y}{U})^2 y') dt + \frac{2 σ_Y}{\sqrt{15}} dW' `$
</p>

<p align="justify">
where $` dy' = U' 0 dt' `$, $` ω_0 l_Y / U = 0.338 · 16/15 ≈ 0.361 `$ and the non-dimensional
quantities $` t' = t U / l_Y `$, $` y' = y / l_Y `$ and the Wiener Process. Without loss of the 
generality, the polar reformulation considers the natural logarithm of the non-dimensional 
velocity $` ν = ln('U) `$. This allows to amplify and better resolve the low-velocity zone, 
where instantaneous velocities are much more strongly correlated compared to high-velocity zone [5].
Therefore, by applying Ito’s formula, the SDE becomes
</p>

<p align="center">
$` dν = - \frac{8}{15} ν dt' + \sqrt{ \frac{2}{5} }σ_Y dW' `$
</p>

## Structure of the codes
<p align="justify">
This section explains how the codes are organized in the different folders. Only relevant codes are included.
The folder “Initial_validation” contains codes used for the validation of the PMVP model:
- “PMVP_model_codes_for_validation_on_statistically_homogeneous_domains” contains the codes used to
plot the temporal evolutions of particle position variances on statistically homogeneous domains (the v-theta
model and the skew-normal RW model).
- “PMVP_model_codes_for_validation_on_statistically_inhomogenous_domains” contains the codes used for
the validation cases on conductivity field with a spatially dependent mean velocity field and with a spatially
dependent log-conductivity variance distributions.

The folder “LOA_(Comparison_between_MC,_LOA_and_PMVP_results)” contains the MATLAB scripts and functions
used to perform the Low Order Approximation. Files “Y_mean.txt” and “Y_variance.txt” are included; they provide the
distributions of the conditional mean and conditional variance of the hydraulic conductivity, for a field conditioned on
a low-value conductivity constraint at the center of the domain. Before running the LOA simulation it is necessary to
produce the file “C_matrix.txt” by means of the script “COVARIANCE_MATRIX.m”. Once the covariance matrix is
generated, the “LOA_simulation.m” can be performed. Plots of the hydraulic head mean and variance, velocity means
and variances, velocity field, conditional mean and variance distributions are generated. It's important to underline that LOA
distributions in the report were obtained with the “SME simulator developed by Dr. Liyong Li of Chevron Energy
Technology Company” which was extensively validated. This SME simulator, unfortunately, could not be provided. The
MATLAB code “LOA_simulation.m” provides distributions in agreement with “SME simulator developed by Dr. Liyong
Li of Chevron Energy Technology Company” as far as the hydraulic head mean, hydraulic head variance, velocity
means are concerned (A picture provides the comparison between hydraulic head mean and variance distributions
obtained with Monte Carlo and with “LOA_simulation.m”). However, distributions of the velocity variances become
more unreasonable in specific log-conductivity scenarios.

The folder “MC_(Comparison_between_MC,_LOA_and_PMVP_results)” contains the python code (MC_simulation.py)
to solve the Darcy Flow problem for multiple realizations of the conductivity field. It also contains the MATLAB file
“Elaborate_MC_results.m” which calculates statistical moments resulting from the Monte Carlo simulation.

The folder “PMVP_model_(Comparison_between_MC,_LOA_and_PMVP_results)” contains the file
“PMVP_particle_generator.py” which allows to generate the particle path line realizations which are then used to
compute statistical moments at each grid point of the computational domain. The MATLAB file
“Numerical_moments.m” takes as input the .txt files generated with “PMVP_particle_generator.py” and computes
the PMVP statistical moments.

The folder “Interpolating_funtions_(Comparison_between_MC,_LOA_and_PMVP_results)” contains the MATLAB
script to generate the distributions of the Cartesian velocity variances with the interpolating functions.
</p>

## References
<p align="justify">
[1] Fanchi, J. R. (2010). 4 - porosity and permeability. In Integrated Reservoir Asset Management, pages 49 – 69. Gulf Professional Publishing, Boston.

[2] Gelhar, L. (1993). Stochastic Subsurface Hydrology. Prentice-Hall.

[3] Meyer, D. W. (2017a). Relating recent random walk models with classical perturbation theory for dispersion predictions in the heterogeneous porous subsurface. Advances in Water Resources, 105:227 – 232.

[4] Gardiner, C. W. (2004). Handbook of stochastic methods for physics, chemistry and the natural sciences, volume 13 of Springer Series in Synergetics. Springer-Verlag, Berlin, third edition.

[5] Le Borgne, T., de Dreuzy, J.-R., Davy, P., and Bour, O. (2007). Characterization of the velocity field organization in heterogeneous media by conditional correlation. Water Resources Research, 43(2).
</p>
