# DS-GHG model equations and implementation map

## 1. Scope and units

The executable depth-state-greenhouse-gas (DS-GHG) model is a reduced-order two-layer water-column and active-sediment system. Time is in days, water depth is in metres, concentrations are in µmol l⁻¹ (numerically mmol m⁻³), areal rates are in mmol m⁻² d⁻¹, and temperature is in °C. `parameters.csv` is the authoritative parameter registry.

An external benchmark for the CH₄ configuration is the qualitative pathway separation reported by [Shen and Raymond (2026)](https://doi.org/10.1029/2026GL124792) under matched lake geometry, forcing and trophic assumptions: shallow-zone ebullition, greater deep-water oxidation and storage before turnover. The reported numerical depth boundaries provide comparison targets. Plant and rhizosphere processes offer an extension following the planted/unplanted and 13C-partitioning logic of [Yan et al. (2022)](https://doi.org/10.1016/j.rhisph.2021.100455).

The state vector is

\[
\mathbf{x}=\{T_e,T_h,T_s,O_{2,e},O_{2,h},NO_{3,h},SO_{4,h},Fe(III),
CH_{4,e},CH_{4,h},CH_{4,s},CO_{2,e},CO_{2,h},CO_{2,s},N_{2}O_e,N_{2}O_h,N_{2}O_s,B,S\},
\]

augmented by cumulative flux and gross-process diagnostics. Subscripts \(e\), \(h\) and \(s\) denote the epilimnion, hypolimnion and active sediment, respectively. \(B\) is sediment bubble-phase CH₄ inventory and \(S\in[0,1]\) is a continuous stratification index.

The reduced model holds the conceptual legacy state \(H\) constant, absorbing its effects into illustrative maximum rates and temperature-response parameters. Site-specific observations would constrain these parameters for application. A later comparison could test fixed parameters against parameters conditioned on independently measured antecedent hydrology, substrate state and microbial capacity.

## 2. Thermal transfer

Water-depth trajectories are prescribed independently of a water balance, and the optical term below modifies heat transfer. Extensions cover light-driven oxygen production, interfacial photochemical methane formation, spatially resolved redox fronts, oxic water-column methane production and lateral boundary fluxes. Figure 2 integrates these pathways in a separately constructed conceptual schematic.

The surface box receives atmospheric forcing and exchanges heat with the lower box:

\[
\tau_e=\frac{\rho_w c_w d_e}{86400\,h_{eff}},\qquad
\frac{dT_e}{dt}=\frac{T_{eq}(t)-T_e}{\tau_e}+\frac{K_z(T_h-T_e)}{d_e},
\]

where \(K_z\) is the state-dependent inter-layer exchange velocity. The factor 86400 converts the heat-transfer timescale from seconds to days. The lower box and active sediment follow

\[
\frac{dT_h}{dt}=\frac{K_z(T_e-T_h)}{d_h}+\frac{T_s-T_h}{\tau_{hs}},\qquad
\frac{dT_s}{dt}=\frac{T_h-T_s}{\tau_s}+\beta_I\exp(-K_d d)I_{day}(t).
\]

This formulation represents depth-dependent amplitude and phase through thermal inertia, while \(K_d\) modifies direct bed heating. It permits thermal isolation when stratification suppresses \(K_z\). A full model would replace these ordinary differential equations with water and sediment heat equations and continuity of interfacial heat flux.

## 3. Oxygen balance

Atmospheric exchange and prescribed production enter the epilimnion. Sediment oxygen demand and benthic oxidation draw from the hypolimnion. The layer balances are

\[
d_e\frac{dO_{2,e}}{dt}=k_{2,O_2}m(O_{2,sat}-O_{2,e})+P_{water}-d_eR_{w,e}-2R_{AeMO,e}-K_z(O_{2,e}-O_{2,h})+Q_{rise},
\]

\[
d_h\frac{dO_{2,h}}{dt}=K_z(O_{2,e}-O_{2,h})-d_hR_{w,h}-SOD-2R_{AeMO,h,s}-2R_{nit}.
\]

`Pwater` is prescribed externally and is zero in the first no-plant scenario set. \(Q_{rise}=(O_{2,sat}-O_{2,e})\max(\dot d,0)\) represents dilution or renewal during rising water. Falling water exports the epilimnetic concentration and therefore does not directly concentrate O₂.

Potential respiration and sediment oxygen demand are multiplied by \(O_2/(K_{O_2,cons}+O_2)\). This smooths the transition to anoxia and prevents the numerical model from consuming unavailable oxygen; the underlying carbon mineralization continues through the anaerobic terms.

## 4. Carbon dioxide

Organic-carbon mineralization is

\[
R_{OM}=V_{OM}Q_{10,OM}^{(T_s-20)/10}f_{labile}.
\]

Sediment CO₂ production includes mineralization, aerobic and anaerobic methane oxidation, and an effective methanogenic stoichiometric term. Sediment-water, inter-layer and air-water exchange are

\[
F_{swi,CO_2}=k_{swi,CO_2}\max(C_s-C_h,0),\quad
F_{z,CO_2}=K_z(C_h-C_e),\quad
F_{air,CO_2}=k_{air,CO_2}(C_e-C_{eq}).
\]

The model keeps gross sediment production, dissolved storage and atmospheric flux separate. A calibrated carbonate-system module using DIC, pH and alkalinity is required before field prediction.

## 5. Methane production and aerobic oxidation

In the Perspective, \(R_{\mathrm{rxn}}\) denotes total biological and abiotic
reactions. This implementation represents biological methane production and
oxidation. A photochemical source module, motivated by illuminated-interface
experiments ([Wu et al.](https://doi.org/10.1016/j.watres.2026.127049)), would
use local irradiance, substrate information and independent formation and
consumption measurements to constrain its rate.

Methanogenesis is

\[
R_{meth}=V_{meth}Q_{10,meth}^{(T_s-20)/10}f_{labile}f_{anox}(O_2)
\prod_{a\in\{NO_3,SO_4,Fe\}}\frac{K_{i,a}}{K_{i,a}+a}.
\]

Aerobic oxidation in medium \(j\in\{e,h,s\}\) is

\[
R_{AeMO,j}=V_{AeMO,j}Q_{10,AeMO}^{(T_j-20)/10}
\frac{CH_4}{K_{CH_4}+CH_4}\frac{O_2}{K_{O_2}+O_2}.
\]

Separating sediment, epilimnetic and hypolimnetic AeMO allows benthic oxygen exposure and state-dependent water-column residence time to act as distinct filters.

## 6. Anaerobic oxidation of methane

For acceptor \(a\in\{NO_3^-,SO_4^{2-},Fe(III)\}\),

\[
R_{AOM,a}=V_{AOM,a}Q_{10,AOM}^{(T_s-20)/10}
\frac{CH_4}{K_{CH_4,AOM}+CH_4}
\frac{a}{K_a+a}\frac{K_{i,O_2}}{K_{i,O_2}+O_2}.
\]

NO₃⁻ and SO₄²⁻ are replenished from measured boundary concentrations through an effective exchange coefficient. Fe(III) is a finite solid inventory with an optional regeneration term. Setting individual AOM pathways to zero supports model comparison. Independent methane and electron-acceptor budgets constrain each pathway's contribution.

This implementation uses one fixed \(Q_{10,AOM}\). Xu et al. (2026; DOI: 10.1021/acs.est.6c06040) reported potential-AOM \(Q_{10}\) values of 2.4 ± 0.2 and 3.1 ± 0.2, and activation energies of 57.2 ± 0.9 and 77.1 ± 4.3 kJ mol⁻¹, under 30 and 45 cm annual mean groundwater-table histories, respectively. The assay used homogenized 0-30 cm peat under CH₄-enriched anoxic conditions. These groundwater-history treatments motivate a fixed-versus-state-dependent kinetic comparison, with surface-water applications drawing parameter estimates from system-specific antecedent treatments and held-out tests.

## 7. Bubble storage and release

The hydrostatic nucleation criterion is represented by

\[
C_{bub}^*(d)=C_{bub,1atm}\left(1+\frac{\rho_wgd}{P_{atm}}\right)+C_{cap}.
\]

Excess dissolved CH₄ enters bubble storage over \(\tau_{nuc}\). Bubble capacity increases with depth, \(B_{cap}=B_0+B_d d\). Release combines capacity exceedance with a drawdown term,

\[
R_{rel}=\frac{\max(B-B_{cap},0)}{\tau_{rel}}+k_{draw}\frac{\max(-\dot d,0)}{d}B.
\]

Atmospheric pressure is a fixed reference parameter; time-varying barometric forcing is a framework extension.

The fraction reaching the atmosphere is \(\exp(-d/L_{diss})\); the remainder dissolves into the hypolimnion and becomes available for AeMO or inter-layer transfer. This formulation represents delayed storage and drawdown pulses. Quantitative prediction would also require bubble size distributions and sediment fracture mechanics.

## 8. Nitrous oxide

The net sediment N₂O rate is

\[
R_{N_2O,net}=y_{nit}(O_2)R_{nit}+y_{denit}R_{denit}-R_{N_2Ored}.
\]

Nitrification is oxic, whereas denitrification and N₂O reduction are inhibited by O₂. A Gaussian multiplier allows nitrification N₂O yield to peak under suboxic conditions. N₂O reduction follows independent Monod and Q₁₀ terms. The separation permits high-O₂ low-denitrification, suboxic N₂O maxima, and anoxic source-to-sink transitions.

## 9. Response measures and operational criteria

Response criteria are specified for the tested system and forcing. For any process \(i\), surface perturbation \(\Delta T_{surface}\), and depth \(d\), calculate

\[
S_i(d)=100\frac{R_i^{perturbed}(d)-R_i^{control}(d)}{R_i^{control}(d)}.
\]

This relative response is undefined when the control rate is zero and unstable near zero; an absolute difference is preferable in that case. Thermal attenuation may be summarized by gain \(G_T=A_{bed}/A_{surface}\) and phase lag for a specified forcing frequency. Gain and phase characterize continuous thermal attenuation. An operational buffering criterion uses a scientifically justified equivalence margin and adequate uncertainty bounds, fixed before testing; evidence of negligible change requires an interval contained within that margin. A regime threshold is identified by independent evidence of a change in exchange, reaction pathway or release mechanism.

## 10. Hydrodynamic state and upgrade path

The implemented reduced model splits the water column into epilimnetic and hypolimnetic boxes. Its continuous stratification index evolves from daytime surface heating, water depth, mixing, convective cooling and prescribed ventilation:

\[
\frac{dS}{dt}=k_{build}f_d(d)I_{day}\frac{\max(T_{surface}-T_h,0)}{T_*}\frac{1-S}{1+m}
-k_{decay}mS-k_{conv}\frac{\max(T_h-T_{surface},0)}{T_*}S-k_{vent}VS.
\]

Inter-layer exchange interpolates between mixed and stratified limits,

\[
K_z=\left[K_{z,strat}+(K_{z,mix}-K_{z,strat})(1-S)^2\right]\max\left(1,\frac{m}{\max(m_0,0.05)}\right).
\]

Here, \(m\) is the regime-adjusted mixing multiplier and \(m_0\) is the scenario mixing factor. The implementation constrains \(m\) to at least 0.05. A ventilation event increases \(m\), erodes \(S\), and exposes stored CH₄, CO₂ and N₂O to oxidation or air-water loss. The phenomenological index requires calibration against temperature, density and turbulence profiles for field application. A third tier resolves one-dimensional vertical heat, oxygen, solute and gas transport. Horizontal advection is added only when inflows or basin circulation invalidate the local-column assumption.

For each perturbation, report forcing-response sensitivity, maximum departure, recovery time, state-dependent hysteresis and cumulative post-forcing legacy flux. These metrics characterize the responses of individual gases and pathways.

The full upgrade sequence is:

1. Calibrate the two-layer geometry and state-transition parameters against observed profiles.
2. Discretize water and sediment into finite volumes and replace effective SWI exchange by molecular and turbulent transport.
3. Resolve millimetre-scale O₂ penetration in the upper sediment.
4. Add DIC-pH-alkalinity speciation and temperature-dependent Henry constants.
5. Split organic matter and methanogenesis into substrate pools only when independent data identify them.
6. Add explicit NO₃⁻ → NO₂⁻ → NO → N₂O → N₂ pools for isotope calibration.
7. Evaluate individual AOM pathways using independent rate and electron-acceptor budgets together with held-out predictive performance. Pathway attribution rests on process measurements and independent budgets alongside predictive skill.
8. Add horizontal transport only where independent observations demonstrate lateral source influence.
