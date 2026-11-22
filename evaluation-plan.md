# Prospective evaluation of the depth-state framework

The accompanying Supplementary Information provides the design narrative and selected source-code excerpts; this document expands the empirical comparison criteria.

The proposed tests connect water-depth history to physical mediators, gross reactions, gas inventories and atmospheric exchange. The accompanying archive supplies an illustrative, uncalibrated implementation. Field and experimental observations would support the forecast comparisons and mechanism tests below.

## Forecast comparisons

All candidates share external forcing, site properties and inundated area. Candidate labels specify the information represented by each model; mathematical form and calibration are determined for the study system.

| Candidate | Added representation | Question |
| --- | --- | --- |
| M0 | Area, external forcing and site properties | How well do these baseline predictors describe the chosen gas response? |
| M1 | Static area-weighted local-depth distribution | Does bathymetric geometry improve prediction? |
| M2 | Local depth history and rate of change | Does the water-level trajectory improve on static geometry? |
| M3 | Transitions in hydrodynamic state | Does representation of mixing, stratification and ventilation improve on depth history alone? |

Use common training records, withheld sites or events, forecast horizons and scoring rules. Fit preprocessing and thresholds using training data. At the forecast origin, each candidate receives the same specification of projected forcing and the information available at that time. M3 propagates future states from those inputs; subsequent density and oxygen observations provide evaluation targets. Withhold complete events or cycles where adjacent observations share thermal or storage histories.

Prespecify primary endpoints and improvement criteria. Endpoints include gas-specific flux trajectories, cumulative exchange through perturbation and recovery, and independently observed gas inventories. Report absolute error, bias, timing error and uncertainty. Retain complexity where it yields reproducible gains on withheld cases or improves an independently evaluated budget. Numerical tolerances follow measurement uncertainty and the intended decision.

## State-conditioned diagnostic

A separate diagnostic, previously denoted M4, compares predictions with and without depth after conditioning on measured bed temperature, oxygen, turbulence, electron acceptors and gas inventories. Its contemporaneous observations address residual information in depth, whereas M0-M3 evaluate forecasts from information available at the forecast origin. Evaluate these two tasks separately.

Loss of depth-associated skill after conditioning is consistent with depth acting as a proxy for measured states. Residual skill can reflect an omitted mechanism, incomplete observations, measurement error or misspecification. Both outcomes guide model refinement. Mechanism attribution additionally draws on interventions and independent budgets. Depth history can still be needed to predict future states when contemporaneous fluxes are well explained by measured conditions.

## Matched process tests

### Thermal transmission

Expose common sediment to identical surface-energy forcing under matched optical properties and a verified well-mixed state. Vary depth and forcing duration, and measure water-column and bed temperature. The conditional prediction is weaker, later bed-temperature variation in deeper water when forcing is brief relative to the thermal response time. An absent or reversed pattern beyond prespecified uncertainty challenges this prediction for the tested conditions. Stratification, inflow and optical attenuation enter the interpretation through their measured effects on heat transfer.

Replay the measured bed-temperature trajectories in separate incubations with resources and oxygen controlled. Gross production and consumption distinguish intrinsic process sensitivity from thermal transmission and reveal compensation beneath a muted net flux.

### Light, oxidants and resources

Manipulate bed illumination and mixing independently at matched temperature. Resolve photosynthetic oxygen supply, vertical renewal and respiratory demand using profiles and independent production and consumption estimates. Reduced illumination of a photosynthetically active bed is expected to lower local oxygen production. Maintained methane oxidation after this reduction challenges oxygen limitation of that sink, with transport and alternative pathways included in its budget. Concurrent measurements of photochemical CH4 formation resolve the source response to illumination.

Resource additions distinguish equal areal loading from equal concentration as water volume changes. For CO2, measure DIC, pH and alkalinity alongside carbon processing and atmospheric exchange. For N2O, resolve formation and reduction separately. For AOM, test the spatial and temporal overlap of methane, alternative electron acceptors and active organisms. Independent carbon and electron-acceptor balances quantify its contribution to whole-system consumption.

For fresh-carbon additions, combine isotope-labelled inputs with matched unamended controls under each depth history. Quantify native-carbon mineralization relative to controls, correcting for residual added carbon and source mixing. Isotope-resolved CO2 and CH4 production, oxidation, storage and release distinguish priming from mineralization of the added material. Cross carbon addition with independently controlled oxidant renewal and matched bed-temperature trajectories to test effects of delivery and redox conditions.

Lake-sediment evidence comes from Bartosiewicz et al. (2024), https://doi.org/10.1111/gcb.17076. Aaltonen et al. (2022), https://doi.org/10.1002/saj2.20419, reported a peat-column null response; Krüger et al. (2025), https://doi.org/10.1016/j.soilbio.2024.109682, compared carbon and oxygen additions. These studies motivate a surface-water-depth test using intact spatial relationships. Model evaluation of priming and isotope partitioning would require the extensions identified in the README.

### Storage and release

Compare filling-drawdown cycles with stable-depth controls in gas-bearing sediment, matching temperature and characterizing production and dissolution. Monitor water depth and atmospheric pressure separately. Estimate release thresholds using observations independent of the flux series used for testing. Pressure-triggered release predicts a pulse after threshold crossing and corresponding inventory depletion after other budget terms are accounted for. Absence of these responses challenges the release mechanism over the tested pressure range. Temperature, sediment structure and conduit history define the comparison domain.

Follow complete cycles and recovery. Retain production, consumption, dissolved and free-gas storage, lateral exchange and atmospheric release as distinct control-volume terms. Later release can offset lower flux during high water. A lasting change in cumulative release corresponds to a change in the integrated budget or final inventory.

## Joint depth-by-forcing tests

The principal comparison is a factorial contrast: perturbed minus unperturbed gas exchange under each depth history. Depth-only controls follow the same trajectory under baseline forcing. Comparing these differences across depth histories identifies how depth modifies the response to warming or resource addition. Prespecify the perturbation, observation window, inventory endpoints and uncertainty needed to resolve the interaction.

Compare gross production and consumption responses while retaining absolute rates and consumption fractions as separate endpoints. Similar net supply can accompany different internal turnover or compensating source and sink changes. Storage and boundary-exchange budgets connect net reaction to atmospheric flux.

Where isotopes or substrate-specific measurements resolve pathway contributions, test whether warming shifts methanogenesis or the balance of aerobic and anaerobic consumption. Replay temperatures with reactants controlled to distinguish intrinsic sensitivity from altered substrate and oxidant delivery. Substrate-addition evidence from Cai et al., DOI: 10.1016/j.watres.2025.124836, motivates this comparison; in situ pathway fractions require independent source partitioning. The present implementation aggregates methanogenesis and prescribes shared AOM temperature sensitivity. Pathway-resolved evaluation would extend those reaction functions using independent observations.

### Thermal attenuation and oxidation loss

With matched sediment, optical properties and forcing, quantify bed temperatures and gross methane production and consumption. Maintained oxidant supply and a similar oxidation fraction define the thermal-attenuation comparison. With escape pathways and imported methane accounted for, a smaller production increase is expected to yield a smaller early emission response. Complete-cycle exchange also depends on integrated production, consumption and final storage.

Independently vary mixing or illumination to test whether oxidant depletion reduces consumption enough to outweigh a smaller production increase in the thermally buffered case. Include alternative electron acceptors and AOM where relevant. Restoring oxidant renewal at the same thermal trajectory is predicted to reduce amplification through renewed oxidation. Gross-process budgets distinguish this response from suppressed production, including illumination-dependent source changes. Maintained consumption together with an independently resolved alternative source of excess release would reject oxidation loss as its explanation. Temperature replay with oxygen and resources controlled separately tests the thermal contribution.

### Retention and consumption

At matched gross production and lateral inputs, longer residence changes the opportunity for oxidation before escape. Attribution of lower cumulative release to consumption requires greater integrated oxidation, with final dissolved and free-gas inventories included in the budget. A flux deficit balanced by accumulation and later depletion indicates delay. Equal cumulative release at matched net input and equal initial and final stores follows from conservation; independent process observations test the ecological mechanisms governing those terms.

### Resource-dependent and cross-gas responses

Prespecify equal areal loading or equal concentration. Evaluate the perturbation-minus-control response under each depth history. Temperature-by-resource studies motivate factorial tests of depth interactions. For CO2, DIC, pH and alkalinity distinguish carbon processing from carbonate partitioning. For N2O, quantify formation and reduction to N2. The same oxidant field can affect each gas differently, so gas-specific budgets define their corresponding response predictions.

## Light-dependent methane formation and consumption

A no-plant comparison can vary illumination while matching sediment temperature and resource supply. Measure incident and bed irradiance, oxygen, methane formation and oxidation. Controls addressing microbial activity and reactive oxygen species evaluate an abiotic source while accounting for treatment-induced substrate and optical changes. These contrasts separate interfacial photochemistry, heating and photosynthetic oxygen supply. Water-column attenuation and methyl-donor availability constrain the depth-dependent source rate. Wu et al., DOI: 10.1016/j.watres.2026.127049, motivates this pathway; measurements determine its contribution under each optical and hydrodynamic state.

## Interpretation of Figure 3

The dimensionless curves are independently constructed schematic hypotheses. Grey shading marks a common warming interval. Temperature and emission anomalies are relative to an unperturbed control under each corresponding depth history. Panel d shows depth change relative to its initial value. The lower panels integrate excess atmospheric methane release.

Panels b and c use the same schematic bed-temperature contrast under different oxidant conditions. In b, maintained oxidation fraction and prompt exchange allow attenuated production to yield a smaller early response. The displayed lower cumulative release additionally assumes lower integrated net supply. In c, prescribed oxidation loss and delayed oxygen recovery produce a larger, longer response despite thermal buffering.

Panel d compares stable depth with filling and drawdown at matched excess net input. The delayed release curve has the same integral as the reference. Its implied excess inventory, obtained by integrating net input minus release, remains non-negative and returns to zero. This construction isolates retention at fixed consumption; pressure-threshold dynamics are evaluated through the storage-and-release test above.

The editable figure source and construction checks are maintained as manuscript artwork assets. They verify the internal consistency of the schematic. Process measurements, inventories and exchange observations provide the separate empirical evaluation. Numerical model files and tests retain their reviewed definitions.
