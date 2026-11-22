# DS-GHG model implementation

## Purpose and scope

The accompanying Supplementary Information presents the framework design, implementation assumptions, empirical evaluation and selected source-code excerpts. Read that note for the scientific design, then use this archive for the complete plain-text Python source and reproducible inputs. Documentation release 1.0.13 preserves the numerical implementation and aligns the detailed equations commentary with its unit conversions and hydrodynamic scaling.

This archive contains an illustrative, uncalibrated two-layer implementation of the depth-state-greenhouse-gas (DS-GHG) framework. It makes the assumptions in equations (2) and (3) explicit and supports the design of empirical tests. Outputs describe the behaviour of the prescribed equations; environmental application requires calibration and independent observations.

The model couples upper and lower water layers to active sediment. It represents depth-, optics- and mixing-dependent heat transfer; oxygen exchange and demand; simplified CO2, CH4 and N2O reactions; dissolved-gas storage and exchange; and pressure-sensitive CH4 ebullition. Aerobic methane oxidation occurs in sediment and water, while nitrate-, sulfate- and iron-dependent pathways represent anaerobic oxidation of methane (AOM). A phenomenological index represents mixed, diel-stratified, persistently stratified and ventilating states.

The following assumptions define the reduced implementation.

| Component | Implemented treatment | Extension for application |
| --- | --- | --- |
| Hydrology and mixing | Prescribed depth trajectories and a phenomenological stratification index | Catchment water balance, density-resolved mixing and horizontal transport |
| Light and productivity | Optical attenuation modifies sediment heating; productivity is prescribed | Light-driven photosynthesis, interfacial photochemical CH4 formation and oxic water-column CH4 production |
| Resources and legacy | Resource supply is prescribed and legacy state H is fixed | Moving redox fronts, mineral-associated carbon mobilization through iron reduction, mineral protection and dynamic legacy effects |
| Carbon cycling | One aggregate methanogenesis rate and simplified CO2 processing | Substrate-resolved methanogenesis, carbonate speciation, fresh-carbon priming and isotope-resolved source partitioning |
| Methane consumption | Separate electron-acceptor limitation terms share Q10_AOM; sediment and water aerobic oxidation share Q10_AeMO | Independently constrained pathway-specific temperature responses |
| Nitrogen and vegetation | Simplified sediment N2O production and reduction in no-plant scenarios | Water-column N2O transformations, plants and rhizosphere processes |
| Pressure | Fixed atmospheric reference pressure and water-level-driven hydrostatic changes | Time-varying barometric forcing |

The framework's total reaction term R_rxn accommodates biological and abiotic pathways. The code represents its biological subset. Independent rate measurements would constrain additional light-dependent sources and state-dependent kinetics. Field profiles and held-out observations would constrain mixing, vertical exchange and reaction parameters.

## Archive contents

- `depth_ghg_model.py`: model implementation and command-line runner.
- `parameters.csv`: illustrative parameter registry.
- `scenarios.csv`: no-plant scenario definitions.
- `test_depth_ghg_model.py`: nine checks of numerical integrity and expected directional behaviour.
- `model-equations.md`: state variables, equations, units and implementation assumptions.
- `evaluation-plan.md`: prospective forecast comparisons, matched process tests and rejection criteria.
- `requirements.txt`: tested Python dependencies.
- `VERSION.txt`: version and execution record.
- `LICENSE`: MIT terms for reuse.

The submission archive contains source, documentation and test inputs. Generated outputs and internal audit records are retained separately.

## Relationship to the figures and evaluation plan

Figure 2 is an independently constructed conceptual schematic integrating reaction zones and the methane budget. It includes light-driven oxygen supply, interfacial photochemical formation, spatial redox zones, conditional water-column production and lateral fluxes. The table above identifies the implementation extensions associated with these pathways.

Figure 3 presents separately constructed response hypotheses. Panel a shows the framework. Panels b and c compare thermal buffering with maintained oxidation against oxidation loss that outweighs a smaller production increase. Panel d isolates delayed release at matched excess net input and equal initial and final excess inventories. The lower plots integrate emission anomalies relative to an unperturbed control under the same depth history. Cumulative differences in b and c are prescribed cases; equal cumulative release in d follows from conservation. Their environmental relevance is tested through the measurements in `evaluation-plan.md`.

The plan distinguishes prospective M0-M3 forecasts, a separate state-conditioned diagnostic and joint depth-by-forcing tests. The code supplies the reduced M3 illustration; the empirical comparisons remain future work.

## Installation

The execution record of 24 October 2026 used Python 3.12.14 and the dependencies pinned in `requirements.txt`. All nine tests and all 23 supplied scenarios completed. These checks assess numerical integrity, execution and expected directional behaviour. Create an isolated environment and install the dependencies:

```bash
python -m venv .venv
```

On Windows:

```powershell
.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

On macOS or Linux:

```bash
source .venv/bin/activate
python -m pip install -r requirements.txt
```

## Run the model and checks

```bash
python depth_ghg_model.py
python -m unittest -v test_depth_ghg_model.py
```

The model creates `outputs/timeseries.csv`, `outputs/summary.csv` and `outputs/illustrative_depth_response.png`.

The checks cover thermal attenuation with depth, intrinsic temperature dependence of methanogenesis, electron-acceptor-dependent AOM, non-negative states and process integrals, conservation of layer depth, stratification under weak mixing, renewed exchange during ventilation and enhanced bubble release during drawdown. Environmental evaluation uses independent process, inventory and exchange observations, with matched controls and rejection criteria specified in `evaluation-plan.md`.

## Availability and reuse

The source is distributed under the MIT License. The versioned repository release at https://github.com/geochem1222/ds-ghg-framework/releases/tag/v1.0.13 retains the reviewed numerical model and inputs, together with the adopted figure scripts.
