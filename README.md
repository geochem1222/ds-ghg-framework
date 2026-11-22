# Depth-state-greenhouse-gas framework

This repository accompanies the Perspective **Water depth reshapes greenhouse-gas responses to global change in freshwater systems**. It contains an illustrative two-layer sediment-water model, its input parameters and numerical checks, and the scripts used to construct the three conceptual figures.

Release **v1.0.13** provides a fixed source version for review and reuse under the MIT License.

## Scientific design

Water-depth history changes the transmission of heat and light, oxidant supply, gas storage and exchange. The framework connects these states to gas-specific production, consumption and atmospheric release. Matched forcing comparisons across depth histories distinguish changes in response magnitude, timing and cumulative exchange.

The numerical implementation couples upper water, lower water and active sediment. It represents heat and oxygen exchange, aggregate CO2, CH4 and N2O reactions, electron-acceptor-dependent anaerobic methane oxidation and pressure-sensitive methane ebullition. Its prescribed parameters and 23 scenarios support mechanistic exploration. Calibration uses independent temperature, chemistry, reaction-rate and gas-inventory observations.

## Contents

| File or directory | Contents |
| --- | --- |
| `depth_ghg_model.py` | Coupled model and scenario runner |
| `test_depth_ghg_model.py` | Nine numerical and directional checks |
| `parameters.csv` | 103 parameter records with units and indicative ranges |
| `scenarios.csv` | 23 prescribed no-plant scenarios |
| `model-equations.md` | State variables, equations, units and implementation assumptions |
| `evaluation-plan.md` | Matched experiments, independent budgets and forecast evaluation |
| `docs/implementation.md` | Detailed implementation guide |
| `figures/` | Conceptual figure scripts and their underlying geometry and curve definitions |
| `docs/figure-provenance.md` | Origin and interpretation of the figure elements |
| `requirements.txt` | Recorded dependency versions |
| `LICENSE` | MIT License |

## Installation and execution

The recorded execution environment used Python 3.12.14. Create and activate a local environment, then install the pinned dependencies:

```bash
python -m venv .venv
# Activate .venv in your shell.
python -m pip install -r requirements.txt
python -m unittest -v test_depth_ghg_model.py
python depth_ghg_model.py
```

The model writes time series, scenario summaries and an illustrative response plot to `outputs/`. The `--no-plot` option selects CSV output.

Generate the manuscript figures from the repository root:

```bash
python figures/generate_current_figures.py
# Or select one figure:
python figures/generate_current_figures.py --figure 3
```

The figure scripts export editable-text SVG, PDF, PNG and TIFF, together with layout and curve checks. Arial, Helvetica or DejaVu Sans are used in that order according to font availability. Numerical arrays retain the supplied definitions across font choices.

## Figure 3 and the numerical model

Figure 3a maps the framework's state and flux relationships. Panels b-d are prescribed, dimensionless response hypotheses. Smooth pulse functions define the temperature and emission trajectories, and trapezoidal integration gives cumulative excess emissions. Panel d normalizes delayed release to the same integrated amount as its reference, with non-negative temporary storage and equal final inventories. These plots express the conditional contrasts to be tested.

The two-layer model has its own scenario calculations and outputs in `outputs/`. Model experiments and conceptual figure generation are separate execution paths. [Figure provenance](docs/figure-provenance.md) explains the construction in detail.

## Version and citation

The model, numerical tests and CSV inputs retain the bytes of the reviewed supplementary archive, version 1.0.13. The release also includes the adopted figure sources. Use the fixed release or its commit identifier when referring to this code:

<https://github.com/geochem1222/ds-ghg-framework/releases/tag/v1.0.13>

`CITATION.cff` supplies the software citation metadata. A DOI archive can preserve this same release for the final publication citation.
