# Figure provenance

## Figure 1

The basin maps and transects use prescribed bathymetric fields. Area-fraction distributions are calculated from those same fields. The water-column panels represent mixed, stratified and ventilating states with qualitative transfer arrows and a separate methane bubble pathway.

## Figure 2

The paired sediment-water sections show conditional light and oxygen environments. The methane-budget diagram connects sediment and water-column sources and sinks, imported methane, dissolved and free-gas inventories, and atmospheric and lateral exchange. Sediment gas cavities and water-column bubbles have distinct schematic forms.

## Figure 3

Panel a diagrams the framework's physical, reaction and transport states and their connection to local and area-integrated atmospheric exchange.

Panels b-d are constructed in `conceptual_series()` in `figures/publication_figure_style.py`. A time vector spans 0 to 1 with 2001 points. Each smooth pulse follows a squared-sine function over a prescribed interval; pulse start, end and height specify the intended relative timing and amplitude.

| Quantity | Prescribed construction |
| --- | --- |
| Shallower-bed temperature | Pulse from 0.15 to 0.43, height 1 |
| Deeper-bed temperature | Pulse from 0.20 to 0.62, height 0.50 |
| Reference CH4 emission anomaly | Pulse from 0.17 to 0.48, height 1 |
| Thermal-buffering emission anomaly | Pulse from 0.22 to 0.64, height 0.48 |
| Oxidant-limitation emission anomaly | Pulse from 0.22 to 0.83, height 1.22 |
| Delayed-release emission anomaly | Pulse from 0.50 to 0.92, normalized to the reference integral |
| Depth history in panel d | Piecewise-linear rise, high-water interval and recession |

Panel b combines attenuated warming with maintained methane oxidation and a prescribed smaller integrated net methane supply. Panel c describes a larger, more persistent emission response when oxidation loss outweighs a smaller production increase. Panel d separates retention from net supply: its matched emission integrals imply equal total excess release, while their time-dependent difference defines temporary storage.

Cumulative curves are trapezoidal integrals of the corresponding emission anomalies. The generator checks that panel d's excess inventory remains non-negative and returns to zero. CSV files retain the underlying numerical arrays. `figures/source_baseline/` contains the reference geometry and pulse definitions used by the regression checks.

The model program computes its own heat, oxygen, reaction and exchange trajectories through coupled state equations. Those outputs are stored separately in `outputs/`. Figure 3's prescribed hypotheses guide matched experiments and model comparisons described in `evaluation-plan.md`.
